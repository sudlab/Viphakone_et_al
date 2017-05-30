'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import pyximport

import collections
from CGAT import GTF, IndexedFasta, IOTools
from iCLIP import TranscriptCoordInterconverter
from CGAT.Genomics import complement

import pysam

pyximport.install()
from _simple_align import simple_align


def convert_alignment(read, converter):
    
    try:
        pos =  converter.genome2transcript((read.reference_start,
                                            read.reference_end))
        if converter.strand == "-":
            pos = (pos[1] + 1 , pos[0] + 1)
        return pos
    
    except ValueError as e:
        if "not in transcript" in e.message:
            return None
        else:
            raise
            

def format_alignment(a, l, unaligned_is_5p):
    ''' a is an alignement of structure (score, [start]), (read.start, read.end)
        l is the length of the unaligned segment that has been mapped'''
    
    alignment, read_pos = a
    a_score, a_start = alignment
    
    results = ((a_start[0], a_start[0] + l), read_pos)
    if not unaligned_is_5p:
        results = results[::-1]
        
    return results


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])
    parser.add_option("-g", "--gtffile", dest="gtffile", type="string",
                      help="GTF file with transcript models to search in")
    parser.add_option("-f", "--fastafile", dest="fastafile", type="string",
                      help="Indexed Fasta file with genome sequence for"
                      "transcripts")
    parser.add_option("-j", "--junctions-file", dest="junctions", type="string",
                      help="Tabix index bed file with junctions to add on"
                      "top of provided GTF file")
    parser.add_option("-m", "--mismatches", dest="mismatches", type="float",
                      default=2,
                      help="Number of mismatches to allow when aligning read"
                      "sections")
    parser.add_option("--minimum-unaligned", dest="min_unaligned", type="int",
                      default=10,
                      help="Minimum number of unaligned bases before"
                      "considering read for realignment")
    parser.add_option("-d", "--minimum-distance", dest="min_dist", type="int",
                      default=5,
                      help="Minimum distance between end of primary alignment"
                      "and novel alignment. This filters out short indels")
    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")
    parser.add_option("-p", "--out-prefix", dest="out_prefix", type="string",
                      default=None,
                      help="Prefix for output files, if not specified bam come out of stdout"
                      "and all others will use input bam name")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args[0])
    fasta = IndexedFasta.IndexedFasta(options.fastafile)
#    fasta.setConverter(IndexedFasta.getConverter("zero-both-open"))
    junctions = pysam.TabixFile(options.junctions, parser=pysam.asBed())

    stats = collections.Counter()
    inner_dist_dist = collections.Counter()
    xl_dist_dist = collections.Counter()
    gene_counts = collections.Counter()

    if options.out_prefix:
        out_bam_name = options.out_prefix + ".bam"
    elif not options.stdout == sys.stdout:
        options.stdout.close()
        out_bam_name = options.stdout.name
    else:
        out_bam_name = "-"
    out_bam = pysam.AlignmentFile(out_bam_name, "wb", template=bamfile)

    if not options.out_prefix:
        if args[0].endswith(".bam"):
            options.out_prefix = args[0][:-4]
        else:
            options.out_prefix = args[0]

    for gene in GTF.gene_iterator(GTF.iterator(IOTools.openFile(
            options.gtffile))):
        
        transcripts = list(gene)
        transcripts = [[exon for exon in transcript if exon.feature == "exon"]
                       for transcript in transcripts]
        start = min(exon.start for transcript in gene for exon in transcript)
        end = max(exon.end for transcript in gene for exon in transcript)
        contig = transcripts[0][0].contig
        strand = transcripts[0][0].strand

        gene_exon = GTF.Entry().fromGTF(transcripts[0][0])
        gene_exon.start = start
        gene_exon.end = end
        gene_exon.transcript_id = "gene"

        transcripts.append([gene_exon])

        reads = list(bamfile.fetch(contig, start, end))
        if len(reads) == 0:
            continue

        converters = [TranscriptCoordInterconverter(transcript)
                      for transcript in transcripts]

        transcript_sequences = [GTF.toSequence(transcript, fasta).upper()
                                for transcript in transcripts]
        transcript_junctions = junctions.fetch(contig, start, end, strand)
        junction_starts = [j.start for j in transcript_junctions]
        junction_ends = [j.end for j in transcript_junctions]

        for read in reads:

            if sum(stats.itervalues()) % 1000 == 0:
                E.debug("Analysed %i reads" % sum(stats.itervalues()))

            # Basic read filters ###
            if not read.is_reverse == (strand == "-"):
                # read is on different strand to gene
                stats["wrong strand"] += 1
                continue

            if not (read.reference_start > start and read.reference_end < end):
                stats["not in gene"] += 1
                continue

            if read.is_secondary or read.is_supplementary:
                stats["not primary"] += 1
                continue

            # Find any unaligned sequence in the read ###
            unaligned_start, unaligned_seq, unaligned_is_5p = \
                get_softclipped(read, options.min_unaligned)

            if unaligned_start is None:
                stats["no softclipped"] += 1
                continue
            
            if not float(options.mismatches).is_integer():
                mm_thresh = options.mismatches * len(unaligned_seq)
            else:
                mm_thresh = options.mismatches

            if (read.reference_end in junction_starts and
                 unaligned_start != 0) or \
               (read.reference_start in junction_ends and
                 unaligned_start == 0):
                stats["junction"] += 1
                continue

            # Get the transcripts and transcript coordinates ###
            aligned_seq_pos, filtered_sequences, filtered_converters = \
                get_read_transcripts(read, transcript_sequences,
                                     converters, start, end, strand)

            # Do the alignments and format results ###
            alignments = [simple_align(s, unaligned_seq.upper())
                          for s in filtered_sequences]

            # add coordinates
            alignments = zip(alignments, aligned_seq_pos)

            # Filter the alignments ###
            min_mm = min(a[0][0] for a in alignments)

            if min_mm > mm_thresh:
                # no good alignments
                stats["No good alignment"] += 1
                continue

            filtered_converters = [c for c, a in zip(filtered_converters,
                                                     alignments)
                                   if a[0][0] == min_mm]
            alignments = [a for a in alignments if a[0][0] == min_mm]

            if any([len(a[0][1]) > 1 for a in alignments]):
                stats["Multiple best hits"] += 1
                continue

            # reformat the alignments in start finish pairs in read order.
            f_alignments = [format_alignment(a, len(unaligned_seq),
                                             unaligned_is_5p)
                            for a in alignments]

            # Find index of alignment with minimum inner distance
            best_inner_dist, best_alignment_index = \
                select_alignement(f_alignments, options.min_dist)

            if best_inner_dist < options.min_dist:
                stats["probable indel"] += 1
                continue

            stats["good"] += 1

            # use alignment with the smallest inner distance as THE alignment
            
            aligned_converter = filtered_converters[best_alignment_index]
            alignment = alignments[best_alignment_index]
            f_alignment = f_alignments[best_alignment_index]

            # This read pair contains two cross links. One is at the
            # 5' end of the 5' fragment (f_alignment[0]) the other is
            # somewhere in the 3' fragment, (f_alignment[1]). best
            # guess - the centre.

            xl_dist = int(f_alignment[1][1] + f_alignment[1][0])/2 - \
                f_alignment[0][0]

            new_read = build_read(alignment,
                                  len(unaligned_seq),
                                  aligned_converter,
                                  unaligned_start,
                                  best_inner_dist,
                                  xl_dist,
                                  read)

            out_bam.write(read)
            out_bam.write(new_read)
            gene_counts[transcript[0].gene_id] += 1
            inner_dist_dist[best_inner_dist] += 1
            xl_dist_dist[xl_dist] += 1
    
    inner_dist_fn = options.out_prefix + ".inner_distance_distribution.tsv.gz"
    xl_dist_fn = options.out_prefix + ".crosslink_distance_distribution.tsv.gz"
    gene_counts_fn = options.out_prefix + ".gene_counts.tsv.gz"

    IOTools.writeLines(inner_dist_fn, sorted(inner_dist_dist.iteritems()),
                       ["Distance", "Count"])
    IOTools.writeLines(xl_dist_fn, sorted(xl_dist_dist.iteritems()),
                       ["Distance", "Count"])
    IOTools.writeLines(gene_counts_fn, sorted(gene_counts.iteritems()),
                       ["gene_id", "Count"])

    stat_lines = [map(str, line) for line in stats.iteritems()]
    stat_lines = sorted(stat_lines)
    stat_lines = ["\t".join(line) for line in stat_lines]
    E.info("\n".join(stat_lines))

    # write footer and output benchmark information.
    E.Stop()


def select_alignement(alignments, min_dist):
    '''Returns the index of the alignment with the minimum
    inner distance. Raises NoGoodAlignment if this is less than
    min_dist'''

    inner_dist = [abs(a[1][0]-a[0][1]) for a in alignments]
    min_inner_dist = min(inner_dist)
    min_dist_i = inner_dist.index(min_inner_dist)

    return min_inner_dist, min_dist_i


def get_softclipped(read, min_unaligned):
    '''Finds any parts of the read where the there are more than
    min_unaligned bases are softclipped and returns the start of the
    unaligned sequence on the read, the sequenced of the unaligned
    part of the read and whether this is on the 5' end of the read
    **in physical read coordinates**.

    unaligned sequence on the 5' end of the read is prioritised.

    Returns None, None, None if sufficient unaligned sequence is not
    found'''

    unaligned_is_5p = None
    unaligned_start = None
    unaligned_seq = None

    if read.is_reverse:
        
        if read.query_length - read.query_alignment_end > min_unaligned:
            unaligned_is_5p = True
            unaligned_start = read.query_alignment_end
            unaligned_seq = complement(read.seq[read.query_alignment_end:])
        elif read.query_alignment_start > min_unaligned:
            unaligned_is_5p = False
            unaligned_start = 0
            unaligned_seq = complement(read.seq[:read.query_alignment_start])
    else:
        if read.query_alignment_start > min_unaligned:
            unaligned_is_5p = True
            unaligned_start = 0
            unaligned_seq = read.seq[:read.query_alignment_start]
        elif read.query_length - read.query_alignment_end > min_unaligned:
            unaligned_is_5p = False
            unaligned_start = read.query_alignment_end
            unaligned_seq = read.seq[read.query_alignment_end:]

    return unaligned_start, unaligned_seq, unaligned_is_5p


def get_read_transcripts(read, transcript_sequences, converters,
                         start, end, strand):
    '''Find which transcripts from the gene the alignment is contained in,
    return a list of the positions in these transcripts, the sequences of
    these transcripts, and the converters for those transcripts'''

    aligned_seq_pos = [convert_alignment(read, converter)
                       for converter in converters]

    assert aligned_seq_pos[-1] is not None

    filtered_sequences = [t for t, c
                          in zip(transcript_sequences, aligned_seq_pos)
                          if c is not None]
    filtered_converters = [c for c, p
                           in zip(converters, aligned_seq_pos)
                           if p is not None]
    aligned_seq_pos = [x for x in aligned_seq_pos if x is not None]

    return aligned_seq_pos, filtered_sequences, filtered_converters


def build_read(alignment, length, converter, unaligned_start,
               inner_dist, xl_dist, old_read):
    '''Build a new read to represent the new supplimentary alignment.'''

    new_read = pysam.AlignedSegment()
    new_read.query_name = old_read.query_name
    new_read.seq = old_read.seq
    new_read.query_qualities = old_read.query_qualities
    new_read.flag = old_read.flag
    new_read.is_supplementary = True
    new_read.reference_id = old_read.reference_id

    if old_read.is_reverse:
        new_read.reference_start = \
                    converter.transcript2genome(alignment[0][1]) - length + 1
    else:
        new_read.reference_start = converter.transcript2genome(alignment[0][1])
    
    cigar = list()
    if unaligned_start > 0:
        cigar.append((4, unaligned_start))

    cigar.append((0, length))

    if (length + unaligned_start) < old_read.query_length:
        cigar.append((4, old_read.query_length - (length + unaligned_start)))

    cigar = tuple(cigar)
    new_read.cigartuples = cigar
    
    new_read.template_length = old_read.template_length

    new_read.tags = (("XL", int(xl_dist)),
                     ("XI", int(inner_dist)),
                     ("NM", int(alignment[0][0])))

    return new_read


if __name__ == "__main__":
    sys.exit(main(sys.argv))

