import pandas as pd
import numpy as np
import CGAT.IOTools as IOTools
from CGATPipelines.Pipeline import cluster_runnable
import CGAT.Database as DUtils
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.FastaIterator as FastaIterator
import iCLIP as iCLIP
import CGAT.Intervals as Intervals
import collections
import pysam
import re
import xml
import os
import itertools

AMBIGUITY_CODES = {'M': 'AC',
                   'R': 'AG',
                   'W': 'AT',
                   'S': 'CG',
                   'Y': 'CT',
                   'K': 'GT',
                   'V': 'ACG',
                   'H': 'ACT',
                   'D': 'AGT',
                   'B': 'CGT',
                   'N': 'CGAT'}

PARAMS = {}

def IUPAC2Regex(sequence):

    for code, regex in AMBIGUITY_CODES.iteritems():
        sequence = re.sub(code, '[%s]' % regex, sequence)

    return sequence


@cluster_runnable
def normalizeIndevidualProfilesToRNASeq(clip_profile_file,
                                        rnaseq_profile_file,
                                        outfile_matrix,
                                        outfile_summary,
                                        pseduo_count=1):

    clip_profile = pd.read_csv(IOTools.openFile(clip_profile_file),
                               sep="\t",
                               index_col=0)
    rnaseq_profile = pd.read_csv(IOTools.openFile(rnaseq_profile_file),
                                 sep="\t",
                                 index_col=0)

    rnaseq_profile = rnaseq_profile + pseduo_count

    normalised_profile = clip_profile/rnaseq_profile
    
    normalised_profile = normalised_profile.apply(lambda x:
                                                  x/x.sum(),
                                                  axis=1)

    average_profile = normalised_profile.sum()

    normalized_average_profile = average_profile/average_profile.sum()

    normalised_profile.to_csv(IOTools.openFile(outfile_matrix, "w"),
                              sep="\t",
                              index_label="Transcript")
    normalized_average_profile.name = "profile"
    normalized_average_profile.to_csv(IOTools.openFile(outfile_summary, "w"),
                                      header=True,
                                      sep="\t",
                                      index_label="position")


@cluster_runnable
def getSingleExonProfiles(clip_profile_file,
                          rnaseq_profile_file,
                          outfile_matrix,
                          outfile_summary,
                          annotations,
                          pseduo_count=1):

    statement = ''' SELECT DISTINCT es.transcript_id as id
                    FROM exon_stats as es
                    INNER JOIN transcript_info as ti
                    ON es.transcript_id = ti.transcript_id
                    GROUP BY ti.gene_id
                    HAVING MAX(nval) = 1 '''

    single_exon_genes = DUtils.fetch_DataFrame(statement, annotations)
    single_exon_genes = single_exon_genes["id"].values

    clip_profile = pd.read_csv(IOTools.openFile(clip_profile_file),
                               sep="\t",
                               index_col=0)

    single_exon_clip_profiles = clip_profile.loc[single_exon_genes]
    
    rnaseq_profiles = pd.read_csv(IOTools.openFile(rnaseq_profile_file),
                                  sep = "\t",
                                  index_col=0)

    single_exon_rnaseq_profiles = rnaseq_profiles.loc[single_exon_genes] + pseduo_count

    normalised_profile = single_exon_clip_profiles/single_exon_rnaseq_profiles
    
    normalised_profile = normalised_profile.apply(lambda x:
                                                  x/x.sum(),
                                                  axis=1)

    average_profile = normalised_profile.sum()

    normalized_average_profile = average_profile/average_profile.sum()

    normalised_profile.to_csv(IOTools.openFile(outfile_matrix, "w"),
                              sep="\t",
                              index_label="Transcript")
    normalized_average_profile.name = "profile"
    normalized_average_profile.to_csv(IOTools.openFile(outfile_summary, "w"),
                                      header=True,
                                      sep="\t",
                                      index_label="position")
    

@cluster_runnable
def averageRegions(infile, resolutions, outfile):

    count_matrix = IOTools.openFile(infile)

    count_matrix.readline()

    outdict = collections.defaultdict(list)

    for line in count_matrix:

        fields = line.split("\t")
 
        gene_id, fields = fields[0], map(float, fields[1:])

        averages = []
        start = 0

        for i in range(len(resolutions)):
            end = start+resolutions[i]
            averages.append(np.mean(fields[start:end]))
            start = end

        outdict["transcript_id"].append(gene_id)
        outdict["upstream"].append(averages[0])
        outdict["5utr"].append(averages[1])
        outdict["CDS"].append(averages[2])
        outdict["3utr"].append(averages[3])
        outdict["downstream"].append(averages[4])

    pd.DataFrame(outdict).to_csv(IOTools.openFile(outfile, "w"),
                                 sep="\t",
                                 index=False)


def scoreRegions(summary_df):

    summary_df.drop("Protein", axis=1, inplace=True)
    summary_df = summary_df.set_index(["Replicate", "transcript_id"])

    summary_df.drop("R1", axis=0, inplace=True)
    summary_df.replace("nan", np.nan, inplace = True)
    df_average = summary_df.groupby(level="transcript_id").mean()
    df_logrank = df_average.rank().apply(np.log10)
    print df_logrank.head()
    df_score = pd.DataFrame({"cds": df_logrank["cds_enrichment"] + df_logrank["cds_count"],
                             "utr3": df_logrank["utr3_enrichment"] + df_logrank["utr3_count"],
                             "utr5": df_logrank["utr5_enrichment"] + df_logrank["utr5_count"]})

    return df_score


@cluster_runnable
def scoreCircularCandidates(outfile):

    statement = ''' SELECT * FROM profile_summaries
                    WHERE Protein='%s' '''

    chtop = DUtils.fetch_DataFrame(statement % "Chtop", "csvdb")
    alyref = DUtils.fetch_DataFrame(statement % "Alyref", "csvdb")
   
    chtop_score = scoreRegions(chtop)
    alyref_score = scoreRegions(alyref)

    clip_score = chtop_score["utr3"] + alyref_score["utr5"]

    clip_score.to_csv(IOTools.openFile(outfile, "w"),
                      sep = "\t",
                      index_label=True)


   
@cluster_runnable
def calculateRegionEnrichments(bamfile, gtffile, outfile):

    import pysam
    import iCLIP

    bam = pysam.Samfile(bamfile)
    outlines = []

    for transcript in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(gtffile))):

        regions = pd.Series()
        exons = GTF.asRanges(transcript, "exon")
        regions["cds"] = GTF.asRanges(transcript, "CDS")

        # skip genes without cds
        if len(regions["cds"]) == 0:
            continue

        utrs = Intervals.truncate(exons, regions["cds"])

        cds_start, cds_end = regions["cds"][0][0], regions["cds"][-1][1]

        if transcript[0].strand == "+":
            regions["utr3"] = [x for x in utrs if x[0] >= cds_end]
            regions["utr5"] = [x for x in utrs if x[0] <= cds_start]
        else:
            regions["utr3"] = [x for x in utrs if x[0] <= cds_start]
            regions["utr5"] = [x for x in utrs if x[0] >= cds_end]

        # check that there is both a 3 and a 5' utr
        if any(regions.apply(len) == 0):
            continue

        # Do the counting:
        region_counts = regions.apply(lambda x:
                                      iCLIP.count_intervals(
                                          bam,
                                          x,
                                          transcript[0].contig,
                                          transcript[0].strand))

        region_counts = region_counts.sum(axis=1)
        region_counts = region_counts.fillna(0)

        transcript_length = sum([x[1] - x[0] for x in exons])

        region_lengths = regions.apply(lambda region: sum(
            [exon[1] - exon[0]
             for exon in region]))

        fractional_lengths = region_lengths/transcript_length
        fractional_counts = region_counts/region_counts.sum()

        enrichments = fractional_counts/fractional_lengths

        region_lengths = region_lengths.sort_index()
        region_counts = region_counts.sort_index()
        enrichments = enrichments.sort_index()

        outline = [transcript[0].transcript_id] + \
                  list(region_lengths.values) + \
                  list(region_counts.values) + \
                  list(enrichments.values) 

        outlines.append(outline)
    
    region_names = sorted(["cds", "utr5", "utr3"])
    header = ["transcript_id"] + \
             [x+"_length" for x in region_names] + \
             [x+"_count" for x in region_names] + \
             [x+"_enrichment" for x in region_names]

    IOTools.writeLines(outfile, outlines, header)


@cluster_runnable
def getTranscriptTypeByGeneID(infile, outfile):

    import CGAT.GTF as GTF
    
    outlines = []
    for transcript in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(infile))):
        outlines.append([transcript[0].gene_id,
                         transcript[0].transcript_id,
                         transcript[0].source])

    IOTools.writeLines(outfile, outlines,
                 header=["gene_id", "transcript_id", "biotype"])

    
@cluster_runnable
def getZagrosRIInputFiles(clusters_file, ri_file, bam_file,
                          target_length,
                          outfiles):
    '''This function takes in a set of clusters and retained introns and
    outputs the input files needed to run zagros, namely:

       Zagros compatible regions: Must all be the same length, and truncated
                                  where possible at the edges of the intron.
                                  
                                  For clusters longer than the specification,
                                  a section of the cluster will be taken so as
                                  to put the new cluster as close to centred
                                  over its maximum height as possible without
                                  including new sequence.

                                  For clusters shorter than the specification,
                                  clusters will be grown on either side until
                                  they hit intron boundaries. If the intron is
                                  smaller than the cluster size, then the 
                                  cluster will be centered on the centre of the
                                  intron.

      Zagros diagnostic events:   for each cluster region, a comma-seperated 
                                  list of mapping heights '''

    clusters = Bed.iterator(IOTools.openFile(clusters_file))
    ris = Bed.readAndIndex(IOTools.openFile(ri_file), with_values=True)
    bam = pysam.Samfile(bam_file)

    clusters_out_fn, des_fn = outfiles

    clusters_out = IOTools.openFile(clusters_out_fn, "w")
    des_out = IOTools.openFile(des_fn, "w")

    for cluster in clusters:

        E.debug("New cluster. Name %s" % cluster.name)
        cluster_depth = iCLIP.count_intervals(bam,
                                              [(cluster.start, cluster.end)],
                                              cluster.contig,
                                              cluster.strand)

        # skip clusters which just have overlapping edges
        if cluster_depth.sum() == 0:
            continue
        overlapping_ris = ris[cluster.contig].find(cluster.start, cluster.end)
       
        for intron in overlapping_ris:

            if not intron[2].strand == cluster.strand:
                continue

            cluster_length = cluster.end-cluster.start
            new_cluster = cluster.copy()
            intron_length = intron[1] - intron[0]
            E.debug("Cluster length = %s, intron length %s"
                    % (cluster_length, intron_length))
            while not cluster_length == target_length:
                if (cluster_length < target_length and
                   cluster_length < intron_length):
                    E.debug("Cluster is too short. Space to expand")
                    E.debug("Cluster is (%s,%s), intron is (%s,%s)" %
                            (new_cluster.start, new_cluster.end,
                             intron[0], intron[1]))
                    difference = target_length - cluster_length

                    right_shift = min((difference+1)/2,
                                      intron[1] - new_cluster.end)
                    E.debug("Shifting right boundary %s bases left"
                            % right_shift)
                    remainder = (difference+1)/2 - right_shift
                    new_cluster.end += right_shift

                    left_shift = min(difference/2 + remainder,
                                     new_cluster.start - intron[0])
                    new_cluster.start = new_cluster.start - left_shift
                    remainder = difference/2 + remainder - left_shift
                    E.debug("shifting left boundary %s bases left"
                            % left_shift)
                    new_cluster.end += min(remainder, intron[1] - new_cluster.end)
                    E.debug("shifting right boundary %s bsaes right" %
                            min(remainder, intron[1] - new_cluster.end))

                elif (cluster_length > target_length and
                      cluster_length < intron_length):
                    E.debug("cluster is too long. Intron is long enough")
                    cluster_peak = cluster_depth.idxmax()
                    intron = (new_cluster.start, new_cluster.end)
                    new_cluster.start = int(cluster_peak)
                    new_cluster.end = int(cluster_peak + 1)

                elif cluster_length >= intron_length:
                    E.debug("cluster is longer than intron")
                    intron_centre = (intron[1] + intron[0])/2
                    E.debug("intron centre is %i" % intron_centre)
                    new_cluster.start = intron_centre - target_length/2
                    new_cluster.end = intron_centre + (target_length+1)/2

                else:
                    raise ValueError(
                        "This shouldn't happen\n"
                        "cluster length is %s. intron length is %s")

                cluster_length = new_cluster.end-new_cluster.start
                E.debug("new cluster length is %s" % cluster_length)

            clusters_out.write("%s\n" % str(new_cluster))
            new_depth = cluster_depth.reindex(np.arange(new_cluster.start,
                                                        new_cluster.end),
                                              fill_value=0)
            des_out.write(','.join(map(str, new_depth.values)) + "\n")

    clusters_out.close()
    des_out.close()


@cluster_runnable
def findRegexMotifs(motifs_file, sequences, outtable, gfffile, len_thresh=0,
                    enrich_thresh=1):
    '''Take a DREME XML motifs file and a FASTA database and find the locations
    of motif matches in the database. Will extract location and gene infomation
    from the sequence file. '''

    motifs = []
    outlines = []
    outgffs = []
    
    tree = xml.etree.ElementTree.ElementTree()
    tree.parse(motifs_file)
    
    model = tree.find("model")
    num_positives = int(model.find("positives").get("count"))
    num_negatives = int(model.find("negatives").get("count"))
    
    for motif in tree.find("motifs").getiterator("motif"):
        p = float(motif.get("p"))
        n = float(motif.get("n"))
        try:
            enrichment = (p/num_positives)/(n/num_negatives)
            if enrichment < enrich_thresh:
                continue   
        except ZeroDivisionError:
            pass

        motif_seq = motif.get("seq")

        if len(motif_seq) < len_thresh:
            continue

        motifs.append(motif_seq)

    E.info("found %i motifs" % len(motifs))
    
    for sequence in FastaIterator.iterate(IOTools.openFile(sequences)):
        # search for coord-like fields in the name

        loc_pattern = "(chr[^\W_]+)[\W_]([0-9]+)[\W_]{1,2}([0-9]+)"
        strand_pattern = "[\W\_]\(?([+-])\)?"
        if re.search(loc_pattern + strand_pattern, sequence.title):
            pattern = loc_pattern + strand_pattern
            chrom, start, end, strand = re.search(
                pattern, sequence.title).groups()
            start, end = (int(start), int(end))
            name = re.sub(pattern, "", sequence.title)

        elif re.search(loc_pattern, sequence.title):
            pattern = loc_pattern
            strand = "+"
            chrom, start, end = re.search(pattern, sequence.title).groups()
            start, end = (int(start), int(end))
            name = re.sub(pattern, "", sequence.title)

        else:
            chrom, start, end, strand = sequence.title, 0, 0, "+"
            name = sequence.title

        for motif in motifs:
            for match in re.finditer(IUPAC2Regex(motif), sequence.sequence):
                
                if strand == "+":
                    match_start = start + match.start()
                    match_end = start + match.end()
                else:
                    match_end = end - (match.start() + 1)
                    match_start = end - (match.end() - 1)

                outlines.append([motif, name,
                                 "%s:%i-%i" % (chrom, match_start, match_end),
                                 strand])
                gff = GTF.Entry()
                gff.contig = chrom
                gff.start = match_start
                gff.end = match_end
                gff.feature = "motif"
                gff.source = "DREME"
                gff.strand = strand
                gff.score = "."
                gff.frame = "."
                gff.addAttribute("name", name)
                gff.addAttribute("motif", motif)
                outgffs.append(gff)

    gff_name = os.path.basename(outtable)
    IOTools.writeLines(outtable, outlines,
                 header=["motif", "sequence", "location", "strand"])
    outgffs = map(str, outgffs)
    IOTools.writeLines(gfffile, [[gffline] for gffline in outgffs],
                 header=["track name='%s' description='%s'" % (gff_name,gff_name)])


@cluster_runnable
def runGOSeq(genes, exprs, outfile, pwf_plot=None):
    '''Each of the params :param:genes and :param:exprs should be pandas
    Series with gene names as the index. :param:`genes` should be a 1 or 0
    indicator variable giving the positive genes and :param:`exprs` gives
    the expression level. Gene_ids should be ENSEMBL. Set pwf_plot to
    save the pwf fit plot as png'''

    import rpy2
    import rpy2.robjects as ro
    from rpy2.robjects import r as R
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()

    goseq = importr("goseq")

    genesv = ro.Vector(genes.values)
    genesv.names = ro.Vector(genes.index.values)

    exprs = exprs.fillna(0)

    exprsv = ro.Vector(exprs.values + 0.1)
    exprsv.names = ro.Vector(exprs.index.values)
    exprsv = R.rank(exprsv, ties_method="first")
 

    pwf = goseq.nullp(genesv, bias_data=exprsv, plot_fit=False)

    if pwf_plot:
        R.png(pwf_plot)
        goseq.plotPWF(pwf)
        R["dev.off"]()

    GO_analysis = goseq.goseq(pwf, "hg19", "ensGene")

    GO_analysis = R["data.frame"](GO_analysis,
                                  over_qvalue=R("p.adjust")(
                                      GO_analysis.rx2("over_represented_pvalue"),
                                      method="BH"))

    GO_analysis = R["data.frame"](GO_analysis,
                                  under_qvalue=R("p.adjust")(
                                      GO_analysis.rx2("under_represented_pvalue"),
                                      method="BH"))

    R["write.table"](GO_analysis, file=outfile, quote=False, sep="\t", row_names=False)


@cluster_runnable
def tRNABaseCounts(gtf_file, bam_file, outfile):

    import pandas

    outs = []
    for tRNA in GTF.gene_iterator(GTF.iterator(IOTools.openFile(gtf_file))):
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        for transcript in tRNA:
            exons = GTF.asRanges(transcript, "exon")
            counts = iCLIP.count_intervals(bamfile,
                                           exons,
                                           strand=transcript[0].strand,
                                           contig=transcript[0].contig)
            converter = iCLIP.TranscriptCoordInterconverter(transcript)
            counts.index = converter.genome2transcript(counts.index.values)
            if len(counts) == 0:
                counts = pandas.Series([0], index=[1])
            counts.name = "count"
            counts.index.name = "base"
            counts = counts.sort_index()
            counts = counts.reset_index()
            counts["tRNA"] = transcript[0].transcript_id
            outs.append(counts)

    outs = pandas.concat(outs)

    outs.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=False)


def get3UTR(gtffile, outfile):

    outlines = []

    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(gtffile))):
        exons = GTF.asRanges(transcript, "exon")
        cds = GTF.asRanges(transcript, "CDS")

        utrs = Intervals.truncate(exons,cds)

        if transcript[0].strand == "+":
            utr3 = [exon for exon in utrs
                    if exon[0] >= cds[-1][1]]
        else:
            utr3 = [exon for exon in utrs
                    if exon[-1] <= cds[0][0]]

        for exon in utr3:
            bed = Bed.Entry()
            bed.contig = transcript[0].contig
            bed.start = exon[0]
            bed.end = exon[1]
            bed.fields = [transcript[0].transcript_id,
                          ".",
                          transcript[0].strand]

            outlines.append(str(bed))

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("\n".join(outlines) + "\n")

        
def extendBedIntervals(infile, outfile, halfwidth):
    '''Get a window of a specified size around the center of each entry'''

    with IOTools.openFile(outfile, "w") as outf:
        for bed in Bed.iterator(IOTools.openFile(infile)):
            center = bed.start + (bed.end-bed.start)/2
            bed.start = center - halfwidth
            bed.end = center + halfwidth
            outf.write(str(bed) + "\n")


def calculateSplicingIndex(bamfile, gtffile, outfile):

    bamfile = pysam.AlignmentFile(bamfile)

    counts = E.Counter()

    for transcript in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(gtffile))):

        exons = GTF.asRanges(transcript, "exon")
        E.debug("Transcript: %s, %i introns" %
                (transcript[0].transcript_id, len(exons)-1))
        ei_juncs = [exon[1] for exon in exons[:-1]]
        ie_juncs = [exon[0] for exon in exons[1:]]

        for junc in ei_juncs:
            reads = bamfile.fetch(
                reference=transcript[0].contig, start=junc, end=junc+1)
            spliced = [read for read in reads if 'N' in read.cigarstring]
            unspliced = [read for read in reads if 'N' not in read.cigarstring]

            if transcript[0].stand == "+":
                direction = "ei"
            else:
                direction = "ie"

            for read in unspliced:
                if (read.reference_end >= junc+3
                   and read.reference_start <= junc-3):
                    counts[direction+"_included"] += 1
                else:
                    counts["uncounted"] += 1

            for read in spliced:
                block_ends = [block[1] for block in read.get_blocks]
                if read.reference_start <= junc-3 and junc in block_ends:
                    counts[direction+"_excluded"] += 1

        for junc in ie_juncs:

            reads = bamfile.fetch(
                reference=transcript[0].contig, start=junc-1, end=junc)
            spliced = [read for read in reads if 'N' in read.cigarstring]
            unspliced = [read for read in reads if 'N' not in read.cigarstring]

            if transcript[0].stand == "-":
                direction = "ei"
            else:
                direction = "ie"

            for read in unspliced:
                if (read.reference_end >= junc+3
                   and read.reference_start <= junc-3):
                    counts[direction+"_included"] += 1

            for read in spliced:
                block_starts = [block[0] for block in read.get_blocks]
                if read.reference_start <= junc-3 and junc in block_starts:
                    counts[direction+"_excluded"] += 1

    header = "\t".join(["exon_intron_included",
                        "exon_intron_excluded",
                        "intron_exon_included",
                        "intron_exon_excluded",
                        "uncounted"])

    with IOTools.openFile(outfile, "w") as outf:

        outf.write(header+"\n")
        outf.write("\t".join(counts["ei_included"],
                             counts["ei_excluded"],
                             counts["ie_included"],
                             counts["ie_excluded"],
                             counts["uncounted"]) + "\n")

@cluster_runnable      
def exonBoundaryProfiles(bamfile, gtffile, outfile):

    bamfile = pysam.AlignmentFile(bamfile)

    counts_collector = []
    nExamined = 0
    nEligable = 0

    E.info("Processing transcripts..")
    for transcript in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(gtffile))):
        
        nExamined += len(transcript) - 1
        E.debug("Processing transcript %s" % transcript[0].transcript_id)

        exons = [(x.start, x.end) for x in transcript if x.feature == "exon"]
        exons.sort()


        E.debug("Processing transcript %s on contig %s, exons %s"
                % (transcript[0].transcript_id,
                   transcript[0].contig,
                   exons))
        counts = iCLIP.count_intervals(bamfile, exons, transcript[0].contig,
                                       transcript[0].strand)

        exons = exons[1:-1]
        coords_translator = iCLIP.TranscriptCoordInterconverter(transcript)
        
        counts.index = coords_translator.genome2transcript(counts.index.values)

        def elen(e):
            return e[1] - e[0]

        exons_starts = []

        for i, exon in enumerate(exons[1:]):
            if elen(exons[i-1]) > 100 and elen(exon) > 100:
                exons_starts.append(exon[0])

        exons_starts = coords_translator.genome2transcript(exons_starts)

        E.debug("%i of %i boundaries eligible, %s total counts" %
                (len(exons_starts), len(transcript)-1, counts.sum()))
        
        nEligable += len(exons_starts)
        
        for boundary in exons_starts:
            boundary_counts = counts[boundary-100:boundary+100]
            boundary_counts.index = boundary_counts.index.values - boundary

            counts_collector.append(boundary_counts)

    E.info("Finished reading transcripts. Constructing matrix")

    final_matrix = pd.concat(counts_collector, axis=1).transpose()
    final_matrix = final_matrix.fillna(0)

    E.info("Normalising and computing profile")
    row_sums = final_matrix.sum(axis=1)
    normed_matrix = final_matrix.div(row_sums.astype("float"), axis=0)
    combined = normed_matrix.sum().reset_index()

    E.info("Examined %i boundaries, of which %i were eligable and contained %i reads" %
           (nExamined, nEligable, row_sums.sum()))
    E.info("Writing Results")
    combined.to_csv(IOTools.openFile(outfile, "w"), sep="\t", index=False, header=["position", "density"])


def findNuclearLocalisation(infile, outfile):

    from rpy2.robjects import r as R, Formula
    from rpy2.robjects.packages import importr

    deseq = importr("DESeq2")

    c = R.c

    counts = R('read.delim("%(infile)s",sep="\t", row.names=1)' % locals())

    counts = R["as.matrix"](counts)
    col_data = R.strsplit(R.colnames(counts), "_")
    col_data = R["do.call"](R.rbind, col_data)
    col_data = R["data.frame"](col_data)
    col_data.rownames = R.colnames(counts)
    col_data.names = c("knockdown", "fraction", "replicate")

    print col_data
    keep = col_data.rx((col_data.rx2("fraction").ro != "Total").ro
                       & (col_data.rx2("knockdown").ro == "Control"),
                       True)

    print keep

    kept_counts = counts.rx(True, R.rownames(keep))
    print R.relevel(keep.rx2("fraction"), "Cytoplasmic")
    keep.rx[True,"fraction"] = R.relevel(keep.rx2("fraction"), "Cytoplasmic")
   
    fraction_ds = deseq.DESeqDataSetFromMatrix(
        kept_counts, keep, Formula("~ replicate + fraction"))
    fraction_ds = deseq.DESeq(fraction_ds)
    fraction_results = deseq.results(fraction_ds, independentFiltering=False,
                                     addMLE=True)

    fraction_results = R("as.data.frame")(fraction_results)
    fraction_results = R("data.frame")(fraction_results,
                                       Gene_id=R.rownames(fraction_results))

    R["write.table"](fraction_results, outfile, quote=False,
                     sep="\t", row_names=False)


def findhnRNPUDependentGenes(connection, outfile):

    from rpy2 import robjects as ro
    py2ri_old = ro.conversion.py2ri
 #   ro.conversion.py2ri = ro.numpy2ri

    import pandas.rpy.common as com

#    ro.numpy2ri.activate()

    counts = DUtils.fetch_DataFrame(''' SELECT DISTINCT gene_id,
                                              Name as Transcript,
                                              WT_NumReads as WT,
                                              hnRNPU1kd_NumReads as kd
                                         FROM fuRNAseq as rna
                                         INNER JOIN biotypes
                                           ON biotypes.transcript_id = rna.Name ''',
                                    connection)

    counts = counts.drop("Transcript", axis=1)
    counts = counts.groupby("gene_id").sum()
    counts.WT = counts.WT.astype(int)
    counts.kd = counts.kd.astype(int)

#    counts_r = ro.r.matrix(counts.as_matrix())
    counts_r = com.convert_to_r_matrix(counts)
    ro.r.assign("counts_r", counts_r)
    ro.numpy2ri.deactivate()
    ro.conversion.py2ri = py2ri_old

    ro.r.assign("rn", ro.StrVector(list(counts.index.values)))
    ro.r.assign("cn", ro.StrVector(list(counts.columns.values)))
    ro.r("rownames(counts_r) <- rn")
    #counts_r.colnames = ro.StrVector(list(counts.index.values))
    ro.r("colnames(counts_r) <- cn")
    
    counts_r = ro.r("counts_r")

    col_data = ro.r("data.frame")(si=ro.r("c(colnames(counts_r))"))
    print ro.r.levels(col_data.rx2("si"))
    col_data[0] = ro.r.relevel(col_data.rx2("si"), "WT")
    print ro.r.relevel(col_data.rx2("si"),"WT")
    print ro.r.levels(col_data.rx2("si"))

    deseq = ro.packages.importr("DESeq2")
    si_ds = deseq.DESeqDataSetFromMatrix(counts_r, col_data, ro.Formula("~si"))
    si_ds = deseq.DESeq(si_ds)
    si_results = deseq.results(si_ds, independentFiltering=False, addMLE=True)
    print (ro.r.head(si_results))
    si_results = ro.r("as.data.frame")(si_results)
    si_results = ro.r("data.frame")(si_results, gene_id=si_results.rownames)

    ro.r("write.table")(si_results, outfile, quote=False, sep="\t",
                        row_names=False)


def mergeExonsAndIntrons(exons, introns):
    ''' Exons are merged with all other exons and introns with all
    other introns. Introns and exons are then combined such that
    each non-overlapping exon and intron is output as a seperate
    GTF exon. Where exons and introns overlap, three seperate
    exons are output one for the unique parts of each exon
    and one for the overlap. '''

    template_entry = exons[0]
    try:
        exons = [(start, end, "exon") for start, end
                 in GTF.asRanges(exons, "exon")]
    except AttributeError:
        exons = [(start, end, "exon") for start, end
                 in Intervals.combine(exons)]

    try:
        introns = [(start, end, "intron") for start, end
                   in GTF.asRanges(introns, "exon")]
    except AttributeError:
        introns = [(start, end, "intron") for start, end
                   in Intervals.combine(introns)]

    combined_exons = exons + introns
    combined_exons.sort()
    
    new_exons = []
    last_from, last_to, last_feature = combined_exons[0]
    
    for this_from, this_to, this_feature in combined_exons[1:]:

        if this_from >= last_to:
            new_exons.append((last_from, last_to, last_feature))
            last_from, last_to, last_feature = this_from, this_to, this_feature
            continue

        new_exons.append((last_from, this_from, last_feature))
            
        if this_to <= last_to:
            new_exons.append((this_from, this_to, this_feature))
            last_from = this_to

        else:
            new_exons.append((this_from, last_to, "exon"))
            last_feature = this_feature
            last_from = last_to
            last_to = this_to
    
    entry = GTF.Entry()

    try:
        entry = entry.copy(template_entry)
    except AttributeError:
        pass

    entry.feature = "exon"
    entry.transcript_id = entry.gene_id
    entry.source = "merged"
    iexon = iintron = 1

    for start, end, feature in new_exons:

        entry.start = start
        entry.end = end
        if feature == "exon":
            entry.attributes["exon_id"] = "E%03d" % iexon
            iexon += 1
        else:
            entry.attributes["exon_id"] = "I%03d" % iintron
            iintron += 1

        yield entry


def getTranscriptsPlusRetainedIntrons(infile, outfile):
    ''' Look for transcripts with retained introns and the
    equivalent transcript without a retained intron. Output
    a merged gene model, where the retained intron is a seperate
    exon. Also merge in the detained introns '''

    outf = IOTools.openFile(outfile, "w")

    for gene in GTF.gene_iterator(
            GTF.iterator(IOTools.openFile(infile))):
        
        gene_out = []
        introns_out = []

        # now find if any of the transcripts are retained intron
        # versions of any of the others
        for first, second in itertools.product(gene, gene):
            
            first = sorted([entry for entry in first
                            if entry.feature == "exon"],
                           key=lambda x: x.start)
            second = sorted([entry for entry in second
                             if entry.feature == "exon"],
                            key=lambda x: x.start)

            first_introns = set(GTF.toIntronIntervals(first))
            second_introns = set(GTF.toIntronIntervals(second))
            
            if len(first_introns-second_introns) > 0 and \
               len(second_introns-first_introns) == 0:
                novel_introns = list(first_introns-second_introns)

                def _filterIntron(intron):
                    return intron[0] > second[0].start and \
                        intron[1] < second[-1].end

                novel_introns = filter(_filterIntron, novel_introns)

                if len(novel_introns) > 0:
                    gene_out.extend(first)

                for intron in novel_introns:
                    introns_out.append(intron)

        if len(gene_out) == 0:
            continue

        for gff in mergeExonsAndIntrons(gene_out, introns_out):
            outf.write("%s\n" % gff)


def insertDetainedIntronsIntoTranscripts(transcripts, introns, outfile):
    ''' This function takes the annotated detained introns of sharp
    and attempts to find the transcript from which they came, and
    returns a flattened gene containing combined exons from the
    transcripts containing the intron, plus an exon containing the
    intron'''

    outf = IOTools.openFile(outfile, "w")

    E.debug("Reading and Indexing bedfile")
    beds = Bed.readAndIndex(IOTools.openFile(introns), with_values=True)

    for gene in GTF.gene_iterator(
            GTF.iterator(IOTools.openFile(transcripts))):

        # Find any detained introns that are in the full interval
        # of the gene
        gene_start = min(exon.start
                         for transcript in gene
                         for exon in transcript)
        gene_end = max(exon.end for transcript in gene for exon in transcript)
        try:
            candidate_beds = set((start, end)
                                 for start, end, bed
                                 in beds[gene[0][0].contig].find(
                                     gene_start, gene_end)
                                 if bed.strand == gene[0][0].strand)
        except KeyError:
            candidate_beds = set()
         
        gene_out = []
        introns_out = []
        # First find if any transcripts have 'detained introns' and
        # add them
        for transcript in gene:
            introns = set(GTF.toIntronIntervals(transcript))
            retained_beds = candidate_beds & introns
        
            if len(retained_beds) > 0:
                gene_out.extend(exon for exon in
                                transcript if exon.feature == "exon")

            for intron in retained_beds:
                introns_out.append((intron[0], intron[1]))

        introns_not_found = set(candidate_beds) - set(introns_out)
        if len(introns_not_found) > 0:
            # exact homes not found. Find approximate homes.

            for transcript in gene:
                exons = GTF.asRanges(transcript, ("exon"))
                starts, ends = zip(*exons)
                for intron in introns_not_found:

                    if Intervals.calculateOverlap(exons, [intron]) == 0 \
                       and (intron[0] in ends or intron[1] in starts):
                        gene_out.extend(transcript)
                        introns_out.append(intron)

        introns_not_found = set(candidate_beds) - set(introns_out)

        if len(introns_not_found) > 0:
            E.warn("Homes not found for %s in gene %s chr=%s start=%i, end=%i"
                   % (introns_not_found, gene[0][0].gene_id, gene[0][0].contig,
                      gene_start, gene_end))
            
        if len(introns_out) == 0:
            continue

        for gff in mergeExonsAndIntrons(gene_out, introns_out):
            outf.write("%s\n" % str(gff))

def mergeDetainedAndRetainedIntrons(infile, outfile):

    outf = IOTools.openFile(outfile, "w")

    for gene in GTF.flat_gene_iterator(GTF.iterator(
            IOTools.openFile(infile))):
        
        exons = [exon for exon in gene
                 if "E" in exon.asDict()["exon_id"]]

        introns = [exon for exon in gene
                   if "I" in exon.asDict()["exon_id"]]

        for gff in mergeExonsAndIntrons(exons, introns):

            outf.write("%s\n" % str(gff))
    
    outf.close()


def convertGenesetToTranscriptomeCoords(infile, outfile):

    with IOTools.openFile(outfile, "w") as outf:
        for transcript in GTF.transcript_iterator(
                GTF.iterator(
                    IOTools.openFile(infile))):

            converter = iCLIP.TranscriptCoordInterconverter(transcript)
            transcript.sort(key=lambda x: x.start)

            exon_intervals = converter.transcript_intervals

            for interval, entry in zip(exon_intervals, transcript):
                entry.contig = entry.transcript_id
                entry.strand = "+"
                entry.start = interval[0]
                entry.end = interval[1]
                outf.write(str(entry) + "\n")


@cluster_runnable
def getOverlapFractions(test_set, category, background, outfile):
    ''' Takes three bed files, and outputs the fraction number of test_set in background,
    conditional on whether there is an overlap with category or not '''

    import pybedtools

    test_set = pybedtools.BedTool(test_set)
    category = pybedtools.BedTool(category)
    background = pybedtools.BedTool(background)

    cat_true = (background + category).count()
    test_cat_true = (background + category + test_set).count()

    cat_false = (background - category).count()
    test_cat_false = (background - category + test_set).count()

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("in_category\toverlap\ttotal\n")
        outf.write("True\t%i\t%i\n" % (test_cat_true,cat_true))
        outf.write("False\t%i\t%i\n" % (test_cat_false, cat_false))


@cluster_runnable
def getUnmappedNucleotideComp(infile, outfile):

    bamfile = pysam.AlignmentFile(infile)

    results = collections.defaultdict(int)
    for read in bamfile.fetch(until_eof=True):
        
        if not read.is_unmapped:
            continue
            
        composition = collections.defaultdict(int)
        l = len(read.query_sequence)
        for na in read.query_sequence:
            composition[na] += 1
            
        for na in composition:
            results[(str(l), na, str(composition[na]))] += 1

    with IOTools.openFile(outfile, "w") as outf:

        outf.write(
            "\t".join(["read_length", "base", "count", "nreads"]) + "\n")

        for (length, base, count), nreads in results.iteritems():
            outf.write("\t".join([length, base, count, str(nreads)]) + "\n")

@cluster_runnable
def intron_profiles(bamfile, gtffile, outfiles):

    long_introns = pd.Series(np.zeros(100), index=range(100))
    small_introns = pd.Series(np.zeros(100), index=range(100))

    bam = pysam.AlignmentFile(bamfile)
    nGenes = 0
    nSmall_introns = 0
    nLarge_introns = 0

    for gene in GTF.flat_gene_iterator(
            GTF.iterator(IOTools.openFile(gtffile))):

        if nGenes % 1000 == 0:
            E.debug("%i genes examined" % nGenes)

        exons = []
        for exon in gene:
            if not exon.feature == "exon":
                continue

            exon = GTF.Entry().copy(exon)
            exon.transcript_id = "1"
            exons.append(exon)
        
        introns = GTF.toIntronIntervals(exons)

        introns = [intron for intron in introns if
                   not (1000 <= (intron[1] - intron[0]) < 10000)
                   and intron[1] - intron[0] > 105]
                 
        if len(introns) == 0:
            nGenes += 1
            continue

        intron_counts = iCLIP.count_intervals(bam, introns,
                                              gene[0].contig,
                                              gene[0].strand)

        if intron_counts.sum() == 0:
            nGenes += 1
            continue

        intron_coords = iCLIP.TranscriptCoordInterconverter(exons,
                                                            introns=True)

        intron_counts.index = intron_coords.genome2transcript(
            intron_counts.index.values)
        
        intron_counts = intron_counts.sort_index()

        for intron in introns:
            
            length = intron[1] - intron[0]

            intron = intron_coords.genome2transcript((intron[0], intron[1]-1))
            intron.sort()

            this_intron_counts = intron_counts.loc[intron[0]:intron[1]-5]
            
            if this_intron_counts.sum() == 0:
                continue

            intron = (intron[0], intron[1] + 1)

            this_intron_counts.index = this_intron_counts.index.values - intron[0]
            bins = np.linspace(0, length, num=101, endpoint=True)

            try:
                this_intron_counts = this_intron_counts.groupby(
                    list(pd.cut(this_intron_counts.index,
                                bins=bins,
                                labels=range(100),
                                include_lowest=True))).sum()
            except:
                E.debug(this_intron_counts)
                raise


            this_intron_counts = this_intron_counts / this_intron_counts.sum()
            
            if intron[1] - intron[0] <= 10000:
                nSmall_introns += 1
                small_introns = small_introns.add(this_intron_counts, fill_value=0)
               # E.debug(small_introns)
            elif intron[1] - intron[0] > 100000:
                nLarge_introns += 1
                long_introns = long_introns.add(this_intron_counts, fill_value=0)
               # E.debug(long_introns)

        nGenes += 1

        
    E.info("Examined %i short and %i long introns" % (nSmall_introns, nLarge_introns))

    reads = 0
    for outfile in outfiles:

        if "long" in outfile:
            introns = long_introns
        elif "short" in outfile:
            introns = small_introns
        else:
            raise ValueError()

        E.info("Normalising and computing profile")

        reads += introns.sum()
        introns = introns / introns.sum()
        combined = introns.reset_index()

        combined.to_csv(IOTools.openFile(outfile, "w"),
                        sep="\t",
                        index=False,
                        header=["position", "density"])

    E.info("Found %i crosslink sites in total" % reads)


def calcAverageDistance(profile1, profile2):
    ''' This function calculates the average distance of all
    pairwise distances in two profiles'''

    def _cartesian(x, y):
        return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])

    positions = _cartesian(profile1.index.values, profile2.index.values)
    counts = _cartesian(profile1.values, profile2.values)
    counts = np.prod(counts, axis=0)
    distances = np.abs(positions[:, 0] - positions[:, 1])

    mean_distance = (distances.astype("float64") * counts) / np.sum(counts)

    return mean_distance


def findMinDistance(profile1, profile2):
    '''Finds mean distance between each read in profile1
    and a read in profile2'''

    locations1 = profile1.index.values.astype("uint16")

    locations2 = profile2.index.values.astype("uint16")

    mat1 = np.repeat(locations1, locations2.size).reshape(
        (locations1.size, locations2.size))
    mat2 = np.tile(locations2, locations1.size).reshape(
        (locations1.size, locations2.size))

    distances = np.abs(mat1-mat2).min(axis=1)

    return distances.mean()


def spread(profile, bases):

    a = pd.Series()

    for i, count in profile.iteritems():
        a.add(pd.Series([count]*(2*bases),
                        range(i-bases, i)+range(i+1, i+bases+1)),
              fill_value=0)

    return profile.add(a, fill_value=0)


def randomiseSites(profile, start, end, keep_dist=True):
    '''Randomise clipped sites within an interval (between start and end)
    if keep_dist is true, then reads on the same base are kept togehter'''

    if keep_dist:

        profile = profile.copy()
        profile.index = np.random.choice(
            np.arange(start, end), profile.size, replace=False)
        profile = profile.sort_index()
        return profile

    else:
        
        randomised = np.random.choice(
            np.arange(start, end), profile.sum(), replace=True)
        randomised = pd.Series(randomised).value_counts().sort_index()
        return randomised

def liftOverFromHg18(infile, outfile):

    import CGATPipelines.Pipeline as P 

    if infile.endswith(".gz"):
        cat_cmd = "zcat"
    else:
        cat_cmd = "cat"

    statement = '''liftOver <( %(cat_cmd)s %(infile)s | grep -P "^chr" )
                            /ifs/mirror/ucsc/hg18/liftOver/hg18ToHg19.over.chain.gz
                            %(outfile)s
                            %(outfile)s.unmapped'''

    P.run()

@cluster_runnable
def getInternalExons(infile, outfile):

    outfile = IOTools.openFile(outfile, "w")

    for transcript in GTF.transcript_iterator(GTF.iterator(
                                   IOTools.openFile(infile))):

        transcript = [exon for exon in transcript if
                      exon.feature == "exon"]

        for exon in transcript[1:-1]:
            outfile.write(str(exon)+"\n")

    outfile.close()

###################################################################
def bamToBigWig(infile, outfile):

    import CGATPipelines.Pipeline as P

    genome_file = os.path.join(PARAMS['annotations_dir'],"contigs.tsv")

    if not infile.endswith(".bam"):
        infile = "-i <( zcat %(infile)s | sort -k1,1 -k2,2n)" % locals()
    else:
        infile = "-ibam %(infile)s" % locals()

    tmp = P.getTempFilename()
    statement = ''' genomeCoverageBed -split -bg %(infile)s
                                      -g %(genome_file)s  2> %(outfile)s.log
                   | sort -k1,1 -k2,2n > %(tmp)s ;
                    
                    checkpoint;

                    bedGraphToBigWig %(tmp)s
                                     %(genome_file)s
                                     %(outfile)s 2>>%(outfile)s.log;

                    checkpoint;

                    rm %(tmp)s'''

    P.run()


###################################################################
def bamToBedGraph(infile, outfile):

    import CGATPipelines.Pipeline as P

    genome_file = os.path.join(PARAMS['annotations_dir'],"contigs.tsv")

    statement = ''' genomeCoverageBed -split -bg -ibam %(infile)s 
                                      -g %(genome_file)s > %(outfile)s 2> %(outfile)s.log;
                '''
    P.run()

###################################################################
def bigWigTrackDB(infiles,
                  long_label_template,
                  group_name,
                  group_long_label,
                  outfile):

    import CGATPipelines.Pipeline as P

    template = '''
       track %(track_name)s
       parent %(group_name)s
       bigDataUrl %(track_data_URL)s
       shortLabel %(short_label)s
       longLabel %(long_label)s
       autoScale on
       visibility full
       alwaysZero on
       maxHeightPixels 32:32:11
       type bigWig
       windowingFunction mean
      '''

    stanzas = []
    for track in infiles:
        track_name = P.snip(os.path.basename(track), ".bigWig")
        track_data_URL = os.path.basename(track)
        short_label = track_name
        long_label = long_label_template % locals()
        
        stanzas.append(template % locals())

    composite_stanaz = '''

    track %(group_name)s
    shortLabel %(group_name)s
    longLabel %(group_long_label)s
    superTrack on
   '''
    
    output = "\n".join([composite_stanaz % locals()] + stanzas)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write(output)


###################################################################
def generateDaParsConfig(condition1_files,
                         condition2_files,
                         utrs,
                         dapars_outfile,
                         config_outfile):
    DaPars_config_template= '''
#The following file is the result of step 1.

Annotated_3UTR=%(utrs)s

#A comma-separated list of BedGraph files of samples from condition 1

Group1_Tophat_aligned_Wig=%(condition1_files)s
#Group1_Tophat_aligned_Wig=Condition_A_chrX_r1.wig,Condition_A_chrX_r2.wig if multiple files in one group

#A comma-separated list of BedGraph files of samples from condition 2

Group2_Tophat_aligned_Wig=%(condition2_files)s

Output_directory=%(outdir)s

Output_result_file=%(dapars_outfile)s

#At least how many samples passing the coverage threshold in two conditions
Num_least_in_group1=%(dapars_num_least_in_group)i

Num_least_in_group2=%(dapars_num_least_in_group)i

Coverage_cutoff=%(dapars_coverage_cutoff)i

#Cutoff for FDR of P-values from Fisher exact test.

FDR_cutoff=%(dapars_fdr_cutoff)s


PDUI_cutoff=%(dapars_pdui_cutoff)s

Fold_change_cutoff=%(dapars_logfc_cutoff)s
'''

    outdir = os.path.dirname(os.path.abspath(dapars_outfile))
    dapars_outfile = os.path.basename(dapars_outfile)
    condition1_files = ",".join([os.path.abspath(f) for f in condition1_files])
    condition2_files = ",".join([os.path.abspath(f) for f in condition2_files])

    local_params = PARAMS.copy()
    local_params.update(locals())

    with IOTools.openFile(config_outfile,"w") as outf:
        outf.write(DaPars_config_template % local_params)
        
    
def get3UTRs(infile, outfile):
    ''' flatten the gene so that all the coding regions are collapsed
    and extract the remaining exon sequence after the end of the last cds
    entry in the GTF - that sequence that is always UTR'''

    with IOTools.openFile(outfile, "w") as outf:
        for gene in GTF.flat_gene_iterator(
                GTF.iterator(IOTools.openFile(infile))):
            cds = GTF.asRanges(gene, "CDS")
            
            if len(cds) == 0:
                continue

            exons = GTF.asRanges(gene, "exon")

            introns = []

            for transcript in GTF.transcript_iterator(gene):
                intron = GTF.toIntronIntervals(transcript)
                introns.extend(intron)

            introns = Intervals.combine(introns)
            exons = Intervals.truncate(exons, introns)
            utrs = Intervals.truncate(exons, cds)

            if len(utrs) == 0:
                continue

            if gene[0].strand == "+":
                utr3 = [utr for utr in utrs if utr[0] >= cds[-1][1]]
            else:
                utr3 = [utr for utr in utrs if utr[1] <= cds[0][0]]

            template = GTF.Entry().fromGTF(gene[0])
            template.feature = "exon"

            for exon in utr3:
                template.start = exon[0]
                template.end = exon[1]
                outf.write(str(template) + "\n")


def find_final_cleavage_sites(geneset_gtf, cleavagesite_bed,
                              outfile):
    '''Iterate through geneset and isolate overlapping cleavage sites
    and select the 3' most of the sites'''

    sites = Bed.readAndIndex(IOTools.openFile(cleavagesite_bed),
                             with_values=True)

    with IOTools.openFile(outfile, "w") as outfile:

        for gene in GTF.flat_gene_iterator(
                GTF.iterator(IOTools.openFile(geneset_gtf))):

            start = min([exon.start for exon in gene])
            end = max([exon.end for exon in gene])

            overlapping_sites = list(sites[gene[0].contig].find(start, end))
            overlapping_sites = [site for _, _, site in overlapping_sites]
            overlapping_sites = [site for site in overlapping_sites
                                 if site.strand == gene[0].strand]

            if len(overlapping_sites) == 0:
                continue

            overlapping_sites = sorted(overlapping_sites,
                                       key=lambda x: x.start)

            if gene[0].strand == "-":
                outfile.write(str(overlapping_sites[0]) + "\n")
            else:
                outfile.write(str(overlapping_sites[-1]) + "\n")


@cluster_runnable
def ExonRatio(bamfile, gtffile, outfile):

    bamfile = pysam.AlignmentFile(bamfile)
    outlines = []

    for gene in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(gtffile))):

        strand = gene[0].strand
        contig = gene[0].contig
        
        exons = GTF.asRanges(gene, "exon")
        introns = GTF.toIntronIntervals(gene)

        if len(exons) == 1:
            continue

        if gene[0].strand == "+":
            first_exon = exons[1]
            first_intron = introns[0]
        else:
            first_exon = exons[-2]
            first_intron = introns[-1]

        exon_counts = iCLIP.count_intervals(bamfile,
                                            [first_exon],
                                            strand=strand,
                                            contig=contig).sum()
        intron_counts = iCLIP.count_intervals(bamfile,
                                              [first_intron],
                                              strand=strand,
                                              contig=contig).sum()
        
        first_exon_length = first_exon[1] - first_exon[0]
        first_intron_length = first_intron[1] - first_intron[0]

        outlines.append([gene[0].transcript_id, exon_counts, first_exon_length,
                         intron_counts, first_intron_length])

    IOTools.writeLines(outfile, outlines,
                       header=["gene_id", "exon_counts", "exon_length",
                               "intron_counts", "intron_length"])


@cluster_runnable
def get_first_exon(infile, outfile):

    outfile = IOTools.openFile(outfile, "w")

    for gene in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(infile))):
        
        gene = [exon for exon in gene if exon.feature == "exon"]
        gene = GTF.CombineOverlaps(gene)
        
        if gene[0].strand == "-":
            gene = sorted(gene, key=lambda x: x.end, reverse=True)
        else:
            gene = sorted(gene, key=lambda x: x.start)

        if len(gene == 1):
            continue

        outfile.write(str(gene[0]) + "\n")


def find_exon_ratios(infile, outfile):
    '''find the ratio of exons in genes in infile'''

    firsts, lasts = [], []

    for transcript in GTF.flat_gene_iterator(
                         GTF.iterator(
                             IOTools.openFile(infile))):
        if not transcript[0].source == "protein_coding":
            continue

        exons = GTF.asRanges(transcript, "exon")
        if len(exons) <= 2:
            continue

        first_exon = exons[0]
        last_exon = exons[-1]

        if transcript[0].strand == "-":
            first_exon, last_exon = last_exon, first_exon

        middle_exons = exons[1:-1]

        first_len = first_exon[0] - first_exon[1]
        last_len = last_exon[0] - last_exon[1]
        middle_len = sum(e[0] - e[1] for e in middle_exons)

        first_ratio = float(first_len)/middle_len
        last_ratio = float(last_len)/middle_len

        firsts.append(first_ratio)
        lasts.append(last_ratio)

    first_median = np.median(firsts)
    lasts_median = np.median(lasts)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("\t".join(map(str, (first_median, 1.0, lasts_median))))


@cluster_runnable
def get_first_last_counts(bamfile, gtffile, outfile):

    getter = iCLIP.make_getter(bamfile=bamfile, centre="GFP" in bamfile)
    outlines = []

    for gene in GTF.flat_gene_iterator(GTF.iterator(
            IOTools.openFile(gtffile))):
        exons = GTF.asRanges(gene, "exon")

        if not gene[0].source == "protein_coding":
            continue

        contig = gene[0].contig
        strand = gene[0].strand

        if len(exons) <= 2:
            if len(exons) == 1:
                single_exon = iCLIP.count_intervals(getter,
                                                    exons,
                                                    contig=contig,
                                                    strand=strand).sum()
                outlines.append(map(str,
                                    [gene[0].gene_id,
                                     "NA", "NA", "NA", "NA", single_exon]))
            continue

        first = iCLIP.count_intervals(getter,
                                      [exons[0]],
                                      contig=contig,
                                      strand=strand).sum()

        middle = iCLIP.count_intervals(getter,
                                       exons[1:-1],
                                       contig=contig,
                                       strand=strand).sum()

        last = iCLIP.count_intervals(getter,
                                     [exons[-1]],
                                     contig=contig,
                                     strand=strand).sum()

        introns = iCLIP.count_intervals(getter,
                                        GTF.toIntronIntervals(gene),
                                        contig=contig,
                                        strand=strand).sum()
        outlines.append(map(str,
                            [gene[0].gene_id, first, middle, last, introns, "NA"]))

    IOTools.writeLines(
        outfile, outlines,
        header=["gene_id", "first_exon",
                "middle_exon", "last_exon", "introns", "single_exon"])


