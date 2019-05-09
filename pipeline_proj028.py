##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline proj028
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline template.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_proj028.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_proj028.tgz
   tar -xvzf pipeline_proj028.tgz
   cd pipeline_proj028
   python <srcdir>/pipeline_proj028.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::
  

Code
====

"""
from ruffus import *
from ruffus.combinatorics import *
import sys
import glob
import os
import re
import sqlite3
import pandas
import collections

import PipelineProj028
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as DUtils

import CGATPipelines.PipelineRnaseq as PipelineRnaseq
import PipelineiCLIP
import CGATPipelines.PipelineMotifs as PipelineMotifs
import CGATPipelines.PipelineGO as PipelineGO

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters(PARAMS["annotations_dir"],
                                      "pipeline_annotations.py",
                                      on_error_raise=True)

PARAMS["project_src"] = os.path.dirname(__file__)
PipelineiCLIP.PARAMS = PARAMS
PipelineMotifs.PARAMS = PARAMS
PipelineGO.PARAMS = PARAMS
PipelineProj028.PARAMS = PARAMS

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob("*.ini" ), "(\S+).ini" )


###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database_name"] )
    statement = '''ATTACH DATABASE '%s' as annotations;
                   ATTACH DATABASE '%s' as chtop_apa''' % (PARAMS["annotations_database"],
                                                           PARAMS["external_chtop_apa_db"])
    cc = dbh.cursor()
    cc.executescript( statement )
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
## worker tasks
###################################################################
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.fasta.gz")
def getGenesetFasta(infile, outfile):
    '''convert geneset annotation in fasta '''

    statement = ''' cgat gff2fasta
                          --is-gtf
                          --feature=exon
                          --genome=%(genome_dir)s/%(genome)s
                          -I %(infile)s
                          -L %(outfile)s.log
                  | cut -d' ' -f1
                  | gzip > %(outfile)s '''
    P.run()


###################################################################
@transform(getGenesetFasta, 
           regex(".+"),
           "salmon_index.dir/txpInfo.bin")
def generateSalmonIndex(infile,outfile):
    ''' generate salmon index for specified geneset '''

    index_dir = os.path.dirname(outfile)
    tmpfile = P.getTempFilename()
    logfile = os.path.basename(outfile) + ".log"
    job_threads=8
    
    statement = ''' zcat %(infile)s > %(tmpfile)s;
                    checkpoint;
                    salmon index -t %(tmpfile)s
                                  -i %(index_dir)s
                    >& %(logfile)s;
                    checkpoint;
                    rm %(tmpfile)s '''
    P.run()


###################################################################
@follows(generateSalmonIndex)
@collate("*fastq*gz",
         regex("(.+).fastq.(?:[12]\.)*gz"),
         r"\1_salmon.dir/quant.sf")
def runSalmon(infiles, outfile):
    ''' Run sailfish on any provided rnaseq files '''

        
    job_memory = "10G"

    track = re.match("(.+)_salmon.dir/quant.sf", outfile).groups()[0]

    if len(infiles) == 1:
        inputs = "-r %s" % infiles[0]
    elif len(infiles) == 2:
        inputs = "-1 <(zcat %s ) -2 <( zcat %s)" % infiles
    else:
        raise ValueError("Don't know how to handle %i input files"
                         % len(infiles))

    statement = ''' salmon quant -i salmon_index.dir
                                   -l '%(sailfish_libtype)s'
                                   %(inputs)s
                                   -o %(track)s_salmon.dir '''

    P.run()


###################################################################
@merge(runSalmon,
       "HEK293_salmon.load")
def loadSalmon(infile, outfile):
    P.concatenateAndLoad(infile, outfile,
                         regex_filename="(.+)_salmon.dir/quant.sf",
                         options="-i transcript_id")


###################################################################
@follows(mkdir("stubbsRNAseq.dir"))
@transform(os.path.join(PARAMS["dir_external"], "stubbsRNAseq/*.bam"),
           regex(".+/(.+)\.merged.+"),
           add_inputs(os.path.join(
               PARAMS["iclip_dir"],
               "mapping.dir/geneset.dir/reference.gtf.gz")),
           r"stubbsRNAseq.dir/\1.feature_counts.tsv.gz")
def countStubbsRNAseq(infiles, outfile):

    bamfile, annotations = infiles
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=0,
        options=PARAMS['featurecounts_options'] + ' -O -M')


###################################################################
@merge(countStubbsRNAseq,
       "stubbsRNAseq.dir/stubbs_counts.tsv.gz")
def mergeStubbsCounts(infiles, outfile):
    '''Merge feature counts data into one table'''

    infiles = " ".join(infiles)

    statement=''' cgat combine_tables
                         -c 1
                         -k 7
                         --regex-filename='(.+).feature_counts.tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(mergeStubbsCounts,
           suffix(".tsv.gz"),
           ".load")
def loadStubbsCounts(infile, outfile):

    P.load(infile, outfile, options="-i gene_id")


###################################################################
@transform(loadStubbsCounts,
           regex("stubbsRNAseq.dir/(.+).load"),
           inputs(r"stubbsRNAseq.dir/\1.tsv.gz"),
           r"stubbsRNAseq.dir/\1.edgeR.tsv")
def runStubbsFractionEdgeR(infile, outfile):
    '''Use edgeR to calculate the distribution between cytoplasmic
    and nucelar fractions in Control and Alyref knockdown, and the
    significance of the difference'''

    statement = '''Rscript %(project_src)s/differential_fraction.R
                     --infiles=%(infile)s --outfiles=%(outfile)s 
                   &> %(outfile)s.log '''

    P.run()


###################################################################
@transform(runStubbsFractionEdgeR, suffix(".tsv"), ".load")
def loadStubbsFractionEdgeR(infile, outfile):

    P.load(infile, outfile, "-i gene_id")



###################################################################
@follows(loadStubbsFractionEdgeR,
         loadSalmon)
def rnaseq():
    pass
###################################################################
# Start/Stop Codon profiles
###################################################################

@transform(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex("(.+)"),
           add_inputs(loadSalmon),
           "expressed_transcripts.gtf.gz")
def filterExpressedTranscripts(infiles, outfile):
    '''Filter the geneset for those genes that are expressed in 
    HEK293 '''


    infile = infiles[0]
    
    query = ''' SELECT Name as transcript_id FROM HEK293_salmon
                GROUP BY transcript_id
                HAVING AVG( TPM ) > 0.5 '''

    transcript_ids = DUtils.fetch(query, connect())
    
    tmp = P.getTempFilename(dir="/ifs/scratch/")

    IOTools.writeLines(tmp, transcript_ids)

    statement = '''cgat gtf2gtf
                          --method=filter --filter-method=transcript
                          --map-tsv-file=%(tmp)s
                          -I %(infile)s
                          -L %(outfile)s.log
                 | gzip -c > %(outfile)s '''

    P.run()

@follows(mkdir("gene_clip_counts.dir"))
@transform([os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
            os.path.join(PARAMS["ejc_iclip_dir"], "deduped.dir/*.bam")],
           formatter(),
           add_inputs(filterExpressedTranscripts),
           r"gene_clip_counts.dir/{basename[0]}.expressed_gene_counts.tsv.gz")
def count_clips_on_genes(infiles, outfile):
    '''Use count_clip_sites to count the number of clips on genes
    useing only expressed transcripts'''

    bamfile, gtffile = infiles

    statement = '''python %(project_src)s/iCLIPlib/scripts/count_clip_sites.py
                   -I %(gtffile)s
                   -f gene
                   %(bamfile)s
                   -L %(outfile)s.log
                   -S %(outfile)s'''

    P.run()


@merge(count_clips_on_genes,
       "gene_clip_counts.dir/expressed_clip_counts.load")
def load_expressed_clip_counts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).expressed_gene_counts.tsv.gz",
                         options = " -i track -i gene_id")

    
###################################################################
@follows(mkdir("gene_profiles.dir"))
@transform([os.path.join(PARAMS["dir_iclip"], "deduped.dir/*-FLAG*.bam"),
            "*.bam"],
           regex("(?:.+/)?(.+).bam"),
           add_inputs(filterExpressedTranscripts),
           r"gene_profiles.dir/\1.CDS_profile.log")
def calculateCDSProfiles(infiles, outfile):
    ''' Some theory suggests that CHTOP might interact with a
    complex that binds near to the stop codon and methylates
    adenosine. Test this by seeing if ChTOP localises near
    the stop codon '''

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles
    outfile = P.snip(outfile, ".CDS_profile.log")
    statement = '''cgat bam2geneprofile
                          --method=utrprofile
                          --bam-file=%(bamfile)s
                          --gtf-file=<(zcat %(gtffile)s | grep protein_codin)
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --scale-flank-length=1
                          --resolution-cds=100
                          --resolution-downstream=100
                          --resolution-upstream=100
                          --resolution-upstream-utr=20
                          --resolution-downstream-utr=70
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.CDS_profile.log 
                          --output-all-profiles '''

    P.run()

    
@follows(mkdir("gene_profiles.dir"))
@transform([os.path.join(PARAMS["dir_iclip"], "deduped.dir/*-FLAG*.bam"),
            "*-*-*.bam"],
           regex("(?:.+/)?(.+).bam"),
           add_inputs(filterExpressedTranscripts),
           r"gene_profiles.dir/\1.iCLIP_utrprofile.tsv.gz")
def calculate_iCLIP_utr_profiles(infiles, outfile):

    bamfile, gtffile = infiles
    
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_transcript_region_metagene.py 
                           -I <(zcat %(gtffile)s | awk '$2=="protein_coding"')
                           %(bamfile)s
                           --regions=5flank,5UTR,CDS,3UTR,3flank
                           --region-length-correction
                           --flanks-length=500
                           -S %(outfile)s'''

    P.run()

@follows(mkdir("gene_profiles.dir"))
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*-FLAG*.bam"),
           regex("(?:.+/)?(.+).bam"),
           add_inputs(filterExpressedTranscripts),
           r"gene_profiles.dir/\1.iCLIP_firstmidlast.tsv.gz")
def calculate_iCLIP_firstmidlast_profiles(infiles, outfile):

    bamfile, gtffile = infiles
    
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_transcript_region_metagene.py 
                           -I <(zcat %(gtffile)s | awk '$2=="protein_coding"')
                           %(bamfile)s
                           --regions=5flank,first_exon,middle_exons,introns,last_exon,3flank
                           --region-length-correction
                           --flanks-length=500
                           -S %(outfile)s'''

    P.run()

    
@follows(mkdir("gene_profiles.dir"))
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*-FLAG*.bam"),
           regex("(?:.+/)?(.+).bam"),
           add_inputs(filterExpressedTranscripts),
           r"gene_profiles.dir/\1.iCLIP_firstmidlast_nonorm.tsv.gz")
def calculate_iCLIP_firstmidlast_nonorm_profiles(infiles, outfile):

    bamfile, gtffile = infiles
    
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_transcript_region_metagene.py 
                           -I <(zcat %(gtffile)s | awk '$2=="protein_coding"')
                           %(bamfile)s
                           --regions=5flank,first_exon,middle_exons,introns,last_exon,3flank
                           --region-length-correction
                           --flanks-length=500
                           --no-gene-norm
                           -S %(outfile)s'''

    P.run()

    
@collate([ calculate_iCLIP_utr_profiles,
           calculate_iCLIP_firstmidlast_profiles,
           calculate_iCLIP_firstmidlast_nonorm_profiles],
         regex("gene_profiles.dir/(.+)\.iCLIP_(.+).tsv.gz"),
        r"gene_profiles.dir/iCLIP_\2.load")
def load_iCLIP_utrprofiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)-(.+)[-\.](.+).iCLIP_.+.tsv.gz",
                         cat="factor,tag,replicate",
                         options="-i factor -i tag -i replicate -i region")

    
###################################################################
@follows(mkdir("gene_profiles.dir"))
@transform( [os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
           "*.bam"],
           regex("(?:.+/)?(.+).bam"),
           add_inputs(filterExpressedTranscripts),
           r"gene_profiles.dir/\1.tssprofile.log")
def calculateSTOPProfiles(infiles, outfile):
    ''' Some theory suggests that CHTOP might interact with a
    complex that binds near to the stop codon and methylates
    adenosine. Test this by seeing if ChTOP localises near
    the stop codon '''

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles
    outfile = P.snip(outfile, ".tssprofile.log")
    statement = '''cgat bam2geneprofile
                          --method=tssprofile
                          --bam-file=%(bamfile)s
                          --gtf-file=<(zcat %(gtffile)s | awk '$3=="CDS"' | sed 's/CDS/exon/')
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --normalize-profile=area
                          --scale-flank-length=1
                          --resolution-downstream=400
                          --resolution-upstream=400
                          --extension-inward=400
                          --extension-outward=400
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.tssprofile.log 
                          --output-all-profiles '''

    P.run()


###################################################################
@subdivide(calculateCDSProfiles, 
           regex("gene_profiles.dir/(.+-FLAG-R[0-9]+).CDS_profile.log"),
           inputs([r"gene_profiles.dir/\1.utrprofile.profiles.tsv.gz",
                   r"gene_profiles.dir/%s.utrprofile.profiles.tsv.gz" %
                   os.path.basename(P.snip(PARAMS["rnaseq_nuclear"], ".bam"))]),
           [r"gene_profiles.dir/\1.rnaseq_normed.matrix.tsv.gz",
            r"gene_profiles.dir/\1.rnaseq_normed.profile.tsv.gz"])
def normalizedToRNASeq(infiles, outfiles):
    ''' Normalise each bin to the corresponding bin in RNAseq data '''

    logfile = P.snip(outfiles[0],"matrix.tsv.gz") + ".log"
    PipelineProj028.normalizeIndevidualProfilesToRNASeq(
        infiles[0],
        infiles[1],
        outfile_matrix=outfiles[0],
        outfile_summary=outfiles[1],
        pseduo_count=1,
        submit=True,
        logfile=logfile,
        job_memory="25G")

    
###################################################################
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
           formatter(),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"gene_profiles.dir/{basename[0]}.introns.profile.tsv")
def get_intron_metagene(infiles, outfile):
    '''Get a metagene profile over concatenated introns'''

    bamfile, gtffile = infiles

    statement = '''cgat gtf2gtf -I %(gtffile)s
                                --method=exons2introns
                                -L %(outfile)s.log
                 | awk -F'\\t' 'BEGIN{OFS=FS} {$3="exon"; print}'
                 | python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                   -f 0
                   %(bamfile)s
                   --exon-bins=100
                   -S %(outfile)s
                   -L %(outfile)s.log; '''

    P.run()
    
###################################################################
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
           formatter(),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"gene_profiles.dir/{basename[0]}.exons.profile.tsv")
def get_exon_metagene(infiles, outfile):
    '''Get a metagene profile over concatenated introns'''

    bamfile, gtffile = infiles

    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                   -f 500
                   --scale-flanks
                   --flank-bins=100
                   --exon-bins=100
                   %(bamfile)s
                   -I %(gtffile)s
                   -S %(outfile)s
                   -L %(outfile)s.log; '''

    P.run()

@merge([get_exon_metagene,
        get_intron_metagene],
       "gene_profiles.dir/basic_iCLIP_metagenes.load")
def load_iCLIP_metagenes(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-.+)[-\.](R.|union).(.+).profile.tsv",
                         cat="pulldown,replicate,interval",
                         options="-i pulldown -i replicate -i interval")

@follows(load_iCLIP_metagenes)
def basic_iclip_metagenes():
    pass

###################################################################
@follows(calculateSTOPProfiles,
         normalizedToRNASeq)
def normalised_profiles():
    pass


###################################################################
# Single Vs Multi Exon genes
###################################################################
@originate("single_exon_genes.tsv")
def getSingleExonGenes(outfile):
    '''Generate a list of genes that have only one exon'''
    
    statement = ''' SELECT ti.gene_id
                    FROM
                       annotations.transcript_stats as es
                    INNER JOIN annotations.transcript_info as ti
                       ON ti.transcript_id = es.transcript_id
                    WHERE contig != 'chrM'
                    GROUP BY ti.gene_id
                    HAVING MAX(nval) = 1 '''

    results = DUtils.fetch(statement, connect())
    IOTools.writeLinesLines(outfile, results)


###################################################################
@transform(filterExpressedTranscripts,
           suffix(".gtf.gz"),
           add_inputs(getSingleExonGenes),
           ".single_exon.gtf.gz")
def getSingleExonGeneModels(infiles, outfile):
    '''Filter out the single exon genes from the list of expressed
    gene models'''

    gtffile, genelist = infiles

    statement = ''' cgat gtf2gtf
                           --method=filter
                           --filter-method=gene
                           --map-tsv-file=%(genelist)s
                           -I %(gtffile)s -L %(outfile)s.log
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(filterExpressedTranscripts,
           suffix(".gtf.gz"),
           add_inputs(getSingleExonGenes),
           ".multi_exon.gtf.gz")
def getMultiExonGeneModels(infiles, outfile):
    '''Filter out the single exon genes from the list of expressed
    gene models'''

    gtffile, genelist = infiles

    statement = ''' cgat gtf2gtf
                           --method=filter
                           --filter-method=gene
                           --map-tsv-file=%(genelist)s
                           --invert-filter
                           -I %(gtffile)s -L %(outfile)s.log
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
@product([os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
          "*.bam"],
         formatter(".+/(?P<TRACK>.+)\.bam"),
         [getSingleExonGeneModels, getMultiExonGeneModels],
         formatter("expressed_transcripts.(?P<GTF>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK[0][0]}.{GTF[1][0]}.geneprofile.tsv.gz")
def SingleVsMultiExonGeneProfiles(infiles, outfile):
    ''' Do metagene profiles for single and multi exon genes '''

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles
    matrix_out = P.snip(outfile, ".tsv.gz") + ".matrix.tsv.gz"
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                           -I %(gtffile)s
                           %(bamfile)s
                           --scale-flanks
                           --output-matrix=%(matrix_out)s
                            -S %(outfile)s
                           -b 25
                           --flank-bins=25
                           --normalised_profile '''

    P.run()


###################################################################
@merge(SingleVsMultiExonGeneProfiles,
       "single_vs_multi_exon_gene_profiles.load")
def loadSingleVsMultiExonGeneProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles,
                         outfile,
                         regex_filename="(.+)[-_\.]([0-9Ru].*)\.(single|multi)_",
                         cat="factor,replicate,exons",
                         options="-i factor -i replicate -i exons -i bin")


###################################################################
@follows(loadSingleVsMultiExonGeneProfiles)
def singleVsMultiExonGenes():
    pass


###################################################################
# Transcripts binned by length
###################################################################
@transform(filterExpressedTranscripts,
           suffix(".gtf.gz"),
           ".stats.tsv.gz")
def getExpressedTranscriptStats(infile, outfile):
    ''' get transcript lengths for expressed transcripts to allow a
    division into different length bins '''

    statement = '''cgat gtf2table
                   -c length -I %(infile)s -L %(outfile)s.log
                   -S %(outfile)s -r transcripts --add-gtf-source'''

    P.run()


###################################################################
@transform(getExpressedTranscriptStats,
           suffix(".tsv.gz"),
           ".load")
def loadExpressedTranscriptStats(infile, outfile):

    P.load(infile, outfile, options="-i transcript_id -i sum")


###################################################################
@follows(loadExpressedTranscriptStats)
@split(filterExpressedTranscripts,
       ["expressed_transcripts.%i_outof_5.gt_?_exons.gtf.gz" % (quantile + 1) for
        quantile in range(5)])
def getTranscriptsBinnedByLength(infiles, outfiles):
    '''divide transcripts into bins based on length '''

    tlens = DUtils.fetch_DataFrame('''SELECT transcript_id, sum as elen, nval
                                   FROM expressed_transcripts_stats
                                   WHERE source = 'protein_coding' ''', 
                                   connect())

    statement_template = ''' cgat gtf2gtf
                                   --method=filter
                                   --filter-method=transcript
                                   -a %(id_list)s
                                   -I %(infiles)s
                                   -L %(outfile)s.log
                                   -S %(outfile)s '''
    statements = []

    for exon_limit in [0,1]:
        temp_df = tlens[tlens.nval > exon_limit]
        temp_df["quantile"] = pandas.qcut(temp_df.elen, 5, labels=False)
        for quantile, group in temp_df.groupby("quantile"):
            id_list = "expressed_transcripts.%i_outof_5.gt_%i_exons.tsv" % \
                      (int(quantile)+1, exon_limit)
            outfile = "expressed_transcripts.%i_outof_5.gt_%i_exons.gtf.gz" % \
                      (int(quantile)+1, exon_limit)
            group.transcript_id.to_csv(id_list, index=False)
            statements.append(statement_template % locals())
    
    P.run()


###################################################################
@product(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
         formatter(".+/(?P<TRACK>.+)\.bam"),
         getTranscriptsBinnedByLength,
         formatter("expressed_transcripts.(?P<GTF>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK[0][0]}.{GTF[1][0]}.quantile.geneprofile.matrix.tsv.gz")
def getBinnedExpressionProfiles(infiles, outfile):
    '''Get gene profiles for each of five bins for gene length '''

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles
    outfile = P.snip(outfile, ".geneprofile.matrix.tsv.gz")
    statement = '''cgat bam2geneprofile
                          --method=geneprofile
                          --bam-file=%(bamfile)s
                          --gtf-file=%(gtffile)s 
                          --normalize-transcript=total-max
                          --use-base-accuracy
                          --normalize-profile=area
                          --scale-flank-length=1
                          --resolution-cds=250
                          --resolution-downstream=250
                          --resolution-upstream=250
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.log 
                          --output-all-profiles '''

    P.run()


###################################################################
@merge(getBinnedExpressionProfiles, "binned_expression_profiles.load")
def loadBinnedExpressionProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-FLAG).([^\.]+).(.)_outof_5\.gt_(.)_exons",
                         cat="factor,rep,quantile,exon_limit",
                         options="-i factor -i rep -i quantile -i bin -i exon_limit")


###################################################################
@follows(loadBinnedExpressionProfiles)
def quantileProfiles():
    pass


###################################################################
# Exon Junction Profiles
###################################################################
@follows(mkdir("transcriptome.dir"))
@transform([os.path.join(PARAMS["dir_transcriptome"],
                         "*-FLAG-R?.bwa.bam"),
            os.path.join(PARAMS["dir_ejc_transcriptome"],
                         "eIF4A3-GFP-R?.bwa.bam")],
           regex(".+/(.+-.+-R[^3])\.bwa\.bam"),
           r"transcriptome.dir/\1.bam")
def dedupAndSplitTranscriptome(infile, outfile):
    '''UMI dedup the transcriptome mapping bams and split multiple 
    alignments in the XA tag into seperate lines'''

    track = P.snip(outfile, ".bam")
    statement = '''samtools view -h %(infile)s
                 | xa2multi.pl
                 | samtools sort -o %(track)s.tmp.bam 
                                 
                   2>> %(track)s.log;

                   checkpoint;

                   samtools index %(track)s.tmp.bam;

                   checkpoint;

                   umi_tools dedup
                          -I %(track)s.tmp.bam
                          --method=directional
                          -L %(track)s.log
                 | samtools sort -o %(outfile)s;

                   checkpoint;

                   samtools index %(outfile)s

                   checkpoint;

                   rm %(track)s.tmp*'''

    P.run()


###################################################################
@collate(dedupAndSplitTranscriptome,
         regex("(.+-[^-]+)-(.+).bam"),
         r"\1-union.bam")
def getTranscritomeUnionBam(infiles, outfile):

    infiles = " ".join(infiles)

    statement = '''samtools merge
                   %(outfile)s
                   %(infiles)s;
 
                   checkpoint;

                   samtools index %(outfile)s'''


    P.run()


###################################################################
@follows(mkdir("transcriptome.dir"))
@transform(os.path.join(PARAMS["iclip_dir"],
                        "mapping.dir/geneset.dir/reference.gtf.gz"),
           regex(".+"),
           r"transcriptome.dir/reference_transcriptome.gtf.gz")
def convertGenesetToTranscriptomeCoords(infile, outfile):

    PipelineProj028.convertGenesetToTranscriptomeCoords(infile, outfile)


###################################################################
@transform([dedupAndSplitTranscriptome,
            getTranscritomeUnionBam],
           regex("(.+).bam"),
           add_inputs(convertGenesetToTranscriptomeCoords),
           r"\1.transcriptome.exon_profile.tsv.gz")
def getTranscriptomeExonBoundaryProfiles(infiles, outfile):

    bam, gtf = infiles

    if "GFP" in outfile:
        center=True
    else:
        center=False
        
    PipelineProj028.exonBoundaryProfiles(bam, gtf, outfile,
                                         center=center,
                                         submit=True,
                                         logfile=outfile + ".log",
                                         job_options="-l mem_free=10G")


###################################################################
@follows(mkdir("exon_boundary_profiles.dir"))
@transform([os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
            os.path.join(PARAMS["ejc_iclip_dir"],
                         "deduped.dir/*.union.bam"),
            "/ifs/projects/proj028/Cheng_ALYREF_iCLIP/deduped.dir/HeLa-Alyref*.bam"],
           regex(".+/(.+-.+-R[0-9]+|.+-.+.union).bam"),
           add_inputs(filterExpressedTranscripts),
           r"exon_boundary_profiles.dir/\1.exon_profile.tsv.gz")
def getExonBoundaryProfiles(infiles, outfile):

    bam, gtf = infiles

    if "GFP" in bam:
        center = True
    else:
        center = False
        
    PipelineProj028.exonBoundaryProfiles(bam, gtf, outfile,
                                         center=center,
                                         submit=True,
                                         logfile=outfile + ".log",
                                         job_memory="15G")


@follows(mkdir("exon_boundary_profiles.dir"))
@transform([os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
            os.path.join(PARAMS["ejc_iclip_dir"],
                         "deduped.dir/eIF4A3-GFP.union.bam")],
           regex(".+/(.+-.+-R[0-9]+|.+-.+.union).bam"),
           add_inputs(filterExpressedTranscripts),
           r"exon_boundary_profiles.dir/\1.centered.exon_profile.tsv.gz")
def getCenteredExonBoundaryProfiles(infiles, outfile):

    bam, gtf = infiles
      
    PipelineProj028.exonBoundaryProfiles(bam, gtf, outfile,
                                         center=True,
                                         submit=True,
                                         logfile=outfile + ".log",
                                         job_memory="15G")

    
###################################################################
@collate([getExonBoundaryProfiles,
          getCenteredExonBoundaryProfiles,
          getTranscriptomeExonBoundaryProfiles],
         regex("(.+)/(.+-.+-R[0-9]|.+-.+.union)\.(transcriptome\.exon|centered\.exon|exon)_profile.tsv.gz"),
         r"\1/\3_boundary_profiles.load")
def loadExonBoundaryProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-FLAG|.+-GFP|HeLa-Alyref).(R[0-9]|union)\..*exon_profile",
                         cat="track,replicate",
                         options="-i position -i replicate")


###################################################################
@follows(loadExonBoundaryProfiles)
def exonBoundaryProfiles():
    pass


###################################################################
# Circularisation Candidates
###################################################################
@follows(mkdir("profile_summaries.dir"))
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
           regex(".+/(.+-FLAG-R[0-9]+).bam"),
           add_inputs(filterExpressedTranscripts),
           r"profile_summaries.dir/\1.region_summary.tsv.gz")
def averageRegionsAndNorm(infiles, outfile):
    '''Average normalised counts over the regions of each transcript
    to identify transcripts strongly bound at the 3' or 5' utr '''

    bamfile, gtffile = infiles
    PipelineProj028.calculateRegionEnrichments(bamfile,
                                               gtffile,
                                               outfile,
                                               submit=True,
                                               job_options="-l mem_free=1")


###################################################################
@merge(averageRegionsAndNorm,
       "profile_summaries.load")
def concatenateAndLoadRegionSummaries(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="profile_summaries.dir/(.+)-FLAG-(R[0-9]+).region_summary.tsv.gz",
                         cat="Protein,Replicate",
                         options="-i transcript_id -i Protein -i Replicate")


###################################################################
@transform(concatenateAndLoadRegionSummaries,
           suffix(".load"),
           ".score.tsv.gz")
def scoreCircularCandidates(infile, outfile):

    PipelineProj028.scoreCircularCandidates(outfile,
                                            submit=True)


###################################################################
@transform(scoreCircularCandidates,
           suffix(".tsv.gz"),
           ".load")
def loadCirCandidateScores(infile, outfile):

    P.load(infile, outfile,
           "--header-names=transcript_id,score -i transcript_id")

###################################################################
@follows(loadCirCandidateScores)
def circular_candidates():
    pass


@follows(mkdir("intron_profiles.dir"))
@subdivide(os.path.join(
                PARAMS["iclip_dir"],
                "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           [r"intron_profiles.dir/\1.long.profile.tsv",
            r"intron_profiles.dir/\1.short.profile.tsv"])
def generateLongAndShortIntronProfiles(infiles, outfiles):
    ''' Divide introns into long and short and generate meta profiles 
    for each seperately'''

    bamfile, gtf_file = infiles
    PipelineProj028.intron_profiles(bamfile, gtf_file, outfiles,
                                    submit=True,
                                    job_memory="10G")


@merge(generateLongAndShortIntronProfiles,
       "intron_profiles.dir/intron_profiles.load")
def loadIntronProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename = ".+/(.+-FLAG).(.+).(short|long)",
                          cat="factor,replicate,length")


###################################################################
# Histone Profiles
###################################################################
@follows(mkdir("histones.dir"))
@transform(filterExpressedTranscripts,
           regex(".+"),
           add_inputs(os.path.join(
               PARAMS["dir_external"],
               "histone_genes.tsv")),
           r"histones.dir/histones.gtf.gz")
def getHistoneAnno(infiles, outfile):
    '''Get histones annotations from annotations file
    using list of gene_ids'''

    anno, genes = infiles

    statement = '''cgat gtf2gtf
                          --method=filter
                          --filter-method=transcript
                          --map-tsv-file <(cut -f2 %(genes)s |
                                           grep -P '^ENST')
                          -I %(anno)s
                          -S %(outfile)s
                          -L %(outfile)s.log'''

    P.run()


###################################################################
@transform(getSingleExonGeneModels,
           regex(".+"),
           add_inputs(os.path.join(
               PARAMS["dir_external"],
               "histone_genes.tsv")),
           r"histones.dir/non-histones.gtf.gz")
def getNonHistoneSingleExon(infiles, outfile):
    '''Get gene models for single exon genes which are not
    histones'''

    models, histones = infiles

    statement = '''cgat gtf2gtf
                          --method=filter
                          --filter-method=transcript
                          --invert-filter
                          --map-tsv-file <(cut -f2 %(histones)s |
                                           grep -P '^ENST' )
                           -I %(models)s
                           -S %(outfile)s
                           -L %(outfile)s.log'''

    P.run()


###################################################################
@product([os.path.join(
                PARAMS["iclip_dir"],
                "deduped.dir/*.bam"),
          os.path.join(
                PARAMS["ejc_iclip_dir"],
                "deduped.dir/*.bam")],
         formatter(".+/(.+).bam"),
         [getHistoneAnno, getNonHistoneSingleExon, getMultiExonGeneModels],
         formatter(".+/(.+).gtf.gz"),
         r"histones.dir/{1[0][0]}_vs_{1[1][0]}_profile.tsv.gz")
def doHistoneMetaGene(infiles, outfile):

    bam, anno = infiles

    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                      %(bam)s
                      -I %(anno)s
                      --exon-bins=10
                      --flanks=100
                      --flank-bins=10
                      --scale-flanks
                      -S %(outfile)s
                      -L %(outfile)s.log'''
    P.run()


###################################################################
@merge(doHistoneMetaGene, "histones.dir/histone_profiles.load")
def loadHistoneMetaGenes(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="histones.dir/(.+-FLAG|.+-GFP).(.+)_vs_(?:expressed_transcripts\.)?(.+)_profile.tsv",
                         cat="factor,replicate,geneset",
                         options="-i factor -i replicate -i geneset")


###################################################################
@files(PARAMS["external_histones"],
       "histone_genes.load")
def load_histone_genes(infile, outfile):

    P.load(infile, outfile)



###################################################################
# Heatmaps
###################################################################
@collate("*.bw",
         regex("(.+-.+)-R([0-9]+).bw"),
         r"\1-agg.bw")
def mergeRNASeqToBigWig(infiles, outfile):
    '''Merge replicate BigWig files so that they can be used for 
    normalisation'''

    if len(infiles) > 0:
        temp_bam = P.getTempFilename() + ".bam"
    
    
@transform(filterExpressedTranscripts,
           suffix(".gtf.gz"),
           ".protein_coding.no_overlapping.gtf.gz")
def filterGenesForHeatmaps(infile, outfile):
    '''Create set of genes suitable for heatmaps: only protein coding
    and remove small genes contained within larger ones'''

    statement=''' cgat gtf2gtf
                     -I %(infile)s
                     -L %(outfile)s.log
                     --method=filter
                     --filter-method=proteincoding
                | cgat gtf2gtf
                     -L %(outfile)s.log
                     -S %(outfile)s
                    --method=filter
                    --filter-method=longest-gene '''
    P.run()


###################################################################
@follows(mkdir("heatmaps"))
@subdivide([os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*union.bam"),
            os.path.join(PARAMS["ejc_iclip_dir"],
                         "deduped.dir/*union.bam")],
           regex(".+/(.+).bam"),
           add_inputs(filterGenesForHeatmaps),
           [r"heatmaps/\1.forward.matrix.tsv.gz",
            r"heatmaps/\1.reverse.matrix.tsv.gz"])
def generate_start_aligned_matricies(infiles, outfiles):
    
    statement_temp = '''
          python %(pipeline_dir)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                   -g %(gtffile)s
                   --no-plot
                    -u 1000
                    -n quantile
                   --outfile-prefix=%(outfile_p)s
                   -L %(outfile_p)s.log
                   %(bamfile)s
                   %(reverse)s'''

    pipeline_dir = PARAMS["project_src"]
    bamfile, gtffile = infiles
    statements = []
    for outfile in outfiles:
        outfile_p = P.snip(outfile, ".matrix.tsv.gz")
        if "reverse" in outfile:
            reverse = "--reverse-strand"
        else:
            reverse = ""

        statements.append(statement_temp % locals())

    P.run()


###################################################################
@transform(PARAMS["external_rnaseq_bw"],
           regex(".+/(.+).bw"),
           add_inputs(filterGenesForHeatmaps),
           r"heatmaps/\1.matrix.tsv.gz")
def generate_start_aligned_RNAseq_matrix(infiles, outfile):

    wigfile, gtffile = infiles
    statement = '''
          python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                   -g %(gtffile)s
                    -r quantile
                    -u 1000
                    -H 200
                    -n quantile
                   --outfile-prefix=%(outfile_p)s
                   -L %(outfile_p)s.log
                   --plus-wig=%(wigfile)s '''

    outfile_p = P.snip(outfile, ".matrix.tsv.gz")
    P.run()


###################################################################
@transform(generate_start_aligned_matricies,
           suffix(".matrix.tsv.gz"),
           add_inputs(filterGenesForHeatmaps),
           ".compressed.png")
def generate_start_aligned_plots(infiles, outfile):

    infile, annotations = infiles
    outfile_p = P.snip(outfile, ".png")
    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                    --use-matrix=%(infile)s
                     -g %(annotations)s
                    -H 200
                    -r quantile
                    
                    --crop=-1000:10000
                    --annotations=start,end
                    --outfile-prefix=%(outfile_p)s'''

    P.run()


###################################################################
@transform(generate_start_aligned_plots,
           regex("(.+-FLAG.+|.+-GFP.+).(forward|reverse).compressed.png"),
           inputs([r"\1.\2.compressed.matrix.tsv.gz",
                   generate_start_aligned_RNAseq_matrix]),  
           r"\1.\2.normed.png")
def RNAseq_norm_start_aligned(infiles, outfile):
    '''Normalized start aligned matricies to RNAseq'''

    iclip_matrix, rna_matrix = infiles
    
    rna_matrix = P.snip(rna_matrix, ".matrix.tsv.gz") + ".compressed.matrix.tsv.gz" 

    outfile_p = P.snip(outfile, ".png")

    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                     -r quantile
                     --crop=-975:10000
                     --outfile-prefix=%(outfile_p)s
                     -L %(outfile_p)s.log
                    --use-matrix=%(iclip_matrix)s
                    --norm-mat=%(rna_matrix)s '''

    P.run()


###################################################################
@follows(mkdir("heatmaps"))
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*union.bam"),
           regex(".+/(.+).bam"),
           add_inputs(filterGenesForHeatmaps),
           r"heatmaps/\1.end_aligned.matrix.tsv.gz")
def generate_end_aligned_matricies(infiles, outfile):
    
    statement= '''
          python %(pipeline_dir)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                   -g %(gtffile)s
                   -a end
                   -n quantile
                   --no-plot
                    -u 1000
                   --outfile-prefix=%(outfile_p)s
                   -L %(outfile_p)s.log
                   %(bamfile)s'''

    pipeline_dir = PARAMS["project_src"]
    bamfile, gtffile = infiles
    outfile_p = P.snip(outfile, ".matrix.tsv.gz")

    P.run()


###################################################################
@transform(PARAMS["external_rnaseq_bw"],
           regex(".+/(.+).bw"),
           add_inputs(filterGenesForHeatmaps),
           r"heatmaps/\1.end_aligned.matrix.tsv.gz")
def generate_end_aligned_RNAseq_matrix(infiles, outfile):

    wigfile, gtffile = infiles
    statement= '''
          python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                   -g %(gtffile)s
                   -a end
                   -H 200
                   -n quantile
                   -r quantile
                    -u 1000
                   --outfile-prefix=%(outfile_p)s
                   -L %(outfile_p)s.log
                   --plus-wig=%(wigfile)s'''

    outfile_p = P.snip(outfile, ".matrix.tsv.gz")
    P.run()


###################################################################
@subdivide(generate_end_aligned_matricies,
           regex("(.+).matrix.tsv.gz"),
           add_inputs(filterGenesForHeatmaps),
           [r"\1.length.compressed.png",
            r"\1.3utr.compressed.png"])
def generate_end_aligned_plots(infiles, outfiles):

    infile, annotations = infiles
   
    statement_t = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                    --use-matrix=%(infile)s
                     -g %(annotations)s
                    -H 200
                    --align-at=end  
                     -r quantile
                    
                    --crop=-10000:1000
                    --annotations=start,end
                    --outfile-prefix=%(outfile_p)s
                    -s %(sort)s'''

    statements = []
    project_src = PARAMS["project_src"]

    for outfile in outfiles:
        outfile_p = P.snip(outfile, ".png")
        sort = re.match(".+(length|3utr).compressed.png", outfile).groups()[0]
        statements.append(statement_t % locals())

    P.run()


###################################################################
@transform(generate_end_aligned_plots,
           regex("(.+-FLAG.+).(length|3utr).compressed.png"),
           inputs([r"\1.\2.compressed.matrix.tsv.gz",
                   generate_end_aligned_RNAseq_matrix]),
           r"\1.\2.normed.png")
def RNAseq_norm_end_aligned(infiles, outfile):
    '''Normalized start aligned matricies to RNAseq'''

    iclip_matrix, rna_matrix = infiles
    
    if "length" in iclip_matrix:
        align = "length"
    else:
        align = "3utr"


    outfile_p = P.snip(outfile, ".png")

    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                     -r quantile
                     --crop=-10000:1000
                     --outfile-prefix=%(outfile_p)s
                     -L %(outfile_p)s.log
                    --use-matrix=%(iclip_matrix)s
                    --norm-mat=%(rna_matrix)s '''

    P.run()


###################################################################
@transform(filterGenesForHeatmaps,
           suffix(".gtf.gz"),
           ".first_exon.gtf.gz")
def get_first_exons(infile, outfile):
    '''Get first exon per flattened gene'''

    PipelineProj028.get_first_exon(infile, outfile)


###################################################################
@transform(PARAMS["external_rnaseq_bw"],
           regex(".+/(.+).bw"),
           add_inputs(get_first_exons),
           r"heatmaps/\1.first_exon.matrix.tsv.gz")
def first_exon_rnaseq_matrix(infiles, outfile):

    inwig, annotations = infiles
    outfile_p = P.snip(outfile, ".matrix.tsv.gz")
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                    --plus-wig=%(inwig)s
                     -g %(annotations)s
                     -H 200
                     -b 5
                     -s length
                    --crop=-500:1500
                    -n quantile
                    -r quantile 
                   --outfile-prefix=%(outfile_p)s
                    '''

    P.run()


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*union.bam"),
           regex(".+/(.+).bam"),
           add_inputs(get_first_exons),
           r"heatmaps/\1.first_exon.matrix.tsv.gz")
def first_exon_heatmaps(infiles, outfile):

    infile, annotations = infiles
    outfile_p = P.snip(outfile, ".matrix.tsv.gz")
    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                    %(infile)s
                    -b 5
                    --crop=-500:1500
                     -g %(annotations)s
                    -H 200
                    -n quantile
                    -r quantile
                    --annotations=start,end
                    --outfile-prefix=%(outfile_p)s
                    -s length
                    -L %(outfile_p)s.log'''

    P.run()


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*union.bam"),
           regex(".+/(.+).bam"),
           add_inputs(get_first_exons),
           r"heatmaps/\1.first_exon.end.matrix.tsv.gz")
def first_exon_end_aligned_heatmaps(infiles, outfile):

    infile, annotations = infiles
    outfile_p = P.snip(outfile, ".matrix.tsv.gz")
    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                    %(infile)s
                    -b 5
                    --crop=-1500:500
                     -g %(annotations)s
                    --align-at=end
                    -H 200
                    -n quantile
                    -r quantile
                    --annotations=start,end
                    --outfile-prefix=%(outfile_p)s
                    -s length
                    -L %(outfile_p)s.log'''

    P.run()


###################################################################
@transform(first_exon_heatmaps,
           suffix(".matrix.tsv.gz"),
           add_inputs(first_exon_rnaseq_matrix),
           ".normed.matrix.tsv.gz")
def norm_first_exon_heatmaps(infiles, outfile):

    iclip_mat, rna_mat = infiles
    outfile_p = P.snip(outfile, ".matrix.tsv.gz")

    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2heatmap.py
                    --use-matrix=%(iclip_mat)s
                    --norm-mat=%(rna_mat)s
                    -H 200
                    -r quantile
                    --annotations=start,end
                    --outfile-prefix=%(outfile_p)s
                    -L %(outfile_p)s.log'''

    P.run()


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(get_first_exons),
           r"gene_profiles.dir/\1.profile.tsv.gz")
def first_exon_profiles(infiles, outfile):

    bamfile, gtffile = infiles
    prefix = P.snip(outfile, ".profile.tsv.gz")
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                          -I %(gtffile)s
                          %(bamfile)s
                          -b 50
                          -m %(prefix)s.matrix.tsv.gz
                  | gzip > %(prefix)s.profile.tsv.gz '''

    P.run()


###################################################################
@merge(first_exon_profiles, "gene_profiles.dir/first_exon_profiles.load")
def load_first_exon_profiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="gene_profiles.dir/(.+)-FLAG.(.+).profile",
                         cat="protein,replicate",
                         options="-i protein -i replicate")

###################################################################
@follows(generate_end_aligned_plots,
         generate_start_aligned_plots,
         RNAseq_norm_end_aligned,
         RNAseq_norm_start_aligned,
         first_exon_heatmaps,
         norm_first_exon_heatmaps)
def heatmaps():
    pass


###################################################################
ENCODE_eCLIP_Samples = {"K562_RBM15-eCLIP-R1": "ENCFF141MFP",
                        "K562_RBM15-eCLIP-R2": "ENCFF532QNF",
                        "K562_RBM15-Control-R1": "ENCFF413OBS",
                        "K562_CPSF6-eCLIP-R1": "ENCFF185XTD",
                        "K562_CPSF6-eCLIP-R2": "ENCFF550ABZ",
                        "K562_CPSF6-Control-R1": "ENCFF747XJC"}


GEO_PARCLIP_Samples = {"HEK293_CPSF6-PARCLIP-R1": "GSM917665_273",
                       "HEK293_CPSF6-PARCLIP-R2": "GSM917664_274"}

@follows(mkdir("encode_eclip.dir"))
@originate(["encode_eclip.dir/%s.bam" % f for f in ENCODE_eCLIP_Samples])
def download_eCLIP(outfile):
    '''Download eCLIP BAM files from ENCODE'''

    track = ENCODE_eCLIP_Samples[P.snip(
        os.path.basename(outfile), ".bam")]

    url = "https://www.encodeproject.org/files/%(track)s/@@download/%(track)s.bam" % locals()
    
    outfile = os.path.basename(outfile)

    statement = '''cd encode_eclip.dir;
                   checkpoint;
                   wget %(url)s --quiet;
                   checkpoint;
                   mv %(track)s.bam %(outfile)s;'''

    P.run()


###################################################################
@transform(download_eCLIP,
           suffix(".bam"),
           ".bam.bai")
def index_encode_eclip(infile, outfile):

    statement = '''samtools index %(infile)s'''

    P.run()


###################################################################
@follows(index_encode_eclip)
@collate(download_eCLIP,
         regex("(.+-.+)-(R.).bam"),
         r"\1.union.bam")
def merge_encode_replicates(infiles, outfile):
    '''Create union tracks for the ENCODE_eCLIP_Samples'''

    if len(infiles) == 1:
        IOTools.cloneFile(infiles[0], outfile)
        return

    infiles = " ".join(infiles)
    statement = ''' samtools merge -f %(outfile)s %(infiles)s '''
    P.run()


###################################################################
@transform(merge_encode_replicates,
           suffix(".bam"),
           ".bam.bai")
def index_encode_union(infile, outfile):

    statement = '''samtools index %(infile)s'''
    P.run()


###################################################################
@originate(["encode_eclip.dir/%s_+.wig.gz" % f for f in GEO_PARCLIP_Samples] +
           ["encode_eclip.dir/%s_-.wig.gz" % f for f in GEO_PARCLIP_Samples])
def download_cpsf_parclip(outfile):

    out_track = re.match(
        "encode_eclip.dir/(.+)_[\+\-].wig.gz", outfile).groups()[0]
    track = GEO_PARCLIP_Samples[out_track]
    geo = track.split("_")[0]
    geo_short = geo[:6]
    infile = re.sub(out_track,
                    track,
                    os.path.basename(outfile))
    infile = re.sub(".wig.gz", "_density.wig.gz", infile)
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/%(geo_short)snnn/%(geo)s/suppl/%(infile)s" % locals()
    outfile = os.path.basename(outfile)

    statement = ''' cd encode_eclip.dir;
                    checkpoint;
                    wget --quiet %(url)s;
                    checkpoint;
                    mv %(infile)s %(outfile)s'''
    P.run()
    
    
###################################################################
@collate(download_cpsf_parclip,
         regex("(.+_.+)-R(.)_(\+|\-).wig.gz"),
         r"\1.union.bw")
def merge_cpsf_wigs(infiles, outfile):
    ''' merge the cpsf wig files into a single bigWig'''

    unzipped = [P.snip(inf, ".gz") for inf in infiles]
    unzip_statement = ["zcat %s > %s" % (inf, outf)
                       for inf, outf in zip(infiles, unzipped)]
    unzip_statement = "; checkpoint; ".join(unzip_statement)

    outfile = P.snip(outfile, ".bw")
    infiles = " ".join(unzipped)

    genome_file = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_contigs"])

    statement = '''%(unzip_statement)s;
                   checkpoint;
                   wiggletools write %(outfile)s.wig sum %(infiles)s;
                   checkpoint;
                   wigToBigWig %(outfile)s.wig %(genome_file)s %(outfile)s.bw;
                   checkpoint;
                   rm %(infiles)s %(outfile)s.wig'''
    P.run()

         
###################################################################
@follows(index_encode_union)
@transform([merge_encode_replicates,
            merge_cpsf_wigs],
           regex("(.+).union.(bam|bw)"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"\1.geneprofilewithintrons.matrix.tsv.gz")
def do_encode_eclip_metagenes(infiles, outfile):
    ''' Run bam to gene profile metagenes for ENCODE eCLIP data'''

    bamfile, gtffile = infiles
    
    if "bam" in bamfile:
        input = "--bam-file=%s" % bamfile
    else:
        input = "--bigwigfile=%s" % bamfile

    outfile = P.snip(outfile, ".geneprofilewithintrons.matrix.tsv.gz")

    statement = '''cgat bam2geneprofile
                           --method=geneprofilewithintrons
                           %(input)s
                           --gtf-file=%(gtffile)s
                           --normalize-transcript=total-sum
                           --normalize-profile=area
                           --log=%(outfile)s.log
                           --output-filename-pattern=%(outfile)s.%%s
                           > %(outfile)s '''

    P.run()


###################################################################
@merge(do_encode_eclip_metagenes,
       "encode_eclip.dir/encode_eclip_metagenes.load")
def load_encode_eclip(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_(.+)-(.+).geneprofile",
                         cat="cell_type,factor,condition",
                         options="-i cell_type -i factor -i condition")


###################################################################
@follows(load_encode_eclip)
def encode_eclip():
    pass


###################################################################
# First exon/last exon profiles
###################################################################
@transform(filterExpressedTranscripts,
           formatter(),
           "exon_stats.tsv")
def find_exon_ratios(infile, outfile):
    '''Find average length of first exons to internal exons
    to last exons'''

    PipelineProj028.find_exon_ratios(infile, outfile)


@transform([os.path.join(PARAMS["ejc_iclip_dir"],"deduped.dir/*.bam"),
            os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
            "*.bam"],
           formatter(),
           add_inputs(filterExpressedTranscripts, find_exon_ratios),
           r"gene_profiles.dir/{basename[0]}.separateexonprofile.log")
def seperate_exon_profile(infiles, outfile):
    '''Build metagene profile with seperate divisions for first, last
    and middle exons. Will do this with bam2geneprofile. Should work
    and making iCLIP_bam2geneprofile.pydo this also. Currently run
    using protein coding genes only'''

    job_memory="15G"
    bamfile, gtffile, ratiofile = infiles
    outfile = P.snip(outfile, ".separateexonprofile.log")

    with IOTools.openFile(ratiofile) as rfile:
        ratios = map(float, rfile.readlines()[0].split("\t"))

    first_res, mid_res, last_res = [int(r*250) for r in ratios]
 
    statement = '''cgat bam2geneprofile
                          --method=separateexonprofile
                          --reporter=gene
                          --bam-file=%(bamfile)s
                          --gtf-file=<(zcat %(gtffile)s | grep protein_codin)
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --normalize-profile=area
                          --scale-flank-length=1
                          --resolution-first-exon=%(first_res)s
                          --resolution-cds=%(mid_res)s
                          --resolution-last-exon=%(last_res)s
                          --resolution-downstream=250
                          --resolution-upstream=250
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.separateexonprofile.log 
                          '''

    
    P.run()


###################################################################
@merge(seperate_exon_profile, "gene_profiles.dir/seperate_exon_profiles.load")
def load_seperate_exon_profiles(infiles, outfile):

    infiles = [re.sub(".log", ".matrix.tsv.gz", inf) for inf in infiles
               if re.match("gene_profiles.dir/(.+).(union|R.)", inf)]
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="gene_profiles.dir/(.+).(union|R.)",
                         cat="factor,replicate",
                         options="-i factor -i replicate -i base")


###################################################################
@follows(mkdir("firstlastexon.dir"))
@transform([os.path.join(PARAMS["ejc_iclip_dir"], "deduped.dir/*.bam"),
            os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam")],
           formatter(),
           add_inputs(filterExpressedTranscripts),
           "firstlastexon.dir/{basename[0]}.tsv.gz")
def get_first_last_exon_counts(infiles, outfile):

    bamfile, gtffile = infiles
    PipelineProj028.get_first_last_counts(bamfile,
                                          gtffile,
                                          outfile,
                                          submit=True)


@follows(mkdir("firstlastexon.dir"))
@transform([os.path.join(PARAMS["ejc_iclip_dir"], "deduped.dir/*.bam"),
            os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam")],
           formatter(),
           add_inputs(getSingleExonGeneModels),
           "firstlastexon.dir/{basename[0]}_single_exons.tsv.gz")
def get_single_exon_counts(infiles, outfile):
    '''Count number of clips in single exon genes using count_clip_sites.
    In theory this is redundant with the featureCounts based counting above
    but this is just more convenient'''

    bamfile, gtffile = infiles
    use_centre = "GFP" in bamfile

    statement = '''python %(project_src)s/iCLIPlib/scripts/count_clip_sites.py
                          -I %(gtffile)s
                          %(use_centre)s
                          %(bamfile)s
                          -L %(outfile)s.log
                          -S %(outfile)s '''

    P.run()


@merge(get_single_exon_counts, "firstlastexon.dir/single_exon_counts.load")
def load_single_exon_counts(infiles, outfile):
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="firstlastexon.dir/(.+)-(.+).(union|R.)",
                         cat = "protein,cell,replicate",
                         options="-i protein -i gene_id -i cell -i replicate")


@merge(get_first_last_exon_counts,
       "firstlastexon.dir/first_last_exon_counts.load")
def load_first_last_exon_counts(infiles, outfile):
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="firstlastexon.dir/(.+)-(.+).(union|R.)",
                         cat="protein,cell,replicate",
                         options="-i protein -i gene_id -i cell -i replicate")


###################################################################
@follows(load_seperate_exon_profiles,
         load_first_last_exon_counts)
def seperate_exon_profiles():
    pass


###################################################################
@follows(singleVsMultiExonGenes,
         normalised_profiles,
         quantileProfiles,
         circular_candidates,
         exonBoundaryProfiles,
         encode_eclip)
def profiles():
    pass

###################################################################
###################################################################
###################################################################
###################################################################

###################################################################
# Clipped Transcripts
###################################################################

@transform(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+"),
           r"all_gene_intervals.gtf.gz")
def mergeGeneStructures(infile, outfile):
    '''Take the complete ensembl geneset and create single intervals'''

    statement = '''cgat gtf2gtf -m merge-transcripts
                         -I %(infile)s -S %(outfile)s -L %(outfile)s.log '''

    P.run()
    
###################################################################
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"\1.count.tsv.gz")
def count_clipped_transcripts(infiles, outfile):
    '''count number of clip tags in each transript. Tags may be in introns,
    or exons, or take complete interval over gene area. Make stranded and
    allow overlaps'''

    bamfile, annotations = infiles

    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=1,
        options=PARAMS['featurecounts_options'] + ' -O')


###################################################################
@transform("*.bam",
           regex("(.+).bam"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"\1.count.tsv.gz")
def countRNAExpression(infiles, outfile):
    ''' count expressoin with feature counts'''

    bamfile, annotations = infiles
        
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=PARAMS['featurecounts_strand'],
        options=PARAMS['featurecounts_options'])


###################################################################
@merge([countRNAExpression, count_clipped_transcripts],
       "track_counts.tsv.gz")
def mergeCounts(infiles, outfile):
    '''Merge feature counts data into one table'''

    infiles = " ".join(infiles)

    statement=''' cgat combine_tables
                         -c 1
                         -k 7
                         --regex-filename='(.+).count.tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(mergeCounts, suffix(".tsv.gz"), ".load")
def loadCounts(infile, outfile):

    P.load(infile, outfile, options="-i geneid")


###################################################################
@follows(mkdir("clusters.dir"))
@collate(os.path.join(PARAMS["iclip_dir"], "clusters.dir/*.bed.gz"),
         regex(".+/(.+-FLAG)-R[0-9].bed.gz"),
         r"clusters.dir/\1.union.bed.gz")
def getUnionClusters(infiles, outfile):

    PipelineiCLIP.callReproducibleClusters(infiles, outfile, 1)


###################################################################
@follows(mkdir("cluster_hits"))
@transform([os.path.join(PARAMS["iclip_dir"], "clusters.dir/*.bed.gz"),
            getUnionClusters],
           regex(".+/(.+).bed.gz"),
           add_inputs(mergeGeneStructures),
           r"cluster_hits/\1.gene_hits.tsv.gz")
def countClustersOverlapingGenes(infiles, outfile):
    
    clusters, genes = infiles 

    statement= '''cgatgff2bed
                         -I %(genes)s
                         -L %(outfile)s
                | bedtools coverage
                  -a %(clusters)s -b stdin
                  -counts -s
                | cut -f4,7
                | gzip > %(outfile)s'''

    P.run()


###################################################################
@merge(countClustersOverlapingGenes,
       "cluster_hits/gene_cluster_counts.tsv.gz")
def combineGeneClusterCounts(infiles, outfile):

    infiles = " ".join(infiles)
    statement = ''' cgat combine_tables
                         -c 1
                         -t
                         --regex-filename='(.+).gene_hits.tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(combineGeneClusterCounts, suffix(".tsv.gz"), ".load")
def loadGeneClusterCounts(infile, outfile):

    P.load(infile, outfile, options = "-i ID")


###################################################################
@follows(mkdir("transcript_chunks.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex("(.+)"),
                 r"transcript_chunks.dir/reference_chunks.gtf.gz")
def getTranscriptChunks(infile, outfile):
    '''Take reference transcriptome and transform so that each part of an
    exon or intron is an "exon" '''


    statement=''' cgat gtf2gtf --method=genes-to-unique-chunks
                  -I %(infile)s
                  -L %(outfile)s.log
                  -S %(outfile)s'''

    P.run()


###################################################################
@transform([os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
            "*.bam",
            os.path.join(PARAMS["dir_external"],
                         "stubbsRNAseq/*.bam"),
            os.path.join(PARAMS["dir_ejc_iclip"], "deduped.dir/*.bam"),
            merge_encode_replicates],
           regex("(?:.+/)?([-\w]+(?:\.union|\.reproducible)?).+"),
           add_inputs(getTranscriptChunks),
           r"transcript_chunks.dir/\1.chunk_counts.tsv.gz")
def countChunks(infiles, outfile):
    ''' Use feature counts to count the number of reads in each 
    chunk, both for RNA and iCLIP'''

    bamfile, gtffile = infiles
    statement = '''python %(project_src)s/iCLIPlib/scripts/count_clip_sites.py
                              %(bamfile)s
                              -I %(gtffile)s
                              --feature=exon
                              -S %(outfile)s'''

    P.run()



###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_exons_gtf"])),
           ".introns.load")
def annotateIntronChunks(infiles, outfile):
    '''Find those transcript chunks that introns. Find consituative introns by looking
    for chunks with intron>=1 and exon==0'''

    chunks, transcripts = infiles

    statement = '''  cgat gtf2gtf
                        -I %(transcripts)s
                        -L %(outfile)s.log
                        --method=set-gene-to-transcript
                 | cgat gtf2gtf
                        -L %(outfile)s.log
                        --method=exons2introns
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\tintron'
                 | %(load_smt)s > %(outfile)s'''
    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options = "-i gene_id -i exon_id")
    
    P.run()

    
###################################################################    
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(filterExpressedTranscripts),
           ".expressed.introns.load")
def annotateExpressedIntronChunks(infiles, outfile):
    '''Find those transcript chunks that introns. Find consituative introns by looking
    for chunks with intron>=1 and exon==0'''

    chunks, transcripts = infiles

    statement = '''  cgat gtf2gtf
                        -I %(transcripts)s
                        -L %(outfile)s.log
                        --method=set-gene-to-transcript
                 | cgat gtf2gtf
                        -L %(outfile)s.log
                        --method=exons2introns
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\tintron'
                 | %(load_smt)s > %(outfile)s'''
    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options = "-i gene_id -i exon_id")
    
    P.run()

    
###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           ".constitive_exons.load")
def annotateConstituativeExonsChunks(infiles, outfile):
    '''Find those transcript chunks that represent consituative exons'''

    chunks, transcripts = infiles

    statement = ''' cgat gtf2gtf
                           -I %(transcripts)s
                           -L %(outfile)s.log
                          --method=intersect-transcripts
                         
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\texon'
                 | %(load_smt)s > %(outfile)s'''
    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options = "-i gene_id -i exon_id")

    P.run()

    
###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           ".exons.load")
def annotateAnyExonsChunks(infiles, outfile):
    '''Find those transcript chunks that represent consituative exons'''

    chunks, transcripts = infiles

    statement = '''zcat %(transcripts)s
                 | awk '$3=="exon"'                        
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\texon'
                 | %(load_smt)s > %(outfile)s'''
    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options = "-i gene_id -i exon_id")

    P.run()

###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(filterExpressedTranscripts),
           ".expressed.exons.load")
def annotateExpressedExonsChunks(infiles, outfile):
    '''Find those transcript chunks that represent consituative exons'''

    chunks, transcripts = infiles

    statement = '''zcat %(transcripts)s
                 | awk '$3=="exon"'                        
                 | bedtools intersect -a %(chunks)s -b - -c -s
                 | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                 | sed '1i gene_id\\texon_id\\texon'
                 | %(load_smt)s > %(outfile)s'''
    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options = "-i gene_id -i exon_id")

    P.run()

    
###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["dir_external"], "sharp_detained_introns.bed.gz")),
           ".detained.load")
def annotateDetainedChunks(infiles, outfile):

    chunks, introns = infiles

    statement = ''' bedtools intersect -a %(chunks)s -b %(introns)s -c
                  | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+([0-9]+)$/\\1\\t\\2\\t\\3/' 
                  | sed '1i gene_id\\texon_id\\texon'
                  | %(load_smt)s > %(outfile)s'''
    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options="-i gene_id -i exon_id")
    

    P.run()

    
###################################################################
@follows(mkdir("transcript_chunks.dir"))
@transform([os.path.join(PARAMS["annotations_dir"],
                         PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
            filterExpressedTranscripts],
           regex("(?:.+/)?(.+).gtf.gz"),
           r"transcript_chunks.dir/\1.last_pc_exons.gtf.gz")
def get_transcript_last_pc_exons(infile, outfile):

    PipelineProj028.get_last_pc_exons(infile, outfile)

    
###################################################################
@follows(mkdir("transcript_chunks.dir"))
@transform([os.path.join(PARAMS["annotations_dir"],
                         PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
            filterExpressedTranscripts],
           regex("(?:.+/)?(.+).gtf.gz"),
           r"transcript_chunks.dir/\1.first_pc_exons.gtf.gz")
def get_transcript_first_pc_exons(infile, outfile):

    PipelineProj028.get_first_pc_exons(infile, outfile)

    
###################################################################
@follows(mkdir("transcript_chunks.dir"))
@transform([os.path.join(PARAMS["annotations_dir"],
                         PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
            filterExpressedTranscripts],
           regex("(?:.+/)?(.+).gtf.gz"),
           r"transcript_chunks.dir/\1.gene_last_pc_exons.gtf.gz")
def get_gene_last_pc_exons(infile, outfile):

    tmpfile = P.getTempFilename(dir=".")
    statement = '''cgat gtf2gtf --method=merge-exons
                                -I %(infile)s
                                -S %(tmpfile)s'''
    P.run()
    
    PipelineProj028.get_last_pc_exons(tmpfile, outfile)

    os.unlink(tmpfile)
    
    
###################################################################
@follows(mkdir("transcript_chunks.dir"))
@transform([os.path.join(PARAMS["annotations_dir"],
                         PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
            filterExpressedTranscripts],
           regex("(?:.+/)?(.+).gtf.gz"),
           r"transcript_chunks.dir/\1.gene_first_pc_exons.gtf.gz")
def get_gene_first_pc_exons(infile, outfile):

    tmpfile = P.getTempFilename(dir=".")
    statement = '''cgat gtf2gtf --method=merge-exons
                                -I %(infile)s
                                -S %(tmpfile)s'''
    P.run()
    
    PipelineProj028.get_first_pc_exons(tmpfile, outfile)

    os.unlink(tmpfile)

    
###################################################################
@transform([get_transcript_last_pc_exons,
            get_gene_last_pc_exons,
            get_transcript_first_pc_exons,
            get_gene_first_pc_exons],
           formatter("(?P<track>expressed_transcripts\..+_pc_exons|(?:gene_)?(?:first|last)_pc_exons).gtf.gz"),
           add_inputs(getTranscriptChunks),
           r"reference_chunks.{track[0]}.load")
def annotateFirstLastExons(infiles, outfile):

    exons, chunks = infiles

    if "first" in exons:
        colname = "first"
    elif "last" in exons:
        colname = "last"
    else:
        raise ValueError(exons)
    
    statement = '''bedtools intersect -a %(chunks)s -b %(exons)s -c
                  | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/' 
                  | sed '1i gene_id\\texon_id\\t%(colname)s'
                  | %(load_smt)s > %(outfile)s'''

    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options="-i gene_id -i exon_id")
    P.run()


################################################################
@follows(mkdir("retained_introns.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_exons_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"retained_introns.dir/\1.retained_introns.bed.gz")
def getRetainedIntrons(infile, outfile):

    statement = '''
                cgat gtf2gtf -I %(infile)s
                                       -m find-retained-introns
                | cgat gff2bed --is-gtf
                | gzip > %(outfile)s '''

    P.run()

    
###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_flat_gtf"])),
           ".nGenes.load")
def find_n_genes_overlapping_chunks(infiles, outfile):
    '''For each chunk, find the number of genes that overlap it'''

    chunks, genes = infiles
    statement = ''' cgat gtf2gtf -I %(genes)s --method=merge-transcripts -L %(outfile)s.log
                 | sort -k1,1 -k2,2n
                 | bedtools intersect -a %(chunks)s -b - -c
                 | sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/'
                 | sed '1i gene_id\\texon_id\\tnGenes'
                 | %(load_smt)s > %(outfile)s'''

    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options="-i gene_id -i exon_id")
    P.run()

    
###################################################################   
@follows(mkdir("transcript_chunks.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_exons_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"transcript_chunks.dir/\1.retained_introns.gtf.gz")
def getMatchedRetainedIntrons(infile, outfile):

    PipelineProj028.getRetainedIntronsFromMatchedTranscripts(infile, outfile)


###################################################################
@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(getMatchedRetainedIntrons),
           ".retained_introns.load")
def annotateRetainedIntronChunks(infiles, outfile):

    chunks, introns = infiles

    statement = '''bedtools intersect -a %(chunks)s -b %(introns)s -s -c
                 |  sed -E 's/.+gene_id \\"(ENS[A-Z]+[0-9]+)\\".+exon_id \\"([0-9]+)\\".+\\t([0-9]+)$/\\1\\t\\2\\t\\3/'
                 | sed '1i gene_id\\texon_id\\tretained'
                 | %(load_smt)s > %(outfile)s'''

    tablename = P.toTable(outfile)
    load_smt = P.build_load_statement(tablename=tablename,
                                      options="-i gene_id -i exon_id")
    P.run()

    connect().executescript('''DROP INDEX IF EXISTS %(tablename)s_joint;
                               CREATE INDEX %(tablename)s_joint ON
                                     %(tablename)s(gene_id, exon_id)'''
                            % locals())

    
###################################################################
@merge(countChunks, "transcript_chunks.dir/chunk_counts.load")
def loadChunkCounts(infiles, outfile):
    
    infiles = " ".join(infiles)

    statement = ''' cgat combine_tables
                           -k 4 -c 1,3
                           --use-file-prefix
                           --regex-filename="(.+).chunk_counts.tsv.gz"
                          %(infiles)s
                  | sed 's/gene_id_exon_id/gene_id\\texon_id/'
                  | sed 's/-/\\t/'
                  | %(load_stmt)s > %(outfile)s'''

    tablename = P.toTable(outfile)
    load_stmt = P.build_load_statement(tablename=tablename,
                                       options="-i gene_id -i exon_id")
    job_memory = "5G"
    P.run()


###################################################################
@jobs_limit(1)
@transform([loadChunkCounts, annotateConstituativeExonsChunks,
            annotateAnyExonsChunks, annotateIntronChunks,
            annotateDetainedChunks, annotateFirstLastExons,
            find_n_genes_overlapping_chunks,
            annotateExpressedExonsChunks,
            annotateExpressedIntronChunks],
           suffix(".load"),
           ".index")
def jointIndexOnChunks(infile, outfile):

    cc = connect()

    table = P.toTable(infile)
    cc.executescript(
            '''DROP INDEX IF EXISTS %(table)s_joint;
               CREATE INDEX %(table)s_joint ON %(table)s(gene_id,exon_id)'''
            % locals())

    P.touch(outfile)


###################################################################
@follows(jointIndexOnChunks, annotateRetainedIntronChunks)
def transcript_chunks():
    pass

##################################################################
# DEXseq on chunks
#################################################################
###################################################################
@follows(mkdir("chunk_splicing.dir"))
@transform(os.path.join(PARAMS["dir_external"],
                         "stubbsRNAseq/*.bam"),
           
           regex(".+/(.+_.+_R.).+bam"),
           add_inputs(getTranscriptChunks),
           r"chunk_splicing.dir/\1.chunk_feature_counts.tsv.gz")
def feature_count_chunks(infiles, outfile):
    '''Use feature counts to count number of RNA seq reads in each chunk'''

    bamfile, annotations = infiles
    
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=0,
        options='-O --primary -M --minOverlap 10 -f')

    
###################################################################
@transform([PARAMS["external_detained_calling_bams"],
            PARAMS["dir_chtop_rnai_bams"]],
           regex(".+/([^\.]+).*bam"),
           add_inputs(getTranscriptChunks),
           r"chunk_splicing.dir/\1.chunk_uniq_counts.tsv.gz")
def feature_count_chunks_uniq(infiles, outfile):
    
    bamfile, annotations = infiles
    
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=0,
        options='-O --minOverlap 10 -f -s 2 -p')

    
 ###################################################################   
@collate([feature_count_chunks,
          feature_count_chunks_uniq],
         regex("chunk_splicing.dir/(.+).chunk_(.+)_counts.tsv.gz"),
         r"chunk_splicing.dir/chunk_\2_counts.tsv.gz")
def merge_feature_count_chunks(infiles, outfile):

    infiles = " ".join(infiles)

    statement=''' cgat combine_tables
                         -c 1,2,3,4,5,6
                         -k 7
                         --regex-filename='(.+).chunk_.+_counts.tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()

    
###################################################################
@transform(getTranscriptChunks,
           regex(".+/(.+).gtf.gz"),
           add_inputs(PARAMS["external_mappability"],
                      annotateIntronChunks),
           r"chunk_splicing.dir/intron_lengths.tsv.gz")
def calculate_introns_lengths(infiles, outfile):

    chunks, mappability, _ = infiles
    statement='''python %(project_src)s/calculate_effective_length.py
                        -I %(chunks)s
                        -l 100
                        -b %(mappability)s
                        -d %(database_name)s
                        -S %(outfile)s'''
    P.run()

    
###################################################################
@transform(calculate_introns_lengths,
           suffix(".tsv.gz"),
           ".load")
def load_intron_lengths(infile, outfile):

    P.load(infile, outfile)

    connect().executescript('''DROP INDEX IF EXISTS intron_length_joint;
                               CREATE INDEX intron_length_joint
                                   ON intron_lengths(gene_id, exon_id)''')

####################################################################
@transform(merge_feature_count_chunks,
           regex("(.+)uniq(.+)"),
           add_inputs(load_intron_lengths),
           r"chunk_splicing.dir/detained_intron_calls.tsv")
def call_detained_introns(infile, outfile):

    basedir = os.path.abspath(PARAMS["workingdir"])
    statement = ''' Rscript -e 'library(knitr); 
                                knitr::opts_knit$set(root.dir="%(basedir)s");
                                knit("%(project_src)s/find_detained_introns.Rmd") ' '''
    P.run()


@transform(call_detained_introns,
           suffix(".tsv"),
           ".load")
def load_detained_intron_calls(infile, outfile):

    P.load(infile, outfile, "-i gene_id -i intron_id -i fraction")

    connect().executescript('''DROP INDEX IF EXISTS detained_intron_joint;
                               CREATE INDEX detained_intron_joint
                                  ON detained_intron_calls(gene_id, intron_id)''')

    
###################################################################
@follows(load_intron_lengths, load_detained_intron_calls, annotateRetainedIntronChunks)
def find_detained_introns():
    pass

###################################################################
@transform("*.ri_design.tsv",
           regex("(.+).ri_design.tsv"),
           add_inputs(merge_feature_count_chunks,
                      getTranscriptChunks),
           r"chunk_splicing.dir/\1.dexseq.tsv")
def runDEXSeqOnChunks(infiles, outfile):
    ''' Run DEXSeq using the stubbs RNAseq and the transcript
    modeles with retained introns '''

    if "chtop" not in infiles[0]:
        infiles = [f for f in infiles if 'uniq' not in f]
    else:
        infiles = [infiles[0]] + [f for f in infiles if "uniq" in f] + [infiles[-1]]
        
    design, counts, models = infiles

    
    infiles = ",".join([models, counts, design])
    outfile = P.snip(outfile, ".tsv")

    job_threads = 6
    job_memory="10G"

    statement = ''' Rscript %(project_src)s/run_dexseq_all.R 
                            --infiles %(infiles)s
                            --outfiles %(outfile)s.tsv,%(outfile)s.gtf.gz,%(outfile)s.RData
                    &> %(outfile)s.log '''

    P.run()

    
###################################################################
@merge(runDEXSeqOnChunks,
       "chunk_splicing.dir/chunk_splicing.load")
def load_dexseq_on_chunks(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="chunk_splicing.dir/(.+).dexseq",
                         options = "-i track -i groupID -i featureID")

    connect().executescript('''DROP INDEX IF EXISTS chunk_splicing_joint;
                               CREATE INDEX chunk_splicing_join ON
                                  chunk_splicing(groupID, featureID)''')


    
###################################################################
@follows(mkdir("exon_intron_ratio.dir"))
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*union.bam"),
           formatter(),
           add_inputs(filterGenesForHeatmaps),
           r"exon_intron_ratio.dir/{basename[0]}.ratios.tsv.gz")
def getExonIntronRatios(infiles, outfile):
    '''Get the length normalised ratio of reads in the first intron
    and the second exon'''

    bamfile, gtffile = infiles

    temp=P.getTempFilename(shared=True)
    statement = '''cgat gtf2gtf
                            -I %(gtffile)s -S %(temp)s
                            --method=sort
                            --sort-order=gene+transcript'''
    P.run()

    PipelineProj028.ExonRatio(bamfile, temp, outfile,
                              submit=True)

    os.unlink(temp)

###################################################################
@merge(getExonIntronRatios,
       "exon_intron_ratio.dir/exon_intron_ratio.load")
def load_exon_intron_ratio(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).union.ratios.tsv.gz")


###################################################################
@follows(load_exon_intron_ratio)
def exon_intron_ratio():
    pass


###################################################################
@follows(loadCounts, loadGeneClusterCounts, exon_intron_ratio)
def transcript_counts():
    pass


###################################################################
# GO
###################################################################
@follows(mkdir("go"), mkdir("go/clipped_genes"), loadSalmon)
@transform(os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(loadCounts),
           r"go/clipped_genes/\1.goseq.tsv")
def runClippedGenesGO(infiles, outfile):
    '''Look for GO enrichments in genes clipped by each factor'''

    track = re.match(".+/(.+).bam", infiles[0]).groups()[0]
    plot_out = P.snip(outfile, ".goseq.tsv") + ".pwf.png"
    column = P.tablequote(track)

    statement = '''SELECT Geneid,
                          %s as clipped,
                          SUM(TPM) as expression
                   FROM track_counts
                    INNER JOIN annotations.transcript_info as ti
                    ON ti.gene_id = track_counts.Geneid
                    INNER JOIN HEK293_salmon as sf
                    ON ti.transcript_id = sf.Name
                   GROUP BY Geneid''' % column

    results = DUtils.fetch_DataFrame(statement,connect())
    results.clipped[results.clipped > 0] = 1
    results = results.set_index("Geneid")

    PipelineProj028.runGOSeq(genes=results.clipped,
                             exprs=results.expression,
                             outfile=outfile,
                             pwf_plot=plot_out,
                             submit=True,
                             logfile=outfile + ".log")


###################################################################
@follows(loadGeneClusterCounts, mkdir("go"), mkdir("go/clusters"))
@transform(countClustersOverlapingGenes,
           regex("cluster_hits/(.+).gene_hits.tsv.gz"),
           add_inputs(loadGeneClusterCounts),
           r"go/clusters/\1.goseq.tsv")
def runClustersGo(infiles, outfile):
    '''Test go enrichments of genes hit by each of the cluster tracks'''

    track = re.match(
        "cluster_hits/(.+).gene_hits.tsv.gz", infiles[0]).groups()[0]
    column = P.tablequote(track)

    statement = '''SELECT ID,
                          clusters.%s as clipped,
                          HEK293_WT_1 as expression
                   FROM
                       gene_cluster_counts as clusters 
                     INNER JOIN track_counts as counts
                     ON clusters.ID = counts.Geneid ''' % column

    results = DUtils.fetch_DataFrame(statement, connect())

    results = results.set_index("ID")

    genes = results.clipped
    genes[genes > 0] = 1

    PipelineProj028.runGOSeq(genes,
                             results.expression,
                             outfile,
                             "go/clusters/%s.pwf.png" % track,
                             submit=True,
                             logfile=outfile + ".log")


##################################################################
@follows(runClippedGenesGO, runClustersGo)
def GO():
    pass


##################################################################
#PIPE-CLIP
###################################################################
@follows(mkdir("pipeclip.dir"))
@transform(os.path.join(PARAMS["iclip_dir"],"deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           r"pipeclip.dir/\1.pipeclip_out.tsv.gz")
def runPipeClip(infile, outfile):
    ''' Run PIPE-CLIP pipeline for comparison purposes '''


    outfile = P.snip(outfile, ".gz")
    statement='''module load bio/PIPE-CLIP;

                 checkpoint;
 
                 pipeclip -i %(infile)s
                          -o %(outfile)s
                          -l 10
                          -m 100
                          -c 3
                          -M 0.05
                          -C 0.05
                          -r 0;
         
                checkpoint;

                gzip %(outfile)s '''

    P.run()


##################################################################
#Cluster optimisation
##################################################################
@follows(mkdir("window_size_scan.dir"))
@subdivide(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*-FLAG-R*.bam"),
       regex(".+/(.+).bam"),
       add_inputs(os.path.join(PARAMS["iclip_dir"],
                               "mapping.dir/geneset.dir/reference.gtf.gz")),
       [r"window_size_scan.dir/\1.%i.bed.gz" % (2**size)
        for size in range(8)])
def scanWindowSizes(infiles, outfiles):
    '''Run find_significant_bases on a range of different cluster
    window sizes to find optimal'''

    bam, gtf = infiles

    tempfile = P.getTempFilename(".")

    statement = ''' zcat %(gtf)s |
                     awk '$1=="chr22"' |
                     gzip > %(tempfile)s.gtf.gz'''
    P.run()

    for outf in outfiles:
        window_size = int(
            re.match(".+/.+-FLAG-R[0-9]+\.([0-9]+).bed.gz",
                     outf).groups()[0]
        )
        bg_file = P.snip(outf, ".bed.gz") + ".bg.gz"
        PipelineiCLIP.callClusters(bam, tempfile + ".gtf.gz", [bg_file, outf],
                                   window_size=window_size)

    os.unlink(tempfile + ".gtf.gz")


##################################################################
@collate(scanWindowSizes,
         regex("(.+-FLAG)-(R[0-9]+)\.([0-9]+).bed.gz"),
         r"\1.reproducible.\3.bed.gz")
def reproducibleScannedWindows(infiles, outfile):

    PipelineiCLIP.callReproducibleClusters(infiles, outfile,
                                           2)


##################################################################
@transform([scanWindowSizes, reproducibleScannedWindows],
           suffix(".bed.gz"),
           ".count")
def countCalledWindows(infile, outfile):

    statement = "zcat %(infile)s | wc -l > %(outfile)s"
    P.run()


##################################################################
@merge(countCalledWindows, "window_scan.load")
def loadWindowCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-FLAG.+)\.([0-9]+)\.count",
                         header="track,size,count",
                         cat="track,size",
                         has_titles=False)


##################################################################
@follows(mkdir("p_threhold_scan.dir"))
@subdivide(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*-FLAG-R*.bam"),
       regex(".+/(.+).bam"),
       add_inputs(os.path.join(PARAMS["iclip_dir"],
                               "mapping.dir/geneset.dir/reference.gtf.gz")),
       [r"p_threhold_scan.dir/\1.%s.bed.gz" % p
        for p in [0.001,0.01,0.05,0.1,0.2]])
def scanPThreshes(infiles, outfiles):
    '''Run find_significant_bases on a range of different cluster
    window sizes to find optimal'''

    bam, gtf = infiles

    tempfile = P.getTempFilename(".")

    statement = ''' zcat %(gtf)s |
                     awk '$1=="chr22"' |
                     gzip > %(tempfile)s.gtf.gz'''
    P.run()

    for outf in outfiles:
        threshold = float(
            re.match(".+/.+-FLAG-R[0-9]+\.(0\.[0-9]+).bed.gz",
                     outf).groups()[0]
        )
        bg_file = P.snip(outf, ".bed.gz") + ".bg.gz"
        PipelineiCLIP.callClusters(bam, tempfile + ".gtf.gz", [bg_file, outf],
                                   pthresh=threshold)

    os.unlink(tempfile + ".gtf.gz")


##################################################################
@collate(scanPThreshes,
         regex("(.+-FLAG)-(R[0-9]+)\.(0\.[0-9]+).bed.gz"),
         r"\1.reproducible.\3.bed.gz")
def reproducibleScannedPs(infiles, outfile):

    PipelineiCLIP.callReproducibleClusters(infiles, outfile,
                                           2)


##################################################################
@collate(scanPThreshes,
         regex("(.+-FLAG)-(R[0-9]+)\.(0\.[0-9]+).bed.gz"),
         r"\1.all.\3.bed.gz")
def allWindowsAtP(infiles, outfile):

    PipelineiCLIP.callReproducibleClusters(infiles, outfile,
                                           1)


##################################################################
@transform([scanPThreshes, reproducibleScannedPs,
            allWindowsAtP],
           suffix(".bed.gz"),
           ".count")
def countCalledPs(infile, outfile):

    statement = "zcat %(infile)s | wc -l > %(outfile)s"
    P.run()


##################################################################
@merge(countCalledPs, "p_threshold_scan.load")
def loadPCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-FLAG)[-|\.](.+)\.(0\.[0-9]+)\.count",
                         header="track,replicate,p,count",
                         cat="track,replicate,p",
                         has_titles=False)


##################################################################
@follows(loadWindowCounts, loadPCounts)
def clusterOptimisaiton():
    pass


##################################################################
@follows(mkdir("rand_clusts.dir"))
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
           formatter(),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"rand_clusts.dir/{basename[0]}.bedGraph.gz")
def call_clusts_by_rand(infiles, outfile):
    '''Call clusters using randomisation - as in previous analyses'''

    bamfile, gtffile = infiles
    job_threads = 6
    job_memory = "0.5G"
    statement = '''python %(project_src)s/iCLIPlib/scripts/significant_bases_by_randomisation.py
                       -I %(gtffile)s
                       -b %(bamfile)s
                       -p %(job_threads)s
                        -S %(outfile)s'''

    P.run()


###################################################################
@transform(call_clusts_by_rand,
           regex("(.+/.+).bedGraph.gz"),
           r"\1.merged_clusters.bed.gz")
def merge_adjacent_clusters(infile, outfile):
    '''Merge bases called as significant if their territories overlap'''

    genome = os.path.join(PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_contigs"])

    statement = '''bedtools slop -b 15 -i %(infile)s -g %(genome)s
                 | sort -k1,1 -k2,2n
                 | bedtools merge -i -
                 | gzip > %(outfile)s'''

    P.run()


@permutations(merge_adjacent_clusters,
              formatter(".+/(?P<track>.+-FLAG-R.).merged_clusters.bed.gz"),
              2,
              r"rand_clusts.dir/{track[0][0]}_v_{track[1][0]}.jaccard.tsv")
def assess_rand_clusters(infiles, outfile):
    '''Use bedtools jaccard to assess the reproducibility of random clusters'''

    first, second = infiles
    statement = '''bedtools jaccard -a %(first)s -b %(second)s > %(outfile)s'''
    P.run()


@merge(assess_rand_clusters, "rand_clusts.dir/rand_clusts_reproducibility.load")
def load_rand_clust_reproducibility(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         "rand_clusts.dir/(.+)-FLAG-(R.)_v_(.+)-FLAG-(R.)",
                         cat="first_factor,first_rep,second_factor,second_rep")


@permutations(os.path.join(PARAMS["iclip_dir"],
                           "clusters.dir/*-FLAG-R?.bed.gz"),
              formatter(".+/(?P<track>.+-FLAG-R.).bed.gz"),
              2,
              r"clusters.dir/{track[0][0]}_v_{track[1][0]}.jaccard.tsv")
def assess_binom_clusters(infiles, outfile):

    first, second = infiles
    statement = '''bedtools jaccard -a <(zcat %(first)s | sort -k1,1 -k2,2n)
                                    -b <(zcat %(second)s |sort -k1,1 -k2,2n)
                                    -split > %(outfile)s '''
    P.run()


@merge(assess_binom_clusters, "clusters.dir/binom_clusts_reproducibility.load")
def load_binom_clust_reproducibility(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         "clusters.dir/(.+-FLAG-R.)_v_(.+-FLAG-R.)",
                         cat="first_sample,second_sample")
##################################################################
# No FlipIn clusters
##################################################################
@transform(os.path.join(PARAMS["iclip_dir"], "clusters.dir",
                        "*.reproducible.bed.gz"),
           regex(".+/clusters.dir/([^(?:FlipIn)].+-FLAG.reproducible).bed.gz"),
           add_inputs(os.path.join(PARAMS["iclip_dir"], "clusters.dir",
                                   "FlipIn-FLAG.reproducible.bed.gz")),
           r"clusters.dir/\1.no_flipin.bed.gz")
def removeFlipInFromClusters(infiles, outfile):
    ''' Remove clusters that overlap with a reproducible
    FlipIn cluster '''

    sample, flipin = infiles
    statement = '''bedtools intersect -a %(sample)s -b %(flipin)s
                                      -v
                 | sort -k1,1 -k2,2n
                 | gzip > %(outfile)s '''

    P.run()


##################################################################
@merge(removeFlipInFromClusters,
       "clusters.dir/no_flipin_cluster_counts.tsv")
def countNoFlipInClusters(infiles, outfile):

    with IOTools.openFile(outfile, "w") as outf:
        outf.write("track\tcount\n")
        for infile in infiles:
            track = re.match(".+/(.+).no_flipin.bed.gz", infile).groups()[0]
            num = IOTools.getNumLines(infile)

            outf.write("%s\t%s\n" % (track, num))


##################################################################
@transform(countNoFlipInClusters, suffix(".tsv"),
           ".load")
def loadNoFlipInClusterCounts(infile, outfile):

    P.load(infile, outfile)


##################################################################
@transform(removeFlipInFromClusters,
           suffix(".bed.gz"),
           add_inputs(os.path.join(PARAMS["iclip_dir"], "reference.context.bed.gz")),
           ".context_stats")
def getNoFlipInContextStats(infiles, outfile):

    clusters, context = infiles
    tmp = P.getTempFilename

    statement = '''cgat bam_vs_bed
                          -a %(clusters)s
                          -b %(context)s
                           -S %(outfile)s
                           -L %(outfile)s.log '''

    P.run()


##################################################################
@merge(getNoFlipInContextStats, "clusters.dir/no_flipin_context_stats.load")
def loadNoFlipInContextStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)-FLAG.reproducible.no_flipin.context_stats",
                         options = "-i category")


##################################################################
@follows(loadNoFlipInContextStats,
         loadNoFlipInClusterCounts)
def no_flipin_clusters():
    pass


##################################################################
# Methylation
##################################################################
@transform(PARAMS["external_methylation_sites"],
           regex("(.+)"),
           "methylation_sites.bed")
def liftOverMethylationSites(infile, outfile):

    PipelineProj028.liftOverFromHg18(infile, outfile)

##################################################################
@transform(liftOverMethylationSites,
           regex("(.+)"),
           add_inputs(filterExpressedTranscripts),
           "methylated_genes.tsv.gz")
def getGenesWithMethylation(infiles, outfile):

    methylation, reference = infiles

    statement = '''zcat %(reference)s
                  
                  | cgat gff2bed --is-gtf
                  | bedtools intersect -a stdin -b %(methylation)s -wa
                  | cut -f4
                  | sort -u
                  | gzip > %(outfile)s'''
    P.run()


###################################################################
@transform(getGenesWithMethylation, suffix(".tsv.gz"),
           ".load")
def loadMethylatedGenes(infile, outfile):

    P.load(infile, outfile, options="-H gene_id -i gene_id")


####################################################################
@transform(PARAMS["external_darnell_methylation"],
           regex("(.+)"),
           "darnell_methylation_sites.bed")
def liftOverDarnellMeth(infile, outfile):

    PipelineProj028.liftOverFromHg18(infile, outfile)


###################################################################
@follows(mkdir("methylation.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex("(.+)"),
           r"methylation.dir/internal_exons.gtf.gz")
def getInternalExons(infile, outfile):

    PipelineProj028.getInternalExons(infile, outfile,
                                     submit=True)


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(liftOverDarnellMeth, getInternalExons),
           r"methyaltion.dir/\1.distances.tsv.gz")
def findMinDistanceMeth(infiles, outfile):

    bamfile, gtffile, bedfile = infiles

    statement = '''zcat %(gtffile)s
                 | python %(project_src)s/iCLIPlib/scripts/reproducibility_by_exon.py
                   --method=min_dist
                   %(bamfile)s %(bedfile)s -L %(outfile)s.log
                   -S %(outfile)s '''

    P.run()


##################################################################
@follows(loadMethylatedGenes)
def methylation():
    pass 

###################################################################
###################################################################
# tRNAs
###################################################################
@follows(mkdir("trnas.dir"))
@subdivide(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*-FLAG-R*.bam"),
       regex(".+/(.+).bam"),
       add_inputs(PARAMS["trna_annotations"]),
       [r"trnas.dir/\1.clusters.bg.gz",
        r"trnas.dir/\1.clusters.bed.gz"])
def calltRNAClusters(infiles, outfiles):
    '''Call significant clusters on tRNAs'''

    bam, gtf = infiles
    PipelineiCLIP.callClusters(bam, gtf, outfiles)


###################################################################
@collate(calltRNAClusters,
         regex("(.+/.+\-FLAG)-(.+)\.bed.gz"),
         r"\1.reproducible.bed.gz")
def callReproducibletRNAClusters(infiles, outfile):
    
    PipelineiCLIP.callReproducibleClusters(infiles, outfile,
                                           PARAMS["clusters_min_reproducible"])
   

###################################################################
@transform(callReproducibletRNAClusters,
           regex(".+/([^(?:FlipIn)].+).reproducible.bed.gz"),
           add_inputs("trnas.dir/FlipIn-FLAG.reproducible.bed.gz"),
           r"trnas.dir/\1.no_flipin.bed.gz")
def removeInputOverlappingClusters(infiles, outfile):
    '''Remove reproducible clusters that overlap with reproducible
    input clusters '''

    sample, control = infiles
    PipelineiCLIP.removeInputOverlappingClusters(sample, control,
                                                 outfile)

    
###################################################################
@merge(removeInputOverlappingClusters,
       "trnas.dir/trnas_with_two_factors.tsv.gz")
def findtRNAsWithTwoFactors(infiles, outfile):
    '''Find clusters that are present in at least two factors'''

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s
                 | sort -k1,1 -k2,2n
                 | cgat bed2bed
                          --method=merge
                          --merge-min-intervals=2
                          --merge-and-resolve-blocks
                          -L %(outfile)s.log
                 | cut -f4,5 
               
                 | sed -r 's/_[0-9]+//'
                 | sed '1i tRNA\\tnum_factors'
                 | gzip > %(outfile)s'''
    P.run()


###################################################################
@transform(findtRNAsWithTwoFactors, suffix(".tsv.gz"),
           ".load")
def loadtRNAsWithTwoFactors(infile, outfile):

    P.load(infile, outfile, "-i tRNA -i num_factors")


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*-FLAG-R*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(PARAMS["trna_annotations"]),
           r"trnas.dir/\1.counts.tsv.gz")
def countSitesOntRNAs(infiles, outfile):
    '''Count clip sites on tRNAs'''

    bamfile, gtf = infiles
    logfile = P.snip(outfile,".tsv.gz") + ".log"
    statement = '''python %(project_src)s/iCLIPlib/scripts/count_clip_sites.py
                          -I %(gtf)s
                          -L %(logfile)s
                          %(bamfile)s
                  | gzip > %(outfile)s '''
    P.run()


###################################################################
@merge(countSitesOntRNAs,
       "trnas.dir/tRNA_counts.load")
def loadtRNACounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="trnas.dir/(.+-FLAG).(.+).counts.tsv.gz",
                         cat = "factor,replicate",
                         options = "-i tRNA -i track")


###################################################################
@transform(os.path.join(PARAMS["trna_annotations"]),
           regex(".+/(.+).gtf.gz"),
           r"trna_stats.tsv.gz")
def gettRNAStats(infile, outfile):

    logfile = P.snip(outfile, ".tsv.gz") + ".log"
    statement = '''cgat gtf2table
                            -I %(infile)s
                            -L %(logfile)s.log
                            --counter=position
                            --counter=length
                | gzip > %(outfile)s '''
    P.run()


###################################################################
@transform(gettRNAStats,
           suffix(".tsv.gz"),
           ".load")
def loadtRNAStats(infile, outfile):

    P.load(infile, outfile)


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*-FLAG*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(PARAMS["trna_annotations"]),
           r"trnas.dir/\1.profile.tsv.gz")
def perBasetRNAProfiles(infiles, outfile):
    '''Calculate a per base profile accross tRNAs '''

    bam_file, gtf_file = infiles

    PipelineProj028.tRNABaseCounts(gtf_file, bam_file,
                                   outfile, submit=True)


###################################################################
@merge(perBasetRNAProfiles,
       "trnas.dir/trna_profiles.load")
def loadtRNAProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename="trnas.dir/(.+).profile.tsv.gz")       


###################################################################
@follows(loadtRNACounts,
         loadtRNAStats,
         loadtRNAsWithTwoFactors,
         loadtRNAProfiles)
def tRNAs():
    pass


##################################################################
# Motifs
##################################################################

###################################################################
# Retained Introns
###################################################################
@transform(os.path.join(PARAMS["iclip_dir"],
                        "mapping.dir/geneset.dir/reference.gtf.gz"),
           regex(".+/(.+).gtf.gz"),
           r"biotypes.tsv.gz")
def getRetainedIntronGeneIDs(infile, outfile):

    PipelineProj028.getTranscriptTypeByGeneID(infile,
                                              outfile,
                                              submit=True)


###################################################################
@transform(getRetainedIntronGeneIDs,
           suffix(".tsv.gz"),
           ".load")
def loadBiotypes(infile, outfile):

    P.load(infile, outfile, options="-i gene_id -i transcript_id -i biotype")


###################################################################
@transform(PARAMS["external_up_regulated_genes"],
           regex(".+/(.+).txt"),
           r"\1.load")
def loadRegulatedGenes(infile,outfile):

    P.load(infile, outfile,
           options = "-i EnsemblID -i EntrezGeneID -i ProbeName")


###################################################################
@transform(PARAMS["external_array_platform"],
           regex(".+/(.+).tsv"),
           r"array_platform.load")
def loadArrayPlatform(infile, outfile):
    ''' Need total coverage for backgruond calculation '''
    P.load(infile, outfile,
           options = "-i EntrezGeneID")


###################################################################
@transform([os.path.join(PARAMS["iclip_dir"],
                        "clusters.dir/*.bed.gz"),
            getUnionClusters],
           regex(".+/(.+).bed.gz"),
           add_inputs(getRetainedIntrons),
           r"retained_introns.dir/\1.clusters.bed.gz")
def getClustersInRetainedIntrons(infiles, outfile):
    ''' find the set of clusters that overlap
    retained introns'''

    clusters, introns = infiles

    statement = '''bedtools intersect -a %(clusters)s -b %(introns)s
                                     -split -s
                  | sort -k1,1 -k2,2n
                  | bedtools merge -i stdin -s -c 4,5,6 -o collapse,sum,distinct
                  | gzip > %(outfile)s '''

    P.run()



###################################################################
@collate(os.path.join(PARAMS["iclip_dir"],
                      "deduped.dir/*-FLAG-R*.bam"),
         regex("(.+)-FLAG-(R.+).bam"),
         r"\1-FLAG.union.bam")
def getUnionBams(infiles, outfile):
    '''Merge all bams for a factor to make a union bam,
    link it to a reproducible bam, which will be identical '''

    outfile = os.path.abspath(outfile)
    factor = P.snip(outfile, ".union.bam")

    infiles = " ".join(infiles)

    statement = ''' samtools merge %(outfile)s %(infiles)s;
                    checkpoint;
                    ln -f -s %(outfile)s %(factor)s.reproducible.bam;
                    checkpoint;
                    samtools index %(outfile)s;
                    checkpoint;
                    samtools index %(factor)s.reproducible.bam '''

    P.run()
      
     
###################################################################
@follows(mkdir("retained_introns.dir/zagros.dir"),
         getUnionBams)
@transform(getClustersInRetainedIntrons,
           regex(".+/(.+).clusters.bed.gz"),
           add_inputs(getRetainedIntrons,
                      os.path.join(PARAMS["iclip_dir"],
                                   r"deduped.dir/\1.bam")),
           r"retained_introns.dir/zagros.dir/\1.bed.gz")
def getZagrosRIClustersAndDes(infiles, outfile):
    ''' Will output a bed file with entry all the same length, and
    truncated as intron boundaries, also outputs a zagros DE file '''

    clusters, introns, bamfile = infiles

    outfiles = [outfile, P.snip(outfile, ".bed.gz") + ".des"]

    PipelineProj028.getZagrosRIInputFiles(
        clusters,
        introns,
        bamfile,
        PARAMS["zagros_cluster_size"],
        outfiles,
        submit=True)


###################################################################
@transform(getZagrosRIClustersAndDes,
           suffix(".bed.gz"),
           ".fa")
def getRIZagrosFasta(infile, outfile):
    
    statement= ''' cgat bed2fasta
                         -g %(genome_dir)s/%(genome)s
                         -m dustmasker
                         -I %(infile)s
                         --use-strand
                         
                         -L %(outfile)s.log
                  | sed 's/[ |\:]/_/g' > %(outfile)s '''

    P.run()


###################################################################
@transform(getRIZagrosFasta,
           suffix(".fa"),
           ".str")
def getRIZagrosStructureFile(infile, outfile):
    
    statement = '''thermo -o %(outfile)s %(infile)s &> %(outfile)s.log'''
    P.run()


###################################################################
@transform(getRIZagrosStructureFile,
           regex("(.+).str"),
           add_inputs(r"\1.fa", r"\1.des"),
           r"\1.zagros")
def runRIZagros(infiles, outfile):

    structure, sequence, des = infiles
    statement = '''zagros
                      -t %(structure)s
                      -d %(des)s
                      %(sequence)s
                      %(zagros_options)s >
                    %(outfile)s '''

    P.run()


###################################################################
@follows(runRIZagros)
def RIZagros():
    pass


@transform(filterExpressedTranscripts,
           regex(".+"),
           r"expressed.3utrs.gtf.gz")
def get3UTRs(infile, outfile):

    PipelineProj028.get3UTRs(infile, outfile)


###################################################################
@follows(mkdir("utr_motifs.dir"))
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(get3UTRs),
           r"utr_motifs.dir/\1.6mers.tsv.gz")
def findUTRenrichedHexamers(infiles, outfile):
    ''' Use clip site randomisation to measure enrichment of pentamers '''

    bamfile, gtffile = infiles
    genome = os.path.join(PARAMS["genome_dir"],
                          PARAMS["genome"] + ".fasta")
    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_kmer_enrichment.py
                      -f %(genome)s
                      -b %(bamfile)s
                      -k 6
                      -p 6
                      -I %(gtffile)s
                      -S %(outfile)s
                      -L %(outfile)s.log'''

    job_threads = 6
    P.run()


###################################################################

###################################################################
@follows(mkdir("utr_motifs.dir"))
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(get3UTRs),
           r"utr_motifs.dir/\1.8mers.tsv.gz")
def findUTRenrichedOctamers(infiles, outfile):
    ''' Use clip site randomisation to measure enrichment of pentamers '''

    bamfile, gtffile = infiles
    genome = os.path.join(PARAMS["genome_dir"],
                          PARAMS["genome"] + ".fasta")
    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_kmer_enrichment.py
                      -f %(genome)s
                      -b %(bamfile)s
                      -k 8
                      -p 6
                      -I %(gtffile)s
                      -S %(outfile)s
                      -L %(outfile)s.log'''

    job_threads = 6
    job_memory = "0.3G"

    P.run()


###################################################################
@follows(mkdir("single_exon_motifs.dir"))
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(getSingleExonGeneModels),
           r"single_exon_motifs.dir/\1.tsv.gz")
def find_single_exon_kmers(infiles, outfile):

    bamfile, gtffile = infiles
    genome = os.path.join(PARAMS["genome_dir"],
                          PARAMS["genome"] + ".fasta")

    statement = ''' python %(project_src)s/iCLIPlib/scripts/iCLIP_kmer_enrichment.py
                      -f %(genome)s
                      -b %(bamfile)s
                      -k 6
                      -p 6
                      -I %(gtffile)s
                      -S %(outfile)s
                      -L %(outfile)s.log'''


    job_threads = 6
    job_memory = "0.3G"

    P.run()


@merge(find_single_exon_kmers, "single_exon_kmers.load")
def load_single_exon_kmers(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename="single_exon_motifs.dir/(.+)-FLAG.(.+).tsv.gz",
                          cat="factor,replicate",
                          options = "-i factor -i replicate -i kmer")

###################################################################
@merge([findUTRenrichedHexamers,
        findUTRenrichedOctamers],
       "utr_motifs.dir/utr_kmers.load")
def loadUTRPentamers(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).([0-9]+)mers.tsv.gz",
                         cat="track,k")


###################################################################
@follows(loadUTRPentamers,
         load_single_exon_kmers)
def pentamers():
    pass


###################################################################
@transform([os.path.join(PARAMS["iclip_dir"], "clusters.dir/*.bed.gz"),
            getUnionClusters],
           regex(".+/(.+).bed.gz"),
           add_inputs(os.path.join(PARAMS["iclip_dir"],
                                   "mapping.dir/geneset.dir/refcoding.introns.gtf.gz")),
           r"clusters.dir/\1_introns.bed.gz")
def getIntronClusters(infiles, outfile):

    clusters, geneset = infiles

    statement = '''
                  bedtools intersect -a %(clusters)s -b %(geneset)s
                      -wa -f 0.5 -split
                 | gzip > %(outfile)s'''

    P.run()


###################################################################
@transform([os.path.join(PARAMS["iclip_dir"], "clusters.dir/*.bed.gz"),
            getUnionClusters],
           regex(".+/(.+).bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"clusters.dir/\1_single_exon.bed.gz")
def getSingleExonClusters(infiles, outfile):
    
    clusters, geneset = infiles

    statement = ''' SELECT DISTINCT es.transcript_id as id
                    FROM annotations.transcript_stats as es
                    INNER JOIN annotations.transcript_info as ti
                    ON es.transcript_id = ti.transcript_id
                    WHERE contig != 'chrM'
                    GROUP BY ti.gene_id
                    HAVING MAX(nval) = 1 '''

    tmp = P.getTempFilename(shared=True)
    IOTools.writeLines(tmp, DUtils.fetch(statement, connect()))

    statement = ''' cgat gtf2gtf
                         -I %(geneset)s 
                         --method=filter
                         --filter-method=transcript
                         -a %(tmp)s
                         -L %(outfile)s.log
                 | sort -k1,1 -k4,4n
                 | bedtools intersect -a %(clusters)s -b stdin
                        -wa -f 0.5 -split
                 | gzip > %(outfile)s'''

    P.run()
    os.unlink(tmp)


###################################################################
@transform(getClustersInRetainedIntrons,
           regex("(.+).clusters.+"),
           r"\1.genes.load")
def loadGenesWithRetainedIntronClusters(infile, outfile):

    tablename = P.toTable(outfile)
    statement = ''' zcat %(infile)s
                 | perl -lane '/(ENS[A-Z]+[0-9]+)/ && print $1'
                 | cgat csv2db 
                   --table=%(tablename)s 
                   --database=%(database)s
                   --header-names=gene_id
                    -i gene_ids >  %(outfile)s '''

    P.run()


@follows(getSingleExonClusters, getIntronClusters,
         getClustersInRetainedIntrons)
def interval_sets():
    pass


###################################################################
###################################################################
@follows(mkdir("retained_introns.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+"),
           r"retained_introns.dir/reference_models_with_ri.gtf.gz")
def getRetainedIntronModelsWithRIasExons(infile, outfile):

    PipelineProj028.getTranscriptsPlusRetainedIntrons(infile,
                                                      outfile)


###################################################################
@follows(mkdir("retained_introns.dir"))
@transform(os.path.join(
               PARAMS["dir_external"], "sharp_detained_introns.bed.gz"),
           regex(".+"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           "retained_introns.dir/detained_introns_in_gene_models.gtf.gz")
def insertDetainedIntronsIntoTranscripts(infiles, outfile):

    introns, transcripts = infiles
    PipelineProj028.insertDetainedIntronsIntoTranscripts(transcripts, introns,
                                                         outfile)


###################################################################
@merge([insertDetainedIntronsIntoTranscripts,
        getRetainedIntronModelsWithRIasExons],
       "retained_introns.dir/merged_de_re_tained_intron_models.gtf.gz")
def mergeDetainedAndRetainedIntrons(infiles, outfile):

    tmp = P.getTempFilename(".")

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s
                 | cgat gtf2gtf --method=sort
                                                    --sort-order=gene+transcript
                                                     -S %(tmp)s'''

    P.run()

    PipelineProj028.mergeDetainedAndRetainedIntrons(tmp, outfile)

    os.unlink(tmp)


###################################################################
@follows(mkdir("retained_introns.dir/exon_counts.dir"))
@transform(os.path.join(PARAMS["dir_external"],
                        "stubbsRNAseq/*.bam"),
           regex(".+/(.+_.+_R[0-9]).+.bam"),
           add_inputs(mergeDetainedAndRetainedIntrons),
           r"retained_introns.dir/exon_counts.dir/\1.counts.tsv.gz")
def countRetainedIntronExons(infiles, outfile):
    '''Count stubbs RNAseq data accross exons of genes with retained 
    introns '''

    bamfile, annotations = infiles

    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=0,
        options=' -O -M  --minOverlap 5 -f')


###################################################################
@merge(countRetainedIntronExons,
       "retained_introns.dir/exon_counts.dir/retained_intron_exon_counts.tsv.gz")
def mergeRetainedIntronExonCounts(infiles, outfile):

    infiles = " ".join(infiles)

    statement=''' cgat combine_tables
                         -c 1,2,3,4,5,6
                         -k 7
                         --regex-filename='(.+).counts.tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()
       

###################################################################
@transform(mergeRetainedIntronExonCounts,
           suffix(".tsv.gz"),
           ".load")
def loadRetainedIntronExonCounts(infile, outfile):

    P.load(infile, outfile,
           "-i Geneid -i Chr -i Start -i End -i strand -i length")


###################################################################
@transform("*.ri_design.tsv",
           regex("^(?!chtop)(.+).ri_design.tsv"),
           add_inputs(mergeRetainedIntronExonCounts,
                      mergeDetainedAndRetainedIntrons),
           r"retained_introns.dir/\1.dexseq.tsv")
def runDEXSeqOnRI(infiles, outfile):
    ''' Run DEXSeq using the stubbs RNAseq and the transcript
    modeles with retained introns '''

    design, counts, models = infiles

    infiles = ",".join([models, counts, design])
    outfile = P.snip(outfile, ".tsv")

    job_threads = 6
    job_memory="10G"

    statement = ''' Rscript %(project_src)s/run_dexseq.R 
                            --infiles %(infiles)s
                            --outfiles %(outfile)s.tsv,%(outfile)s.gtf.gz,%(outfile)s.RData
                    &> %(outfile)s.log '''

    P.run()


###################################################################
@jobs_limit(1)
@transform(runDEXSeqOnRI, suffix(".tsv"),
           ".load")
def loadDEXSeq(infiles, outfile):

    P.load(infiles, outfile, options = "-i groupID -i featureID -i padj")

    tablename = P.toTable(outfile)
    connect().executescript(''' DROP index if exists %(tablename)s_join_index;
                                CREATE INDEX %(tablename)s_join_index
                                   ON %(tablename)s(groupID,featureID)'''
                            % locals())


###################################################################
@transform(mergeDetainedAndRetainedIntrons,
           regex(".+"),
           "retained_introns.dir/re_de_detained_introns.gtf.gz")
def extractRetainedIntrons(infile, outfile):

    statement = ''' zcat %(infile)s
                 | grep -P 'exon_id \\"I[0-9]+\\";'
                 | gzip > %(outfile)s'''

    P.run()


###################################################################
@transform(os.path.join(PARAMS["iclip_dir"],
                        "deduped.dir/*.bam"),
           regex(".+/(.+-FLAG.(?:R[0-9]+|union)).bam"),
           add_inputs(extractRetainedIntrons),
           r"retained_introns.dir/\1.retained_intron_counts.tsv.gz")
def quantifyRetainedIntronCLIPTags(infiles, outfile):
    '''Count the number of CLIP tags in each reatined intron'''

    bamfile, gtffile = infiles
    statement = '''python %(project_src)s/iCLIPlib/scripts/count_clip_sites.py
                              %(bamfile)s
                              -I %(gtffile)s
                              --feature=exon
                              -S %(outfile)s'''

    P.run()


###################################################################
@merge(quantifyRetainedIntronCLIPTags,
       "retained_introns.dir/retained_intron_clip_tag_counts.load")
def loadRetainedIntronCLIPTagCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).retained_intron",
                         options="-i gene_id -i transcript_id -i exon_id")

    statement = '''DROP INDEX IF EXISTS retained_intron_clip_tag_counts_joint_index;
                   CREATE INDEX retained_intron_clip_tag_counts_joint_index 
                          ON retained_intron_clip_tag_counts(gene_id, exon_id);'''

    connect().executescript(statement)


###################################################################
@transform(loadDEXSeq,
           regex("retained_introns.dir/(.+).dexseq.load"),
           inputs([r"retained_introns.dir/\1.dexseq.gtf.gz",
                   extractRetainedIntrons,
                   os.path.join(
                       PARAMS["dir_external"], "sharp_detained_introns.bed.gz")]),
           r"retained_introns.dir/\1.detained_fraction.tsv")
def getDetainedIntronFraction(infiles, outfile):

    test_set, background, category = infiles
    PipelineProj028.getOverlapFractions(test_set=test_set,
                                        category=category,
                                        background=background,
                                        outfile=outfile)


###################################################################
@merge(getDetainedIntronFraction,
       "retained_introns.dir/detained_intron_fractions.load")
def loadDetainedIntronFractions(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).detained_fraction.tsv")


###################################################################
@follows(loadDEXSeq,
         loadRetainedIntronCLIPTagCounts)
def differential_intron_usage():
    pass


###################################################################
@follows(loadRegulatedGenes, loadArrayPlatform,
         loadBiotypes,  differential_intron_usage,
         loadDetainedIntronFractions)
def retained_introns():
    pass



###################################################################
@follows(interval_sets,
         pentamers)
        
def motifs():
    pass


###################################################################
# NSUN6
###################################################################
@follows(mkdir("nsun6.dir"))
@transform([os.path.join(PARAMS["iclip_dir"],
                         "clusters.dir",
                         "*.bed.gz"),
            "clusters.dir/*.union.bed.gz"],
           regex(".+/([^/]+).bed.gz"),
           r"nsun6.dir/\1.region.bed.gz")
def getclusterRegions(infile, outfile):
    '''Get a bed region a fixed amount around the centre of 
    each cluster'''

    PipelineProj028.extendBedIntervals(infile, outfile, 100)


###################################################################
@transform(os.path.join(PARAMS["nsun6_dir"], "clusters.dir",
                        "NSun6-Endo.reproducible.bed.gz"),
           regex(".+/(.+).bed.gz"),
           r"\1.bed.gz")
def recompressAndIndexNSun6(infile, outfile):
    ''' Recompress NSun6 clusteres with bgzip and index with tabix'''


    statement = ''' zcat %(infile)s
                  | sort -k1,1 -k2,2n
                  | bgzip > %(outfile)s;

                  checkpoint;
 
                  tabix -p bed %(outfile)s '''

    P.run()


###################################################################
@transform(getclusterRegions,
           suffix(".region.bed.gz"),
           add_inputs(recompressAndIndexNSun6),
           ".intervalprofile.matrix.tsv.gz")
def getNSUN6IntervalProfile(infiles, outfile):

    regions, clusters = infiles
    track = P.snip(outfile, ".intervalprofile.matrix.tsv.gz")
    statement = ''' cgat bam2geneprofile
                                         --method=intervalprofile
                                         --bedfile=%(clusters)s
                                         --gtf-file=<(cgat bed2gff -I %(regions)s)
                                         --resolution-cds=200
                                         --output-filename-pattern=%(track)s.%%s
                                     &> %(outfile)s.log '''

    P.run()


###################################################################
@merge(getNSUN6IntervalProfile, "nsun6.dir/nsun6_interval_profiles.load")
def loadNSUN6Profiles(infiles, outfile):
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="nsun6.dir/(.+).intervalprofile.matrix.tsv.gz",
                         options = "-i bin -i region -i track")


###################################################################
#@transform(getclusterRegions,
#           suffix(".bed.gz"),
#           ".fasta")
#defGetClusterRegionFasta(infile, outfile):

@follows(loadNSUN6Profiles)
def NSUN6():
    pass


###################################################################
# HNRNPU1
###################################################################
@follows(mkdir("hnrnpu1.dir"))
@transform(os.path.join(PARAMS["external_hnrnpu1"],
                        "*.bed.gz"),
           regex(".+/GSM[0-9]+_(.+)"),
           r"hnrnpu1.dir/\1")
def liftOverhnRNPUReads(infile, outfile):
    ''' Sort and dedup reads, and then liftOver to hg19 '''

    tmp = P.getTempFilename()

    statement = '''liftOver %(infile)s
                            /ifs/mirror/ucsc/hg18/liftOver/hg18ToHg19.over.chain.gz
                            %(tmp)s
                            /dev/null;
                   checkpoint;

                  sed -E  's/rnpu2?[:_][0-9]+/\\./' %(tmp)s
                 | sort -k1,1 -k2,2n -u
                 | bgzip > %(outfile)s;
 
                   checkpoint;

                   tabix -p bed %(outfile)s;
                   
                   checkpoint;

                   rm %(tmp)s'''

    P.run()


###################################################################
@transform(liftOverhnRNPUReads,
           suffix(".bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           ".geneprofile.matrix.tsv.gz")
def hnRNPUMetagene(infiles, outfile):
    '''Plot hnRNPU1 metagene profile'''

    bedfile, gtffile = infiles
    outfile = P.snip(outfile, ".geneprofile.matrix.tsv.gz")
    statement = '''cgat bam2geneprofile
                          --method=geneprofile
                          --bedfile=%(bedfile)s
                          --gtf-file=%(gtffile)s
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --normalize-profile=area
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.log
                          '''
    P.run()


###################################################################
@merge(hnRNPUMetagene,
       "hnrnpu1.dir/hnrnpu1_geneprofiles.load")
def loadhnRNPUMetagene(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+(rep[0-9])\.geneprofile",
                         options="-i bin -i region -i region_bin")


###################################################################
@transform(liftOverhnRNPUReads,
           suffix(".bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           "wholetranscript.geneprofile.matrix.tsv.gz")
def hnRNPUWholeTranscriptMeta(infiles, outfile):
    '''Plot hnRNPU1 metagene profile'''

    bedfile, gtffile = infiles
    outfile = P.snip(outfile, ".geneprofile.matrix.tsv.gz")
    statement = '''cgat bam2geneprofile
                          --method=geneprofile
                          --bedfile=%(bedfile)s
                          --gtf-file=<( zcat %(gtffile)s | grep exon
                                       | cgat gtf2gtf
                                        --method=join-exons 
                                        -L /dev/null
                                       | sed 's/\\ttranscript\\t/\\texon\\t/')
                          --normalize-transcript=total-sum
                          --use-base-accuracy
                          --normalize-profile=area
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.log
                          '''
    P.run()


###################################################################
@merge(hnRNPUWholeTranscriptMeta,
       "hnrnpu1.dir/hnrnpu1_wholetranscript_geneprofiles.load")
def loadRNPUWholeTranscriptMeta(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+(rep[0-9])whole",
                         options="-i bin -i region -i region_bin")


###################################################################
@transform(liftOverhnRNPUReads,
           suffix(".bed.gz"),
           add_inputs(os.path.join(PARAMS["iclip_dir"],
                                   "reference.context.bed.gz")),
           ".context_stats.tsv")
def gethnRNPUContext(infiles, outfile):

    reads, context = infiles
    statement = '''cgat bam_vs_bed
                          -a %(reads)s
                          -b %(context)s
                           -S %(outfile)s
                           -L %(outfile)s.log '''
    P.run()


###################################################################
@merge(gethnRNPUContext,
       "hnrnpu1_context_stats.load")
def loadhnRNPUConext(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+(rep[0-9])\.",
                         options="-i category")

    
###################################################################
@follows(mkdir("fractions.dir"))
@transform("*RiboZ*bam",
           suffix(".bam"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"fractions.dir/\1.feature_counts.tsv.gz")
def countNuclearCytoRNAseq(infiles, outfile):

    bamfile, annotations = infiles
    PipelineRnaseq.runFeatureCounts(
        annotations,
        bamfile,
        outfile,
        job_threads=PARAMS['featurecounts_threads'],
        strand=0,
        options=PARAMS['featurecounts_options'] + ' -O -M')

    
###################################################################
@merge(countNuclearCytoRNAseq,
       "fractions.dir/nuclear_cyto_counts.tsv.gz")
def mergeNuclearCounts(infiles, outfile):
    '''Merge feature counts data into one table'''

    infiles = " ".join(infiles)

    statement=''' cgat combine_tables
                         -c 1
                         -k 7
                         --regex-filename='(.+).feature_counts.tsv.gz'
                         --use-file-prefix
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


@transform(mergeNuclearCounts,
           suffix(".tsv.gz"),
           ".load")
def loadNuclearCytoCounts(infile, outfile):

    P.load(infile, outfile, "-i geneid")


###################################################################
@transform(mergeNuclearCounts,
           regex(".+"),
           "fractions.dir/fraction_diff.deseq.tsv")
def getNuclearLocalisation(infile, outfile):
    '''Use DESeq to compute a relative measure of nucelar localistation
    as the ratio of nuclear to cytoplasmic normalised counts from
    the stubbs data '''

    PipelineProj028.findNuclearLocalisation(infile, outfile)


###################################################################
@transform(getNuclearLocalisation,
           suffix(".tsv"),
           ".load")
def loadNuclearLocalisation(infile, outfile):

    P.load(infile, outfile, "-i gene_id")


###################################################################
@follows(generateSalmonIndex)
@collate(os.path.join(PARAMS["dir_external"],
                      "Fu_RNAseq/*fastq.?.gz"),
         regex(".+/(.+).fastq.(?:[12]\.)*gz"),
         r"hnrnpu1.dir/\1_salmon.dir/quant.sf")
def runSalmonFu(infiles, outfile):
    ''' Run salmon on any provided rnaseq files '''

    track = re.match("(.+)_salmon.dir/quant.sf", outfile).groups()[0]
    if len(infiles) == 1:
        inputs = "-r %s" % infiles[0]
    elif len(infiles) == 2:
        inputs = "-1 <(zcat %s ) -2 <( zcat %s)" % infiles
    else:
        raise ValueError("Don't know how to handle %i input files"
                         % len(infiles))

    statement = ''' salmon quant -i salmon_index.dir
                                   -l '%(sailfish_libtype)s'
                                   %(inputs)s
                                   -o %(track)s_salmon.dir '''

    P.run()


###################################################################
@merge(runSalmonFu,
       r"hnrnpu1.dir/fuRNAseq.tsv.gz")
def mergeSalmonFu(infiles, outfile):

    headers = ",".join([re.search("Hela_(.+)_R1", inf).groups()[0]
                        for inf in infiles])
    #infiles = " ".join(["<( tail -n +10 %s | sed 's/# //' )" % infile
    #                    for infile in infiles])
    infiles = " ".join(infiles)
    statement = '''cgat combine_tables
                    -c 1
                    -k 4,5
                    --prefixes=%(headers)s
                    --regex-filename="hnrnpu1_Hela_(.+)_R1"
                    -L %(outfile)s.log
                    -S %(outfile)s
                    %(infiles)s
           '''

    P.run()


###################################################################
@transform(mergeSalmonFu, suffix(".tsv.gz"), ".load")
def loadSalmonFu(infile, outfile):

    P.load(infile, outfile, "-i gene_id")


###################################################################
@transform(loadSalmonFu,
           suffix(".load"),
           "_diff.deseq.tsv")
def findhnRNPUDependentGenes(infile, outfile):

    PipelineProj028.findhnRNPUDependentGenes(connect(), outfile)


###################################################################
@transform(findhnRNPUDependentGenes,
           suffix(".tsv"),
           ".load")
def loadhnRNPDependentGenes(infile, outfile):

    P.load(infile, outfile)


###################################################################
@follows(loadhnRNPUMetagene,
         loadRNPUWholeTranscriptMeta,
         loadhnRNPUConext,
         loadStubbsCounts,
         loadNuclearLocalisation,
         loadhnRNPDependentGenes)
def hnRNPU1():
    pass


###################################################################
# ChTOP PolyA binding
###################################################################
@follows(mkdir("unmapped_polyA.dir"))
@transform(os.path.join(PARAMS["iclip_dir"],
                        "mapping.dir/star.dir/merged*bam"),
           regex(".+/merged_(.+)\..+\.bam"),
           r"unmapped_polyA.dir/\1.composition.tsv.gz")
def getUnmappedNucleotideComp(infile, outfile):
    '''If ChTop binds to the polyA tail, then it might pull back reads
    that are polyA. Thus look at the nucleotide composition of unmapped
    reads. Other factors for comparision '''

    PipelineProj028.getUnmappedNucleotideComp(infile, outfile,
                                              submit=True)


###################################################################
@merge(getUnmappedNucleotideComp,
       "unmapped_polyA.dir/unmapped_composition.load")
def loadUnmappedComposition(infile, outfile):

    P.concatenateAndLoad(infile, outfile,
                         regex_filename=".+/(.+).composition.tsv.gz",
                         options = "-i base -i track")


###################################################################
@follows(loadUnmappedComposition)
def chtop_polya():
    pass


###################################################################
# Alyref and 3' end processing
###################################################################
@follows(mkdir("stubbs_profiles.dir"))
@transform(os.path.join(PARAMS["dir_external"], "stubbsRNAseq/*.bam"),
           regex(".+/(.+)\.merged.+"),
           add_inputs(os.path.join(
               PARAMS["iclip_dir"],
               "mapping.dir/geneset.dir/reference.gtf.gz")),
           r"stubbs_profiles.dir/\1.geneprofile.matrix.tsv.gz")
def stubbsProfiles(infiles, outfile):
    ''' create metagene profiles over transcripts for control and
    Alyref RNAi to look for effects on run through '''

    bam, genes = infiles
    outprefix = P.snip(outfile, ".geneprofile.matrix.tsv.gz")
    statement = '''cgat bam2geneprofile
                         --method=geneprofile
                         --bam-file=%(bam)s
                         --gtf-file=%(genes)s
                         --reporter=transcript
                         --use-base-accuracy
                         --normalize-transcript=total-sum
                         --normalize-profile=area
                         --scale-flank-length=1
                         --output-filename-pattern=%(outprefix)s.%%s'''

    P.run()


###################################################################
@merge(stubbsProfiles, "stubbs_profiles.load")
def loadStubbsProfiles(infiles, outfile):
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_(.+)_(R[0-9]+).geneprofile.matrix.tsv.gz",
                         cat="condition,fraction,replicate",
                         options="-i fraction,condition,replicate")

###################################################################
# DaPars
###################################################################
@follows(mkdir("dapars.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+"),
           "dapars.dir/geneset.bed")
def getGenesetBed12(infile, outfile):
    '''Convert geneset to BED12 format'''

    statement = '''cgat gff2bed 
                           --bed12-from-transcripts
                            -I %(infile)s
                            -S %(outfile)s
                            -L %(outfile)s.log '''

    P.run()


###################################################################
@follows(mkdir("dapars.dir"))
@originate("dapars.dir/transcripts_to_genes.txt")
def generateDaParsTranscriptsToGenes(outfile):

    statement = '''SELECT DISTINCT transcript_id, gene_id
                   FROM biotypes'''

    data = DUtils.fetch_DataFrame(statement, connect())

    data.to_csv(outfile, header = False, index = False, sep = "\t")


###################################################################
@merge([getGenesetBed12, generateDaParsTranscriptsToGenes],
       "dapars.dir/geneset_extracted.bed")
def getDaParsGeneset(infiles, outfile):
    ''' Process geneset to generate the input file for DaPars '''

    geneset, symbols = infiles

    statement=''' DaPars_Extract_Anno.py -b %(geneset)s
                                         -s %(symbols)s
                                         -o %(outfile)s '''

    P.run()

###################################################################
# Get bedGraph files
@follows(mkdir("dapars.dir"))
@transform(os.path.join(PARAMS["dir_rnaseq"],
                        "*-si*.bam"),
           regex(".+/(GB2-si.+-R.)\..+\.bam"),
           r"dapars.dir/\1.bedGraph")
def ChtopRNAi2BedGraph(infile, outfile):

    PipelineProj028.bamToBedGraph(infile, outfile)


###################################################################
@follows(mkdir("dapars.dir"))
@transform(os.path.join(PARAMS["dir_external"],
                        "stubbsRNAseq/*.bam"),
           regex(".+/(.+)_(.+)_(R[0-9]).+"),
           r"dapars.dir/HEK293\2-\1-\3.bedGraph")
def StubbsToBedGraph(infile, outfile):

    PipelineProj028.bamToBedGraph(infile, outfile)


###################################################################
@collate([ChtopRNAi2BedGraph, StubbsToBedGraph],
         regex(".+/(.+)-(.+)-(R[0-9]).bedGraph"),
         add_inputs(getDaParsGeneset),
         r"dapars.dir/\1.dapars_config.txt")
def generateDaParsConfig(infiles, outfile):
    '''Generate config file for DaPars. Files all have same first name
    part. Conditions are divided on second part, and reps on third'''

    conditions = collections.defaultdict(list)
    track = P.snip(outfile, ".dapars_config.txt")
    utrs = infiles[0][1]
    for infile in infiles:
        conditions[re.match("(.+)-R[0-9]+", infile[0]).groups()[0]].append(infile[0])

    condition1_files = conditions[conditions.keys()[0]]
    condition2_files = conditions[conditions.keys()[1]]

    PipelineProj028.generateDaParsConfig(condition1_files,
                                         condition2_files,
                                         utrs,
                                         os.path.join(track,
                                                      "dapars_out.tsv"),
                                         outfile)


###################################################################
@transform(generateDaParsConfig,
           regex(".+/(.+).dapars_config.txt"),
           r"dapars.dir/\1/dapars_out.tsv")
def runDaPars(infile, outfile):

    statement = '''DaPars_main.py %(infile)s > %(infile)s.log'''
    P.run()


###################################################################
###################################################################
# Processing Index
###################################################################
@follows(mkdir("processing_index.dir"))
@originate("processing_index.dir/HEK293_cleavage_sites.bed.gz")
def download_cleavage_sites(outfile):
    '''Download A-seq cleavage site data from GSM909242'''

    statement = '''wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM909nnn/GSM909242/suppl/GSM909242%%5F178%%2Ebed%%2Egz;
                mv GSM909242_178.bed.gz %(outfile)s'''

    P.run()


###################################################################
@transform(download_cleavage_sites,
           suffix(".bed.gz"),
           add_inputs(filterExpressedTranscripts),
           ".3prime.bed.gz")
def get_3prime_cleavage_sites(infiles, outfile):
    '''For each expressed gene find the 3'-most cleavage site'''

    cleavage_sites, geneset = infiles
    PipelineProj028.find_final_cleavage_sites(geneset, cleavage_sites,
                                              outfile)
 
   
###################################################################
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
           regex(".+/(.+).bam"),
           add_inputs(get_3prime_cleavage_sites),
           r"processing_index.dir/\1.tsv")
def get_processing_index(infiles, outfile):
    '''Calculate processing_index for each bam file using
    3' most cleavage sites'''


    bamfile, sites = infiles

    statement = '''python %(project_src)s/iCLIPlib/scripts/processing_index.py
                             -I %(sites)s
                             -L %(outfile)s.log
                             -S %(outfile)s
                             %(bamfile)s'''

    P.run()


###################################################################
@merge(get_processing_index, "processing_index.dir/processing_index.load")
def load_processing_index(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)-FLAG.(.+).tsv",
                         cat="factor,replicate",
                         has_titles=False,
                         header="factor,replicate,stat,filename,processing_index",
                         options="-i factor -i replicate")


###################################################################
@follows(load_processing_index)
def processing_index():
    pass

###################################################################
# Chimeric Reads
###################################################################
@follows(mkdir("chimeric_reads.dir"))
@transform([os.path.join(PARAMS["iclip_dir"],
                        "mapping.dir/star.dir/merged*.bam"),
            os.path.join(PARAMS["ejc_iclip_dir"],
                         "mapping.dir/star.dir/merged*.bam")],
           regex(".+/merged_(.+).star.bam"),
           r"chimeric_reads.dir/\1.unmapped.fastq.gz")
def get_unmapped_fastq(infile, outfile):
    '''Filter out the unmapped reads from the iCLIP pipeline and output
    them to pairs of fastq files ready for mapping'''

    # read unmapped is flag 4
   
    statement = ''' samtools view -hf4 %(infile)s 
                  | samtools fastq -
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(get_unmapped_fastq,
           formatter(),
           "chimeric_reads.dir/{basename[0]}.bam")
def map_unmapped_reads(infile, outfile):
    ''' Take the reads that STAR failed to map and map them with BWA'''

    from CGATPipelines import PipelineMapping
    bwa_threads = 4
    bwa_index_dir="/ifs/mirror/genomes/bwa-0.7.5a"
    bwa_mem_options = "-M"

    m = PipelineMapping.BWAMEM(remove_non_unique="1",
                               strip_sequence=None,
                               set_nh=None)
    statement = m.build((infile,),outfile)
    P.run()


###################################################################
@collate([map_unmapped_reads,
          os.path.join(PARAMS["iclip_dir"],
                       "mapping.dir/star.dir/merged_*-FLAG-R*bam"),
         os.path.join(PARAMS["ejc_iclip_dir"],
                       "mapping.dir/star.dir/merged_*-GFP-R*bam")],
         regex(".+/(?:merged_)?(.+-FLAG-R.|.+-GFP-R.).+bam"),
         r"chimeric_reads.dir/\1.merged.bam")
def merge_mapped_and_remapped(infiles, outfile):

    infiles = " ".join(infiles)

    statement = '''samtools merge %(outfile)s %(infiles)s;
                   checkpoint;
                   samtools index %(outfile)s'''

    P.run()

    
###################################################################
@transform(merge_mapped_and_remapped,
           regex("(.+).merged.bam"),
           r"\1.deduped.bam")
def dedup_remapped_reads(infile, outfile):
    '''Do UMI based deduplication of the reads that BWA mapped but
    STAR doesn't'''

    tmpfile = P.getTempFilename()

    statement = '''umi_tools dedup -I %(infile)s -S %(tmpfile)s.bam -L %(outfile)s.log;
                   checkpoint;
                   samtools view -bq 20 %(tmpfile)s.bam
                 | samtools sort -o %(outfile)s;
                   checkpoint;
                   rm %(tmpfile)s.bam;
                   checkpoint;
                   samtools index %(outfile)s'''

    P.run()


###################################################################
@transform(dedup_remapped_reads,
           suffix(".deduped.bam"),
           add_inputs(os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"] + ".fasta"),
                      filterExpressedTranscripts,
                      PARAMS["external_junctions_db"]),
           ".chimeric.bam")
def get_chimeric_reads(infiles, outfile):
         
    bam, fasta, gtf, junctions = infiles
    tmp = P.getTempFilename()
    out_prefix = P.snip(outfile, ".bam")
    log = out_prefix + ".log"
    statement = '''python %(project_src)s/find_chimeric_reads.py
                      -g %(gtf)s
                      -f %(fasta)s
                      -j %(junctions)s
                      -m 0.2
                      --minimum-unaligned=5
                      --out-prefix=%(out_prefix)s
                      -L %(log)s
                      %(bam)s;
                checkpoint;
                samtools sort %(outfile)s -o %(tmp)s.bam;
                checkpoint;
                mv %(tmp)s.bam %(outfile)s;
                checkpoint;
                samtools index %(outfile)s'''

    P.run()


###################################################################
@collate(get_chimeric_reads,
         regex("(.+-FLAG|.+-GFP).+"),
         r"\1-union.chimeric.bam")
def merge_chimeric_reads(infiles, outfile):

    infiles = " ".join(infiles)
    
    statement = '''samtools merge - %(infiles)s
                 | samtools sort - -o %(outfile)s;

                 checkpoint;
 
                 samtools index %(outfile)s '''

    P.run()

    
###################################################################
@collate(get_chimeric_reads,
         regex("(.+).bam"),
         inputs(r"\1.log"),
         "chimeric_log.load")
def load_chimeric_read_log(infiles, outfile):
    '''Process and load the general log file for the chimeric reads'''

    loglines = []

    for infile in infiles:
        track = P.snip(os.path.basename(infile), ".chimeric.log")
        F=track.split("-")
        protein = "-".join(F[:-1])
        rep = F[-1]
        
        for line in IOTools.openFile(infile):
            if line.startswith('#'):
                continue
            F = line.strip().split("\t")
            loglines.append([protein, rep, F[0].split("INFO")[-1], F[1]])

    tmpfn = P.getTempFilename(shared=True)
    IOTools.writeLines(tmpfn, loglines, header=["protein", "replicate", "category",  "count"])

    P.load(tmpfn, outfile)

    os.unlink(tmpfn)

    
@follows(get_chimeric_reads)
@collate("chimeric_reads.dir/*-*-R*.chimeric.*.tsv.gz",
         regex("(.+).chimeric.(.+).tsv.gz"),
         r"chimeric_reads.dir/chimeric_\2.load")
def load_chimeric_stats(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-FLAG|.+-GFP)-(R.).chimeric",
                         cat="protein,rep",
                         options = "-i gene_id -i Distance -i protein -i rep")


@follows(load_chimeric_read_log, load_chimeric_stats, merge_chimeric_reads)
def chimeric_reads():
    pass

###################################################################
# PAPERCLIP AND ChTop
###################################################################
@follows(mkdir("paperclip_chtop.dir"))
@originate("paperclip_chtop.dir/HEK293_paperclip.bed.gz")
def download_paperclip(outfile):
    '''download the paperclip data from GEO as a bed file'''

    address = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66092/suppl/GSE66092%5FHEK293%5FBC2%5F16414%5FGeneAssigned%2Ebed%2Egz"
    statement = '''wget 
                   %(address)s
                   -O %(outfile)s'''
    P.run()

    
###################################################################
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           formatter(),
           add_inputs(download_paperclip),
           "paperclip_chtop.dir/single_polyA_windows.gtf.gz")
def get_single_polyA_windows(infiles, outfile):
    '''Extracts 100bp around polyA sites where a transcript only contains
    a single site'''

    gtffile, bedfile = infiles
    PipelineProj028.get_single_polyA_windows(gtffile, bedfile, outfile)

    
###################################################################
@transform([os.path.join(PARAMS["iclip_dir"],
                         "deduped.dir/*union.bam"),
            merge_cpsf_wigs],
           regex(".+/(.+).union.(bam|bw)"),
           add_inputs(get_single_polyA_windows),
           r"paperclip_chtop.dir/\1.single_polyA_metagene.tsv",
           r"\2")
def get_metagene_over_single_polyA_windows(infiles, outfile, informat):
    '''Generate metagene plots around polyA signals where a transcript
    has just a single polyA site in the paperclip'''

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles

    if informat == "bam":
        instat = bamfile
    elif informat == "bw":
        instat = "--plus-wig=%s" % bamfile
        
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                           -I %(gtffile)s
                           %(instat)s
                           --flanks=0
                           -S %(outfile)s
                           -b 60
                           -L %(outfile)s.log'''
    P.run()
    

@merge(get_metagene_over_single_polyA_windows,
       "paperclip_chtop.dir/single_polyA_window_metagenes.load")
def load_single_polyA_metagenes(infiles, outfile):

    P.concatenateAndLoad(infiles,
                          outfile,
                          regex_filename=".+/(.+).single_polyA_metagene.tsv")


@subdivide(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           formatter(),
           add_inputs(download_paperclip),
           ["paperclip_chtop.dir/apa_5prime.gtf.gz",
            "paperclip_chtop.dir/apa_3prime.gtf.gz"])
def get_apa_windows(infiles, outfiles):

    gtffile, bedfile = infiles
    out_5prime, out_3prime = outfiles

    PipelineProj028.get_alternate_polyAs(gtffile, bedfile,
                                         out_5prime, out_3prime)


@product([os.path.join(PARAMS["iclip_dir"],
                         "deduped.dir/*union.bam"),
            merge_cpsf_wigs],
         formatter(".+/(?P<TRACK>.+).union.(?P<ftype>.+)"),
         get_apa_windows,
         formatter(".+/apa_(?P<PA>.+).gtf.gz"),
         "paperclip_chtop.dir/{TRACK[0][0]}.apa_{PA[1][0]}_normed.tsv",
         "{ftype[0][0]}")
def get_normed_metagene_over_apa_windows(infiles, outfile, informat):

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles

    if informat == "bam":
        instat = bamfile
    elif informat == "bw":
        instat = "--plus-wig=%s" % bamfile
        
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                           -I %(gtffile)s
                           %(instat)s
                           --flanks=0
                           -S %(outfile)s
                           -b 30
                           -L %(outfile)s.log'''
    P.run()

    
@product([os.path.join(PARAMS["iclip_dir"],
                         "deduped.dir/*union.bam"),
            merge_cpsf_wigs],
         formatter(".+/(?P<TRACK>.+).union.(?P<ftype>.+)"),
         get_apa_windows,
         formatter(".+/apa_(?P<PA>.+).gtf.gz"),
         "paperclip_chtop.dir/{TRACK[0][0]}.apa_{PA[1][0]}_unnormed.tsv",
         "{ftype[0][0]}")
def get_unnormed_metagene_over_apa_windows(infiles, outfile, informat):

    job_options = "-l mem_free=15G"
    bamfile, gtffile = infiles

    if informat == "bam":
        instat = bamfile
    elif informat == "bw":
        instat = "--plus-wig=%s" % bamfile
        
    statement = '''python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                           -I %(gtffile)s
                           %(instat)s
                           --flanks=0
                           -S %(outfile)s
                           -b 30
                           --no-gene-norm
                           -L %(outfile)s.log'''
    P.run()

    
@merge([get_normed_metagene_over_apa_windows,
        get_unnormed_metagene_over_apa_windows],
         "apa_metagenes.load")
def load_apa_metagenes(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).apa_(.+)_(normed|unnormed).tsv",
                         cat="track,pA,normed",
                         options="-i track -i pA -i normed")



@follows(load_apa_metagenes)
def paperclip():
    pass

@transform(merge_cpsf_wigs,
           suffix(".bw"),
           ".clusters.bedGraph.gz")
def get_cpsf_clusters(infile, outfile):
    '''Use randomisation based cluster calling to find cpsf clusters'''

    tmpfile = P.getTempFilename()
    statement = '''bigWigToBedGraph %(infile)s %(tmpfile)s;

                   checkpoint;
 
                   bedtools merge -i %(tmpfile)s -d 10 -c 4 -o sum
                 | bgzip > %(outfile)s;


                 checkpoint;

                 tabix -p bed %(outfile)s;

                 checkpoint;

                 rm %(tmpfile)s '''

    P.run()

@follows(mkdir("paperclip_chtop.dir"))
@transform( merge_adjacent_clusters,
            regex(".+/(.+).merged_clusters.bed.gz"),
            add_inputs(merge_cpsf_wigs,
                       os.path.join(
                           PARAMS["annotations_dir"],
                           PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
            r"paperclip_chtop.dir/\1.cpsf_around_clusters.tsv")
def get_cpsf_metagenes_around_clusters(infiles, outfile):
    '''Get CPSF metagenes around clusters. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf = infiles
    
    PipelineProj028.get_profile_around_clusters(
        clusters,
        to_profile,
        gtf,
        outfile,
        submit=True)

@follows(mkdir("paperclip_chtop.dir"))
@transform( merge_adjacent_clusters,
            regex(".+/(.+).merged_clusters.bed.gz"),
            add_inputs(merge_cpsf_wigs,
                       os.path.join(
                           PARAMS["annotations_dir"],
                           PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
             PARAMS["external_chtop_apa_db"]),
            r"paperclip_chtop.dir/\1.cpsf_around_apa_clusters.tsv")
def get_cpsf_metagenes_around_apa_clusters(infiles, outfile):
    '''Get CPSF metagenes around clusters. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf, apa_db = infiles
    
    apa_transcripts = DUtils.fetch_DataFrame(
        '''SELECT transcript_id from dapars
           WHERE PDUI_Group_diff < -0.25 AND
                 adjusted_P_Val < 0.05''',
        connect())

    use_transcripts = list(apa_transcripts.transcript_id)
    PipelineProj028.get_profile_around_clusters(
        clusters,
        to_profile,
        gtf,
        outfile,
        use_transcripts=use_transcripts,
        submit=True)

    
@follows(mkdir("paperclip_chtop.dir"))
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*union.bam"),
           regex(".+/(.+.union).bam"),
           add_inputs(merge_cpsf_wigs,
                      os.path.join(
                          PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"paperclip_chtop.dir/\1.cpsf_around_sites.tsv")
def get_cpsf_metagenes_around_clipsites(infiles, outfile):
    '''Get CPSF metagenes around sites. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf = infiles
    
    PipelineProj028.get_profile_around_clipsites(
        clusters,
        to_profile,
        gtf,
        outfile,
        submit=True)

    

@follows(mkdir("paperclip_chtop.dir"))
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*union.bam"),
           regex(".+/(.+.union).bam"),
           add_inputs(merge_cpsf_wigs,
                      os.path.join(
                          PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
                      PARAMS["external_chtop_apa_db"]),
           r"paperclip_chtop.dir/\1.cpsf_around_apa_sites.tsv")
def get_cpsf_metagenes_around_clipsites_in_apa(infiles, outfile):
    '''Get CPSF metagenes around sites. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf, apa_db = infiles

    apa_transcripts = DUtils.fetch_DataFrame(
        '''SELECT transcript_id from dapars
           WHERE PDUI_Group_diff < -0.25 AND
                 adjusted_P_Val < 0.05''',
        connect())

    use_transcripts = list(apa_transcripts.transcript_id)
    
    PipelineProj028.get_profile_around_clipsites(
        clusters,
        to_profile,
        gtf,
        outfile,
        use_transcripts=use_transcripts,
        submit=True)




@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*union.bam"),
           regex(".+/(.+.union).bam"),
           add_inputs(merge_cpsf_wigs,
                      os.path.join(
                          PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"paperclip_chtop.dir/\1.cpsf_around_sites_normed.tsv")
def get_normed_cpsf_metagenes_around_clipsites(infiles, outfile):
    '''Get CPSF metagenes around sites. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf = infiles
    
    PipelineProj028.get_profile_around_clipsites(
        clusters,
        to_profile,
        gtf,
        outfile,
        norm=True,
        submit=True)


@follows(mkdir("paperclip_chtop.dir"))
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*union.bam"),
           regex(".+/(.+.union).bam"),
           add_inputs(merge_cpsf_wigs,
                      os.path.join(
                          PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
                      PARAMS["external_chtop_apa_db"]),
           r"paperclip_chtop.dir/\1.cpsf_around_apa_sites_normed.tsv")
def get_normed_cpsf_metagenes_around_clipsites_in_apa(infiles, outfile):
    '''Get CPSF metagenes around sites. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf, apa_db = infiles

    apa_transcripts = DUtils.fetch_DataFrame(
        '''SELECT transcript_id from dapars
           WHERE PDUI_Group_diff < -0.25 AND
                 adjusted_P_Val < 0.05''',
        connect())

    use_transcripts = list(apa_transcripts.transcript_id)
    
    PipelineProj028.get_profile_around_clipsites(
        clusters,
        to_profile,
        gtf,
        outfile,
        use_transcripts=use_transcripts,
        norm=True,
        submit=True)



@follows(mkdir("paperclip_chtop.dir"))
@transform( merge_adjacent_clusters,
            regex(".+/(.+).merged_clusters.bed.gz"),
            add_inputs(merge_cpsf_wigs,
                       os.path.join(
                           PARAMS["annotations_dir"],
                           PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
            r"paperclip_chtop.dir/\1.cpsf_around_clusters_normed.tsv")
def get_cpsf_metagenes_around_clusters_normed(infiles, outfile):
    '''Get CPSF metagenes around clusters. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf = infiles
    
    PipelineProj028.get_profile_around_clusters(
        clusters,
        to_profile,
        gtf,
        outfile,
        norm=True,
        submit=True)

    
@follows(mkdir("paperclip_chtop.dir"))
@transform( merge_adjacent_clusters,
            regex(".+/(.+).merged_clusters.bed.gz"),
            add_inputs(merge_cpsf_wigs,
                       os.path.join(
                           PARAMS["annotations_dir"],
                           PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
             PARAMS["external_chtop_apa_db"]),
            r"paperclip_chtop.dir/\1.cpsf_around_apa_clusters_normed.tsv")
def get_cpsf_metagenes_around_apa_clusters_normed(infiles, outfile):
    '''Get CPSF metagenes around clusters. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf, apa_db = infiles
    
    apa_transcripts = DUtils.fetch_DataFrame(
        '''SELECT transcript_id from dapars
           WHERE PDUI_Group_diff < -0.25 AND
                 adjusted_P_Val < 0.05''',
        connect())

    use_transcripts = list(apa_transcripts.transcript_id)
    PipelineProj028.get_profile_around_clusters(
        clusters,
        to_profile,
        gtf,
        outfile,
        use_transcripts=use_transcripts,
        norm=True,
        submit=True)


@follows(mkdir("paperclip_chtop.dir"))
@transform( merge_adjacent_clusters,
            regex(".+/(.+union).merged_clusters.bed.gz"),
            add_inputs(merge_cpsf_wigs,
                       os.path.join(
                           PARAMS["annotations_dir"],
                           PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
            r"paperclip_chtop.dir/\1.cpsf_around_clusters_normed_rand.tsv")
def get_cpsf_rand_metagenes_around_clusters_normed(infiles, outfile):
    '''Get CPSF metagenes around clusters. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf = infiles
    
    PipelineProj028.get_profile_around_clusters(
        clusters,
        to_profile,
        gtf,
        outfile,
        norm=True,
        rands=100,
        submit=True)

    
@transform(os.path.join(PARAMS["iclip_dir"], "deduped.dir/*union.bam"),
           regex(".+/(.+.union).bam"),
           add_inputs(merge_cpsf_wigs,
                      os.path.join(
                          PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"paperclip_chtop.dir/\1.cpsf_around_sites_normed_rand.tsv")
def get_normed_cpsf_metagenes_around_clipsites_rand(infiles, outfile):
    '''Get CPSF metagenes around sites. Only works with 
    last exons of transcripts, and doesn't bother to look outside
    those'''

    clusters, to_profile, gtf = infiles
    
    PipelineProj028.get_profile_around_clipsites(
        clusters,
        to_profile,
        gtf,
        outfile, 
        norm=True,
        rands=20,
        submit=True)


    
@collate([get_cpsf_metagenes_around_clusters,
          get_cpsf_metagenes_around_apa_clusters,
          get_cpsf_metagenes_around_clipsites,
          get_normed_cpsf_metagenes_around_clipsites_rand,
          get_cpsf_metagenes_around_clipsites_in_apa,
          get_normed_cpsf_metagenes_around_clipsites,
          get_normed_cpsf_metagenes_around_clipsites_in_apa,
          get_cpsf_metagenes_around_clusters_normed,
          get_cpsf_metagenes_around_apa_clusters_normed,
          get_cpsf_rand_metagenes_around_clusters_normed],
        regex(".+/(.+-FLAG).(R.|union).+cpsf_around_(.+).tsv"),
        r"paperclip_chtop.dir/cpsf_metagenes_around_\3.load")
def load_cpsf_metagenes_around_clusters(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+-FLAG).(R.|union).cpsf_around.+",
                         cat="track,replicate")

    
###################################################################
@follows(mkdir("chtop_apa.dir"))
@originate(["chtop_apa.dir/HEK293-Aseq_Cntrl-R1.bed.gz",
            "chtop_apa.dir/HEK293-Aseq_siCPSF6-R1.bed.gz",
            "chtop_apa.dir/HEK293-Aseq_nosi-R1.bed.gz",
            "chtop_apa.dir/HEK293-Aseq_CSTF64-R1.bed.gz"])
def get_Aseq_data(outfile):
    '''Download Aseq data from GSE37401'''

    addresses = {"HEK293-Aseq_Cntrl-R1.bed.gz":
                 "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM909nnn/GSM909243/suppl/GSM909243%5F179%2Ebed%2Egz",
                 "HEK293-Aseq_siCPSF6-R1.bed.gz":
                 "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM909nnn/GSM909245/suppl/GSM909245%5F183%2Ebed%2Egz",
                 "HEK293-Aseq_nosi-R1.bed.gz":
                 "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM909nnn/GSM909242/suppl/GSM909242%5F178%2Ebed%2Egz",
                 "HEK293-Aseq_CSTF64-R1.bed.gz":
                 "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM909nnn/GSM909244/suppl/GSM909244%5F181%2Ebed%2Egz"}

    infile = addresses[os.path.basename(outfile)]

    statement = '''wget %(infile)s -O %(outfile)s'''

    P.run()


@transform(get_Aseq_data,
           suffix(".bed.gz"),
           ".tpm.bed.gz")
def convert_aseq_to_tpm(infile, outfile):

    PipelineProj028.bed_score_to_TPM(infile, outfile)


@merge(get_Aseq_data,
       "aseq_clusters.bed.gz")
def merge_aseq_to_clusters(infiles, outfile):

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s
                 | sort -k1,1 -k2,2n
                 | bedtools merge -i - -c 4,5,6 -o first,sum,first -d 10
                 | gzip > %(outfile)s'''
    P.run()

@merge([convert_aseq_to_tpm, merge_aseq_to_clusters],
       "chtop_apa.dir/aseq_clusters.filtered.bed.gz")
def quantify_and_filter_clusters(infiles, outfile):

    PipelineProj028.quantify_and_filter_clusters(infiles[-1],
                                                 infiles[:-1],
                                                 outfile)
    
@transform( quantify_and_filter_clusters,
           suffix(".filtered.bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           ".prox_distil.tsv")
def get_Aseq_ratios(infiles, outfile):
    ''' Find the ratio of proximal to distil usage for transcripts'''

    bedfile, gffile = infiles
    
    PipelineProj028.get_alternate_polyA_distil_usage(gffile, bedfile,
                                                      outfile)


@transform(get_Aseq_ratios,
           suffix(".tsv"),
           ".load")
def load_aseq_ratio(infile, outfile):

    P.load(infile, outfile, options = "-i gene_id -i transcript_id")


@transform(getTranscriptChunks,
           suffix(".gtf.gz"),
           add_inputs(PARAMS["external_chtop_apa_db"]),
           ".chtop_apa.tsv")
def annotate_chtop_apa_chunks(infiles, outfile):

    chunks = infiles[0]
    table = DUtils.fetch_DataFrame('''SELECT chrom as chr,
                                             Predicted_Proximal_APA as position
                                      FROM chtop_apa.dapars
                                      WHERE adjusted_P_Val < 0.05
                                       AND PDUI_Group_diff < -0.1''',
                                   connect())
    PipelineProj028.get_apa_chunks(table, chunks, outfile)


@transform(annotate_chtop_apa_chunks,
           suffix(".tsv"),
           ".load")
def load_chtop_apa_chunks(infile, outfile):

    P.load(infile, outfile, options = "-i gene_id -i exon_id")

    connect().executescript('''DROP INDEX IF EXISTS chtop_apa_joint;
                               CREATE INDEX chtop_apa_joint 
                                ON reference_chunks_chtop_apa(gene_id, exon_id)''')


@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           formatter(),
           r"transcripts_without_last_exon.gtf.gz")
def get_transcripts_exon_last_exons(infile, outfile):

    PipelineProj028.geneset_without_last_exons(infile, outfile)
           
    
###################################################################
@follows(mkdir("splicing_ratios.dir"))
@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"splicing_ratios.dir/\1.firt_two_exons.gtf.gz")
def get_first_two_exons(infile, outfile):

    import CGAT.GTF as GTF
    with IOTools.openFile(outfile, "w") as outf:
        for transcript in GTF.transcript_iterator(
                GTF.iterator(IOTools.openFile(infile))):

            transcript = [e for e in transcript if e.feature == "exon"]

            if transcript[0].strand == "+":
                transcript = sorted(transcript, key=lambda x: x.start)
            else:
                transcript = sorted(transcript, key=lambda x: x.end,
                                    reverse=True)

            for e in transcript[:2]:
                outf.write(str(e) + "\n")


###################################################################
@transform([os.path.join(PARAMS["iclip_dir"], "deduped.dir/*.bam"),
            os.path.join(PARAMS["ejc_iclip_dir"],
                         "deduped.dir/*.bam")],
           regex(".+/(.+).bam"),
           add_inputs(get_first_two_exons),
           r"splicing_ratios.dir/\1.splicing_ratio.tsv")
def get_first_exon_splicing_ratio(infiles, outfile):

    bamfile, gtffile = infiles
    PipelineProj028.calculateSplicingIndex(bamfile, gtffile, outfile,
                                           submit=True)


@merge( get_first_exon_splicing_ratio, "first_exon_splicing_index.load")
def load_first_exon_si(infiles, outfile):

    P.concatenateAndLoad(infiles,
                         outfile,
                         regex_filename=".+/(.+)-(.+).(R.|union).splicing_ratio.tsv",
                         cat="factor,tag,replicate")

    
###################################################################
@follows(load_first_exon_si)
def first_exon_si():
    pass
###################################################################
# Export
###################################################################
@follows(mkdir("export/hg19"))
@transform(callReproducibletRNAClusters,
           regex(".+/(.+).bed.gz"),
           r"export/hg19/\1_tRNA.bigBed")
def exporttRNAClusters(infile, outfile):
    
    PipelineiCLIP.clustersToBigBed(infile, outfile)


###################################################################
@merge(exporttRNAClusters,
       "export/hg19/trna_trackDb.txt")
def maketRNAClusterUCSC(infiles, outfile):

    PipelineiCLIP.makeClustersUCSC(infiles, outfile, "tRNACluster",
                                   "Clusters from tRNAs")


###################################################################
@follows(mkdir("export/hg19"))
@transform(os.path.join(PARAMS["dir_external"],
                        "stubbsRNAseq/*.bam"),
           regex(".+/(.+_.+_R[0-9]).+"),
           r"export/hg19/\1.bigWig")
def stubbsRNAseqToBigWig(infile, outfile):

    PipelineProj028.bamToBigWig(infile, outfile)

###################################################################
@follows(mkdir("export/hg19"))
@merge(stubbsRNAseqToBigWig,
       "export/hg19/stubbs_trackDb.txt")
def generateStubbsTrackDb(infiles, outfile):

    PipelineProj028.bigWigTrackDB(infiles,
                                  "RNASeq data from stubbs et al track %(track)s",
                                  "StubbsRNAseq",
                                  "RNAseq data from Stubbs et al",
                                  outfile)


###################################################################
@follows(mkdir("export/hg19"))
@transform(os.path.join(PARAMS["dir_rnaseq"],
                        "*-si*-*.bam"),
           regex(".+/(.+\-si.+-.+)\..+\.bam"),
           r"export/hg19/\1.bigWig")
def chtopRNASeqToBigWig(infile, outfile):

    PipelineProj028.bamToBigWig(infile, outfile)


###################################################################
@merge(chtopRNASeqToBigWig,
       "export/hg19/chtop_trackDb.txt")
def generateChTopTrackDb(infiles, outfile):

    PipelineProj028.bigWigTrackDB(infiles,
                                  "RNASeq from ChTop knockdown track %(track_name)s",
                                  "ChTop_RNAseq",
                                  "RNASeq data from ChTop known experiments",
                                  outfile)

@transform(PARAMS["dir_chtop_rnai_bams"],
           regex(".+/([^\.]+)\..+\.bam"),
           r"export/hg19/\1.bigWig")
def our_chtop_rnaseq_to_bigwig(infile, outfile):

    PipelineProj028.bamToBigWig(infile, outfile)

@merge(our_chtop_rnaseq_to_bigwig,
       "export/hg19/out_chtop_trackDb.txt")
def generate_our_chtop_rnaseq_trackDb(infiles, outfile):

    PipelineProj028.bigWigTrackDB(infiles,
                                  "Our RNASeq from ChTop knockdown track %(track_name)s",
                                  "Our_Chtop",
                                  "Our RNAseq from ChTop knockdown",
                                  outfile)
    


@transform(get_Aseq_data,
           regex(".+/(HEK293-Aseq_(?:Cntrl|siCPSF6)-R.).bed.gz"),
           r"export/hg19/\1.bigWig")
def get_CPSF6_Aseq_bigWig(infile, outfile):

    tmpfile = P.getTempFilename()
    genome_file = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_contigs_tsv"])
    
    statement = '''zcat %(infile)s
                 | sort -k1,1 -k2,2n
                 | bedtools merge -i - -c 5 -o sum -d -1
                 > %(tmpfile)s;

                 checkpoint;

                 bedGraphToBigWig %(tmpfile)s %(genome_file)s %(outfile)s;

                 checkpoint;

                 rm %(tmpfile)s'''

    P.run()


@merge(get_CPSF6_Aseq_bigWig,
       "export/hg19/CPSF6_RNAi_Aseq_trackDb.txt")
def generate_CPSF6_RNAi_Aseq_trackDb(infiles, outfile):

    PipelineProj028.bigWigTrackDB(infiles,
                                  "Aseq from Martin et al, track %(track_name)s",
                                  "CFSP6_Aseq",
                                  "Aseq from Martin et al",
                                  outfile)

###################################################################
@transform(os.path.join(PARAMS["dir_external"], "Fu_RNAseq/*.bed.gz"),
           regex(".+/GSM[0-9]+_Hela_hnRNPU_RNAseq_(.+).bed.gz"),
           r"hnrnpu1.dir/HeLa_RNASeq_\1.bed.gz")
def liftOverFuRNASeq(infile, outfile):

    outfile = P.snip(outfile, ".gz")
    statement = '''liftOver <( zcat %(infile)s | grep -P "^chr" )
                            /ifs/mirror/ucsc/hg18/liftOver/hg18ToHg19.over.chain.gz
                            %(outfile)s
                            %(outfile)s.unmapped;

                   checkpoint;

                   gzip %(outfile)s'''

    P.run()

    

###################################################################
@transform([liftOverFuRNASeq,
            liftOverhnRNPUReads],
           regex(".+/(.+)\.bed.gz"),
           r"export/hg19/\1.bigWig")
def FuRNAToBigWig(infile, outfile):

    PipelineProj028.bamToBigWig(infile, outfile)

###################################################################
@collate(FuRNAToBigWig, 
         regex(".+/HeLa_(.+)_(.+)\.bigWig"),
         r"export/hg19/Fu_\1_trackDb.txt")
def generateFuTrackDb(infiles, outfile):

    if "iCLIP" in infiles[0]:
        type = "iCLIP"
    else:
        type = "RNASeq"

    long_label_template = "%s data from Fu et al track %%(track)s" % type

    PipelineProj028.bigWigTrackDB(infiles,
                                  long_label_template,
                                  "Fu_et_al_%s" % type,
                                  "%s data from Fu et al" % type,
                                  outfile)


###################################################################
@merge([os.path.join(PARAMS["iclip_dir"], "export/hg19/trackDb.txt"),
        os.path.join(PARAMS["ejc_iclip_dir"], "export/hg19/trackDb.txt"),
        maketRNAClusterUCSC,
        generateStubbsTrackDb,
        generateFuTrackDb,
        generateChTopTrackDb,
        generate_our_chtop_rnaseq_trackDb,
        generate_CPSF6_RNAi_Aseq_trackDb],
       "export/hg19/trackDb.txt")
def mergeTrackDbs(infiles, outfile):
    
    infiles = " ".join(infiles)

    statement = "cat %(infiles)s > %(outfile)s"

    P.run()


###################################################################
@follows(mkdir("export/hg19"))
@transform([os.path.join(PARAMS["iclip_dir"], "export/hg19/*"),
            os.path.join(PARAMS["ejc_iclip_dir"], "export/hg19/*")],
           regex(".+/(export/hg19/.+)"),
           r"\1")
def linkiCLIPtracks(infile, outfile):

    if "trackDb" not in infile:
        try:
            os.symlink(os.path.abspath(infile), os.path.abspath(outfile))
        except OSError:
            os.unlink(outfile)
            os.symlink(os.path.abspath(infile), os.path.abspath(outfile))


###################################################################
@follows(mkdir("export"))
@originate(["export/hub.txt",
            "export/genomes.txt"])
def makeHubFiles(outfiles):

    hub_file = '''
    hub proj028
    shortLabel CGAT project 28
    longLabel All browser tracks associated with CGAT project 28
    genomesFile genomes.txt
    email i.sudbery@sheffield.ac.uk'''

    with IOTools.openFile("export/hub.txt", "w") as outf:
        outf.write(hub_file)

    genomes_file = '''
    genome hg19
    trackDb hg19/trackDb.txt'''

    with IOTools.openFile("export/genomes.txt", "w") as outf:
        outf.write(genomes_file)


###################################################################
@follows(mergeTrackDbs,
         linkiCLIPtracks,
         makeHubFiles)
def export():
    pass




###################################################################
# primary targets
###################################################################
@follows(profiles,
         transcript_counts,
         GO,
         clusterOptimisaiton,
         no_flipin_clusters,
         methylation,
         tRNAs,
         retained_introns,
         motifs,
         NSUN6,
         hnRNPU1,
         transcript_chunks,
         chtop_polya,
         processing_index,
         export,
         loadStubbsProfiles)
def full():
    pass

@follows(mkdir("notebooks"),
         mkdir("notebooks/imgs"))
@transform(os.path.join(os.path.dirname(__file__), "notebooks/*.Rmd"),
           formatter(),
           "notebooks/{basename[0]}.html")
def build_R_notebooks(infile, outfile):

    basedir = os.path.abspath(PARAMS["workingdir"])
    statement = ''' Rscript -e 'library(knitr); 
                                knitr::opts_knit$set(root.dir="%(basedir)s/notebooks");
                                knit("%(infile)s", output="%(outfile)s") ' '''
    P.run()

         
@follows( mkdir( "report" ),
          build_R_notebooks)
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )
    statement = ''' cp -a notebooks report/html/notebooks'''
    P.run()
    
@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

    statement = ''' cp -au notebooks report/html/notebooks '''
    P.run()

@follows( update_report ,
          build_R_notebooks)
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
