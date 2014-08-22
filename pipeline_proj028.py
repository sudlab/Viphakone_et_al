################################################################################
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

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGATPipelines.PipelineUtilities as PUtils


###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

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

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
## worker tasks
###################################################################
@transform(PARAMS["annotations_geneset"],
           suffix(".gtf.gz"),
           ".fasta.gz")
def getGenesetFasta(infile, outfile):
    '''convert geneset annotation in fasta '''

    statement = ''' python %(scriptsdir)s/gff2fasta.py
                          --is-gtf
                          --genome=%(genome_dir)s/%(genome)s
                          -I %(infile)s
                          -L %(outfile)s.log
                  | cut -d' ' -f1
                  | gzip > %(outfile)s '''
    P.run()


###################################################################
@transform(getGenesetFasta, 
           regex(".+"),
           "sailfish_index.dir/transcriptome.sfi")
def generateSailfishIndex(infile,outfile):
    ''' generate sailfish index for specified geneset '''

    tmpfile = P.getTempFilename()
    logfile = os.path.basename(outfile) + ".log"
    statement = ''' zcat %(infile)s > %(tmpfile)s;
                    checkpoint;
                    sailfish index -t %(tmpfile)s
                                  -o sailfish_index.dir
                                  -k %(sailfish_kmer)s 
                    >& %(logfile)s;
                    checkpoint;
                    rm %(tmpfile)s '''
    P.run()


###################################################################
@collate("*fastq*gz",
         regex("(.+).fastq.(?:[12]\.)*gz"),
         r"\1_sailfish.dir/quant.sf")
def runSailFish(infiles,outfile):
    ''' Run sailfish on any provided rnaseq files '''


    track = re.match("(.+)_sailfish.dir/quant.sf", outfile).groups()[0]
    if len(infiles) == 1:
        inputs = "-r %s" % infiles[0]
    elif len(infiles) == 2:
        inputs = "-1 <(zcat %s ) -2 <( zcat %s)" % infiles
    else:
        raise ValueError("Don't know how to handle %i input files"
                         % len(infiles))

    statement = ''' sailfish quant -i sailfish_index.dir 
                                   -l '%(sailfish_libtype)s'
                                   %(inputs)s
                                   -o %(track)s_sailfish.dir '''

    P.run()


###################################################################
@transform(runSailFish, 
           regex("(.+).dir/quant.sf"),
           r"\1.load")
def loadSailfish(infile, outfile):
    P.load(infile, outfile, 
           options="--header=transcript_id,length,TPM,RPKM,KPKM,nKmers,nReads -i transcript_id")


###################################################################
# Start/Stop Codon profiles
###################################################################

@transform(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex("(.+)"),
           add_inputs(loadSailfish),
           "expressed_transcripts.gtf.gz")
def filterExpressedTranscripts(infiles, outfile):
    '''Filter the geneset for those genes that are expressed in 
    HEK293 '''


    infile = infiles[0]
    
    query=''' SELECT transcript_id FROM HEK293_sailfish
              WHERE TPM > 1 '''

    transcript_ids=PUtils.fetch(query)
    
    tmp = P.getTempFilename(dir="/ifs/scratch/")

    PUtils.write(tmp, transcript_ids)

    statement = '''python %(scriptsdir)s/gtf2gtf.py
                          --filter=transcript
                          --apply=%(tmp)s
                          -I %(infile)s
                          -L %(outfile)s.log
                 | gzip -c > %(outfile)s '''

    P.run()

@follows(mkdir("gene_profiles.dir"))
@transform([os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
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
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=utrprofile
                          --bamfile=%(bamfile)s
                          --gtffile=<(zcat %(gtffile)s | grep protein_codin)
                          --normalization=total-max
                          --base-accuracy
                          --normalize-profile=area
                          --scale-flank-length=1
                          --resolution-downstream=1000
                          --resolution-upstream=1000
                          --resolution-upstream-utr=200
                          --resolution-downstream-utr=700
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.CDS_profile.log 
                          --output-all-profiles '''

    P.run()

###################################################################
@follows(mkdir("gene_profiles.dir"))
@transform([os.path.join(PARAMS["dir_iclip"], "deduped.dir/*.bam"),
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
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=tssprofile
                          --bamfile=%(bamfile)s
                          --gtffile=<(zcat %(gtffile)s | awk '$3=="CDS"' | sed 's/CDS/exon/')
                          --normalization=total-sum
                          --base-accuracy
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

@follows(calculateSTOPProfiles, calculateCDSProfiles)
def profiles():
    pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(  )
def full(): pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
