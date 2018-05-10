'''
calculate_effective_length.py - calculate mappability adjusted lengths
====================================================

:Author:
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
from CGAT import GTF
from CGAT import Database
import pyBigWig
import numpy

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-b", "--mappability-file", dest="mappability", type="string",
                      help="Bigwig file with mappability")
    parser.add_option("-d", "--database", dest="db", type="string",
                      help="Database containing intron chunck table")
    parser.add_option("-l", "--read-length", dest="rlen", type="int",
                      default=50,
                      help="Read length")
    parser.add_option("-o", "--overlap-length", dest="olen", type="int",
                      default=10,
                      help="Min overlap before read is counted")
    parser.add_option("-M","--multimap", dest="mm", action="store_true",
                      default=False,
                      help="Allow multimapping")
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    introns = Database.fetch("""SELECT gene_id, exon_id 
                                 FROM reference_chunks_introns 
                                 WHERE intron>0""",
                             Database.connect(options.db))

    introns = [tuple(x) for x in introns]
    
    mappability = pyBigWig.open(options.mappability)
    options.stdout.write("\t".join(["gene_id", "exon_id", "efflen"])+"\n")
    for exon in GTF.iterator(options.stdin):

        if (unicode(exon.gene_id), int(exon.exon_id)) not in introns:
            continue
        
        vals = mappability.values(exon.contig,
                                  int(exon.start)-(options.rlen - options.olen),
                                  int(exon.end) - options.olen)
        if options.mm:
            eff_len = sum(vals)

        else:
            eff_len = int(exon.end) - int(exon.start) \
                  + options.rlen \
                  - 2 * options.olen \
                  - len([x for x in vals if x < 1])

        options.stdout.write(
            "\t".join([exon.gene_id, exon.exon_id, str(eff_len)]) + "\n")
    
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
