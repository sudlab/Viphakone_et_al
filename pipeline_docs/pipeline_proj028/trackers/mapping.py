from ProjectTracker import *
import pandas
import re

class FinalStatsTable(ProjectTracker):

    def __call__(self, track):

        statement = ''' SELECT dd.track as track, 
                              
                               dd.counts as 'Unique Reads',
                               clb.count as 'Genic cross-linked bases'
                               
                        FROM 
                              iclip.deduped_bam_stats as dd
                        INNER JOIN 
                              iclip.cross_linked_bases as clb ON
                           dd.track = clb.track
                        WHERE dd.category = 'reads_mapped' '''

        counts = self.getDataFrame(statement)
        counts.set_index("track", inplace=True)
        statement = ''' SELECT * FROM track_counts '''

        gene_counts = self.getDataFrame(statement)

        gene_counts = pandas.melt(gene_counts, id_vars="Geneid",
                                  var_name="track", value_name="counts")

        genes_with_clip = gene_counts.groupby(
            ["track"]).apply(lambda x: (x["counts"] > 0).sum())

        genes_with_clip.index = [re.sub("_","-", x) for x in genes_with_clip.index.values]
        
        print genes_with_clip
        counts["Genes with clipped bases"] = genes_with_clip

        counts.reset_index(inplace=True)
        return counts

