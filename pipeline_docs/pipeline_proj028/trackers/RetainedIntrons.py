from ProjectTracker import *

#import CGATReport.Tracker.TrackerMultipleLists as TrackerMultipleLists


class UpregulatedRetainedIntrons(TrackerMultipleLists, ProjectTracker):

    tracks = ["up", "down"]
    slices = ["total", "cyto"]
    statements = {"Up regulated genes":
                  '''SELECT DISTINCT e2e.gene_id
                     FROM %(slice)s_UAP56_DDX39_%(track)s as array
                       INNER JOIN annotations.ensembl_to_entrez as e2e
                         ON e2e.entrez_id = array.EntrezGeneId
                       INNER JOIN biotypes
                         ON biotypes.gene_id = e2e.gene_id
                       INNER JOIN array_platform
                         ON e2e.entrez_id = array_platform.EntrezGeneID''',

                  "Genes with retained introns":
                  ''' SELECT DISTINCT biotypes.gene_id
                      FROM biotypes
                       INNER JOIN annotations.ensembl_to_entrez as e2e
                          ON biotypes.gene_id = e2e.gene_id
                       INNER JOIN array_platform
                          ON array_platform.EntrezGeneID = e2e.entrez_id
                      WHERE biotype = "retained_intron" ''',
                  
                  "background":''' SELECT DISTINCT biotypes.gene_id
                                   FROM biotypes
                                   INNER JOIN annotations.ensembl_to_entrez as e2e
                                        ON biotypes.gene_id = e2e.gene_id
                                   INNER JOIN array_platform
                                     ON e2e.entrez_id = array_platform.EntrezGeneID
                                '''
              }       


class ClustersInRetainedIntrons(ProjectTracker):

    pattern = "(.+)_genes"

    def __call__(self, track):


        statement = '''SELECT DISTINCT gi.gene_id as "Gene ID",
                              gi.gene_name as "symbol"
                       FROM %(track)s_genes as genes
                         INNER JOIN
                            annotations.gene_info as gi
                         ON gi.gene_id = genes.gene_id
                        '''

        
        results = self.getDataFrame(statement)
        

        return results


class RetainedIntronsExpressionVsClips(ProjectTracker):

    pattern = "(.+)_dexseq"

    def getSlices(self):

        slices = self.getValues("SELECT DISTINCT track from retained_intron_clip_tag_counts")
        print "slices are %s" % slices
        return slices

    def __call__(self, track, slice):

        print "called"
        statement = ''' SELECT dex.*,
                               exon_count as clip_tags
                        FROM retained_intron_clip_tag_counts as tag_counts
                          INNER JOIN %(track)s_dexseq as dex
                          ON dex.groupID = tag_counts.gene_id AND
                             dex.featureID = tag_counts.exon_id
                        WHERE track = '%(slice)s' AND
                              genomicData_width >= 10 '''

        return self.getDataFrame(statement)


statement_template = ''' SELECT groupID || featureID
                                    FROM %s_dexseq
                                    WHERE padj < 0.05 
                                       AND %%(track)s'''


class ChangedRIVenn(ProjectTracker, TrackerMultipleLists):

    statements = {"nuclear": statement_template % "nuclear",
                  "cytoplasmic": statement_template % "cytoplasmic",
                  "total": statement_template % "total",
                  "background": '''SELECT groupID || featureID
                                   FROM total_dexseq'''}

    tracks = ["log2fold_Control_Alyref > 0", "log2fold_Control_Alyref < 0"]


class RICandidateRanking(ProjectTracker):

    pattern = "(.+)_dexseq"

    def __call__(self, track):

        statement = '''SELECT groupID as gene_id,
                          featureID as exon,
                          padj,
                          log2fold_Control_Alyref as lfc,
                          (Control + 0.15)/genomicData_width as rna_density,
                          (exon_count + 0.0)/genomicData_width as clip_density,
                          (Control * exon_count)/(genomicData_width * genomicData_width) as score
                   FROM %(track)s_dexseq as dexseq
                        INNER JOIN retained_intron_clip_tag_counts as tag_counts
                         ON dexseq.groupID = tag_counts.gene_id AND
                            dexseq.featureID = tag_counts.exon_id
                   WHERE track = 'Alyref-FLAG.union'
                    AND padj < 0.05 AND log2fold_Control_Alyref > 0.58
                    AND genomicData_width > 20
                   ORDER BY clip_density DESC'''

        return self.getDataFrame(statement)


class FractionsOfDiffIntronsDetained(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT track, in_category, overlap, total-overlap as no_overlap
                    FROM detained_intron_fractions'''

    fields = ("track", "in_category")


class DifferntialIntronsVsCellFraction(ProjectTracker):

    pattern = "(.+)_dexseq"

    def __call__(self, track):

        statement = '''SELECT groupID as gene_id,
                           featureID as exon,
                           exons.padj as exon_padj,
                           log2fold_Control_Alyref as exon_lfc,
                           fractions.ALYREF_Fraction as ALYREF_Nuclear,
                           fractions.Control_Fraction as Control_Nuclear,
                           fractions.padj as Fraction_padj
                       FROM %(track)s_dexseq as exons
                       INNER JOIN stubbs_counts_edgeR as fractions
                        ON fractions.gene_id = exons.groupID     '''

        result =  self.getDataFrame(statement)
        result = result.dropna()
        result["exon_sig"] = ((result.exon_padj < 0.05) & (abs(result.exon_lfc) > 0.58)) * (result.exon_lfc/result.exon_lfc.abs())
        return result
