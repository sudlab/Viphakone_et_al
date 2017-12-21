from ProjectTracker import *
from scipy.stats import hypergeom
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

class ChunkSplicing(ProjectTracker):

    tracks = [" ", "not"]

    slices = ["nuclear", "cytoplasmic", "total", "chtop"]

    def __call__(self, track, slice):

        if slice=="chtop":
            comp="Chtop"
        else:
            comp="Alyref"
            
        statement = ''' SELECT DISTINCT gene_name as symbol,
                               genomicData_seqnames as contig,
                               genomicData_start as intron_start,
                               genomicData_end as intron_end,
                               genomicData_width,
                               featureID, groupID,
                               padj,
                               log2Fold_%(comp)s_Control as l2fold,
                               Alyref_FLAG_union as clip_tags
                        FROM chunk_counts
                          INNER JOIN chunk_splicing as dex
                          ON dex.groupID = chunk_counts.gene_id AND
                             dex.featureID = chunk_counts.exon_id
                          INNER JOIN annotations.gene_info as gi
                          ON gi.gene_id = dex.groupID
                        WHERE dex.track = '%(slice)s' AND
                              genomicData_width >= 10 '''

        results = self.getDataFrame(statement)
        
        if track=="combined":
            track = "retained or detained"
        detained_calls = self.getDataFrame('''SELECT * FROM detained_intron_calls ''')

        detained_calls["detained"] = (detained_calls.padj < 0.01) & (detained_calls.log2FoldChange > 2)
        detained_calls = detained_calls.groupby(["gene_id", "intron_id"]).detained.max()
        retained_calls = self.getDataFrame('''SELECT * FROM reference_chunks_retained_introns''')
        detained_calls = retained_calls.join(detained_calls, on=["gene_id", "exon_id"], how="left")
        detained_calls.detained = detained_calls.detained.fillna(False)
        detained_calls.retained = detained_calls.retained > 0
        detained_calls = detained_calls.query("%s (%s)" % (track, self.intron_set), engine="python")
        detained_calls.set_index(["gene_id", "exon_id"], inplace=True)
        results = results.join(detained_calls, on=["groupID", "featureID"], how="inner")
        results.l2fold = results.l2fold.astype("float64")
        results.to_csv("%s_%s.csv" %(track, slice), sep="\t")
        return results
        

class DetainedChunkSplicing(ChunkSplicing):

    intron_set="detained"


class RetainedChunkSplicing(ChunkSplicing):

    intron_set="retained"


class CombinedChunkSplicing(ChunkSplicing):

    intron_set="detained or retained"


class ChunkEnrichment(ChunkSplicing):

    tracks = ["all"]
    
    def __call__(self, track, slice):

        diff = ChunkSplicing.__call__(self, " ", slice)
        notdiff = ChunkSplicing.__call__(self, "not", slice)

        results = pandas.concat([diff, notdiff], keys=["tained", "notained"], names=["intron"]).reset_index()
                                
        results["clip_bin"] = pandas.cut(pandas.np.log2((results["clip_tags"]/results["genomicData_width"]) + 1e-5), bins=20,
                                         labels=False)

        def _frac(x):
            sig = sum(pandas.notnull(x["padj"]) & (x["padj"] < 0.01) & (x["l2fold"] < -0.58))
            total = x.shape[0]
            return pandas.Series([sig,total], index=["sig","total"])
        
        fractions = results.groupby(["clip_bin", "intron"]).apply(_frac)
        fractions = fractions.reset_index()

        def _enrich(x):

            x = x.set_index("intron")


            try:
                N = x.loc["tained"]["total"]
                k = x.loc["tained"]["sig"]
            except KeyError:
                N = 0
                k = 0
                
            M = N + x.loc["notained"]["total"]
            n = k + x.loc["notained"]["sig"]
           

            p = hypergeom.sf(k-1, M, n, N)

            try:
                enrich = (float(k)/N)/(float(n)/M)
            except ZeroDivisionError:
                enrich = pandas.np.nan
            except KeyError:
                enrich = 0
                
            return pandas.Series([enrich, p], index=["enrichment", "p"])
            

        enrichment = fractions.groupby("clip_bin").apply(_enrich)

        enrichment = enrichment.reset_index()

        enrichment = enrichment.sort_values("clip_bin")
        
        return enrichment

    
class DetainedEnrichment(ChunkEnrichment):

    intron_set="detained"

       
class RetainedEnrichment(ChunkEnrichment):

    intron_set="retained"

class CombinedEnrichment(ChunkEnrichment):

    intron_set="detained or retained"


class CountDiffChunks(ProjectTracker):

    slices  = ["> 0.58", "< -0.58"]
    tracks = ["nuclear", "cytoplasmic", "total"]

    def __call__(self, track, slice):
        statement = ''' SELECT COUNT (DISTINCT groupID) as genes
                    FROM
                    chunk_splicing
                    WHERE padj < 0.05 AND log2fold_Alyref_Control %(slice)s
                          AND track='%(track)s' '''

        return self.getValue(statement % locals())
