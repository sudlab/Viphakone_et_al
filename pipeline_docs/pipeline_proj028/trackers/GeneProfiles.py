from CGATReport.Tracker import *
from ProjectTracker import *
import numpy as np

class UTRProfiles(TrackerDataframes):
    ''' Read in normalized profile matricies as dataframes for plotting'''

    pass


class GeneProfiles(TrackerDataframes):
    pass


class SingleVsMultiExonProfiles(TrackerSQL):

    def getSlices(self):

        return self.getValues("SELECT DISTINCT replicate FROM single_vs_multi_exon_gene_profiles")

    def getTracks(self):
        return self.getValues("SELECT DISTINCT factor FROM single_vs_multi_exon_gene_profiles")

    def __call__(self, track, slice):

        statement = '''SELECT bin, density, region, exons
                       FROM single_vs_multi_exon_gene_profiles
                       WHERE replicate='%(slice)s' AND factor='%(track)s' '''

        df = self.getDataFrame(statement)
        
        return df



class AverageSingleVsMultiExonProfiles(TrackerSQL):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT factor FROM single_vs_multi_exon_gene_profiles")

    def __call__(self, track):
        
        statement = '''SELECT bin, density, region, exons, replicate
                       FROM single_vs_multi_exon_gene_profiles
                       WHERE factor='%(track)s' '''

        df = self.getDataFrame(statement)

        df = df.groupby(["exons","bin"]).mean()
        df = df.reset_index()
 
        return df


class BinnedExpressionProfiles(ProjectTracker):
    
    def getTracks(self):
        return self.getValues(
            "SELECT DISTINCT factor FROM binned_expression_profiles")

    def getSlices(self):
        return self.getValues(
            "SELECT DISTINCT rep FROM binned_expression_profiles WHERE rep != 'reproducible'")

    def __call__(self, track, slice):

        statement = ''' SELECT bin, area, region, quantile, exon_limit
                        FROM binned_expression_profiles
                        WHERE factor='%(track)s'
                        AND rep='%(slice)s' '''

        df = self.getDataFrame(statement)
        df.exon_limit.replace([0,1],["All transcripts","Multi-exon"], inplace=True)
        return df


class ExpressedTranscriptStats(ProjectTracker):

    def __call__(self, track):

        statement = '''SELECT sf.transcript_id, TPM, sum as "Exon Length"
                              
                       from expressed_transcripts_stats as ets
                       INNER JOIN
                       HEK293_sailfish as sf
                       ON sf.transcript_id = ets.transcript_id
              
                       WHERE source='protein_coding'  '''

        df = self.getDataFrame(statement)

        df["quantile"] = pandas.qcut(df["Exon Length"], 5, labels=False)
        return df


class ExonBoundaryProfiles(ProjectTracker):

    table = "exon_boundary_profiles"

    def getTracks(self):
        return self.getValues("SELECT DISTINCT track FROM %(table)s")

    def getSlices(self):
        return self.getValues("SELECT DISTINCT replicate FROM %(table)s")

    def __call__(self, track, slice):

        statement = ''' SELECT track as factor, position,density FROM %(table)s
                        WHERE track = '%(track)s' AND replicate = '%(slice)s' '''

        results = self.getDataFrame(statement)

        results["density"] = results["density"]/sum(results["density"])
        results["density"] = pandas.rolling_mean(results["density"], 3)

        return results


class ExonBoundaryProfilesWith4A3(ExonBoundaryProfiles):

    def __call__(self, track, slice):

        results =  ExonBoundaryProfiles.__call__(self, track, slice)
        eif4a3 = ExonBoundaryProfiles.__call__(self, "eIF4A3-GFP", "union")

        results = pandas.concat([results, eif4a3])

        return results


class CenteredExonBoundaryProfiles(ExonBoundaryProfilesWith4A3):

    table = "centered_exon_boundary_profiles"

    
class FirstExonBoundaryProfiles(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT protein, replicate, base, density 
                   FROM first_exon_profiles'''
    fields = ('protein',)


class FirstExonBoundaryProfiles2(TrackerDataframes):

    glob = "heatmaps/*.first_exon.end.matrix.tsv.gz"
    regex = "heatmaps/(.+).first_exon.end.matrix.tsv.gz"

    def __call__(self, track):

        df = pandas.read_csv("heatmaps/%s.first_exon.end.matrix.tsv.gz" % track,
                             sep='\t',
                             index_col=0)

        df.columns = df.columns.astype("float")
        rowsums= df.sum(axis=1)
        df_normed = df.div(rowsums.astype("float"), axis=0)
        df_summed = df_normed.sum()
        df_summed.name="density"
        return df_summed.reset_index().query("index > -200 & index < 50", engine="python")

    
class TranscriptomeExonBoundaryProfiles(ExonBoundaryProfiles):

    table = "transcriptome_exon_boundary_profiles"


class TranscriptomeExonBoundaryProfilesWith4A3(ExonBoundaryProfilesWith4A3):

    table = "transcriptome_exon_boundary_profiles"

    
class StubbsProfiles(ProjectTracker, SQLStatementTracker):

    fields = ("fraction", "condition")

    statement = '''SELECT fraction, condition, region, replicate, bin, area
                   FROM stubbs_profiles'''
 
    def __call__(self):

        results = self.getDataFrame(self.statement)

        results = results.groupby(["fraction", "condition", "region","bin"])["area"].mean().reset_index()

        results["area"] = results.groupby(["fraction","condition","region"]
                        )["area"].transform(pandas.rolling_mean, window=3
                        )

        return results


class IntronProfiles(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT factor, replicate, length, position, density
                   FROM intron_profiles'''

    fields = ("factor","replicate")


class HighlyExpressed3BiasGenes(ProjectTracker):

    def __call__(self, track):

        statement = ''' SELECT gi.gene_id, gi.gene_name, ti.transcript_id as transcript_id,
                               utr3_enrichment, utr3_count,
                               500000000.0 * (Control_Total_R1 + Control_Total_R2) /
                                ((SELECT sum(Control_Total_R1 + Control_Total_R2) 
                                  FROM stubbs_counts) * gs.sum ) as expression
                        FROM profile_summaries
                           INNER JOIN annotations.transcript_info as ti
                             ON ti.transcript_id = profile_summaries.transcript_id 
                           INNER JOIN annotations.gene_info as gi
                             ON gi.gene_id = ti.gene_id
                           INNER JOIN stubbs_counts as counts
                             ON gi.gene_id = counts.geneid
                           INNER JOIN annotations.gene_stats as gs
                             ON gs.gene_id = gi.gene_id
                        WHERE
                             protein = 'Chtop' AND
                             expression > 10 '''

        
        results = self.getDataFrame(statement)
        

        results.replace("nan", np.nan, inplace=True)
        results.dropna(inplace=True)
        results.set_index("transcript_id", inplace=True)
        info = results[["gene_id","gene_name"]]
        data = results.drop(["gene_id", "gene_name"], axis=1)
        data = data.astype("float64")
        
        data = data.groupby(
            level="transcript_id").mean()

        log_ranks = data.rank().apply(np.log10)
        log_ranks = log_ranks.drop("expression", axis = 1)
        #log_ranks["expression"] = log_ranks["expression"]/4
        log_ranks = log_ranks.sum(axis=1)
        log_ranks.name = "score"
        data = data.applymap(lambda x: "%.2f" % x)
        
        results = info.join(log_ranks).join(data).drop_duplicates()
        results.sort("score", inplace=True, ascending=False)
        results["score"] = results["score"].apply(lambda x: "%.3f" % x)
        
        return results.iloc[0:100].reset_index()


class HistoneMetaGenes(ProjectTracker):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT factor FROM histone_profiles")

    def getSlices(self):
        return self.getValues("SELECT DISTINCT replicate FROM histone_profiles")

    def __call__(self, track, slice=None):
        print "called with %s, %s" % (track,slice)
        statement = '''SELECT *  FROM histone_profiles WHERE factor = '%(track)s' AND replicate = '%(slice)s' AND geneset != 'non-histones' '''
        data = self.getDataFrame(statement % locals())
    #    data["density"] = data.groupby("geneset")["density"].apply(lambda x: x/x.sum())
        return data


