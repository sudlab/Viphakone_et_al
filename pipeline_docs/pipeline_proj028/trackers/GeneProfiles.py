from CGATReport.Tracker import *
from ProjectTracker import *


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

        statement = '''SELECT bin, area, region, exons
                       FROM single_vs_multi_exon_gene_profiles
                       WHERE replicate='%(slice)s' AND factor='%(track)s' '''

        df = self.getDataFrame(statement)
        
        df["area"] = pandas.rolling_mean(df["area"], window=10)

        return df



class AverageSingleVsMultiExonProfiles(TrackerSQL):

    def getTracks(self):
        return self.getValues("SELECT DISTINCT factor FROM single_vs_multi_exon_gene_profiles")

    def __call__(self, track):
        
        statement = '''SELECT bin, area, region, exons, replicate
                       FROM single_vs_multi_exon_gene_profiles
                       WHERE factor='%(track)s' '''

        df = self.getDataFrame(statement)

        df = df.groupby(["exons","bin"]).mean()
        df = df.reset_index()
        df = df.set_index(["bin","exons"]).groupby(
            level="exons").transform(lambda x: pandas.rolling_mean(x, window=10)).reset_index()
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

        statement = ''' SELECT position,density FROM %(table)s
                        WHERE track = '%(track)s' AND replicate = '%(slice)s' '''

        results = self.getDataFrame(statement)

        results["density"] = results["density"]/sum(results["density"])
        results["density"] = pandas.rolling_mean(results["density"], 3)

        return results


    
class TranscriptomeExonBoundaryProfiles(ExonBoundaryProfiles):

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
