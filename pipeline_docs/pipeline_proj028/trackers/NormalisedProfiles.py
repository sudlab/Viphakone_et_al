from CGATReport.Tracker import TrackerDataframes
from ProjectTracker import ProjectTracker
import numpy as np
import pandas as pd
import CGAT.IOTools as IOTools

class NormalizedUTRMatrix(TrackerDataframes):
    ''' Read in normalized profile matricies as dataframes for plotting'''

    pass


    
class NormalisdSummaryUTRMatrix(TrackerDataframes):
    ''' Read in normalized profile matricies as dataframes for plotting'''

    def __call__(self, track, **kwargs):

        profile = TrackerDataframes.__call__(self, track, **kwargs)

        rnaseq = pd.read_csv(IOTools.openFile(
            "gene_profiles.dir/HEK293-WT-1.utrprofile.matrix.tsv.gz"),
                             sep="\t",
                             header=0)

        profile["normalized"] = profile["area"]/rnaseq["area"]

        profile["normalized"] = profile["normalized"].replace(np.inf, 0)
        profile["normalized"] = profile["normalized"].fillna(0)

        return profile


class NormalizedTSSMatrix(TrackerDataframes):

    def __call__(self, track, **kwargs):

        profile = TrackerDataframes.__call__(self, track, **kwargs)

        rnaseq = pd.read_csv(IOTools.openFile(
            "gene_profiles.dir/HEK293-WT-1.tssprofile.matrix.tsv.gz"),
                             sep="\t",
                             header=0)

        profile["normalized"] = profile["area"]/rnaseq["area"]

        profile["normalized"] = profile["normalized"].replace(np.inf, 0)
        profile["normalized"] = profile["normalized"].fillna(0)

        return profile


class FirstLastExons(ProjectTracker):

    table = "seperate_exon_profiles"
    
    def getTracks(self):
        
        return self.getValues("SELECT DISTINCT factor FROM %(table)s")

    def getSlices(self):

        return self.getValues("SELECT DISTINCT replicate FROM %(table)s")


    def __call__(self, track, slice=None):


        statement = '''SELECT bin, area, region_bin
                         FROM %(table)s
                        WHERE factor = '%(track)s'
                              AND replicate = '%(replicate)s' '''

        return self.getDataFrame(statement)


class NormalisedFirstLastExons(FirstLastExons):

    norm_tracks = "Nuclear-RiboZ"

    def __call__(self, track, slice):

        norm_data = self.getDataFrame('''SELECT bin, area
                                           FROM %(table)s
                                         WHERE factor='%(norm_tracks)s' ''')

        norm_data.groupby("bin").mean()

        track_data = self.getDataFrame('''SELECT bin, area, region_bin
                                           FROM %(table)s
                                           WHERE factor = '%(track)s'
                                               AND replicate = '%(slice)s' ''')

        if track_data.shape[0] == 0:
            return pd.DataFrame()

        track_data.set_index("bin", inplace=True)

        track_data["normed"] = track_data["area"]/norm_data["area"]

        return track_data

class GFPNormalisedFirstLastExons(FirstLastExons):

    norm_tracks = "FlipIn-GFP"

    def __call__(self, track, slice):

        norm_data = self.getDataFrame('''SELECT bin, area
                                           FROM %(table)s
                                         WHERE factor='%(norm_tracks)s' ''')

        norm_data.groupby("bin").mean()

        track_data = self.getDataFrame('''SELECT bin, area, region_bin
                                           FROM %(table)s
                                           WHERE factor = '%(track)s'
                                               AND replicate = '%(slice)s' ''')

        if track_data.shape[0] == 0:
            return pd.DataFrame()

        track_data.set_index("bin", inplace=True)

        track_data["normed"] = track_data["area"]/norm_data["area"]

        return track_data


class FirstLastExonCount(ProjectTracker):

    def getTracks(self):

        return list(self.getValues(
            "SELECT DISTINCT cell FROM first_last_exon_counts"))

    def getSlices(self):

        return list(self.getValues(
            "SELECT DISTINCT replicate FROM first_last_exon_counts"))

    def __call__(self, track, slice):

        statement = '''SELECT protein, first_exon, middle_exon, last_exon, introns
                         FROM first_last_exon_counts
                        WHERE cell = '%(track)s'
                              AND replicate = '%(slice)s' '''

        data = self.getDataFrame(statement)
        if data.shape[0] == 0:
            return pd.DataFrame({"protein":[], "exon":[], "normed_count": []})
        data.set_index("protein", inplace=True)
        row_sums = data.sum(axis=1)
        gene_normed = data.loc[row_sums > 0]
        
        row_sums = gene_normed.sum(axis=1)
        print data.shape
        print row_sums.shape

        gene_normed = gene_normed.div(row_sums, axis=0)

        summed_data = gene_normed.groupby(level="protein").mean()

        normed = summed_data/summed_data.loc["FlipIn"]
        normed = normed.drop("FlipIn")


        melted = pd.melt(normed.reset_index(), id_vars="protein",
                         var_name="exon", value_name="normed_count")

        summed_data = data.groupby(level="protein").mean()

        normed = summed_data/summed_data.loc["FlipIn"]
        normed = normed.drop("FlipIn")

        melted2 = pd.melt(normed.reset_index(), id_vars="protein",
                          var_name="exon", value_name="unnormed_count")

        melted = pd.merge(melted, melted2, how='inner', on=['protein','exon'])

        return melted



