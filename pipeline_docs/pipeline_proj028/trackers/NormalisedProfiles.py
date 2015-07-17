from CGATReport.Tracker import *
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
