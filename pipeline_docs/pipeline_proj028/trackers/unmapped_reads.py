from ProjectTracker import ProjectTracker
import pandas
import numpy
import re


class UnmappedComposition(ProjectTracker):

    def getTracks(self):
        paths = self.getValues(
            "SELECT DISTINCT track FROM unmapped_composition")
        paths = [re.sub("-FLAG-", "-", track) for track in paths]
        return paths

    slices = ["A", "T", "G", "C"]

    def __call__(self, track, slice):

        track = re.sub("-","-FLAG-", track)
 
        statement = ''' SELECT read_length, "count", "nreads"
                        FROM unmapped_composition
                        WHERE track='%(track)s' AND base='%(slice)s' '''

        data = self.getDataFrame(statement)

        data["density"] = data["count"].astype("float")/data["read_length"]

        data["bin"] = pandas.cut(data["density"], numpy.arange(0, 1.1, 0.1),
                                 labels=numpy.arange(0, 1.0, 0.1))

        sums = data.groupby("bin")["count"].sum().reset_index()

        return sums
