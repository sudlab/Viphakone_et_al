from ProjectTracker import *

class FractionReproducible(ProjectTracker):

    slices=["p_threshold_scan"]

    def getTracks(self):

        table = self.slices[0]
        statement = "SELECT DISTINCT track FROM %(table)s"
        return self.getValues(statement)

    def __call__(self, track, slice):

        statement = ''' SELECT * FROM %(slice)s
                        WHERE track='%(track)s' AND
                              (replicate='reproducible' OR
                              replicate='all') '''

        results = self.getDataFrame(statement)
    
        results = results.pivot(index="p",
                                columns="replicate",
                                values="count")

        results = results.reset_index()

        return results

class ClustersCalled(ProjectTracker):

    def getTracks(self):
        table = "p_threshold_scan"
        statement = "SELECT DISTINCT track FROM %(table)s"
        return self.getValues(statement)

    def __call__(self, track):

        statement = ''' SELECT replicate,p,count FROM p_threshold_scan
                        WHERE track='%(track)s' AND 
                              NOT 
                              replicate='all' '''

        return self.getDataFrame(statement)
