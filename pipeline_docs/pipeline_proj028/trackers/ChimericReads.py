from ProjectTracker import *

class ChimericReadProportions(TrackerSQL):

    table="chimeric_log"

    def __call__(self, track):

        statement = '''SELECT * FROM %(table)s'''

        results = self.getDataFrame(statement)

        results["fraction"] = results.groupby(["protein", "replicate"])["count"].apply(lambda x: x/x.sum())

        def _add_se(x):
            x["se"] = numpy.sqrt(x["fraction"]*(1-x["fraction"])/x["count"].sum())
            return x
        results = results.groupby(["protein", "replicate"]).apply(_add_se)
        results = results[results.category == "good"]
        return results
