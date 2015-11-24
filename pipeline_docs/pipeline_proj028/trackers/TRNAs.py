from ProjectTracker import *

import numpy as np
class TRNAScores(ProjectTracker):

    def __call__(self, track):

        statement = ''' SELECT score.tRNA as tRNA,
                               score.num_factors as "Number of Factors",
                               counts.factor as Pulldown,
                               AVG(counts.exon_count) as "Average tags",
                               ts.contig as contig,
                               ts.start as start,
                               ts.end as end
                        FROM
                          trnas_with_two_factors as score
                        INNER JOIN tRNA_counts as counts
                          ON score.tRNA = counts.gene_id
                        INNER JOIN trna_stats as ts
                          ON score.tRNA = ts.gene_id
                        GROUP BY counts.gene_id, counts.factor 
                        ORDER BY score.num_factors'''

        results = self.getDataFrame(statement)
     

        results = results.set_index(["tRNA", "Number of Factors",
                                     "Pulldown", "contig", "start", "end"])["Average tags"]
        results = results.unstack(level="Pulldown")
        results = results.reset_index()
        
#        project_id = P.getProjectId()
#        prefix = PARAMS["report_prefix"]
    
#        tracks_url = "https://www.cgat.org/downloads/%(project_id)s/%(prefix)sexport/tRNAs/UCSC.txt" % locals()

#        url_template = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=h19&position=%s&hg_customText=%s"

#        positions = results.apply(
#            lambda x: "%s:%s-%s" % (x.contig, x.start-500, x.end+500), axis=1)

#        link_urls = positions.apply(lambda x: url_template % (x, tracks_url))
        
#        links = link_urls.apply(lambda x: "`LINK <%s>`_" % x)
        
#        results["Link"] = links
        results = results.drop(["contig", "start", "end"], axis=1)
        results = results.sort("Number of Factors", ascending=False )
        return results


class TRNAHeatMap(ProjectTracker):

    def getTracks(self):
        tracks = self.getValues("SELECT DISTINCT track FROM trna_profiles")
        return [track for track in tracks 
                if "reproducible" not in track]

    pseudo_counts = 1 
    def getProfile(self, track):

        statement = ''' SELECT count, base, tRNA, sum
                         FROM trna_profiles as profiles
                           INNER JOIN trna_stats as stats
                           ON profiles.trna = stats.gene_id
                         WHERE
                           track = '%(track)s'
                         ORDER BY sum '''

        results = self.getDataFrame (statement)
        end = results["sum"].max()
        
        def _expandandnorm(x):
            x = x.set_index("base")
            x = x["count"]/(sum(x["count"])+ self.pseudo_counts)
            x = x.reindex(range(int(end) + 1))
            return x
        results = \
            results.groupby(["tRNA", "sum"]).apply(_expandandnorm)
        results.reset_index(level="tRNA", drop=True, inplace=True)
        results = results.sort_index()
        results = results.reset_index()
        results = results.reset_index()
        results = pandas.melt(results, id_vars=["sum","index"])

        results = results.fillna(0)
        results = results.replace(np.inf, 0)
        return results


    def __call__(self, track):
        return self.getProfile(track)

class TRNAHeatMapWithStruc(TRNAHeatMap):

    slices = [73,74,82]

    slice2struc = {73: ">>>>>>>..>>>>........<<<<.>>>>>.......<<<<<....>>>>>.......<<<<<<<<<<<<.",
                   74: ">>>>.>>..>>>>.........<<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<<.<<<<.",
                   82: ">>>>>>>..>>>...........<<<.>>>>>.......<<<<<.>>>>..<<<<..>>>>>.......<<<<<<<<<<<<."}

    def __call__(self, track, slice):

        results = self.getProfile(track)

        results = results[results["sum"] == slice]

        struc = [x for x in self.slice2struc[slice]]

        struc = pandas.DataFrame({"struc": struc})

        results = pandas.merge(results, struc, left_on="base", right_index=True)

        return results


class tRNAlengthDist(ProjectTracker):

    def __call__(self, track):

        statement = "SELECT sum FROM trna_stats"

        return self.getDataFrame(statement)


class FlipInSubtractedtRNAHeatMaps(TRNAHeatMap):
    

    def getTracks(self):
        tracks = self.getValues("SELECT DISTINCT track FROM trna_profiles")
        tracks= [track for track in tracks 
                 if "union" in track]
        print tracks
        tracks.remove("FlipIn-FLAG.union")

        return tracks

    pseduo_counts = 0

    def __call__(self, track):

        track_results = self.getProfile(track)

        FlipIn_results = self.getProfile("FlipIn-FLAG.union")

        track_results["value"] = track_results["value"] - FlipIn_results["value"]
        
        return track_results



class FlipInSubtractedAverageStruc(TRNAHeatMapWithStruc):

    def getTracks(self):
        tracks = self.getValues("SELECT DISTINCT track FROM trna_profiles")
        tracks= [track for track in tracks 
                 if "union" in track]
        print tracks
        tracks.remove("FlipIn-FLAG.union")

        return tracks

    pseudo_counts = 0

    def __call__(self, track, slice):

        track_result = TRNAHeatMapWithStruc.__call__(self, track, slice)
        flipin_result = TRNAHeatMapWithStruc.__call__(self, "FlipIn-FLAG.union", slice)

        track_result["FlipIn"] = flipin_result["value"]
        track_result["subtracted"] = track_result["value"] - track_result["FlipIn"]
        
        return track_result.groupby(["base","sum","struc"]).mean().reset_index()
