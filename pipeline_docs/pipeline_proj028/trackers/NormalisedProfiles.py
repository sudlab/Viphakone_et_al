from CGATReport.Tracker import TrackerDataframes
from ProjectTracker import ProjectTracker
import numpy as np
import pandas as pd
import CGAT.IOTools as IOTools
import iCLIP.random
import re

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



class FirstLastExonBoot(ProjectTracker):

    data = None
    
    statement = '''SELECT DISTINCT cc.gene_id, cc.exon_id,
                                   %(columns)s,
                                   first,
                                   last,
                                   exon,
                                   intron
                     FROM chunk_counts as cc
                      INNER JOIN annotations.gene_stats as gs
                       ON gs.gene_id = cc.gene_id
                      INNER JOIN annotations.gene_info as gi
                       ON cc.gene_id = gi.gene_id
                      INNER JOIN reference_chunks_gene_first_pc_exons as fe
                       ON fe.gene_id = cc.gene_id AND
                          fe.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_gene_last_pc_exons as le
                       ON le.gene_id = cc.gene_id AND
                          le.exon_id = cc.exon_id 
                      INNER JOIN reference_chunks_exons as ce
                       ON ce.gene_id = cc.gene_id AND
                          ce.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_introns as ie
                       ON ie.gene_id = cc.gene_id AND
                          ie.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_ngenes as ng
                       ON ng.gene_id = cc.gene_id AND
                          ng.exon_id = cc.exon_id
                      WHERE ng.ngenes == 1
                        AND gene_biotype=='protein_coding'
                        AND contig != 'chrM'
                     '''

    track_pattern = "(.+_(?:FLAG|GFP))_union"
    
    def track2control(self, track):
        return "FlipIn_" + re.match(".+_(.+_[Ru].+)", self.track2col[track]).groups()[0]
    
    def grouping(self, row):
        if row["intron"]>0 and row["exon"]==0:
            return "intron"
        elif row["first"]>=1 :
            return "first"
        elif row["last"]>=1 :
            return "last"
        elif  row["intron"]==0:
            return "CDS"
        else:
            return None
    

    def getTracks(self):

        cols = self.getColumns("chunk_counts")
        cols = [track for track in cols if re.search(self.track_pattern, track)]
        tracks = [re.search(self.track_pattern, track).groups()[0]
                  if len(re.search(self.track_pattern, track).groups()) > 0
                  else track
                  for track in cols]
                  
        self.track2col = {track: col for track, col in zip(tracks, cols)}

        return tracks
    
    def fill(self):

        print "filling data...",
        self.getTracks()
        
        columns = self.track2col.values()
        
        columns = ",".join(columns)
        
        self.data = self.getDataFrame(self.statement)
        print self.statement % locals()
        self.data["grouping"] = self.data.apply(self.grouping, axis=1)
        self.data = self.data[self.data.grouping.notnull()]
        print "done"
        
    def __call__(self, track):

        if self.data is None:
            self.fill()
            
        test = self.track2col[track]
        control = self.track2control(track)
        return self.data.groupby("grouping").apply(lambda x: iCLIP.random.ratio_and_ci(x[test], x[control]))
        

class DetainedBoot(FirstLastExonBoot):

    statement = '''SELECT DISTINCT cc.gene_id, cc.exon_id,
                                   %(columns)s,
                                   retained,
                                   exon,
                                   intron,
                                   sig
                     FROM chunk_counts as cc
                      INNER JOIN annotations.gene_stats as gs
                       ON gs.gene_id = cc.gene_id
                      INNER JOIN annotations.gene_info as gi
                       ON cc.gene_id = gi.gene_id
                      LEFT JOIN detained_intron_calls di
                       ON di.gene_id = cc.gene_id AND
                          di.intron_id = cc.exon_id
                      INNER JOIN reference_chunks_retained_introns as ri
                       ON  ri.gene_id = cc.gene_id AND
                           ri.exon_id = cc.exon_id
                       INNER JOIN reference_chunks_exons as ce
                       ON ce.gene_id = cc.gene_id AND
                          ce.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_introns as ie
                       ON ie.gene_id = cc.gene_id AND
                          ie.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_ngenes as ng
                       ON ng.gene_id = cc.gene_id AND
                          ng.exon_id = cc.exon_id
                      WHERE ng.ngenes == 1
                        AND gene_biotype=='protein_coding'
                        AND contig != 'chrM'
                     '''

    def grouping(self, row):

        if row["sig"] == 'TRUE':
            return "Detained"
        elif row["retained"]==1:
            return "Retained"
        elif row["intron"]>0 and row["exon"]==0:
            return "Intron"
        elif  row["intron"] == 0:
            return "Exon"
        else:
            return None

class CombinedDeReBoot(DetainedBoot):

    def grouping(self, row):

        if row["sig"] == 'TRUE':
            return "Retained"
        elif row["retained"]==1:
            return "Retained"
        elif row["intron"]>0 and row["exon"]==0:
            return "Intron"
        elif  row["intron"] == 0:
            return "Exon"
        else:
            return None
    
class ChTopAPABoot(FirstLastExonBoot):

     statement = '''SELECT DISTINCT cc.gene_id, cc.exon_id,
                                   %(columns)s,
                                   apa_site
                     FROM chunk_counts as cc
                      INNER JOIN annotations.gene_stats as gs
                       ON gs.gene_id = cc.gene_id
                      INNER JOIN annotations.gene_info as gi
                       ON cc.gene_id = gi.gene_id
                      LEFT JOIN reference_chunks_gene_last_pc_exons as le
                       ON le.gene_id = cc.gene_id AND
                          le.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_chtop_apa as ca
                       ON  ca.gene_id = cc.gene_id AND
                           ca.exon_id = cc.exon_id
                      INNER JOIN reference_chunks_ngenes as ng
                       ON ng.gene_id = cc.gene_id AND
                          ng.exon_id = cc.exon_id
                      WHERE ng.ngenes == 1 AND
                       le.last > 0
                        AND gene_biotype=='protein_coding'
                        AND contig != 'chrM'
                     '''

     def grouping(self, row):

         if row["apa_site"] == 1:
             return "APA"
         else:
             return "No APA"
         
