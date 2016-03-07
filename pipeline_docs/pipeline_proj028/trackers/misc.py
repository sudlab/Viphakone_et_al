from ProjectTracker import *
from Sample_QC import ContextStats

class CustomisedContextStats(ProjectTracker):

    def getPaths(self):
        tracks = self.getValues("SELECT DISTINCT track from iclip.deduped_context_stats")

        paths = [re.match("(.+)-FLAG-(.+)",x).groups() for x in tracks]
        paths = zip(*paths)

        return paths
    
    categories = {"antisense": "Long ncRNA",
                  "intron": "Intron",
                  "retained_intron": "Retained intron",
                  "lincRNA": "Long ncRNA",
                  "miRNA": "Small ncRNA",
                  "miscRNA": "Other",
                  "none": "Genomic",
                  "nonsense_mediated_decay": "Other",
                  "polymorphic_pseudogene": "Long ncRNA",
                  "processed_pseudogene": "Long ncRNA",
                  "processed_transcript": "Other",
                  "protein_coding": "Protein coding exon",
                  "pseudogene": "Long ncRNA",
                  "rRNA": "ribosomal RNA",
                  "snoRNA": "Small ncRNA",
                  "snRNA": "Small ncRNA",
                  "unitary_pseudogene": "Long ncRNA",
                  "unprocessed_pseudogene": "Long ncRNA"}

    def __call__(self, track, slice=None):

       # categories = "','".join(self.categories)
        statement = ''' SELECT category, alignments
                      FROM deduped_context_stats
                      WHERE track = '%(track)s-FLAG-%(slice)s' AND category != 'total' '''

        results = self.getDataFrame(statement)
        results.category = [self.categories[x] if x in self.categories else "Other"
                            for x in results.category]
        results = results.groupby("category").sum().reset_index()
                            
        return results

class NoFlipInContext(ProjectTracker, ContextStats):

    method="no_flipin"


class GeneProfiles3(ProjectTracker):

    def getSlices(self):

        return self.getValues("SELECT DISTINCT rep FROM gene_profiles")

    def getTracks(self):
        return self.getValues("SELECT DISTINCT factor FROM gene_profiles")

    def __call__(self, track, slice):

        statement = '''SELECT bin, area, region
                       FROM gene_profiles
                       WHERE rep='%(slice)s' AND factor='%(track)s' AND
                       region != "introns" '''

        df = self.getDataFrame(statement)
        
        df["area"] = pandas.rolling_mean(df["area"], window=5)
        df.bin = numpy.arange(1,df.shape[0]+1)
        return df

class GeneFractions(ProjectTracker):

    def getTracks(self):
        return set([re.match("(.+)_FLAG",x).groups()[0]
                    for x in self.getColumns("track_counts")
                    if re.match(".+_FLAG_[Ru]",x)])

    slices = ["R1","R2","R3","R4", "union"]

    def __call__(self, track, slice):

        statement = '''SELECT DISTINCT counts.Geneid, 
                              %(track)s_FLAG_%(slice)s as count
               FROM track_counts as counts
                INNER JOIN annotations.gene_info as gi
                ON gi.gene_id = counts.Geneid
               WHERE gene_biotype = 'protein_coding'
                  AND counts.HEK293_WT_1 > 0 '''
               
        tags_per_gene = self.getDataFrame(statement).set_index("Geneid")
               
        return (tags_per_gene['count'] > 0).sum()/float(tags_per_gene.shape[0])

statement_template = ''' SELECT DISTINCT counts.Geneid
                             FROM track_counts as counts
                               INNER JOIN annotations.gene_info as gi
                               on gi.gene_id = counts.Geneid
                             WHERE      counts.HEK293_WT_1 > 0  AND
                                  %s_FLAG_%%(slice)s > 1 '''
 

class ClippedFractionOverlaps(ProjectTracker,TrackerMultipleLists):

    slices = ["R1","R2","R3","R4", "union"]


        
    ListA = statement_template % "Alyref"
    ListB = statement_template % "Chtop"
    ListC = statement_template % "Nxf1"

    labels = ["Alyref", "Chtop", "Nxf1"]


class UnclippedFraction(ProjectTracker, SQLStatementTracker):


    fields = ("Gene ID",)
    statement = ''' SELECT DISTINCT counts.Geneid as "Gene ID",
                           gene_name as Symbol,
                         (HEK293_banks_2_star + HEK293_banks_3_star + HEK293_banks_6_star +0.0)/3 as "RNA counts"
                    FROM track_counts as counts
                        INNER JOIN annotations.gene_info as gi
                        ON gi.gene_id = counts.Geneid
                    WHERE Nxf1_FLAG_union == 0 AND
                          Alyref_FLAG_union  == 0 AND
                          Chtop_FLAG_union == 0 AND
                    "RNA counts" >= 1
                    ORDER BY "RNA Counts" DESC '''



