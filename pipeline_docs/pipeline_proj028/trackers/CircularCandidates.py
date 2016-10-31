from ProjectTracker import *

from rpy2.robjects.packages import importr
from rpy2.robjects import r as R



class GenomePlot(ProjectTracker):
    ''' This tracker returns a R Gviz plot of the genome and some datatracks.

    One user rendered plot is returned for each track.

    Subclasses should implement a *getTracks* method that regions either a list
    of gene ids, transcript ids or regions (chr:start-end). You should specify
    which you are returning by setting the attribute track_type to either
    "gene", "transcript" or "region".

    The *DataTracks* attribute is a dictionary that specifies what data tracks
    to plot. The top level key will be the title for the track, and the
    dictionary, should contain at least the file name (in the "file" entry) and
    the Gviz plot type (in the "type" entry). You may also specify an
    "import options" entry that passes options to the Gviz DataTrack importer.
    An example is:

    DataTracks = {"First Bam": {"file": "1.bam",
                                "type": "h"},
                  "Second Bam": {"file": "2.bam",
                                 "type": "h"}}

    Specifies that you want to DataTracks called First Bam and Second Bam, that
    they should be loaded from the 1.bam and 2.bam files respectively and that
    they should be plotted with histograms. You may also specify a default plot
    type with DataTracks_type.

    If you only wish some of the DataTracks to be shown for each track, then
    write a getDataTracks method, which should take a track name and a gene
    info dataframe and return a list of DataTracks objects. This is made easier
    by the fact that when getDataTracks is called, self.DataTracks should
    contain a "obj" slot with a DataTrack object in it. For example you might
    want to filter based on strand.

    Various other attributes control the plot:

    :attribute: gene_track_options options dictionary for creating the gene
                                   track
    :attribute: plot_options       options dictionary to the plotTracks command
    :attribute: extend             bp to extend either side of the specified
                                   region
    :attribute: height, width      height and width of plot in pixels '''

    DataTracks = {}
    DataTracks_type = "h"
    track_type = "gene"
    gene_track_options = {"showId": True,
                          "transcriptAnnotation": "transcript"}
    plot_options = {}
    extend = 100
    height = 7
    width = 14

    def getTracks(self):
        raise NotImplementedError(
            "This tracker is not fully implemented. No getTracks method")
    
    def populateDataTracks(self):
        return self.DataTracks

    def getDataTracks(self, track):

        return [data_track["obj"] for data_track in self.DataTracks]
    
    def __init__(self, *args, **kwargs):
        
        ProjectTracker.__init__(self, *args, **kwargs)

        global Gviz, GenomicFeatures, AnnotationDbi

        Gviz = importr("Gviz")
        GenomicFeatures = importr("GenomicFeatures")
        AnnotationDbi = importr("AnnotationDbi")
        try:
            db_name = P.snip(self.annotations, ".gtf.gz") + ".sqlite"
        except AttributeError:
            raise NotImplementedError(
                "No annotations file set")

        if not os.path.exists(db_name):
            print "Making Transcription Database"
            GenomicFeatures = importr("GenomicFeatures")
            txdb = GenomicFeatures.makeTxDbFromGFF(self.annotations,
                                                   format = "gtf")
            self.txdb = txdb
            print "Caching Annotation Database"
            AnnotationDbi.saveDb(txdb, db_name)
            print "Done"
        else:
            print "loading cached database"
            self.txdb = AnnotationDbi.loadDb(db_name)
            print "done"

        self.DataTracks = self.populateDataTracks()

        for track in self.DataTracks:
            self.DataTracks[track]["obj"] = Gviz.DataTrack(
                self.DataTracks[track]["file"],
                type=self.DataTracks[track].get("type", self.DataTracks_type),
                name=track,
                **self.DataTracks[track].get("import options", {}))
        
    def __call__(self, track):

        region_statement = '''
            SELECT MIN(exons.start) as start,
                   MAX(exons.end) as end,
                   exons.contig as contig
            FROM   annotations.exon_stats as exons
              INNER JOIN
                   annotations.transcript_info as ti
                ON exons.transcript_id = ti.transcript_id
              INNER JOIN
                   annotations.gene_info as gi
                ON gi.gene_id = ti.gene_id
            WHERE
                  %s '''

        if self.track_type == "region":
            try: 
                chrom, start, end = re.match(
                    "(.+):([0-9]+)-([0-9]+)", track).groups()
            except AttributeError:
                raise ValueError(
                    "%s is not a valid region specification" % track)
        else:
            if self.track_type == "gene id":
                where = "gi.gene_id = '%s'" % track

            elif self.track_type == "gene name":
                where = "gi.gene_name = '%s' " % track
            elif self.track_type == "transcript":
                where = "ti.transcript_id == '%s'" % track
            else:
                raise NotImplementedError(
                    "Track type: %s not implemented" % self.track_type)

            start, end, chrom = self.getFirstRow(region_statement % where)

        gene_track = Gviz.GeneRegionTrack(self.txdb,
                                          chromosome=chrom,
                                          start=start,
                                          end=end,
                                          **self.gene_track_options)

        data_tracks = self.getDataTracks(track)        
        axisTrack = Gviz.GenomeAxisTrack()

        all_tracks = [axisTrack,gene_track] + data_tracks

        if not os.path.exists("export/GenomePlots"):
            os.makedirs("export/GenomePlots")

        # Hack to get around problem with user render not being able
        # to find font "sans"
        filename = os.path.join("export/GenomePlots",
                                self.__class__.__name__ + track + ".png")

        R.png(filename,
              units="in",
              res=200,
              height=self.height, 
              width=self.width)

        Gviz.plotTracks(all_tracks, main=track, **self.plot_options)
        R["dev.off"]()

        return odict((('name', track), ('filename', filename)))


class CircularCandidates(GenomePlot):

    track_type = "gene name"
    annotations = "expressed_transcripts.gtf.gz"

    def getTracks(self):

        statement = '''SELECT DISTINCT gi.gene_name
                       FROM   profile_summaries_score as score
                        INNER JOIN annotations.transcript_info as ti
                           ON ti.transcript_id = score.transcript_id
                        INNER JOIN annotations.gene_info as gi
                           ON gi.gene_id = ti.gene_id
                       GROUP BY gi.gene_name
                       ORDER BY max(score.score) DESC
                       LIMIT 10'''

        return self.getValues(statement)

    def populateDataTracks(self):
        
        alyref = glob.glob(os.path.join(PARAMS["iclip_dir"], "bigWig/Alyref*.bw"))
        chtop = glob.glob(os.path.join(PARAMS["iclip_dir"], "bigWig/Chtop*.bw"))
        nxf1 = glob.glob(os.path.join(PARAMS["iclip_dir"], "bigwig/Nxf1*.bw"))

        files = alyref + chtop + nxf1
        files = [f for f in files if "R1" not in f]

        tracks = [re.match("(?:.+/)?(.+)-FLAG.([Ru].+)_(plus|minus).bw", f).groups()
                  for f in files]
        print tracks
        track_names = [factor + " " + replicate if strand == "plus"
                       else factor + "  " + replicate
                       for factor, replicate, strand in tracks]

        strands = ["+" if strand == "plus"
                   else "-"
                   for factor, replicate, strand in tracks]

        DataTracks = {track_name: {"file": file, "strand": strand}
                      for track_name, file, strand in
                      zip(track_names, files, strands)}

        DataTracks["Expression"] =  {"file": "HEK293-WT-1.bam"}
        return DataTracks


    def getDataTracks(self, track):

        gene_strand = self.getValue('''SELECT DISTINCT strand 
                                         FROM annotations.exon_stats as es
                                         INNER JOIN annotations.transcript_info as ti
                                            ON ti.transcript_id = es.transcript_id
                                         INNER JOIN annotations.gene_info as gi
                                            ON gi.gene_id = ti.gene_id
                                         WHERE gi.gene_name  == '%s' ''' % track)

        tracks = [(t,self.DataTracks[t]["obj"]) for
                  t in self.DataTracks if
                  "strand" not in self.DataTracks[t] or 
                  self.DataTracks[t]["strand"] == gene_strand]

        tracks = [t for name, t in sorted(tracks)]

        return tracks


class CircularCandidateTable(ProjectTracker):

    def __call__(self, track):

        statement = ''' SELECT DISTINCT gi.gene_name as Name,
                                        score.transcript_id as "Transcript Id",
                                        score.score as Score
                       FROM   profile_summaries_score as score
                        INNER JOIN annotations.transcript_info as ti
                           ON ti.transcript_id = score.transcript_id
                        INNER JOIN annotations.gene_info as gi
                           ON gi.gene_id = ti.gene_id 
                       ORDER BY score.score DESC LIMIT 20
                       '''

        return self.getAll(statement)

class SingleExonProfiles(ProjectTracker):

    def __call__(self, track):

        statement = ''' SELECT DISTINCT profiles.*, es.nval
                     FROM profile_summaries as profiles
                      INNER JOIN annotations.exon_stats as es
                      ON profiles.transcript_id = es.transcript_id '''

        return self.getDataFrame(statement)
