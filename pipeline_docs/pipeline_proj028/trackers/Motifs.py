from ProjectTracker import *
import xml.etree.ElementTree

class DremeResults(ProjectTracker):

    def getTracks(self):
        tracks = glob.glob(self.glob_files)
        print "glob returned %i tracks" % len(tracks)
        tracks = [re.match(self.glob_pattern, track).groups()[0] for track in tracks]

        return tracks

    def __call__(self, track, slice=None):

        print "called with track %s" % track

        pattern = re.sub("\(\.\+\)", "%s", self.glob_pattern)
        resultsdir = os.path.abspath(
            os.path.join(pattern % track))

        if not os.path.exists(resultsdir):
            print "resultsdir %s does not exist"
            return []

        tree = xml.etree.ElementTree.ElementTree()
        tree.parse(os.path.join(resultsdir, "dreme.xml"))
        model = tree.find("model")
        num_positives = int(model.find("positives").get("count"))
        num_negatives = int(model.find("negatives").get("count"))

        motifs = tree.find("motifs")
        nmotif = 0
        result = pandas.DataFrame(columns = ["sequence","evalue", "positives",
                                             "negatives","enrichment", "link",
                                             "img"])

        for motif in motifs.getiterator("motif"):
            nmotif += 1
            id = motif.get("id")
            seq = motif.get("seq")

            motif_img = "%(resultsdir)s/%(id)snc_%(seq)s.png" % locals()
            img, rc_img = "na", "na"
            if os.path.exists(motif_img):
                img = '''.. image:: %s
   :scale: 25%% ''' % motif_img

            p = float(motif.get("p"))
            n = float(motif.get("n"))
            
            try:
                enrichment = (p/num_positives)/(n/num_negatives)
                if enrichment < self.enrichment_threshold:
                    continue
                enrichment = "{:.0%}".format(
                    enrichment)
            except ZeroDivisionError:
                enrichment = "inf"

            result = result.append(dict((
                ("sequence", seq),
                ("evalue", motif.get("evalue")),
                ("positives", "{:d}/{} ({:.0%})".format(int(p), num_positives,
                                                        p/num_positives)),
                ("negatives", "{:d}/{} ({:.0%})".format(int(n), num_negatives,
                                                        n/num_negatives)),
                ("enrichment", enrichment),
                ("link", "`dreme_%s <%s/dreme.html>`_" %
                 (track, resultsdir)),
                ("img", img),
            )), ignore_index=True)

        if result.shape[0] > 0:
            return result


class AllAgainstFlipIn(DremeResults):
    enrichment_threshold = 1.5
    glob_files = "export/dreme/*_vs_FlipIn"
    glob_pattern = "export/dreme/(.+)_vs_FlipIn"
    

class RetainedIntronDreme(DremeResults):
    enrichment_threshold = 1.5
    glob_files = "export/dreme/*.clusters_retained_introns"
    glob_pattern = "export/dreme/(.+).clusters_retained_introns"


class DREMELocations(ProjectTracker):
    table = "motif_locations"
    
    def getTracks(self):
        return self.getValues("SELECT DISTINCT track||'-'||motif from %(table)s")

    def buildUCSCLink(self, location, config):
        return '''http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=%(location)s&hgct_customText=%(config)s''' %locals()

    def addLink(self, name, address):
        return '''`%(name)s <%(address)s>`_''' % locals()

    def __call__(self, track):

        statement = ''' SELECT sequence as name, motif, location, strand FROM motif_locations
                        WHERE track||'-'||motif='%(track)s' '''

        fn = "-".join(track.split("-")[:-1])

        query_results = self.getDataFrame(statement)
        if query_results.shape[0] == 0:
            return

        project_id = P.getProjectId()
        prefix = PARAMS["report_prefix"]
        track_config="https://www.cgat.org/downloads/%(project_id)s/%(prefix)sexport/dreme/%(fn)s_vs_FlipIn/UCSC.txt" % locals()
        
        query_results.location = query_results.location.apply(
            lambda x: self.addLink(x, self.buildUCSCLink(x, track_config))
        )

        query_results.strand = query_results.strand.apply(lambda x: x)

        return query_results
