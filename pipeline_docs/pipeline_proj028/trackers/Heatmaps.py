from ProjectTracker import *

class TrackerHeatmaps(TrackerImages):

    glob = "heatmaps/*.compressed.png"

    def __init__(self, *args, **kwargs ):
        Tracker.__init__(self, *args, **kwargs )
 
    def getTracks(self, subset = None ):
        
        filenames = glob.glob( self.glob )
        pattern = re.sub("\*","(.+)",self.glob)
        tracks = [re.match(pattern,fn).groups()[0] for fn in filenames]
        return tracks

    def __call__(self,track, **kwargs):

        fn = re.sub("\*",track,self.glob)

        track = " ".join(track.split("."))

        return odict( ( ('name', track), ( 'filename', fn) ) )


class NormedHeatmaps(TrackerHeatmaps):

     glob = "heatmaps/*.normed.png"


class FirstExonHeatmaps(TrackerHeatmaps):

    glob = "heatmaps/*first_exon.png"

class NormedFirstExonHeatmaps(TrackerHeatmaps):

    glob = "heatmaps/*first_exon.normed.png"

