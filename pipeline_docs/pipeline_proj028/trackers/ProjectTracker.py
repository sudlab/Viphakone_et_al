from CGATReport.Tracker import *
import CGATPipelines.Pipeline as P

P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS


class ProjectTracker(TrackerSQL):

    def __init__(self,  backend=None, *args, **kwargs):

        attach = [(os.path.join(PARAMS["iclip_dir"], PARAMS["iclip_database"]),
                   "iclip"),
                  (PARAMS["annotations_database"], "annotations"),
                  (os.path.join("/ifs/projects/proj028/EJC_iCLIP_new", "csvdb"), "ejc"),
                  (PARAMS["external_chtop_apa_db"], "chtop_apa")]

        TrackerSQL.__init__(self, backend, attach=attach, *args, **kwargs)
