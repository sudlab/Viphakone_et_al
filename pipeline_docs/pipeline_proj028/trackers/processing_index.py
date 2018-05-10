from ProjectTracker import *

class ProcessingIndex(ProjectTracker, SQLStatementTracker):

    statement = ''' SELECT factor, replicate as track, processing_index
                    FROM processing_index'''

    fields = ["track"]
