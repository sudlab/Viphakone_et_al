from ProjectTracker import ProjectTracker
from CGATReport.Tracker import SQLStatementTracker
import pandas as pd

class SinglePolyAMetagenes(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM single_polyA_window_metagenes'''
    

class NormedAPAMetagenes(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM apa_metagenes where normed=='normed' '''
    
class UnNormedAPAMetagenes(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM apa_metagenes where normed=="unnormed" '''

class AseqVsDapars(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT DISTINCT 'all' as track,
                          gene_name,
                          dapars.gene_id,
                          dapars.transcript_id as transcript_id, 
                          Cntrl, 
                          siCPSF6,
                          dapars.adjusted_P_Val as padj,
                          dapars.PDUI_Group_diff
                    FROM aseq_clusters_prox_distil as aseq
                      INNER JOIN chtop_apa.dapars as dapars
                        ON dapars.transcript_id == aseq.transcript_id
                      INNER JOIN annotations.gene_info as gi
                        ON dapars.gene_id == gi.gene_id '''

    
class AseqVsDaparsFishers(AseqVsDapars):

    def __call__(self):
        df = pd.DataFrame(AseqVsDapars.__call__(self))

        tab = pd.crosstab((df.padj < 0.05) & (df.PDUI_Group_diff < -0.25),
                          (df.Cntrl - df.siCPSF6) < 0)

        from scipy.stats import fisher_exact

        OR, pvalue = fisher_exact(tab)
        return {'track': 'all', 'Odds Ratio': OR, 'pvalue': pvalue}
    

class CPSF_around_clusters(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_cluster
                   WHERE replicate = 'union'  '''
    
class CPSF_around_sites(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_sites
                   WHERE replicate = 'union'  '''


class CPSF_around_apa_sites(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_apa_sites
                   WHERE replicate = 'union'  '''

class CPSF_around_apa_clusters(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_apa_clusters
                   WHERE replicate = 'union'  '''

class CPSF_around_apa_sites_normed(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_apa_sites_normed
                   WHERE replicate = 'union'  '''

class CPSF_around_sites_normed(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_sites_normed
                   WHERE replicate = 'union'  '''

class CPSF_around_clusters_normed(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_clusters_normed
                   WHERE replicate = 'union'  '''

class CPSF_around_apa_clusters_normed(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT * FROM cpsf_metagenes_around_apa_clusters_normed
                   WHERE replicate = 'union'  '''

class CPSF_around_clusters_normed_rand(ProjectTracker, SQLStatementTracker):

    def __call__(self, track=None):

        statement1 = '''SELECT *
                        FROM cpsf_metagenes_around_clusters_normed
                        WHERE replicate = 'union' '''

        df1 = self.getDataFrame(statement1)
        df1["rand"] = "Observed"

        statement2 = '''SELECT * 
                        FROM cpsf_metagenes_around_clusters_normed_rand
                        WHERE replicate = 'union' '''
        df2 = self.getDataFrame(statement2)
        df2["rand"] = "rand"

        return pd.concat([df1,df2])
    

