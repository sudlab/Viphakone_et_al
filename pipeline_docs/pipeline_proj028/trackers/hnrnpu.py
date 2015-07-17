from ProjectTracker import *
from Sample_QC import ContextStats

class HnRNPUMeta(ProjectTracker, SQLStatementTracker):

    statement = ''' SELECT track, region_bin, area
                    FROM hnrnpu1_geneprofiles
                    WHERE region='exons' '''


class HnRNPUContext(ProjectTracker, ContextStats):

    method='hnrnpu1'



class HnRNPUWholeTranscriptMeta(ProjectTracker, SQLStatementTracker):

    statement = ''' SELECT track, region_bin, area
                    FROM hnrnpu1_wholetranscript_geneprofiles
                    WHERE region='exons'
                    '''


class KDVsNuclearLinc(ProjectTracker, SQLStatementTracker):

    fields = ('gene_id','localisation')
    statement = ''' SELECT DISTINCT
                           biotypes.gene_id as gene_id,
                           gi.gene_name as symbol,
                           nuc.log2FoldChange as nuclear,
                           nuc.padj as nuclear_padj,
                           (sum(kd.hnrnpu1kd_EstimatedNumReads) + 0.1)/(sum(kd.WT_EstimatedNumReads)+0.1) as KD,
                           biotype,
                           CASE WHEN nuc.log2FoldChange >= 2 THEN 'Nuclear'
                                WHEN nuc.log2FoldChange < 2 THEN 'Cytoplasmic'
                                END as localisation
                    FROM
                           fraction_diff_deseq as nuc
                           INNER JOIN
                           biotypes ON nuc.gene_id = biotypes.gene_id
                           INNER JOIN
                           fuRNAseq as kd
                           on kd.Transcript = biotypes.transcript_id
                           INNER JOIN annotations.gene_info as gi
                           ON gi.gene_id = biotypes.gene_id 
                   WHERE nuc.baseMean > 50 AND
                         (kd.hnrnpu1kd_EstimatedNumReads + kd.WT_EstimatedNumReads) > 100 AND
                         biotype = "lincRNA"
                   GROUP BY nuc.gene_id'''


class KDVsNuclear(ProjectTracker, SQLStatementTracker):

    fields = ('gene_id',)
    statement = ''' SELECT DISTINCT
                           biotypes.gene_id as gene_id,
                           nuc.log2FoldChange as nuclear,
                           nuc.padj as nuclear_padj,
                           (sum(kd.hnrnpu1kd_EstimatedNumReads) + 0.1)/(sum(kd.WT_EstimatedNumReads)+0.1) as KD,
                           biotype
                    FROM
                           fraction_diff_deseq as nuc
                           INNER JOIN
                           biotypes ON nuc.gene_id = biotypes.gene_id
                           INNER JOIN
                           fuRNAseq as kd
                           on kd.Transcript = biotypes.transcript_id 
                   WHERE nuc.baseMean > 50 AND
                         (kd.hnrnpu1kd_EstimatedNumReads + kd.WT_EstimatedNumReads) > 100 
                      
                   GROUP BY nuc.gene_id'''
