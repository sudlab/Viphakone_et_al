from ProjectTracker import *


class ExpressionCor(ProjectTracker):

    table = "track_counts"
    col_pattern = ".+_FLAG_R.+"
    expression_col = "HEK293_WT_1"

    def getTracks(self):

        tracks = self.getColumns(self.table)
        tracks = [track for track in tracks
                  if re.match(self.col_pattern, track)]
        return tracks


    def __call__(self, track):

        factor, replicate = re.match("(.+)_.+_(.+)", track).groups()
        statement=''' SELECT Geneid as gene_id,
                             %(expression_col)s as expression,
                             %(track)s as clip,
                             '%(factor)s' as protein,
                             '%(replicate)s' as replicate
                             
                      FROM %(table)s '''

        return self.getDataFrame(statement)


    
class ExpressionCorStats(ExpressionCor):

    def __call__(self,track):

        result = ExpressionCor.__call__(self, track)

        result_grouped = result.reset_index().groupby(["protein","replicate"])

        def col_corr(x):
            return x.loc[:,["expression","clip"]].corr(method="spearman").loc[["expression"],["clip"]]

        cors = result_grouped.apply(col_corr).reset_index()

        return cors

class ProbOfClipByExpression(ExpressionCor):

    col_pattern = ".+_FLAG_[R.+|union]"

    def __call__(self,track):

        results = ExpressionCor.__call__(self, track)
        results["clip"] = results["clip"] > 0
        results_gb = results["clip"].groupby(by=pandas.qcut(results.expression.rank(method="first"), 10, labels=range(10)))
        fractions = results_gb.aggregate(lambda x: float(sum(x))/len(x))
        fractions = fractions.reset_index()
        return fractions


class ProbOfClipTPM(ProbOfClipByExpression):

    def getTracks(self):
        tracks = self.getColumns("track_counts")
        tracks = [track for track in tracks
                  if re.match(self.col_pattern, track)]
        return tracks

    table = '''track_counts as counts
             INNER JOIN annotations.transcript_info as gi
             ON gi.gene_id = counts.Geneid
             INNER JOIN HEK293_sailfish as sf
             ON gi.transcript_id = sf.transcript_id
             GROUP BY Geneid'''

    expression_col = "SUM(TPM)"


class TranscriptChunks(ProjectTracker):


    table = "chunk_counts"
    data = None

    def fillData(self, *args, **kwargs):
	print "query started"
	self.data= self.getDataFrame('''SELECT counts.*,
                              intron > 0 as Constitutive_intron,
                              exon >0 as Constitutive_exon
                       FROM %(table)s as counts 
                           INNER JOIN
                            reference_chunks_introns as introns
                           ON introns.gene_id = counts.gene_id 
                           AND introns.exon_id = counts.exon_id
                           INNER JOIN
                            reference_chunks_constitive_exons as exons
                           ON exons.gene_id = counts.gene_id  
                           AND exons.exon_id = counts.exon_id
                       WHERE (exon > 0 OR intron >0) AND NOT (exon >0 AND intron > 0)
                       ''')
        print "query finished"

        print "normalise rna"
        all_rna_cols = [col for col in self.data.columns if "HEK293" in col
                        or "Control_Nuclear" in col]

        def _getQant(s):
            return (s[s>0]).quantile(0.75)

        def _norm(s, lib_size):
            return s*(lib_size/_getQant(s))

        lib_size = self.data[all_rna_cols].apply(_getQant).mean()

        self.data[all_rna_cols] = self.data[all_rna_cols].apply(_norm, lib_size = lib_size)
        

    def getTracks(self):

        columns = self.getColumns(self.table)
        return [col for col in columns if "FLAG" in col]

    slices = ["total", "chromotin", "nuclear"]

    def __call__(self, track, slice):

        if self.data is None:
            self.fillData()

        slice2cols = {"total": [col for col in self.getColumns(self.table)
                                if "banks" in col],
                      "chromotin": ["HEK293_WT_1"],
                      "nuclear": [col for col in self.getColumns(self.table)
                                  if "Control_Nuclear" in col]}

        extra_cols = ["gene_id", "Constitutive_intron", "Constitutive_exon"]

        rna_cols = slice2cols[slice]
        all_cols = [track] + rna_cols + extra_cols

        results = self.data.ix[:, all_cols]

        results["RNA"] = results[rna_cols].mean(axis=1)

        results.rename(columns={track: "iCLIP"}, inplace=True)
        results.drop(rna_cols, axis=1, inplace=True)

        results = results.groupby(
            ["gene_id", "Constitutive_intron", "Constitutive_exon"]
        ).sum().reset_index()

        results = results[(results.iCLIP + results.RNA) > 0]
        return results


class DetainedChunks(TranscriptChunks):

    def fillData(self, *args, **kwargs):

        self.data = self.getDataFrame('''SELECT counts.*,
                              intron > 0 as Constitutive_intron,
                              di.exon > 0 as Constitutive_exon
                       FROM %(table)s as counts
                           INNER JOIN
                            reference_chunks_introns as introns
                           ON introns.gene_id = counts.gene_id
                           AND introns.exon_id = counts.exon_id
                           INNER JOIN
                            reference_chunks_constitive_exons as exons
                           ON exons.gene_id = counts.gene_id
                           AND exons.exon_id = counts.exon_id
                           INNER JOIN
                            reference_chunks_detained as di
                           ON di.gene_id = counts.gene_id 
                           AND di.exon_id = counts.exon_id
                       WHERE intron >0 AND exons.exon = 0
                       ''')
        print "query finished"

        print "normalise rna"
        all_rna_cols = [col for col in self.data.columns if "HEK293" in col]

        def _getQant(s):
            return (s[s>0]).quantile(0.75)

        def _norm(s, lib_size):
            return s*(lib_size/_getQant(s))

        lib_size = self.data[all_rna_cols].apply(_getQant).mean()

        self.data[all_rna_cols] = self.data[all_rna_cols].apply(_norm, lib_size = lib_size)


class LongExpressedGenes(ProjectTracker, SQLStatementTracker):


    statement = '''SELECT DISTINCT gs.gene_id as gene_id,
                          gene_name,
                          (gs.end - gs.start)/1000 as length,
                          500000000 *(Control_Total_R1 + Control_Total_R2) / 
                                 ((SELECT sum(Control_Total_R1 + Control_Total_R2) FROM stubbs_counts)
                                  * gs.sum) as RPKM
                          FROM
                            stubbs_counts as sc
                           INNER JOIN annotations.gene_info as gi on gi.gene_id = sc.geneid
                           INNER JOIN annotations.gene_stats as gs on gi.gene_id = gs.gene_id
                          WHERE length > 50 AND RPKM > 10
                          ORDER BY length DESC
                         
                          LIMIT 100'''

    fields = ('gene_id',)

class LongCDSExpressedGenes(ProjectTracker, SQLStatementTracker):

    statement = '''SELECT DISTINCT gs.gene_id as gene_id,
                          gene_name,
                          gs.sum/1000.0 as length,
                          500000000 *(Control_Total_R1 + Control_Total_R2) / 
                                 ((SELECT sum(Control_Total_R1 + Control_Total_R2) FROM stubbs_counts)
                                  * gs.sum) as RPKM
                          FROM
                            stubbs_counts as sc
                           INNER JOIN annotations.gene_info as gi on gi.gene_id = sc.geneid
                           INNER JOIN annotations.gene_stats as gs on gi.gene_id = gs.gene_id
                          WHERE length > 3.5 AND RPKM > 10
                          ORDER BY RPKM DESC
                         
                          LIMIT 100'''

    fields = ('gene_id',)
