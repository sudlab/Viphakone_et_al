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
