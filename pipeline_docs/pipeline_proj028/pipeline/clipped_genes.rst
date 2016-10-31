Which genes have clip tags Clipped Genes
==========================================


List of unclipped genes
------------------------

.. report:: misc.UnclippedFraction
   :render: table
   :max-cols: 10

   List of "expressed" genes with no clip tags for Alyref, ChTop or Nxf1 in any replicate


Expression of unclipped genes
-------------------------------

.. report:: Expression.ExpressionVsClipped
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=factor(clipped), y=RNA.counts) + geom_boxplot() + scale_y_log10() + theme_bw() + xlab("Gene Clipped?") + ylab("Average RNA counts") + scale_x_discrete(breaks=c("0","1"), labels=c("Clipped", "Not Clipped"))

   Number of reads accross either clipped or unclipped genes


.. report:: Expression.ExpressionVsClipped
   :render: r-ggplot
   :groupby: all
   :statement: aes(y=1-clipped, x=RNA.counts) + geom_smooth() + theme_bw() + scale_x_log10() + xlab("Average gene expression (counts)") + ylab("Probability of Clip tag")

   Probability of a gene being clipped is dependent on its expression level






