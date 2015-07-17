Expression
===========


Does expression predict the probability of being CLIP'ed?

Expression was divided into bins and the percentage of genes in each expression bin with at least one clip tag was calculated.

Expression was first calculated as the RNA-seq counts: this is a convolution of expression and length:

.. report:: Expression.ProbOfClipByExpression
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=expression+1, y=clip) + geom_point() + facet_wrap(~track) + xlab("Expression decile") + scale_x_continuous(breaks=1:10)

   Effect of RNA-seq count rank on probability of CLIP

Then expression was calculated as sailfish TPMs: this excludes length from the influance

.. report:: Expression.ProbOfClipTPM
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=expression, y=clip) + geom_point() + facet_wrap(~track)+ xlab("Expression decile") + scale_x_continuous(breaks=1:10)

   Effect of RNA-seq TPM rank on probability of CLIP.


