.. _rnaseq:

Correlation with RNAseq data
=============================


As we have discussed previously, it doesn't make much sense to ask if reads are enriched in various different types of sequence (protein coding, lincRNA etc) based on the genomic size of these regions. However, simply asking for the fraction of reads that map to each type isn't very good either. 

One way to ask if the distribution of reads across samples is to corrolation the counts of iCLIP reads with counts from a gene expression experiment. Luckily we have just such and experiment, and can look at the corrolations:


.. report:: Expression.ExpressionCor
   :render: r-ggplot
   :groupby: all
   :statement: aes(expression+1, clip+1) + geom_point(alpha=0.1, size=2) + stat_smooth(method="lm") + scale_x_log10() + scale_y_log10() + facet_grid(protein~replicate) + xlab("Gene expression") +  ylab("Clip counts") + theme_bw()

   Relationship between expression and clip tags

For Alyref and Chtop there is a strong relationship at on the log scale between expression and number of clip tags. The relationship is less strong for Nxf1 and for FLAG only. This looks to be due to a higher number of transcript at a high expression level that have few or no clip tags. We can quantify this by looking at the correlation coefficients between expression and iCLIP density:

.. report:: Expression.ExpressionCorStats
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=protein, y=clip, fill=protein, group=replicate) + geom_bar(col = "white", position="dodge", stat="identity") + xlab("Factor") +  ylab(expression(paste("Spearmans ", rho))) + theme_bw(base_size=18) + theme(legend.position="none")

   Correlation between gene expression and CLIP density


In general the correlation coefficients are better for Alyref and Chtop. We are now in a position to ask if some transcript types are bound more or less than would be expected given their expression level.

..TODO::
    * Are the slopes of the correlation different for different categories of transcript
    * If the corrolation with length normalized expression better or worst.






