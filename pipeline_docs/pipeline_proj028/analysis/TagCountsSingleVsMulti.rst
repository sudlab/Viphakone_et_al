Tag counts on single vs multiexon genes
========================================


.. report:: Expression.SingleVsMultiExon
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=category, y=log2(ratio), col=slice, group=slice) + geom_point(position=position_dodge(width=0.3)) + facet_wrap(~protein)

   Tag counts per RNA 


There is a big difference depending on the source of the RNA counts. With the total RNA counts there is a large excess of tags on histone genes. This extends accross all the samples, including FlipIn, suggesting that this is an artefact of the RNA counts. There is also a depletion of counts on single exon in Alyref only. This isn't present in the other samples, or at least not so clearly. 

If the normalisation is done by nuclear RNA, then in alyref there appears to be ablsolutely no effect of the category. In Chtop there appears to be an excess of of tags on the Single exon genes and a depletion on the histones. this is mirrored in the FlipIn samples to an extent, although less clearly than for ChTop. In Nxf1 there is a and excess of tags for both Histones and Single exon.
