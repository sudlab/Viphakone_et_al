.. _geneprofiles:

Gene Profiles
===============


In the pilot data we saw that Alyref bound the 5' end of the transcript and Chtop at the 3' end. Does this still hold in the new data?

.. report:: Profiles.GeneProfiles2
   :render: r-ggplot
   :groupby: all
   :statement: aes(bin,area, col=track) + geom_line(alpha=0.8) + geom_vline(xintercept=c(1000,2000,3000), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Exons","Introns","Downstream"), breaks=c(500,1500,2500,3500)) + theme_bw() + facet_grid(slice~.) + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL)

   Metagene profiles


An alternate way of viewing this is as a heatmap


.. report:: Profiles.GeneProfiles2
   :render: r-ggplot
   :groupby: all
   :width: 600
   :height: 600
   :statement: aes(bin,paste(track," ",slice), fill=area) + geom_raster() + geom_vline(xintercept=c(1000,2000,3000), lwd=0.5, lty=2, col="white") +  scale_x_continuous(labels=c("Upstream","Exons","Introns","Downstream"), breaks=c(500,1500,2500,3500)) + theme_bw()  + theme( aspect.ratio = 0.5, legend.position = "none")  + xlab("") + ylab("")  + scale_fill_gradientn(colours=c("black","#56B1F7"))

   Heatmap of metagene profiles


We can see that the patterns that we saw in the pilot data are still present. However, we it should be noted that there is a slightly worrying enrichment of the x-links at the 5' end of genes for replicate 1 of FLAG only.

Breaking transcripts down by length
++++++++++++++++++++++++++++++++++++

It has been observed by others that the bias of Alyref for the 5' end is length dependent. We wondered if this was due to single exon genes being shorter than multi-exon genes. We can't test this by comparing single exon genes to multi-exon genes, because there are so few single exon genes. But what we can do is see if the length bias is still present when only looking at multi-exon genes. To test this I broken down expressed protein coding transcripts into 5 "quintiles" based on length, each the same number of genes and then computed metagene profiles across these sets of transcripts for both multi-exon transcripts and all  transcripts.

For Alyref, it is clear that the 5' bias gets substainally stronger with increasing length of transcript, and that this trend is indepentent of the number of transcripts:

.. report:: GeneProfiles.BinnedExpressionProfiles
   :render: r-ggplot
   :groupby: track
   :tracks: Alyref-FLAG
   :statement: aes(x=bin, y=area, col=factor(quantile, levels=sort(quantile, decreasing=T))) + geom_line() + facet_grid(slice~exon_limit, scale="free_y") + scale_y_continuous(breaks=NULL) + ylab("Relative coverage") + geom_vline(mapping=aes(xintercept=c(250,500)), lty=2,lwd=0.5) + scale_x_continuous(breaks=c(125, 375, 625), labels = c("Upstream","CDS", "Downstream")) + xlab("") + scale_color_manual(values=colorRampPalette(c("#132B43","#56B1F7"))(5), name = "Length\nQuintile") + theme_bw()

   Metagene profiles for Alyref in transcripts binned by length. 


The pattern is even more interesting for Chtop. Here we can see a 3' bias for the first 4 Quintiles. However, the quintile containing the longest genes suddenly inverts this, with these genes showing a strong 5' bais, as well as a smaller 3' bias (note that int he plot below, I've swaped replicates and quintiles, so that panels are qunintiles and color is replicate)


.. report:: GeneProfiles.BinnedExpressionProfiles
   :render: r-ggplot
   :groupby: track
   :tracks: Chtop-FLAG
   :statement: aes(x=bin, y=area, col=slice) + geom_line() + facet_grid(quantile~exon_limit, scale="free_y") + scale_y_continuous(breaks=NULL) + ylab("Relative coverage") + geom_vline(mapping=aes(xintercept=c(250,500)), lty=2,lwd=0.5) + scale_x_continuous(breaks=c(125, 375, 625), labels = c("Upstream","CDS", "Downstream")) + xlab("") + scale_color_manual(values=colorRampPalette(c("#132B43","#56B1F7"))(5), name = "Replicate") + theme_bw()

   Metagene profiles for Alyref in transcripts binned by length. 


Neither Nxf1 nor FlipIn showed any trend with length or number of exons. The complete results, along with plots showing the distribution of transcript lengths and expression levels can be found in the :ref:`profiles-by-quantile` pipeline section. 

