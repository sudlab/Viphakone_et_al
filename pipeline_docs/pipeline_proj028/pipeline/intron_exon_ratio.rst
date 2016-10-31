Intron Exon Ratios
======================

Ratios are the ratio of counts in the first intron to the counts in the second exon, calculated over all expressed (TPM > 1) transcripts

.. report:: intron_exon_ratio.RatioComparison
   :render: r-ggplot
   :groupby: all
   :statement: aes(y=(intron_counts+0.1)*exon_length/((exon_counts+0.1)*intron_length), x=track) + scale_y_log10() + geom_boxplot(data=subset(rframe, intron_counts+exon_counts>0)) + xlab("Factor") + ylab("Ratio of densities") + theme_bw()

   Distributions of the ratios between the density of tags in firt introns and second exons


.. report:: intron_exon_ratio.RatioComparison
   :render: r-ggplot
   :groupby: all
   :statement: aes(y=intron_counts+0.1, x=(exon_counts*intron_length/exon_length) + 0.1, col = track) + stat_smooth() +  scale_x_log10() + scale_y_log10() + xlab("Predicted tag count (+0.1)") + ylab("Observed tag count (+0.1)")

   Predicted tag count vs observed assuming same density in intron as exon. 


.. report:: intron_exon_ratio.RatioComparison
   :render: r-ggplot
   :groupby: all
   :statement: aes(y=as.numeric(intron_counts>0), x=(exon_counts*intron_length/exon_length), col = track) + stat_smooth() +  scale_x_log10() + xlab("Predicted tag density") + ylab("Probability of observing at least one tag")

   Predicted tag density vs probability of observing at least one tag. 
