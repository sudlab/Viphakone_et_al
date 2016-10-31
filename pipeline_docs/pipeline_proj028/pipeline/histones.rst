Replication dependent histone genes
====================================

Metagene of histone genes
--------------------------

.. report:: GeneProfiles.HistoneMetaGenes
   :render: r-ggplot
   :slices: union
   :groupby: all
   :transform: pandas
   :tf-statement: query('geneset == "histones"')
   :statement: aes(bin+0.5, density) + geom_line() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=16) + ylab("Relative Read Density")+facet_grid(track~., scale="free_y")

   Meta gene plot of the replicate dependent histone genes


Comparison to non-histone single exon genes
-------------------------------------------

.. report:: GeneProfiles.HistoneMetaGenes
   :render: r-ggplot
   :groupby: all
   :slices: union
   :statement: aes(bin, density, col=geneset) + geom_step() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=16) + ylab("Relative Read Density")+facet_grid(track~geneset, scale="free_y") + theme(legend.position="none", aspect.ratio=0.5)

   Metagenes of histone and non-histone single exon genes for union tracks


.. report:: GeneProfiles.HistoneMetaGenes
   :render: r-ggplot
   :groupby: all
   :slices: r(R)
   :statement: aes(bin, density, col=replicate) + geom_step() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=16) + ylab("Relative Read Density")+facet_grid(track~geneset, scale="free_y") + theme(legend.position="none", aspect.ratio=0.5)

   Metagenes of histone and non-histone single exon genes for indevidual replicates



