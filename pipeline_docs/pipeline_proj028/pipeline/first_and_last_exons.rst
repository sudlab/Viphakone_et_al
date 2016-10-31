First and Last exon profiles
=============================

Gene models expressed in HEK293 cells were used. Genes were flattened. total-sum per transcript normalisation was used. 


All the unnormalised profiles
------------------------------

.. report:: NormalisedProfiles.NormalisedFirstLastExons
   :render: r-ggplot
   :groupby: track
   :layout: column-3
   :statement: aes(x=bin, y=area, col=slice) +  facet_grid(slice~.) + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   Unnormalised metagene profiles seperating out the first and last exons of the gene



Normalised to Nuclear RNA
--------------------------

.. report:: NormalisedProfiles.NormalisedFirstLastExons
   :render: r-ggplot
   :groupby: track
   :layout: column-3
   :statement: aes(x=bin, y=normed, col=slice) +  facet_grid(slice~.) + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   Unnormalised metagene profiles seperating out the first and last exons of the gene


EJC components normalised to their FlipIn track
--------------------------------------------------

.. report:: NormalisedProfiles.GFPNormalisedFirstLastExons
   :render: r-ggplot
   :groupby: track
   :tracks: r(GFP)
   :layout: column-3
   :statement: aes(x=bin, y=normed, col=slice) +  facet_grid(slice~.) + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))


Union tracks by cell type
--------------------------


.. report:: NormalisedProfiles.NormalisedFirstLastExons
   :render: r-ggplot
   :groupby: slice
   :slices: union
   :tracks: r(FLAG)
   :layout: column-3
   :statement: aes(x=bin, y=area, col=track) +  facet_grid(track~.) + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   TREX components unnormalised metagene profiles seperating out the first and last exons of the gene


.. report:: NormalisedProfiles.NormalisedFirstLastExons
   :render: r-ggplot
   :groupby: slice
   :slices: union
   :tracks: r(GFP)
   :layout: column-3
   :statement: aes(x=bin, y=area, col=track) +  facet_grid(track~.) + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   EJC components unnormalised metagene profiles seperating out the first and last exons of the gene


.. report:: NormalisedProfiles.NormalisedFirstLastExons
   :render: r-ggplot
   :groupby: slice
   :slices: union
   :tracks: r(FLAG)
   :layout: column-3
   :statement: aes(x=bin, y=normed, col=track) +  facet_grid(track~.) + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   TREX components normalised metagene profiles seperating out the first and last exons of the gene


.. report:: NormalisedProfiles.NormalisedFirstLastExons
   :render: r-ggplot
   :groupby: slice
   :slices: union
   :tracks: r(GFP)
   :layout: column-3
   :statement: aes(x=bin, y=normed, col=track) +  facet_grid(track~., scale="free_y") + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   EJC components normalised metagene profiles seperating out the first and last exons of the gene


.. report:: NormalisedProfiles.GFPNormalisedFirstLastExons
   :render: r-ggplot
   :groupby: slice
   :slices: union
   :tracks: r(GFP)
   :layout: column-3
   :statement: aes(x=bin, y=normed, col=track) +  facet_grid(track~., scale="free_y") + theme_bw() + guides(color=F) + geom_vline( xintercept=c(250, 308, 558, 744), lty=2, col="grey75") + geom_line() + scale_x_continuous(breaks=c(125,279,433,651,869), labels=c("Upstream", "First", "CDS", "Last", "Downstream"))

   EJC components metagene profiles normalised to FlipIn seperating out the first and last exons of the gene


Aggregated, normalised counts
------------------------------


These counts are summed across all genes and then normalised by the counts for FlipIn across all genes. They are then normalised to the counts for the middle exon. 

.. report:: NormalisedProfiles.FirstLastExonCount
   :render: r-ggplot
   :slices: r(R.)
   :layout: row
   :statement: aes(x=exon, y=log2(normed_count), col=slice) + geom_point() + facet_wrap(~protein) + theme_bw() + scale_x_discrete(limits=c("first_exon", "middle_exon", "last_exon"), labels=c("First", "CDS", "Last")) + theme(aspect.ratio=1) + ylab("LogFC compared to middle exons")

   Normalised counts across first, last and middle exons.
