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
   :tracks: Alyref-FLAG,Chtop-FLAG,FlipIn-FLAG
   :plot-width: 3
   :plot-height: 2
   :statement: aes(bin, density, col=track) + geom_step() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=9) +facet_grid(geneset ~track, scale="free_y", labeller=labeller(track=function(x) {x<-gsub("-FLAG","",x); return(gsub("FlipIn", "Cntrl", x))}, geneset=function(x) gsub("_", " ", x))) + theme(legend.position="none", aspect.ratio=0.5, axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous(name="Relative tag density", breaks=NULL) + theme(strip.background=element_blank(), panel.spacing=grid::unit(0,"lines")) + scale_color_manual(values=c("Alyref-FLAG"="#D55E00", "Chtop-FLAG"="#009E73", "eIF4A3-GFP"="#E69F00", "BTZ-GFP"="#0072B2","FlipIn-FLAG"="#56B4E9" ))

   Metagenes of histone and non-histone single exon genes for union tracks


.. report:: GeneProfiles.HistoneMetaGenes
   :render: r-ggplot
   :groupby: all
   :slices: union
   :tracks: eIF4A3-GFP,BTZ-GFP,FlipIn-GFP
   :plot-width: 3
   :plot-height: 2
   :statement: aes(bin, density, col=track) + geom_step() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=9) +facet_grid(geneset ~track, scale="free_y", labeller=labeller(track=function(x) {x<-gsub("-GFP","",x); return(gsub("FlipIn", "Cntrl", x))}, geneset=function(x) gsub("_", " ", x))) + theme(legend.position="none", aspect.ratio=0.5, axis.text.x=element_text(angle=45, hjust=1)) + scale_y_continuous(name="Relative tag density", breaks=NULL) + theme(strip.background=element_blank(), panel.spacing=grid::unit(0,"lines")) + scale_color_manual(values=c("Alyref-FLAG"="#D55E00", "Chtop-FLAG"="#009E73",  "eIF4A3-GFP"="#E69F00",  "BTZ-GFP"="#0072B2", "FlipIn-GFP"="#56B4E9"))

   EJC Metagenes of histone and non-histone single exon genes for union tracks

   
.. report:: GeneProfiles.HistoneMetaGenes
   :render: r-ggplot
   :groupby: all
   :slices: r(R)
   :statement: aes(bin, density, col=replicate) + geom_step() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=16) + ylab("Relative Read Density")+facet_grid(track~geneset, scale="free_y") + theme(legend.position="none", aspect.ratio=0.5)

   Metagenes of histone and non-histone single exon genes for indevidual replicates



