Standard plots requested on EJC
==================================


Context maps
------------

.. report:: misc.EJCContext
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=track,  y=alignments, fill=category) + geom_bar(position="fill", stat="identity") + coord_flip() + theme_bw() + scale_fill_brewer(type="qual", palette="Paired", name="") + theme(aspect.ratio=0.5) + scale_y_continuous(labels=function(x) paste(format(x*100, digits=2),"%", sep=""), name = "Fraction of Reads") + xlab("")

   Mapping contexts of the EJC components


.. report:: misc.EJCContext
   :render: r-ggplot
   :tracks: eIF4A3,BTZ
   :groupby: all
   :plot-width: 4
   :plot-height: 2
   :statement: aes(x=track,  y=alignments, fill=category) + geom_bar(position="fill", stat="identity") + coord_flip() + theme_bw(base_size=9) + scale_fill_brewer(type="qual", palette="Paired", name="") + theme(aspect.ratio=0.2) + scale_y_continuous(labels=function(x) paste(format(x*100, digits=2),"%", sep=""), name = "Fraction of Reads") + xlab("")

    As above but eIF4A3 only

    
Metagene Profiles
------------------

.. report:: misc.EJCGeneProfiles
   :render: r-ggplot
   :groupby: all
   :statement: aes(bin,area, col=track) + geom_line(alpha=0.8, lwd=1) + geom_vline(xintercept=c(1000,2000), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Exons","Downstream"), breaks=c(500,1500,2500)) + theme_bw() + facet_grid(slice~.) + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL) + scale_color_brewer(type="qual", palette="Paired", name = "")

   Metagene profiles over EJC components

 
.. report:: misc.EJCGeneProfiles
   :render: r-ggplot
   :groupby: all
   :tracks: eIF4A3
   :statement: aes(bin,area, col=track) + geom_line(alpha=0.8, lwd=1) + geom_vline(xintercept=c(1000,2000), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Exons","Downstream"), breaks=c(500,1500,2500)) + theme_bw() + facet_grid(slice~.) + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL) + scale_color_brewer(type="qual", palette="Paired", name = "")

   As above but eIF4A3 only


Single exon gene meta plots
---------------------------


.. report:: GeneProfiles.HistoneMetaGenes
   :render: r-ggplot
   :tracks: Alyref-FLAG,eIF4A3-GFP,Chtop-FLAG,Nxf1-FLAG,FlipIn-FLAG
   :groupby: all
   :slices: union
   :statement: aes(bin, density, col=geneset) + geom_step() + geom_vline(xintercept=c(10,20), lty=2) + scale_x_continuous(breaks=c(5,15,25), labels = c("upstream", "exons", "downstream"), name="") + theme_bw(base_size=16) + ylab("Relative Read Density")+facet_grid(track~geneset, scale="free_y") + theme(legend.position="none", aspect.ratio=0.5)

   Metagene plots over histone and non-histone single exon genes


Interesting here is that eIF4A3 and Alyref don't look much like each other, either on histone or non-histone single exon genes.



