Figures for publication
========================

Figure 1
----------

.. report:: misc.NoFlipInContext
   :render: r-ggplot
   :groupby: all
   :width: 600
   :height: 600
   :statement: aes(x=track, y=alignments, fill=category) + geom_bar(position="fill", stat="identity") + coord_flip() + theme_bw() + scale_fill_brewer(type="qual", palette="Paired", name="") + theme(aspect.ratio=0.5) + scale_y_continuous(labels=function(x) paste(format(x*100, digits=2),"%", sep=""), name = "Fraction of Reads") + xlab("") + scale_x_discrete(limits=c("Nxf1","Chtop","Alyref"))

   Figure 1a: Mapping contexts of clusters found in two samples but not in FlipIn


.. report:: misc.GeneFractions
   :render: r-ggplot
   :transform: melt
   :tracks: Alyref,Chtop,Nxf1
   :groupby: all
   :statement: aes(paste(track,variable),value, fill = track) + geom_bar(stat="identity") + scale_y_continuous(labels=function(x) paste(format(x*100, digits=2),"%", sep=""), limit = c(0,1), name = "Percentage of Genes with crosslink") + theme_bw() + xlab("") + theme(axis.text.x=element_text(angle=90), legend.position="none") + scale_fill_brewer("Paried", type = "qual")

   Figure 1b: Fraction of expressed protein coding transcriptome with at least one crosslinked base.


.. report:: misc.GeneFractions
   :render: table
   :transform: melt
   :tracks: Alyref,Chtop,Nxf1
   :groupby: all

   Above as table


.. report:: misc.ClippedFractionOverlaps
   :render: venn-plot
   :transform: venn

     
   Figure 1b: Fraction of expressed protein coding transcriptome with at least one crosslinked base.

 
.. report:: misc.ClippedFractionOverlaps
   :render: table
   :transform: venn

     
   Figure 1b: Fraction of expressed protein coding transcriptome with at least one crosslinked base.
 

.. report:: misc.GeneProfiles3
   :render: r-ggplot
   :tracks: Nxf1,Chtop,Alyref
   :groupby: all
   :statement: aes(bin,area, col=track) + geom_line(alpha=0.8, lwd=1) + geom_vline(xintercept=c(1000,2000), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Exons","Downstream"), breaks=c(500,1500,2500)) + theme_bw() + facet_grid(slice~.) + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL) + scale_color_brewer(type="qual", palette="Paired", name = "")

   Figure 1b: Meta-gene plot over protein coding genes



.. report:: misc.GeneProfiles3
   :render: r-ggplot
   :tracks: Nxf1,Chtop,Alyref
   :groupby: all
   :statement: aes(bin,paste(track," ",slice), fill=area) + geom_raster() + geom_vline(xintercept=c(1000,2000), lwd=0.5, lty=2, col="white") +  scale_x_continuous(labels=c("Upstream","Exons","Downstream"), breaks=c(500,1500,2500)) + theme_bw()  + theme( aspect.ratio = 0.5, legend.position = "none")  + xlab("") + ylab("")  + scale_fill_gradientn(colours=c("black","#56B1F7"))

   Figure 1b: Heatmap plot over protein coding genes


Figure S1
-----------

.. report:: misc.CustomisedContextStats
   :render: r-ggplot
   :tracks: Alyref,Chtop,Nxf1
   :groupby: all
   :width: 600
   :height: 600
   :statement: aes(x=paste(track,slice),  y=alignments, fill=category) + geom_bar(position="fill", stat="identity") + coord_flip() + theme_bw() + scale_fill_brewer(type="qual", palette="Paired", name="") + theme(aspect.ratio=0.5) + scale_y_continuous(labels=function(x) paste(format(x*100, digits=2),"%", sep=""), name = "Fraction of Reads") + xlab("")

   Figure 1a: Mapping contexts of raw reads 
 
