Bootstraped enrichments for first/middle/last exons and introns
=================================================================

.. report:: NormalisedProfiles.FirstLastExonBoot
   :render: r-ggplot
   :tracks: Alyref_FLAG,Chtop_FLAG,Nxf1_FLAG,eIF4A3_GFP,BTZ_GFP
   :groupby: all
   :plot-width: 8
   :plot-height: 1.5
   :statement: aes(x=grouping, y=log2(ratio), ymax=log2(q0.975), ymin=log2(q0.025), col=track) + geom_errorbar(width=0.5) + geom_point() + facet_wrap(~track, scale="free_y", nrow=1) + theme_bw(base_size=9) + scale_x_discrete(limits=c("first", "CDS", "intron", "last"), labels=c("Frist Exon", "Middle Exons", "Introns", "Last Exons"), name=NULL) + scale_y_continuous(limits=c(0,NA), name=expression(paste(Log[2], " enrichment"))) + scale_color_manual(values=c("Alyref_FLAG"="#D55E00", "Chtop_FLAG"="#009E73", "Nxf1_FLAG"="#CC79A7", "eIF4A3_GFP"="#E69F00", "PTB_GFP"="#999999", "BTZ_GFP"="#0072B2"), guide=FALSE) + theme(aspect.ratio=1, axis.text.x=element_text(angle=45,hjust=1), strip.text.x=element_text(margin=margin(1,1,2,1)), strip.background = element_blank(), panel.spacing=grid::unit(0,"lines"))

   Enrichment of proteins over matched control for first/middle/last exons and introns.

Confirming what we more or less already knew from the metagene plots we can see that Alyref is significantly enrichmed over the first exon compared to the middle/last exons and the opposite is true for ChTop. Intresetingly it seems that Nxf1 has an ever so slight enrichment in the 3' UTR, same with eIF4A3. In fact almost everything does, which is a bit worrying. When you compare this to the metagenes, which don't look enriched at the 3', it look as if this must be due to lower counts of the control, rather than higher counts of the protein of interest.

.. report:: NormalisedProfiles.FirstLastExonBoot
   :render: r-ggplot
   :tracks: eIF4A3_GFP,BTZ_GFP
   :groupby: all
   :plot-width: 3
   :plot-height: 1.5
   :statement: aes(x=grouping, y=log2(ratio), ymax=log2(q0.975), ymin=log2(q0.025), col=track) + geom_errorbar(width=0.5) + geom_point() + facet_wrap(~track, nrow=1) + theme_bw(base_size=9) + scale_x_discrete(limits=c("first", "CDS", "intron", "last"), labels=c("Frist Exon", "Middle Exons", "Introns", "Last Exons"), name=NULL) + scale_y_continuous(limits=c(0,NA), name=expression(paste(Log[2], " enrichment"))) + scale_color_manual(values=c("Alyref_FLAG"="#D55E00", "Chtop_FLAG"="#009E73", "Nxf1_FLAG"="#CC79A7", "eIF4A3_GFP"="#E69F00", "PTB_GFP"="#999999", "BTZ_GFP"="#0072B2"), guide=FALSE) + theme(aspect.ratio=1, axis.text.x=element_text(angle=45,hjust=1), strip.text.x=element_text(margin=margin(1,1,2,1)), strip.background = element_blank(), panel.spacing=grid::unit(0,"lines"))

   Enrichment of proteins over matched control for first/middle/last exons and introns.


Normalised to RNA rather than FlipIn.
-------------------------------------


.. report:: NormalisedProfiles.FirstLastExonBootTotalRNA
   :render: r-ggplot
   :tracks: Alyref_FLAG
   :plot-width: 2
   :plot-height: 2
   :statement: aes(x=grouping, y=ratio, ymax=q0.975, ymin=q0.025) + geom_errorbar(width=0.5) + geom_point() + theme_bw(base_size=9) + scale_x_discrete(limits=c("first", "CDS", "intron", "last"), labels=c("Frist Exon", "Middle Exons", "Introns", "Last Exons"), name=NULL) + scale_y_continuous(limits=c(0,NA), name=expression(paste(Log[2], " enrichment"))) + theme(aspect.ratio=1, axis.text.x=element_text(angle=45,hjust=1), strip.text.x=element_text(margin=margin(1,1,2,1)), strip.background = element_blank(), panel.spacing=grid::unit(0,"lines"))

   Normalised to Total RNA


.. report:: NormalisedProfiles.FirstLastExonBootNucRNA
   :render: r-ggplot
   :tracks: Alyref_FLAG
   :plot-width: 2
   :plot-height: 2
   :statement: aes(x=grouping, y=ratio, ymax=q0.975, ymin=q0.025) + geom_errorbar(width=0.5) + geom_point() + theme_bw(base_size=9) + scale_x_discrete(limits=c("first", "CDS", "intron", "last"), labels=c("Frist Exon", "Middle Exons", "Introns", "Last Exons"), name=NULL) + scale_y_continuous(limits=c(0,NA), name=expression(paste(Log[2], " enrichment"))) + theme(aspect.ratio=1, axis.text.x=element_text(angle=45,hjust=1), strip.text.x=element_text(margin=margin(1,1,2,1)), strip.background = element_blank(), panel.spacing=grid::unit(0,"lines"))

   Normalised to Nuclear RNA


.. report:: NormalisedProfiles.FirstLastExonBootChromRNA
   :render: r-ggplot
   :tracks: Alyref_FLAG
   :plot-width: 2
   :plot-height: 2
   :statement: aes(x=grouping, y=ratio, ymax=q0.975, ymin=q0.025) + geom_errorbar(width=0.5) + geom_point() + theme_bw(base_size=9) + scale_x_discrete(limits=c("first", "CDS", "intron", "last"), labels=c("Frist Exon", "Middle Exons", "Introns", "Last Exons"), name=NULL) + scale_y_continuous(limits=c(0,NA), name=expression(paste(Log[2], " enrichment"))) + theme(aspect.ratio=1, axis.text.x=element_text(angle=45,hjust=1), strip.text.x=element_text(margin=margin(1,1,2,1)), strip.background = element_blank(), panel.spacing=grid::unit(0,"lines"))

   Normalised to Chromatin RNA
