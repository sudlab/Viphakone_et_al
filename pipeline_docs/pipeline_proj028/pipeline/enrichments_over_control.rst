Bootstrapped enrichments over control
========================================

First/middle/last exons
------------------------

.. report:: NormalisedProfiles.FirstLastExonBoot
   :render: r-ggplot
   :tracks: r(^[^F].+)
   :groupby: all
   :statement: aes(x=grouping, y=log2(ratio), ymax=log2(q0.975), ymin=log2(q0.025)) + geom_errorbar(width=0.5) + geom_point() + facet_wrap(~track, scale="free_y") + theme_bw() + scale_x_discrete(limits=c("first", "CDS", "intron", "last"), labels=c("first", "middle", "intron", "last"), name="Region") + ylab("Enrichment over control")

   Enrichment of proteins over matched control for first/middle/last exons and introns.


Detained/retained
-----------------

.. report:: NormalisedProfiles.DetainedBoot
   :render: r-ggplot
   :tracks: r(^[^F])
   :groupby: all
   :plot-width: 4
   :plot-height: 6
   :statement: aes(x=grouping, y=log2(ratio), ymax=log2(q0.975), ymin=log2(q0.025), col=track) + geom_errorbar(width=0.5) + geom_point() + facet_wrap(~track, scale="free_y", labeller=labeller(track=function(x) gsub("_.+","",x)), ncol=3) +theme_bw(base_size=9)  + ylab(expression(paste(Log[2]," enrichment over control"))) + scale_x_discrete(limits=c("Intron", "Detained", "Retained", "Exon"), name="") + scale_color_manual(values=c("Alyref_FLAG"="#D55E00", "Chtop_FLAG"="#009E73", "Nxf1_FLAG"="#CC79A7", "eIF4A3_GFP"="#E69F00", "PTB_GFP"="#999999", "BTZ_GFP"="#0072B2", "UPF3B_GFP"="#0072B2", "RNPS1_GFP"="#F0E442"), guide=FALSE) + theme(aspect.ratio=1, strip.background=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

   Ratios of factors to their controls for introns, exons and retained/detained introns.


If we combine Retained and Detained introns:

.. report:: NormalisedProfiles.CombinedDeReBoot
   :render: r-ggplot
   :tracks: r(^[^F])
   :groupby: all
   :statement: aes(x=grouping, y=log2(ratio), ymax=log2(q0.975), ymin=log2(q0.025)) + geom_errorbar(width=0.5) + geom_point() + facet_wrap(~track, scale="free_x") + theme_bw(base_size=9)  + ylab(expression(paste(Log[2]," enrichment over control"))) + scale_x_discrete(limits=c("Intron", "Retained", "Exon"), name="")

   Ratios of factors to their controls for introns, exons and retained/detained introns.
