Metagene profiles at Exon boundaries
=====================================

Only boundaries with 100bp either side of the boundary were included in the analysis. This removed a large number
of boundaries from the analysis. Profiles were normalized by sum for each boundary and then the whole profile normalized
by sum after averaging. The final displayed profiles have been smoothed with a 3bp rolling window.

.. report:: GeneProfiles.ExonBoundaryProfiles
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=position, y=density, col=slice) + geom_line() + facet_wrap(~track) + theme_bw()

   Metagene profiles at exon boundarys

.. report:: GeneProfiles.ExonBoundaryProfiles
   :render: r-ggplot
   :groupby: all
   :slices: union
   :tracks: r(FLAG)
   :plot-width: 6
   :plot-height: 1.5
   :statement: aes(x=position, y=density, col=track) + geom_hline(yintercept=0.005, lty=2, alpha=0.5) + geom_vline(xintercept=0, lty=2, alpha=0.5) + geom_line() + facet_wrap(~track, labeller=labeller(track=function(x) {x<-gsub("-.+","",x); return(gsub("FlipIn", "Cntrl", x))}), nrow=1) + theme_bw(base_size=9) + theme(strip.background=element_blank(), aspect.ratio=1) + scale_y_continuous(breaks=NULL, name="Relative read density") + coord_cartesian(ylim=c(0.003,0.0075)) + scale_x_continuous(breaks=c(-100,-50,0,50,100), name="Distance from Exon-Exon junction") + scale_color_manual(values=c("Alyref-FLAG"="#D55E00", "Chtop-FLAG"="#009E73", "Nxf1-FLAG"="#CC79A7", "eIF4A3-GFP"="#E69F00", "PTB-GFP"="#999999", "BTZ-GFP"="#0072B2", "UPF3B-GFP"="#0072B2", "RNPS1-GFP"="#F0E442", "FlipIn-FLAG"="#56B4E9"), guide=FALSE) 

   Metagene profiles at exon boundarys


.. report:: GeneProfiles.ExonBoundaryProfiles
   :render: r-ggplot
   :groupby: all
   :slices: union
   :tracks: eIF4A3-GFP,BTZ-GFP
   :plot-width: 3
   :plot-height: 2
   :statement: aes(x=position, y=density, col=track) +  geom_hline(yintercept=0.005, lty=2, alpha=0.5) + geom_vline(xintercept=0, lty=2, alpha=0.5) + geom_line() + facet_wrap(~track, labeller=labeller(track=function(x) {x<-gsub("-.+","",x); return(gsub("FlipIn", "Cntrl", x))}), nrow=1) + theme_bw(base_size=9) + theme(strip.background=element_blank(), aspect.ratio=1) + scale_y_continuous(breaks=NULL, name="Relative read density") + scale_x_continuous(breaks=c(-100,-50,0,50,100), name="Distance from Exon-Exon junction") + scale_color_manual(values=c("Alyref-FLAG"="#D55E00", "Chtop-FLAG"="#009E73", "Nxf1-FLAG"="#CC79A7", "eIF4A3-GFP"="#E69F00", "PTB-GFP"="#999999", "BTZ-GFP"="#0072B2", "UPF3B-GFP"="#0072B2", "RNPS1-GFP"="#F0E442", "FlipIn-FLAG"="#56B4E9"), guide=FALSE)

   Metagene profiles at exon boundarys

   
transcriptome mapping
---------------------------

.. report:: GeneProfiles.TranscriptomeExonBoundaryProfiles
   :render: r-ggplot
   :groupby: all
   :slices: union
   :statement: aes(x=position, y=density) + stat_summary(fun.y="mean", geom = "line") + facet_wrap(~track, scale="free_y") + theme_bw() + geom_vline(xintercept=c(0,-24), lty=2, lwd=0.5)

   As above but mapped to transcriptome.


First exons
------------------


.. report:: GeneProfiles.FirstExonBoundaryProfiles2
   :render: r-ggplot
   :regex: heatmaps/(.+).first_exon.end.matrix.tsv.gz
   :glob: heatmaps/*.first_exon.end.matrix.tsv.gz
   :groupby: all
   :statement: aes(index,density,col=track) + geom_line() + facet_wrap(~track)

   Meta genes at boundaries of first exons.

Boundary profiels with eIF4A3
------------------------------

.. report:: GeneProfiles.ExonBoundaryProfilesWith4A3
   :render: r-ggplot
   :groupby: all
   :tracks: Alyref-FLAG,Chtop-FLAG,Nxf1-FLAG,FlipIn-FLAG
   :slices: union
   :statement: aes(x=position, y=density, fill=relevel(factor, ref="eIF4A3-GFP"), col=relevel(factor, ref="eIF4A3-GFP")) + geom_area(alpha=0.5, position="identity") + geom_line() + facet_wrap(~track, scale="free_y") + theme_bw() + scale_fill_discrete(name="Protein") + guides(color=FALSE)

   Metagene profiles at exon boundarys with 4A3 superimposed.


Boundary profiles using the center of the read
-----------------------------------------------

.. report:: GeneProfiles.CenteredExonBoundaryProfiles
   :render: r-ggplot
   :groupby: all
   :tracks: Alyref-FLAG,Chtop-FLAG,Nxf1-FLAG,FlipIn-FLAG
   :slices: union
   :statement: aes(x=position, y=density, col=factor, fill=factor) + geom_area(alpha=0.5, position="identity") + facet_wrap(~track, scale="free_y") + theme_bw() + scale_fill_discrete(name="Protein") + guides(color=FALSE)

   Metagene profiles at exon boundarys with 4A3 and XL sites at center of read


Transcriptome mapped profiles with 4A3
---------------------------------------


.. report:: GeneProfiles.TranscriptomeExonBoundaryProfilesWith4A3
   :render: r-ggplot
   :groupby: all
   :tracks: Alyref-FLAG,Chtop-FLAG,Nxf1-FLAG,FlipIn-FLAG
   :slices: union
   :statement: aes(x=position, y=density, fill=relevel(factor, ref="eIF4A3-GFP"), col=relevel(factor, ref="eIF4A3-GFP")) + geom_area(alpha=0.5, position="identity") + geom_line() + facet_wrap(~track, scale="free_y") + theme_bw() + scale_fill_discrete(name="Protein") + guides(color=FALSE)	     
	     
   As above but mapped to transcriptome.
