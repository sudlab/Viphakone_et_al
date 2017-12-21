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
