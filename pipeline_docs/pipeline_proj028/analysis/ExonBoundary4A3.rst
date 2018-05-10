Exon Boundary metagenes with eIF4A3
====================================

Adding eIF4A3 to the metagenes has been instructive. You can see the results here:


.. report:: GeneProfiles.ExonBoundaryProfilesWith4A3
   :render: r-ggplot
   :groupby: all
   :tracks: Alyref-FLAG,Chtop-FLAG,Nxf1-FLAG,FlipIn-FLAG
   :slices: union
   :statement: aes(x=position, y=density, fill=relevel(factor, ref="eIF4A3-GFP"), col=relevel(factor, ref="eIF4A3-GFP")) + geom_area(alpha=0.5, position="identity") + geom_line() + facet_wrap(~track, scale="free_y") + theme_bw() + scale_fill_discrete(name="Protein") + guides(color=FALSE)

   Metagene profiles at exon boundarys with 4A3 superimposed.

(I've also exprimented with a different way of ploting it. You can see it plotted the old way here:

.. report:: GeneProfiles.ExonBoundaryProfilesWith4A3
   :render: r-ggplot
   :groupby: all
   :tracks: Alyref-FLAG,Chtop-FLAG,Nxf1-FLAG,FlipIn-FLAG
   :slices: union
   :statement: aes(x=position, y=density, col=relevel(factor, ref="eIF4A3-GFP")) + geom_line() + facet_wrap(~track, scale="free_y") + theme_bw() + scale_color_discrete(name="Protein")

   Metagene profiles at exon boundarys with 4A3 superimposed.


Looking at it this way does emphasise how strong and sharp the eIF4A3 peak is compared to the Alyref peak, which appears both less high and wider than 4A3. It is pleasing to see small dips in both Alyref and Chtop biding at exactly the point of maximum eIF4A3 binding.

