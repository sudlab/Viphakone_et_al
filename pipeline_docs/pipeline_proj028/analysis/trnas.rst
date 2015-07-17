
tRNAs
======


In order examine weather or not TREX complex components bound at particular positions on tRNA genes I looked at the tag densities accross each tRNA in each sample. Densities were ordered by length and plotted as heatmaps. Densities were normalised per tRNA, with a pseduo count to allow for unbound tRNAs, and down weight tRNAs were there was only a single tag. i.e.

.. math:: ntags_{i,j} = \frac{tags_{i,j}}{\sum\limits_{j=1}^{l_i} tags_{i,j} + 1}

The resultuant plots look like the below example for replicate 2:

.. report:: TRNAs.TRNAHeatMap
   :render: r-ggplot
   :layout: column-4
   :tracks: r(R2)
   :display: png,png,50
   :statement: aes(x=base, y=index, fill=value) + geom_raster() + theme_minimal() + scale_y_reverse(expand=c(0,0), breaks=NULL, name=NULL) + ylab("tRNAs (ordered by length)") + xlab("Position") + coord_cartesian(xlim=c(0,120)) + scale_x_continuous(expand=c(0,0)) + geom_line(mapping=aes(x=sum, y=index),size = 0.25, col="white") + scale_fill_continuous(limits = c(0,1))

   Binding pattern of tRNAs sorted by length (white line indicates the length of the tRNA)


The most striking feature of this plot is the strong line on the 3' end of the tRNA about 20bp away from the end. Unfortunately this is perhaps due to the fact that the tagged base is the base 5' of the mapped read, and very short reads are discarded. Thus there is a minimal distance from the 3' end it is possible to find a tag. Secondly, this pattern is also present, albeit less strongly, in the FlipIn samples.

There are several other features of these plots. Indeed there appear to be maybe three strong lines in the plot. It is unclear if all of these lines are present in the FlipIn sample. It is also not clear how these lines match with the single and double stranded portions of the tRNA. 

The full set of plots is available at :ref:`trnaprofiles`.

