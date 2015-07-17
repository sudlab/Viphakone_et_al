Analysis of tRNA clip sites
==============================

The following tRNAs are those that have a significant cluster in two replicates of two factors, but not in the FLAG only pull down. They are ranks by the number of reproducible clusters they contain. Also shown is the average number of tags for each factor.

.. report:: TRNAs.TRNAScores
   :render: table
   :force:

   tRNAs bound by two or more factors

.. _trnaprofiles:

Binding patterns on tRNAs
---------------------------

.. report:: TRNAs.TRNAHeatMap
   :render: r-ggplot
   :layout: column-5
   :display: png,png,40
   :statement: aes(x=base, y=index, fill=value) + geom_raster() + theme_minimal() + scale_y_reverse(expand=c(0,0), breaks=NULL, name=NULL) + ylab("tRNAs (ordered by length)") + xlab("Position") + coord_cartesian(xlim=c(0,120)) + scale_x_continuous(expand=c(0,0)) + geom_line(mapping=aes(x=sum, y=index),size = 0.25, col="white")

   Binding pattern of tRNAs sorted by length


With structures
++++++++++++++++

In order add structure we look at the tRNAs all of a similar length. In this case 73bp:

.. report:: TRNAs.TRNAHeatMapWithStruc
   :render: r-ggplot
   :layout: column-5
   :slices: 73
   :display: png,png,50
   :statement: aes(x=base, y=index, fill=value) + geom_raster() + theme_minimal() + scale_y_reverse(expand=c(0,0), breaks=NULL, name=NULL) + ylab("tRNAs (ordered by length)") + xlab("Position") +  scale_x_continuous(expand=c(0,0), labels=rframe$struc[rframe$index==rframe$index[1]], breaks=rframe$base[rframe$index==rframe$index[1]])

   Binding pattern on 73bp tRNAs with structure


Because they are all now the same length, we can look at averages:

.. report:: TRNAs.TRNAHeatMapWithStruc
   :render: r-ggplot
   :groupby: none
   :layout: column-3
   :display: png,png,50
   :statement: aes(x=base, y=value) + stat_summary(fun.y="mean", geom="line") + theme_bw() + xlab("Position") +  ylab("Density") + scale_x_continuous(expand=c(0,0), labels=rframe$struc[rframe$index==rframe$index[1]], breaks=rframe$base[rframe$index==rframe$index[1]])

   Average density of clip tags accross 73bp tRNAs.


FlipIn subtracted
+++++++++++++++++++

.. report:: TRNAs.FlipInSubtractedAverageStruc
   :render: r-ggplot
   :layout: column-3
   :groupby: none
   :display: png,png,50
   :statement: aes(x=base, y=subtracted) + geom_line() + theme_bw() + scale_x_discrete(labels = rframe$struc)

   Average density of clip tags accross tRNAs with FlipIn density subtracted


