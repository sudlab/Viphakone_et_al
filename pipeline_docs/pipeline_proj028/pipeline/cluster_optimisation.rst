Optimisaiton of cluster calling
===============================

P value thresholds
------------------

Number of clusters called:

.. report:: Clusters.ClustersCalled
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=p, y=count, col = replicate) + facet_wrap(~track, scale="free_y") + geom_line() + geom_point() + theme_bw() + xlab("P value cutoff") + ylab("Number of clusters")

   Number of clusters called at a variety of pvalue thresholds


Fraction of the called clusters that were reproducilbe

.. report:: Clusters.FractionReproducible
   :render: r-ggplot
   :groupby: slice
   :statement: aes(x=p, y=reproducible/all, col = track) + geom_line() + geom_point() + theme_bw() + xlab("P value cutoff") + ylab("Fraction of significant windows reproducilbe")
   
   Fraction of clusters called in any replicate that are reproducilbe in at least two replicates at various p value thresholds



