Cluster Calling
================


.. report:: clusters.ClusterStats
   :render: r-ggplot
   :regex: .+/(.+).bg.gz
   :glob: clusters.dir/*.bg.gz
   :groupby: all
   :statement: aes(track,Fraction_Significant, fill=substr(track,1,2)) + geom_bar(stat="identity") + theme_bw(base_size=16) + theme(axis.text.x=element_text(angle=90), legend.position="none") + scale_y_continuous(labels=function(x) sprintf("%.0f%%",x*100), limits=c(0,1)) + ylab("Percent Significant") + xlab("")

    Percent of calls that are significant at a 1% FDR
