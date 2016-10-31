Comparison with ENCODE CPSF6 and RBM15
=======================================


.. report:: misc.EclipProfiles
   :render: r-ggplot
   :groupby: all
   :statement: aes(bin,area, col=condition) + geom_line(alpha=0.8, lwd=1) + geom_vline(xintercept=c(1000,2000,3000), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Exons","Introns", "Downstream"), breaks=c(500,1500,2500,3500)) + theme_bw() + facet_grid(track + slice~.) + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL) + scale_color_brewer(type="qual", palette="Paired", name = "")

   Metagenes for RBM15 and CPSF6



.. report:: misc.NormedEclipProfiles
   :render: r-ggplot
   :groupby: all
   :statement: aes(bin,area, col=condition) + geom_line(alpha=0.8, lwd=1) + geom_vline(xintercept=c(1000,2000,3000), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Exons","Introns", "Downstream"), breaks=c(500,1500,2500,3500)) + theme_bw() + facet_grid(track + slice~., scale="free_y") + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL) + scale_color_brewer(type="qual", palette="Paired", name = "")

   eCLIP signals normalised to control. 

