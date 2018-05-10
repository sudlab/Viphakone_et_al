Differential splicing effcicency of detained and retained introns
=================================================================


A previous analysis of the differential splicing of detained and retained introns suggested that knock out of Alyref could lead to an increase in splicing efficiency of detained introns. However, that analysis used a set of detained introns defined by Phil Sharp and Co in other cell types.

To investigate this further I have defined two sets of introns:

**Detained Introns**: These are defined using a method similar to that used in the Sharp paper, but using HEK293T RNAseq data (polyA). This is repurposed from the controls for the methylation knockdown project and has both nuclear and cytoplasmic samples. Sharp used a variety of different cutoffs when defining his detained introns. I've just a pretty conservative cut off here (FDR 1% and four times as much signal as the expected). I did this in both nuclear and cytoplasmic samples and the detained set is any intron that is signficant in either of these data sets.

**Retained Introns:** These are introns that for which there both spliced and unspliced isofroms in the reference dataset. They may or may not be expressed in HEK293T cells.

There are far more retained introns than there are detained introns.

I then did differential splicing analysis on all exons and introns (not just the detained/retained sets), and overlapped these results with our introns of interest afterwards.

.. report:: RetainedIntrons.CountDiffChunks
   :render: table
   :force:

   Differentially spliced genes at a 5% FDR


I then divided the retained and retained introns sets into bins based on the density of Alyref CLIP tags in them (tags/bp) and overlapped them with the detained/retained set.

Detained Introns
-----------------


Looking at the detained introns first, we can see that the pattern of increasing fractions of the introns being differentially spliced with increased Alyref usage, peaking in the middle, hold for the new detained intron set.

.. report:: RetainedIntrons.DetainedChunkSplicing
   :render: r-ggplot
   :groupby: all
   :statement: aes(cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=20)) + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold<0), col=T),fun.y=mean, geom="point") + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold < -0.58), col=F),fun.y=mean, geom="point")+ facet_grid(slice~track, scale="free_y") +theme_bw() + scale_color_manual(values=c("FALSE"="red", "TRUE"="black"), labels=c("FALSE"="FC>1.5", "TRUE"="FC>0"), name="Fold Change") + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Fraction significant increase in splicing efficiency")

    Fraction of introns in detained category that show increased splicingon Alyref Knockdown


In the above we have the detained introns on the right, and they are compared to all other transcript regions on the right, and then shown for regions that are differentially spliced in the cytoplasmic RNAseq (top), nuclear RNAseq (middle) and total RNAseq (bottom). This plot only shows those introns that are at lower levels in the knockout - i.e. those where splicing efficiency has increased. The biggest effects are in the nuclear RNAseq where for some levels of Alyref binding, over 25% of all regions are marked as differentially spliced. This is a much greater fraction than for non-detained intron regions.

We can see this more effectively by looking at enrichment of detained introns over other regions rather than just the raw numbers:


 .. report:: RetainedIntrons.DetainedEnrichment
   :render: r-ggplot
   :statement: aes(x=clip_bin, y=enrichment, col=p<0.05/20) + geom_vline(xintercept=0, lwd=8, col="grey75", alpha=0.5) + geom_point() + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Enrichment") + scale_color_manual(name="", labels=c("TRUE"="Significant", "FALSE"="Not Significant"), values=c("TRUE"="red", "FALSE"="black")) + facet_grid(slice~.) + theme_bw()

   Enrichment of detained introns in the differentially spliced set.


In the above plot red points mark those bins where the enrichment is significant at p<0.05. Note that the relationship between enrichment and amount of Alyref is not simple. In particular, the grey highlighted region is where there is no Alyref on the introns, yet there is a large and significant enrichment in this bin. It is difficult to know how to interpret this, because the 0 Alyref clip tags could be due to a low affinity for Alyref, or due to the introns being short. Apart from this though, it appears that even if the fraction of regions changed is fewer at higher Alyref levels, the enrichment is still larger.


Retained Introns
-----------------


The same patterns hold if we consider Retained introns rather than detained.

.. report:: RetainedIntrons.RetainedChunkSplicing
   :render: r-ggplot
   :groupby: all
   :slices: nuclear
   :statement: aes(cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=20)) + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold<0), col=T),fun.y=mean, geom="point") + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold < -0.58), col=F),fun.y=mean, geom="point")+ facet_grid(slice~track, scale="free_y") +theme_bw() + scale_color_manual(values=c("FALSE"="red", "TRUE"="black"), labels=c("FALSE"="FC>1.5", "TRUE"="FC>0"), name="Fold Change") + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Fraction significant increase in splicing efficiency") + theme(aspect.ratio=0.5)

   Fraction of introns in retained category that show increased splicing in the nucleus on Alyref Knockdown


and the enrichment:

 .. report:: RetainedIntrons.RetainedEnrichment
   :render: r-ggplot
   :slices: nuclear
   :statement: aes(x=clip_bin, y=enrichment, col=p<0.05/20) + geom_vline(xintercept=0, lwd=8, col="grey75", alpha=0.5) + geom_point() + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Enrichment") + scale_color_manual(name="", labels=c("TRUE"="Significant", "FALSE"="Not Significant"), values=c("TRUE"="red", "FALSE"="black")) + facet_grid(slice~.) + theme_bw() + theme(aspect.ratio=0.5)

    Enrichment of Retained introns in the differentially spliced set on Alyref knockdown.


I suspect that the reason the patterns are less pronounced is that the retained introns are from a reference set that wasn't defined in HEK293 cells exclusively, so there will be some sequence marked as "Retained" that is always spliced out already in HEK293 cells.

Displaying for paper
--------------------

In any paper figure, we need to decide wether to show seperate or combined plots. I'm not sure if the idea of showing introns from the middle section makes sense any more as we can see that the enrichment increases with increase Alyref amount. Personally, I think perhaps the following plot:


 .. report:: RetainedIntrons.CombinedEnrichment
   :render: r-ggplot
   :slices: nuclear
   :statement: aes(x=clip_bin, y=enrichment, col=p<0.05/20) + geom_vline(xintercept=0, lwd=8, col="grey75", alpha=0.5) + geom_point() + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Enrichment") + scale_color_manual(name="", labels=c("TRUE"="Significant", "FALSE"="Not Significant"), values=c("TRUE"="red", "FALSE"="black")) + facet_grid(slice~.) + theme_bw() + theme(aspect.ratio=0.5)

    Enrichment of Retained and Detained introns in the differentially spliced set on Alyref knockdown.

But I am open to other suggestions.


Below are lists of the differentially detained and retained introns in the nucleus and the cytoplasm.

.. report:: RetainedIntrons.DetainedChunkSplicing
   :render: table
   :transform: pandas
   :groupby: all
   :slices: nuclear
   :large: xls
   :tf-statement: dropna().query("(track == \" \")  and (padj < 0.05) and (l2fold < -0.58)", engine="python")[["symbol", "groupID", "contig", "intron_start", "intron_end", "padj", "l2fold"]].sort_values("l2fold")

   Differentially spliced detained introns (nuclear)


 
.. report:: RetainedIntrons.RetainedChunkSplicing
   :render: table
   :transform: pandas
   :groupby: all
   :slices: nuclear
   :large: xls
   :tf-statement: dropna().query("(track == \" \")  and (padj < 0.05) and (l2fold < -0.58)", engine="python")[["symbol", "groupID", "contig", "intron_start", "intron_end", "padj", "l2fold"]].sort_values("l2fold")

    Differentially spliced retained introns (nuclear)


Plots with no breakdown by Alyref quantity
------------------------------------------

It seems that the whole relationship with Alyref density is less clear now. The following plots remove the subdivision by Alyref density, and just show the *number* of altered introns, comparing introns annotated as detained or retained against introns/exons without this annotation. 

.. report:: RetainedIntrons.CombinedChunkSplicing
   :render: r-ggplot
   :groupby: slice
   :layout: row
   :tracks: nuclear,chtop	    
   :plot-width: 3
   :plot-height: 1.5
   :statement: aes(x=retained | detained, y=as.numeric(is.finite(padj) & is.finite(l2fold) & padj < 0.1 & l2fold < -0.58)) + stat_summary(fun.y=mean, geom="bar") + coord_flip() + theme_bw(base_size=10) + theme(aspect.ratio=0.5) + scale_x_discrete(labels=c("Other", "Retained"), name="") + scale_y_continuous(labels=scales::percent, name="% with increased splicing\non knockdown")

   Enrichment of (combined detained) and retained introns


As we can see, there is clearly still an enrichment. Or we could breakdown into only three groups...

.. report:: RetainedIntrons.CombinedChunkSplicing
   :render: r-ggplot
   :groupby: slice
   :layout: row
   :slices: nuclear,chtop
   :plot-width: 3
   :plot-height: 2
   :statement: aes(x=cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=3, labels=c("High", "Medium", "Low")), fill=retained | detained, y=as.numeric(is.finite(padj) & is.finite(l2fold) & padj < 0.1 & l2fold < -0.58)) + stat_summary(fun.y=mean, geom="bar", position="dodge") + theme_bw(base_size=10)  + coord_flip() + scale_y_continuous(labels=scales::percent, name="% regions with increased splicing\non knockdown") + xlab("Alyref Tag Density") + theme(legend.position="bottom", aspect.ratio=0.5, legend.key.size=unit(0.8, "lines"), legend.margin=margin(0,0,0,0,"cm")) + scale_fill_manual(name=NULL, values=c("grey35", rgb(0,0.45,0.7)), labels=c("Other", "Retained"))

   Enrichment of (combined detained and) retained introns in changes on TREX knockdown.


Or only bound and not bound.

.. report:: RetainedIntrons.CombinedChunkSplicing
   :render: r-ggplot
   :groupby: slice
   :layout: row
   :slices: nuclear,chtop
   :plot-width: 3
   :plot-height: 2
   :statement: aes(x=clip_tags>0, fill=retained | detained, y=as.numeric(is.finite(padj) & is.finite(l2fold) & padj < 0.1 & l2fold < -0.58)) + stat_summary(fun.y=mean, geom="bar", position="dodge") + theme_bw(base_size=10)  + coord_flip() + scale_y_continuous(labels=scales::percent, name="% regions with increased splicing\non knockdown") + xlab("Alyref Tag Density") + theme(legend.position="bottom", aspect.ratio=0.5, legend.key.size=unit(0.8, "lines"), legend.margin=margin(0,0,0,0,"cm")) + scale_fill_manual(name=NULL, values=c("grey35", rgb(0,0.45,0.7)), labels=c("Other", "Retained")) + scale_x_discrete(labels=c("Not Bound","Alyref Bound"), name=NULL)

   Enrichment of (combined detained and) retained introns in changes on TREX knockdown.
   

splitting out the retained and detained....

.. report:: RetainedIntrons.DetainedChunkSplicing
   :render: r-ggplot
   :groupby: slice
   :layout: row
   :slices: nuclear,chtop
   :plot-width: 3
   :plot-height: 2
   :statement: aes(x=cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=3, labels=c("Low", "Medium", "High")), fill=detained, y=as.numeric(is.finite(padj) & is.finite(l2fold) & padj < 0.1 & l2fold < -0.58)) + stat_summary(fun.y=mean, geom="bar", position="dodge") + theme_bw(base_size=10)  + coord_flip() + scale_y_continuous(labels=scales::percent, name="% regions with increased splicing\non knockdown") + xlab("Alyref Tag Density") + theme(legend.position="bottom", aspect.ratio=0.5, legend.key.size=unit(0.8, "lines"), legend.margin=margin(0,0,0,0,"cm")) + scale_fill_manual(name=NULL, values=c("grey35", rgb(0,0.45,0.7)), labels=c("Other", "Detained"))
 
   Enrichment of detained introns in changes on TREX knockdown.


Or only bound and not bound.

.. report:: RetainedIntrons.RetainedChunkSplicing
   :render: r-ggplot
   :groupby: slice
   :layout: row
   :slices: nuclear,chtop
   :plot-width: 3
   :plot-height: 2
   :statement: aes(x=cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=3, labels=c("Low", "Medium", "High")), fill=retained, y=as.numeric(is.finite(padj) & is.finite(l2fold) & padj < 0.1 & l2fold < -0.58)) + stat_summary(fun.y=mean, geom="bar", position=position_dodge()) + theme_bw(base_size=10)  + coord_flip() + scale_y_continuous(labels=scales::percent, name="% regions with increased splicing\non knockdown") + xlab("Alyref Tag Density") + theme(legend.position="bottom", aspect.ratio=0.5, legend.key.size=unit(0.8, "lines"), legend.margin=margin(0,0,0,0,"cm")) + scale_fill_manual(name=NULL, values=c("grey35", rgb(0,0.45,0.7)), labels=c("Other", "Retained"))

   Enrichment of retained introns in changes on TREX knockdown.
   
