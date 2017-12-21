Retained Introns
=================


.. report:: RetainedIntrons.UpregulatedRetainedIntrons
   :render: venn-plot
   :transform: venn
   

   Overlap between genes regulated with UAP56 and DDX39 are knocked down and genes with retained introns



.. report:: RetainedIntrons.UpregulatedRetainedIntrons
   :render: table
   :transform: hypergeometric
   

   Enrichment of regulated genes in genes with retained introns



.. _retainedintonclusters:
Genes with retain intron clusters
-----------------------------------

.. report:: RetainedIntrons.ClustersInRetainedIntrons
   :render: table
   :tracks: Nxf1_FLAG_reproducible
   :force:

   Genes with Nxf1 clusters in retained introns



Differential Exon Usage
==========================

MA plots
---------

.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :groupby: slice
   :slices: Alyref-FLAG-R1
   :statement: aes(exonBaseMean + 1.5e-1, log2fold_Control_Alyref, col = padj<0.05 & is.finite(padj), alpha=padj<0.05 & is.finite(padj)) + geom_point(size=1) + scale_x_log10() + scale_alpha_manual(limits=c(T,F), values=c(1,0.25), name="significant") + scale_color_manual(limits=c(T,F), values=c("red","black"), name="significant") + geom_smooth(alpha=0.3) + coord_cartesian(ylim=c(-5,5)) + facet_grid(track~.)

   Classic MA plots


.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :groupby: slice
   :slices: Alyref-FLAG-R1
   :statement: aes(Control + 1.5e-1, log2fold_Control_Alyref, col = padj<0.05 & is.finite(padj), alpha=padj<0.05 & is.finite(padj)) + geom_point(size=1) + scale_x_log10() + scale_alpha_manual(limits=c(T,F), values=c(1,0.25)) + scale_color_manual(limits=c(T,F), values=c("red","black")) + geom_smooth(alpha=0.3) + coord_cartesian(ylim=c(-5,5)) + facet_grid(track~.)

   Control expression rather than average expression


.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :groupby: slice
   :slices: Alyref-FLAG-R1
   :statement: aes((exonBaseMean+0.1)/(genomicData_width+0.1), log2fold_Control_Alyref, col = padj<0.05 & is.finite(padj), alpha=padj<0.05 & is.finite(padj)) + geom_point(size=1) + scale_x_log10() + scale_alpha_manual(limits=c(T,F), values=c(1,0.25),name="significant") + scale_color_manual(limits=c(T,F), values=c("red","black"),name="significant" ) + geom_smooth(alpha=0.3) + coord_cartesian(ylim=c(-5,5)) + facet_grid(track~.)

   Read density vs fold change



.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :groupby: slice
   :slices: r(union)
   :layout: column-2
   :statement: aes(clip_tags+1, log2fold_Control_Alyref, col = padj<0.05 & is.finite(padj), alpha=padj<0.05 & is.finite(padj)) + geom_point(size=1) + scale_x_log10() + scale_alpha_manual(limits=c(T,F), values=c(1,0.25), name="significant") + scale_color_manual(limits=c(T,F), values=c("red","black"), name="significant") + geom_smooth(alpha=0.3) + coord_cartesian(ylim=c(-5,5)) + facet_grid(track~.)

   CLIP tags vs fold change


.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :groupby: slice
   :slices: r(union)
   :statement: aes((clip_tags/genomicData_width) + 1e-5, log2fold_Control_Alyref, col = padj<0.05 & is.finite(padj), alpha=padj<0.05 & is.finite(padj)) + geom_point(size=1) + scale_x_log10() + scale_alpha_manual(limits=c(T,F), values=c(1,0.25), name="significant") + scale_color_manual(limits=c(T,F), values=c("red","black"), name="significant") + geom_smooth(alpha=0.3) + coord_cartesian(ylim=c(-5,5)) + facet_grid(track~.)

   Density of clip tags vs fold change



Overlaps
-----------


.. report:: RetainedIntrons.ChangedRIVenn
   :render: venn-plot
   :transform: venn
   

   Overlaps between cell fractions


.. report:: RetainedIntrons.ChangedRIVenn
   :render: table
   :transform: hypergeometric

   
   stats on the above




RNAseq density vs CLIP density vs significance
-----------------------------------------------

.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :tracks: nuclear
   :slices: r(union)
   :statement: aes((Control/genomicData_width) + 1e-5, (clip_tags/genomicData_width) + 1e-5, col=is.finite(padj) & padj<0.05 & log2fold_Control_Alyref > 0.58, alpha=is.finite(padj) & padj<0.05 & log2fold_Control_Alyref > 0.58) + geom_point(size=1) + facet_wrap(~slice, scale="free_y") + scale_x_log10() + scale_y_log10() + scale_alpha_manual(limits=c(T,F), values=c(1,0.25), name="significant") + geom_smooth(alpha=0.3) + scale_color_manual(limits=c(T,F), values =c("red","black"), name="significant") + xlab("RNAseq Density") + ylab("CLIP density")

   Density of RNAseq reads vs density of clip tags on retained introns, introns significantly downregulated in KD are highlighted



Detained vs. Retained introns
------------------------------

.. report:: RetainedIntrons.FractionsOfDiffIntronsDetained
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=in_category, y = overlap/(overlap+no_overlap)) + geom_bar(stat="identity") + facet_wrap(~track, scale="free_y") + theme_bw() + theme(aspect.ratio=1) + xlab("Detained Intron?") + ylab("Fraction differential")

   Fraction of introns in the detained category that are differential on ALyref knockdown


HEK293 detained and Retained
----------------------------


.. report:: RetainedIntrons.DetainedChunkSplicing
   :render: r-ggplot
   :groupby: all
   :statement: aes(cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=20)) + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold<0), col=T),fun.y=mean, geom="point") + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold < -0.58), col=F),fun.y=mean, geom="point")+ facet_grid(slice~track, scale="free_y") +theme_bw() + scale_color_manual(values=c("FALSE"="red", "TRUE"="black"), labels=c("FALSE"="FC>1.5", "TRUE"="FC>0"), name="Fold Change") + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Fraction significant increase in splicing efficiency")

   Fraction of introns in detained category that show increased splicingon Alyref Knockdown

Detained introns defined used FDR < 0.01 lfc > 2


.. report:: RetainedIntrons.RetainedChunkSplicing
   :render: r-ggplot
   :groupby: all
   :statement: aes(cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=20)) + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold<0), col=T),fun.y=mean, geom="point") + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold < -0.58), col=F),fun.y=mean, geom="point")+ facet_grid(slice~track, scale="free_y") +theme_bw() + scale_color_manual(values=c("FALSE"="red", "TRUE"="black"), labels=c("FALSE"="FC>1.5", "TRUE"="FC>0"), name="Fold Change") + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Fraction significant increase in splicing efficiency")

   Fraction of introns in retained category that show increased splicingon Alyref Knockdown

Detained introns defined used FDR < 0.01 lfc > 2


.. report:: RetainedIntrons.CombinedChunkSplicing
   :render: r-ggplot
   :groupby: all
   :statement: aes(cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=20)) + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold<0), col=T),fun.y=mean, geom="point") + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & l2fold < -0.58), col=F),fun.y=mean, geom="point")+ facet_grid(slice~track, scale="free_y") +theme_bw() + scale_color_manual(values=c("FALSE"="red", "TRUE"="black"), labels=c("FALSE"="FC>1.5", "TRUE"="FC>0"), name="Fold Change") + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Fraction significant increase in splicing efficiency")

   Fraction of introns in retained category that show increased splicingon Alyref Knockdown

Detained introns defined used FDR < 0.01 lfc > 2


Enrichment
+++++++++++

The following show the *enrichment* of De/Retained Introns in transcript regions that are differentially reduced on Alyref knockdown.

.. report:: RetainedIntrons.DetainedEnrichment
   :render: r-ggplot
   :statement: aes(x=clip_bin, y=enrichment, col=p<0.05/20) + geom_point() + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Enrichment") + scale_color_manual(name="", labels=c("TRUE"="Significant", "FALSE"="Not Significant"), values=c("TRUE"="red", "FALSE"="black")) + facet_grid(slice~.) + theme_bw()

   Fold Enrichment of detained introns in chunks more efficiently spliced on ALyref knockdown


.. report:: RetainedIntrons.RetainedEnrichment
   :render: r-ggplot
   :statement: aes(x=clip_bin, y=enrichment, col=p<0.05/20) + geom_point() + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Enrichment") + scale_color_manual(name="", labels=c("TRUE"="Significant", "FALSE"="Not Significant"), values=c("TRUE"="red", "FALSE"="black")) + facet_grid(slice~., scale="free_y") + theme_bw()

   Fold Enrichment of retained introns in chunks more efficiently spliced on ALyref knockdown


.. report:: RetainedIntrons.CombinedEnrichment
   :render: r-ggplot
   :statement: aes(x=clip_bin, y=enrichment, col=p<0.05/20) + geom_point() + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Enrichment") + scale_color_manual(name="", labels=c("TRUE"="Significant", "FALSE"="Not Significant"), values=c("TRUE"="red", "FALSE"="black")) + facet_grid(slice~., scale="free_y") + theme_bw()

   Fold Enrichment of detained or retained introns in chunks more efficiently spliced on ALyref knockdown
