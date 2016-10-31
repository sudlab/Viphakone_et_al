Gene profiles normalised to RNAseq data
=========================================


CDS/UTR matricies
------------------


.. report:: NormalisedProfiles.NormalizedUTRMatrix
   :render: r-ggplot
   :regex: gene_profiles.dir/(.+).rnaseq_normed
   :glob: gene_profiles.dir/*.rnaseq_normed.profile.tsv.gz
   :groupby: all
   :statement: aes(position, track, fill=profile) + geom_raster() + scale_x_continuous(labels=c("upstream", "5 UTR", "CDS", "3 UTR", "downstream"), breaks=c(500,1100,1700,2550,3400)) + theme_bw() + theme( aspect.ratio = 0.5, legend.position = "none") + xlab("") + ylab("") + geom_vline(xintercept=c(1000,1200,2200,2900,3900), col = "white", lwd=0.25, lty=2) + scale_fill_gradientn(colours=c("black","#56B1F7","#56B1F7"), values = c(0,quantile(rframe$profile,0.995)/quantile(rframe$profile,1),1))

   UTR/CDS profiles, normalised per gene to RNAseq data



.. report:: NormalisedProfiles.NormalisdSummaryUTRMatrix
   :render: r-ggplot
   :regex: gene_profiles.dir/(.+).utr
   :glob: gene_profiles.dir/*FLAG*.utrprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(bin, track, fill=normalized) + geom_raster() + scale_x_continuous(labels=c("upstream", "5 UTR", "CDS", "3 UTR", "downstreaam"), breaks=c(500,1100,1700,2550,3400)) + theme_bw() + theme( aspect.ratio = 0.5, legend.position = "none") + xlab("") + ylab("") + geom_vline(xintercept=c(1000,1200,2200,2900,3900), col = "white", lwd=0.25, lty=2) + scale_fill_gradientn(colours=c("black","#56B1F7","#56B1F7"), values = c(0,quantile(rframe$profile,0.995)/quantile(rframe$profile,1),1))

   UTR/CDS profiles, normalised using average RNASeq profile


Same plots as line plots:

.. report:: NormalisedProfiles.NormalizedUTRMatrix 
   :render: r-ggplot
   :regex: gene_profiles.dir/(.+).rnaseq_normed
   :glob: gene_profiles.dir/*.rnaseq_normed.profile.tsv.gz
   :groupby: all
   :statement: aes(position,profile) + facet_wrap(~track) + geom_line() + geom_vline(xintercept=c(1000,1200,2200,2900,3900), lwd=0.25, lty=2) + scale_x_continuous(labels=c("upstream", "5 UTR", "CDS", "3 UTR", "downstream"), breaks=c(500,1100,1700,2550,3400)) + theme_bw()

    UTR/CDS profiles, normalised per gene to RNAseq data


.. report:: NormalisedProfiles.NormalisdSummaryUTRMatrix
   :render: r-ggplot
   :regex: gene_profiles.dir/(.+).utr
   :glob: gene_profiles.dir/*.utrprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(bin,normalized) + facet_wrap(~track) + geom_line() + geom_vline(xintercept=c(1000,1200,2200,2900,3900), lwd=0.25, lty=2) + scale_x_continuous(labels=c("upstream", "5 UTR", "CDS", "3 UTR", "downstream"), breaks=c(500,1100,1700,2550,3400)) + theme_bw()

    UTR/CDS profiles, normalised to average RNAseq data

TSS profiles
-------------

.. report:: GeneProfiles.UTRProfiles
   :render: r-ggplot
   :regex: gene_profiles.dir/Chtop-FLAG-(.+).tssprofile
   :glob: gene_profiles.dir/Chtop-FLAG-*tssprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(region_bin-400, area, col=track) + geom_line() + facet_wrap(~region) + geom_vline(xintercept=0, lty=2, lwd=0.5) + theme_bw(base_size=18) + xlab("Relative Position") + ylab("") + scale_y_continuous(breaks=NULL) 

   Unnormalized tss profiles



.. report:: NormalisedProfiles.NormalizedTSSMatrix
   :render: r-ggplot
   :regex: gene_profiles.dir/Chtop-FLAG-(.+).tssprofile
   :glob: gene_profiles.dir/Chtop-FLAG-*tssprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(region_bin-400, normalized, col=track) + geom_line() + facet_wrap(~region) + geom_vline(xintercept=0, lty=2, lwd=0.5) + theme_bw(base_size=18) + xlab("Relative Position") + ylab("") + scale_y_continuous(breaks=NULL)

   Normalised to RNASeq



Single vs Multi Exon profiles
--------------------------------

.. report:: GeneProfiles.SingleVsMultiExonProfiles
   :render: r-ggplot
   :statement: aes(bin, density, colour=exons) + geom_line() + facet_grid(slice~., scale="free_y") + theme_bw() + geom_vline(xintercept=c(25,50), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Gene","Downstream"), breaks=c(12,37,62)) + theme_bw() + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL)

   Gene Profiles divided into single and multi exon genes

.. report:: GeneProfiles.AverageSingleVsMultiExonProfiles
   :render: r-ggplot
   :groupby: all
   :statement: aes(bin, density, colour=exons) + geom_line() + facet_grid(track~., scale="free_y") + theme_bw() + geom_vline(xintercept=c(25,50), lwd=0.5, lty=2) + scale_x_continuous(labels=c("Upstream","Gene","Downstream"), breaks=c(12,37,62)) + theme_bw() + xlab("")+ ylab("Relative Read depth") + scale_y_continuous(breaks=NULL)

   Average profiles

.. _profiles-by-quantile:

Gene Profiles by gene length quantile
--------------------------------------

Protein coding genes were divided into 5 bins with equal numbers of transcripts based on their length, and coverage profile calculated for each bin.

.. report:: GeneProfiles.BinnedExpressionProfiles
   :render: r-ggplot
   :groupby: track
   :statement: aes(x=bin, y=area, col=factor(quantile, levels=sort(quantile, decreasing=T))) + geom_line() + facet_grid(slice~exon_limit, scale="free_y") + scale_y_continuous(breaks=NULL) + ylab("Relative coverage") + geom_vline(xintercept=c(250,500), lty=2,lwd=0.5) + scale_x_continuous(breaks=c(125, 375, 625), labels = c("Upstream","CDS", "Downstream")) + xlab("") + scale_color_manual(values=colorRampPalette(c("#132B43","#56B1F7"))(5), name = "Length\nQuantile") + theme_bw()

   Metagene profiles for genes binned by length



Some statistics on the sets of transcripts used above:

.. report:: GeneProfiles.ExpressedTranscriptStats
   :render: r-ggplot
   :statement: aes(x=factor(quantile+1), y=Exon.Length) + geom_boxplot() + scale_y_log10(breaks=c(10,100,1000,10000,100000), labels = c("10 bp","100 bp","1 kb", "10 kb", "100 kb"), limits=c(10,100000)) + ylab("Exonic Length") + xlab("Length Quintile") + theme_bw()

   Distribution of transcript lengths in each quintile


.. report:: GeneProfiles.ExpressedTranscriptStats
   :render: r-ggplot
   :statement: aes(x=factor(quantile), y=TPM) + geom_boxplot() + scale_y_log10(breaks = c(1,10,100,1000,10000)) + theme_bw() + ylab("Expression Level (TPM)") + xlab("Length Quintile")

   Distribution of transcript expression levels in each quintile

