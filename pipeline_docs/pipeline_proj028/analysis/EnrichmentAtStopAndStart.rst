.. _utrprofiles:

Enrichment around the start and stop codons
===========================================

Recent work has shown that RNA is often methylated around the start and stop codons of the coding sequence of a transcript. The factor responsible for this is WTAP, and WTAP interacts with Chtop. Chtop is enriched at the end of genes. Could this enrichment at the end be an enrichment at the stop codon?

The methylation paper produced a profile showing the enrichment for 400bp either side of the start and stop codons. Here we reproduce that plot from Chtop, note that although TSS and TTS are labeled in the below, these are actaully the start and stop sites of the CDS, not of the transcript. 


.. report:: GeneProfiles.UTRProfiles
   :render: r-ggplot
   :regex: gene_profiles.dir/Chtop-FLAG-(.+).tssprofile
   :glob: gene_profiles.dir/Chtop-FLAG-*tssprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(region_bin-400, area, col=track) + geom_line() + facet_wrap(~region) + geom_vline(xintercept=0, lty=2, lwd=0.5) + theme_bw(base_size=18) + xlab("Relative Position") + ylab("") + scale_y_continuous(breaks=NULL) 

   Unnormalized tss profiles


This looks fairly similar to the plot in the methylation paper. However, the drop in density after the stop could be due to some UTRs being shorter than 400, so that the averages further away are including some non-utr sequence. Note that you can also see this in the original mehtylation paper as a drop in the level of there input sample. There are two possible ways to fix this. One would be to normalize by RNA expression - if the drop off is because of short UTRs, the RNA signal should drop off as well:

.. report:: NormalisedProfiles.NormalizedTSSMatrix
   :render: r-ggplot
   :regex: gene_profiles.dir/Chtop-FLAG-(.+).tssprofile
   :glob: gene_profiles.dir/Chtop-FLAG-*tssprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(region_bin-400, normalized, col=track) + geom_line() + facet_wrap(~region) + geom_vline(xintercept=0, lty=2, lwd=0.5) + theme_bw(base_size=14) +coord_cartesian(xlim=c(-350,350)) + xlab("Relative Position") + ylab("") + scale_y_continuous(breaks=NULL)

   Normalised to RNASeq


This indeed shows that when you normalise for RNA expression, the peak at the start stays, but at the stop codon, the density becomes high in at the stop and then stays high afterwards.

A different way to approach this would be to normalize for each transcripts UTR length, as shown below:


.. report:: NormalisedProfiles.NormalisdSummaryUTRMatrix
   :render: r-ggplot
   :regex: gene_profiles.dir/(.+).utr
   :glob: gene_profiles.dir/*FLAG*.utrprofile.matrix.tsv.gz
   :groupby: all
   :statement: aes(bin, track, fill=normalized) + geom_raster() + scale_x_continuous(labels=c("upstream", "5 UTR", "CDS", "3 UTR", "downstream"), breaks=c(500,1100,1700,2550,3400)) + theme_bw() + theme( aspect.ratio = 0.5, legend.position = "none") + xlab("") + ylab("") + geom_vline(mapping=aes(xintercept=c(1000,1200,2200,2900,3900)), col = "white", lwd=0.25, lty=2) + scale_fill_gradientn(colours=c("black","#56B1F7"))

   UTR/CDS profiles, normalised using average RNASeq profile


In both these cases we see that the enrichment for ChTop appears to be across the whole UTR, and not just at the STOP codon. However, it remains to be seen if the finding in the original paper is suffering from the same artifact.

