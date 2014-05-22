.. _context:

Contexts of mapped reads
-------------------------

If the reads we are sequencing come primarily from RNA pulled down the iCLIP process, we would expect that there most reads would map to areas of the genome annotated as RNA. If the reads were effectively random across the genome, we would expect that reads would be found in areas of the genome annotated as RNA in proportion to the fraction of the genome such contexts cover.

Below I show for each "context", the log ratio of the observed fraction of reads mapping to that context compared to the expected based on the fraction of the genome that that context covers. The 'none' category reffers to areas of the genome not annotated as being covered by any type of RNA.


.. report:: Sample_QC.ContextRepresentation
   :render: r-ggplot
   :statement: aes(category, log2(precent_alignments/percent_bases)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=90,hjust=1)) + ylab("log2 enrichment")
   :layout: column-3
   :groupby: track

   Enrichments of contexted over expectation


The most striking observation is that reads annotated as mapping to ribosomal RNA sequences are massively enriched compared to expectation in all samples (greater than :math:`2^10` fold) , even the control pull downs. The same is true to a slightly lesser extent for snRNAs and snoRNA. These are highly expressed RNA sequences, and their presence could be due to contamination of RNA not associated with the antibody in the pull down. A second alternative is that the FLAG antibody is showing some affinity for an RNA binding factor that binds to them, and they are being pulled down in all samples. 

The second point to note is that reads are depleted in mapping to areas not annotated as RNA in the genome in all samples. This suggests that we are sequencing RNA, rather than genomic contamination.

In Alyref and Chtop there is a enrichment in exonic sequence (protein_coding) and a depletion in intronic sequencing, showing that these factors have an affinity for mature mRNA rather than primary transcript. These patterns are replicated, but less stringly in Nxf1. In the control pull downs the enrichment for exonic sequence is reduced, and the depletion for intronic sequence abolised or even inverted. 

In all the samples the enrichment for lincRNA is similar to that for protein coding genes. 

My conclusions from this are that Alyref and Chtop are binding to mature mRNA and similar transcripts, while the FLAG pulldown is most likely returning RNA in proportion to its abundance in the cell. 
