Presence of poly-A reads
-------------------------

If the factors bind to the poly-A tail, then perhaps there reads are 
pulled down and sequenced, but not mapped (because they don't match anything
in the genome). To look for this I measured the nucleotide composition of unmapped reads

.. report:: unmapped_reads.UnmappedComposition
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=bin, y=count, fill=slice) + geom_bar(stat="identity") + facet_grid(track~slice, scale = "free_y") + xlab("Fraction of read base N") +  theme_bw(base_size=10) + theme(legend.position="none")

   Nucleotide composition of unmapped reads


ALthough there is a large number of reads that are mostly A, this is true of all factors including the FlipIn samples. 
