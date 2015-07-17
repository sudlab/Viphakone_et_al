.. _mapping:

Mapping and quality control
===========================


The new data set contains many more reads than the pilot data. 

The total number of reads for each sample is as follows:

.. report:: Sample_QC.TotalMergedReads
   :render: r-ggplot
   :groupby: all
   :force:
   :statement: aes(track,Reads/1000000) + geom_bar(stat="identity") + theme_bw() + ylab("Reads (millions)") + theme(axis.text.x=element_text(angle=90))

   Total number of (unmapped) reads per sample


The increase isn't so drastic once we map and remove PCR duplicates, but there is still a significant increase in the amount of data.

.. report:: Sample_QC.FinalReads
   :render: r-ggplot
   :statement: aes(y=reads_mapped, x=sample) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x,...) format(x,...,big.mark=",", scientific= F, trim = T)) + ylab("Total unique mapped reads")

   Final total mapped reads


Before we extended the pilot to the production run, we had analysed the saturation of the libraries. We concluded for the majority of libraries, saturation hadn't been reached, but wasn't far off. Have we now reached saturation?

.. report:: Sample_QC.AlignmentSaturation
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=subset, y=counts, color = factor, shape = factor) + geom_point() + geom_line() + facet_wrap(~replicate) + theme_bw() + theme(aspect.ratio = 1) + xlab("Fraction Sampled") + ylab("Unique Read Count")

   Subsampling of alignments

So no, we havn't. What is more we don't seem to be any closer than we were with the pilot data. It seems that the data begins to saturate, but the changes from being an exponentially saturating curve to being a linear one. One possible explaination for this is that we are sequencing each read so many times that we are picking up base errors in the UMI sequence. This could be causing us to be call reads as different when they are not. This requires investigating further, but could be a worrying developement. 


Finally we can build a table of summary statisitics:

.. report:: mapping.FinalStatsTable
   :render: table
   

   Mapping statistics

