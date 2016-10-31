Exploring profiles using heatmaps
==================================

Headlines:

*  5' Enrichment in Alyref is a proportaion of gene length rather than a fixed amount
*  The above might be negated by normalisation to RNA, but the results are inconclusive
*  The 3' enrichment on ChTop looks like it is probably UTRs rather that 3' end. 
*  There is an enrichment of tags on PROMPTs

If there is a plot that is missing, you will probably find it at :ref:`pipeline-heatmaps` pipeline page.


Details - Standard heatmaps
----------------------------

In standard heatmaps of this type, we look at the tag counts in 25bp bins across the genes and then sort the genes according their length. Profiles are normalised to the 99th centile and then averaged veritcally to make them fit on the image. They are then quantile normalised again. This means that like metagene plots, the results we see here are relative.  I've trimmed the plots to 10kb to make them fit. 


.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(forward)

   Tag heatmaps for reads on forward strand of genes sorted by length

Looking at Alyref, we can clearly see the 5' bias. This bias appears to lengthen as the genes get long. That is it looks like the first x% of the gene is enriched, rather than the first y base pairs. We can also see somthing that looks a bit like the start of a 5' among the longest genes for ChTop, but a much stronger pattern for ChTop is the bias to the 3' end. It also looks a bit like there might be something of a 3' bias for Alyref. We get a clearer look at this if we align genes at their ends rather than at their starts:


.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: Alyref-FLAG.union.end_aligned.length, Chtop-FLAG.union.end_aligned.length

   Tag heatmaps for reads on forward strand of genes sorted by length

The ChTop 3' enrichment is very clear on these plots, and in constrast to the Alyref 5' end enrichment, does look to be a fixed distance from the end of the gene (expect for the very shortest genes). The AlyRef 3' enrichment looks very weak here, and restricted to the shorter genes. 

Sorting by UTRs
----------------

One of the questions that we hoped to answer with this analysis is weather the 3' enrichment is for the 3' end of the gene in general or for the 3' UTR particularly. One way to attack this is to look at sort the genes by the length of their UTR, rather than by the length of gene. If the enrichment is to the UTR, we might see the enrichment lengthen as the UTR gets longer, and stop suddenly rather than fading out. The following plots do this:

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :tracks: Chtop-FLAG.union.end_aligned.3utr

   Tag heatmap for ChTop, end aligned and sorted by 3' UTR length

This was not the expected result - It looks for all the world like the enrichment is at the stop codon rather than at the end of the gene or accross the 3' UTR. However there are some reasons to be sceptical of this enrichment. The frist is that it is present for the other factors, including, very faintly, FlipIn

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: Alyref-FLAG.union.end_aligned.3utr, FlipIn-FLAG.union.end_aligned.3utr

   Tag heatmaps for Alref and FlipIn, end aligned and sorted by 3' UTR length

The second reason to be sceptical of this effect is that the analysis is performed on merged genes - where all transcripts for a gene are merged into a single annotations for the gene. It is not unheard of for genes to have isofroms that have very long UTRs. If these long UTRs were spurious or not expressed in HEK293 cells (or expressed at a lower level - the gene set is filtered for very lowly expressed transcripts), then this might lead to the pattern we see. 

One way to look at that is to look at the RNA expression levels. If we generate similar plots using total nuclear RNA, they look like so:

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: Nuclear-RiboZ-R1.star.end_aligned.3utr, Nuclear-RiboZ-R1.star

   Coverage heatmap for nuclear RNA

Here we can see the same enrichment at the stop codon we saw in the iCLIP heatmaps, suggesting that indeed, this effect is not real. 

Normalisation by RNA-seq
-------------------------

We have two questions that we might conceivably answer by normalising to the amount of RNA. 

1) Is the enrichment of AlyRef at the 5' end purely due to the higher quantity of 5' RNA in transcribing genes, or is there some cap mediated enrichment on top of this?
2) Is the 3' enrichment in the 3' UTR for ChTop purely due to the problems with defining the 3' end of genes?

When we are looking at the normalised profiles we should be aware that it is not there are going to be fewer RNAseq reads right at the start and right at the end of genes - this is because RNA seq fragments are 200-300bp long, and so there are fewer possible reads to be found at the very start and end of genes. We can see this is we look at the start of the start aligned heatmap above. Care must be taken because normalising by this will give an enrichment at the 5' end of genes even if the iCLIP trace is level - that is it might seem we have an enrichment at the end which could be caused by a reduction in the level of RNA, rather than an incrase in the level of iCLIP. 

With that caveat in mind, the RNA normalised heatmaps:

.. report:: Heatmaps.NormedHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(forward)

   RNAseq normalised tag heatmaps

First we can note that the enrichment that was proportional to the length of the gene has gone from the AlyRef trace. There is however, a strong enrichment in a fixed length window at the start of the gene. However, if is unclear if this is real because of the reasons outlined above. Indeed we see similar enrichments at the 5' for all factors, including FlipIn. The enrichment is particularly strong in Alyref, but it is unclear how to rule out the normalisation artifact. 

For the 3' UTR, again the result is not entirely clear:

.. report:: Heatmaps.NormedHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(end_aligned.3utr)

   RNAseq normalied, end_aligned, 3' UTR length sorted tag heatmaps


We can see that a large amount of the enrichment at the stop codon is removed by this normalisation. However there are still some artefacts - first is the light band around the stop codon in Alyref. Second, the enrichement is still missing between the stop codon and the termination site for the longest of genes in ChTop. Finally "something" is happening in FlipIn and Nxf1. My take on this is that the data favours a UTR enrichment, rather than a 3' end enrichment, but is not conclusive. 


PROMPTs
-----------

Finally we can use this analysis to look for evidence that PROMPTs are bound by performing the same analysis, but this time on reads from the minus strand.


.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(reverse)

   Tag heatmaps for reads on the reverse strand to the gene

There is a clear enrichment just upstream of the TSS in Alyref. It is weaker, but present in ChTop and absent in FlipIn. I'd say this is pretty good evidence that PROMPTs are bound by at least Alyref and ChTop. 


