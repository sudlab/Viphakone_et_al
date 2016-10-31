Processing Index
=================


Method
-------

Defining the processing index
++++++++++++++++++++++++++++++

The processing index is basically a ratio of the reads in the unprocessed
pre-mRNA to the processed mature RNA. **Higher ratios** mean that a higher
proportion of the reads fall in the **unprocessed pre-mRNA.**

It was originally introduced by Baejen et al Mol Cell 5(55):745-757
:PMID:`25192364`, where they describe it as:


  We assume that read counts Ndown downstream of a pA site can only occur from pre-mRNAs,
  Ndown = Nprem, whereas read counts Nup upstream of a pA site are a mixture of mature mRNA
  counts Nmat and pre-mRNA counts Nprem. . Therefore, Nup = Nmat + Nprem. For increased
  robustness with regard to different transcript isoforms and uncertainties in the exact location of
  pA sites, we computed Nup and Ndown as average of the read counts within 50 nt upstream and
  downstream NPM of the pA site over all G gene transcripts, respectively. Similar to the ‘splicing
  index’, we define the ‘processing index’ (PI) as follows (Equation 2). 

and they give equation 2 as:

.. math::
   PI==log_2(\frac{1}{G} \frac{\sum_{i=1}^G N_i^{PM}}{\sum_{i=1}^G N_i^M})


this definition didn't make a lot of sense to me, as I couldn't figure what the
1/G was doing. So I emailed the authors, and indeed, this definition is
incorrect. What they actaully did was to calculate the average depth per base in
the up and down stream windows (remember this is PAR-CLIP, not iCLIP and so tags
cover multiple bases) for each gene, calculate the average for up and down
stream seperately, and then calculate the following:

.. math::
   PI == log_2( \frac{N^{down}}{min (1, N^{up}-N^{down})}

where :math:`N^{up}` and :math:`N^{down}` are genome wide averages. This is
unlikely to be suitable for iCLIP because tags are single bases, and so the
averge :math:`N^M` is unlikely to ever be as a large as 1. Thus I actaully
calculated

.. math::
   PI==log_2( \frac{\sum_{i=1}^G N_i^{down}}{min(1,\sum_{i=1}^G N_i^{up}-N_i^{down})})

I used 50 bp windows up and downstream of the defined cleavage sites.


Defining the set of cleavage sites
+++++++++++++++++++++++++++++++++++

Cleavage sites were obtained from Martin *et al* (2012) :PMID:`22813749`. In
Baeken *et al*, they use the most highly expressed transcript and find
"annotated each gene’s TSS and pA". However, I don't think this is appropriate,
because if a there are other transcripts, less well expressed, that are longer
than the highest expressed, they will have a large effect on the estimation of
the amount of pre-mRNA, infalting it artefactually.

Instead, I take the 3' most (or 5' most on the - strand) cleavage site that
overlapps any expressed transcript from a gene (TPM > 1), as any reads more 3'
than than this are bound to be pre-mRNA.

Results
---------

.. report:: processing_index.ProcessingIndex
   :render: r-ggplot
   :transform: pandas
   :tracks: union
   :statement: aes(x=factor, y=processing_index) + geom_bar(stat="identity") + theme_bw(base_size=18) + coord_flip() + ylim(-7,7) + xlab("Sample") + ylab("Processing Index") + geom_hline(yintercept=0) + theme(aspect.ratio=0.5)
   :tf-statement: loc["union"]

   Processing index from union tracks


As we can see the processing index is positive in all cases. This means that for
all of our factors, there is more bound downstream of the cleavage point than
upstream. i.e. the factors all bind to the pre-mRNA, not the mature.

The second point to notice is that ChTop and Nxf1 have a higher index
than Alyref or FlipIn, but notice that FlipIn *does* have a positive index
(see below). 

Limitations
------------

There are several limitations to this analysis. The two serious are :

1. Mapping bias at the end of molecules
2. Undue influance of highly expressed transcripts.

Mapping bias
++++++++++++
 
There will be a mapping bias towards reads from pre-mRNAs. Consider:
we cannot map reads closer than 15bp to the end of a molecule, and the further
we are away from the end, the more likely it is that we will be able to map. Thus
we are likely to miss reads that come from the mature mRNA, but will capture those 
that come from the pre-mRNA.


Highly expressed transcripts
+++++++++++++++++++++++++++++

Consider the following situation: 

Transcript 1 has 101 tags upstream of the cleavage site and 100 downstream
Transcript 2 has 11 tags upstream of the cleavage site, and 1 downstream.
Transcript 3 has 3 tags upstream of the cleavage site, and 0 downstream.

The total PI would be 115/101, which is more or less 1, but in two of the
three cases there is clearly a bias towards there being more tags upstream
rather than downstream. I.e., if we had ignored transcript 1, the result
would have been 14/1 = 14, a very different processing index. 





