Results of Motif finding
===========================


There are several ways to find motifs in sequences that have different strengths and weaknesses. The difference between them mostly depends on what we mean by "expect" when we say that the sequence motif is more common than we expect. So far in this project we have used two different approaches:

1. The EM approach as implemented by the MEME algorithmn. In this approach the alogrithmn builts up a model of what the probability of the next base being X based on the previous base being Y and then computes the probability of a seeing a particular set of sequences based on this. This approach is very flexible, allows us to check for many different types and widths of motifs, from very long to very short, probabilitic motifs (eg. first base is probably T but could be C) etc. However it is a) very slow meaning we can't use all of our sequences as input and b) very dependent on the model used to generate the expected probabilities. As such, if for example your sequecnes are from 3' UTRs, this model will pick up the ATAA site just because it is enriched within 3' UTRs, not because the protein is enriched at that seqeuence within 3' UTRs.

2. The Regex and hypergeometric approach as implemented by the DREME alogrigthm. In this approach we simply ask is a particular string of letters is more common in set of sequences A, than in set of sequences B. Seqeunces B can either come from a control factor, or by shuffling the sequences in A randomly (to maintain base composition). This is much faster than method 1, and also accounts properly for the negative control, but only works with a fixed set of short motif sizes and the motifs are not probabilistic (i.e. it can cope with base 1 is T or C, but not that T is more likely).

We have implemented both these methods on the data.


Results of MEME analysis
--------------------------

MEME was run on clusters from all tracks as well as the set of reproducible clusters. It was run on control clusters seperately to test clusters and so we are looking to see if the same clusters come up in both. I will summerize these results below.

The MEME results are judge on several criteria: 

* The evalue tells you the statistical significance: it is the number of motifs of that width and site counts that you would expect to see with the same or better log likelihood ratio

* The information content: How well defined the motif is.

* The number of sites contributing to the motif.


TAATTTTTGTATTTTT
+++++++++++++++++++

The most prominent motif found in the reproducilbe peaks from both Alyref and ChTop is the sequence GGCTAATTTTTGTATTTTTAGTAGAGATGG and varients thereof:

.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: all
   :large-html-class: sortable
   :tracks: r(reproducible)
   :slices: r(TAATTTTTGTATTTTT)
   :force:

   Motifs containing TAATTTTTGTATTTTT in reproducible clusters


It is also present in many of the indevidual replicates. Unfortunately it is also present in indevidual replicates of the FlipIn control, meaning that it is probably either an artifact of the method, or a piece of RNA that sticks to FLAG. It is probably not found in the reproducible FLAG replicates simply because there are so few. 

.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: all
   :large-html-class: sortable
   :tracks: r(FLAG-R[0-9])
   :slices: r(TAATTTTTGTATTTTT)
   :force:

   Motifs containing TAATTTTTGTATTTTT in indevidual replicates.


GAGNNNGAG
++++++++++

The next most attension grabbing motif is GAGNNNGAG, which is seen in the reproducible tracks for both Alyref and ChTop. Although the E-value is very high for ChTop and would usually be considered significant I don't think. 

.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: all
   :tracks: r(reproducible)
   :slices: r(GAG...GAG)
   :force:

   Motifs containing GAGNNNGAG in reproducible clusters


Interestingly, while this motif is not found in indevidual replicates of FLAG tag only, it is found in 3 of the 4 replicates of Alyref, but only one of the replicates of ChTop. Even more interesingly, slighly degenerate forms are seen in two of the Nxf1 replicates.

.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: all 
   :large-html-class: sortable
   :tracks: r(FLAG-R[0-9])
   :slices: r(GAGA..GAG)
   :force:

   Motifs containing GAGNNNG in indevidual replicates


Motifs containing GAGNNNGAG in reproducible clusters


TCTCCA containing motifs
++++++++++++++++++++++++++++++


Motifs containig this sub motif are interesting because they match the motif found for m6C bases from Van Haute et al. Such a motif was located in reproducible clusters from Nxf1. Unfortunately, a very similar motif was also found in the reproducible clusters from the FLAG only samples. 


.. report:: iCLIPMotifs.MemeResults
   :render: table
   :groupby: all
   :large-html-class: sortable
   :slices: r(TCTCCA)
   :force:

   Motifs containing TCTCCA


Other motifs
+++++++++++++++

I consider all the other motifs found in the reproducible clusters to either be too low on one of the criteria outline above, but feel free to look through the full results at `this page <https://www.cgat.org/downloads/N6Cduavf7p/iCLIP_fullrun2/report/pipeline/Motifs.html>`_

Results of Randomised DREME analysis
----------------------------------------

In these analyses the negative set for DREME is built by randomising the positive set, and the freqeuncy of various motifs compared between the positive and negative sets. Using a randomised set instead of the negative set gets around the problem of only having a very small negative set if we use the FLAG only as a control, but suffers from the fact that any motifs caused by higher order patterns are not accounted for.

All DREME motifs are fiexed at 7 nt long for this analysis. 

There are many other Dreme motifs enriched in these samples. Find the full results `here <https://www.cgat.org/downloads/N6Cduavf7p/iCLIP_fullrun2/report/pipeline/Motifs.html#dreme-motifs-within-clusters>`_


AATTAG
++++++

The most stiking feature of these tests is the motifs containing AAATTAG. This is present in both the ChTop and the Alyref reproducible clusters. 


.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :tracks: r(reproducible)
   :large-html-class: sortable
   :groupby: all
   :slices: r(AATTAG)
   :force:

   Motifs containing AAATTAG in reproducible clusters


Unfortunately, as well as being present in indevidual replicates for these factors as well as for Nxf1, such motifs are also present in replicates 2 and 4 of the FLAG only samples. 


.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :tracks: r(FLAG-R[0-9])
   :large-html-class: sortable
   :groupby: all
   :slices: r(AATTAG)
   :force:

   Motifs containing AAATTAG in reproducible clusters

CGCCATG
+++++++++

This motifs shows up in the reproducible Alyref results, but not in any of the indevidual replicates for any factor. 

.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :large-html-class: sortable
   :groupby: all
   :slices: r(CGCCATG)
   :force:

   Motifs containing CGCCATG


CYCCR
++++++

The strongest motif in the Chtop results is the motif ACYCCRTC. This is of course of interest because if could be a match for the m5C motif from  Van Haute et al. If we look for similar motifs in other samples we find very simlar matches in replicate 2 of Chtop, and unfortunately in replicate 3 of FLAG only:

.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :groupby: all
   :slices: r(C[YWKT]C[CYSM][ARWM])
   :force:

   Close matches to the CYCCR motif


GWTACAGA
+++++++++

The only motif enriched in reproducilbe Nxf1 clusters, this nothing like this motif is found in any other sample

.. report:: iCLIPMotifs.SimpleDremeResults
   :render: table
   :groupby: all
   :slices: r(TACAG)
   :force:

   The enriched motif in reproducible Nxf1 clusters




Todo
------


* Write about descriminative Dreme clusters
* Motifs in subsets of RNA: retained introns, UTRs etc. 
* zagros motifs

