Comparison of first and last exons to middle exons
==================================================

The aim of this analysis is to detect if there is a difference between
the first and last exons of a gene and the middle exons. The principle
reason why such differences might occur are that there is not EJC on
the final exon. Thus we will also study the eIF4A3 component of the
EJC that we have implicated in Alyref loading.

We have decided to use the frist and last exons, rather than the UTRs
and coding regions because from a nuclear stand point there should be
no difference between UTR and coding sequence per-se. Current models
would suggest that only ribosommes can disginuish between coding and
non-coding sequence.

Method
------

Gene-set
++++++++

Expressed transcripts (TPM>1) have been "merged" so that all
overlapping exons are combined. We have then filtered out only protein
coding genes with more than 3 exons in their merged form. Note that
this might lead to minor isoforms with longer first or last exons
changing what is generally regarded as the canonical first of last
exons. Next on the todo list is to isolate a set of genes with only
those genes with a single dominant isoform.

Counting and normalisation
++++++++++++++++++++++++++

We use the FlipIn counts as our primary source of normalisation data.
* After first and last exons had been isolated, tags were counted in
  first, last and middle exons, as well as sequence that is never
  inluded in an exon in any isoform ("introns").
* For TREX components the base preceeding the first base of the read
  is used, for EJC components the base at the centre of the read is
  used.
* The data will be presented in two ways.
   1) These counts are then normalised to the total for the gene,
      producing a **fraction that is present in that category for that
      gene**.
   2) Left as raw counts, representing **absolute amount of binding**
      that region. Note that in these data there will be a lot of skew
      introduced by the most highly expressed genes.
* Counts accross all genes are summed.
* The summed counts for each protein is then normalsed to the control
  (HEK293 FlipIn-FLAG or HeLa FlipIn-GFP).

Although in most spatial analyses the flipin track is too noisey, here
we are summing over sufficiently large areas that we hope this isn't a
problem. Normalising by FlipIn should solve several problems: It
should control for different expression levels of different genes and
different parts of genes. It should control for the relative real
lengths of the exons (i.e. what is actaully expressed rather than just
what we are computing). It will do this better than RNAseq for two
reasons:

1) We still havn't resolved whether chromatin associated, nuclear or
    total RNA is the best source for normalisation, and we getting
    different results depending on what we use
2) The EJC work was done in HeLa cells, thus the RNAseq appropriate
    for TREX is not appropriate for EJC. Even if we find HeLa RNAseq,
    we can't guarenttee it is equivalent to the HEK293 RNAseq we have
    in all ways.

Results
-------


First we look at eIF4A3:



.. report:: NormalisedProfiles.FirstLastExonCount
   :render: r-ggplot
   :slices: r(R.)
   :tracks: GFP
   :transform: pandas
   :tf-statement: query('protein=="eIF4A3"', engine='python')
   :statement: aes(x=exon, y=log2(normed_count), col=slice, group=slice) + geom_point() + geom_line() + facet_wrap(~protein, scale="free_y") + theme_bw() + scale_x_discrete(limits=c("first_exon", "middle_exon", "last_exon", "introns"), labels=c("First", "CDS", "Last", "intron")) + theme(aspect.ratio=1) + ylab("LogFC compared to FlipIn")

   Fraction of eIF4A3 tags in each region, normalized to FlipIn.


Note first of all that there is quite a lot of difference between the
replicates. This suggests that the enrichment over FlipIn varies from
replicate to replicate. This mean that FlipIn might not be a
particulalry good normalizer to captrue absolute amounts of
enrichment.

However, the pattern is *mostly* reflected accross the
replicates. Thus all areas areas of coding sequence are enriched
compared to introns. Last exons, i.e. supposedly after the last EJC
location, are enrichment compared to introns, but have less eIF4A3
than middle exons. Infact that have a similar fraction of the reads as
the controls does. But introns have even less. Interestingly, the
first exon, which should have an EJC on it, looks similar to the last
one. If we look at number, rather than fraction we see that there 
are far more eIF4A3 reads in final exons than in introns compared
to FlipIn.

.. report:: NormalisedProfiles.FirstLastExonCount
   :render: r-ggplot
   :slices: r(R.)
   :tracks: GFP
   :transform: pandas
   :tf-statement: query('protein=="eIF4A3"', engine='python')
   :statement: aes(x=exon, y=log2(unnormed_count), col=slice, group=slice) + geom_point() + geom_line() + facet_wrap(~protein, scale="free_y") + theme_bw() + scale_x_discrete(limits=c("first_exon", "middle_exon", "last_exon", "introns"), labels=c("First", "CDS", "Last", "intron")) + theme(aspect.ratio=1) + ylab("LogFC compared to FlipIn")

   Enrichment of eIF4A3 tags in each region, normalized compared to FlipIn.

Moving on to look at TREX components:


.. report:: NormalisedProfiles.FirstLastExonCount
   :render: r-ggplot
   :slices: r(R.)
   :tracks: FLAG
   :statement: aes(x=exon, y=log2(normed_count), col=slice, group=slice) + geom_point() + geom_line() + facet_wrap(~protein, scale="free_y") + theme_bw() + scale_x_discrete(limits=c("first_exon", "middle_exon", "last_exon", "introns"), labels=c("First", "CDS", "Last", "intron")) + theme(aspect.ratio=1) + ylab("LogFC compared to FlipIn")

   Fraction of TREX tags in each region, normalized to FlipIn.

First Alyref. There is no consensus between the replicates on the
First/Middle/Last exon split, with one replicate having the CDS up,
one having it down, and the final rep having them equal. In all cases
the exons are higher than the introns. 

For ChTop, each replicate shows lower enrichment in the CDS than in 
either the first or last exon. The enrichment is also almost zero in
the introns. The same is mostly true in Nxf1.

None of this means that the introns are only background. We are
looking at the *fraction* of tags in the region, not the number. If we
look at number, as in below, we see that, for all replicates for all
proteins, except one replicate of Nxf1, there are more tags in the
introns for the pull down than for the control.


.. report:: NormalisedProfiles.FirstLastExonCount
   :render: r-ggplot
   :slices: r(R.)
   :tracks: FLAG
   :statement: aes(x=exon, y=log2(unnormed_count), col=slice, group=slice) + geom_point() + geom_line() + facet_wrap(~protein, scale="free_y") + theme_bw() + scale_x_discrete(limits=c("first_exon", "middle_exon", "last_exon", "introns"), labels=c("First", "CDS", "Last", "intron")) + theme(aspect.ratio=1) + ylab("LogFC compared to FlipIn")

   Enrichment of TREX tags in each region, compared to FlipIn.


Conclusions
-----------

* There are fewer eIF4A3 tags in the last exon than in the middle
  exons, but not none.
* There is not much to be said about the split of Alyref across coding
  exons, only that there is more than introns.
* ChTop has some enrichment in first and last exons compared to middle
  exons.

Why is there not more Alref in first exons? One possibility is that
there is an enrichment at the 5' end of the middle exons, but that
this is not present in the first exon. This would explain the metagene
plots. Alternatively our annotations of first exons are not good.

Investigations continue. 
