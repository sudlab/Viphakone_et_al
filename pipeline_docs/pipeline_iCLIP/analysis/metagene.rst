.. _metagene:

Metagene profiles
===================

Meta-gene profile show how reads are distributed across an "average" gene. Thus they can show use whether reads tend to accumulate at the start or end of genes, upstream or downstream, or tend to be found in the introns. 

In a standard metagene profile transformations are applied to equalise the length of gene, so the x-axis represents the percent of the gene length rather than a gene length.

Profiles over whole gene models
--------------------------------

.. report:: Profiles.GeneProfiles
   :render: gallery-plot
   :glob: gene_profiles.dir/*geneprofilewithintrons.detail.png
   :width: 200
   :layout: column-3
   :groupby: none

   Gene Profiles


There are several things to be noted here. Firstly, for Alyref and Chtop, reads tend to map to the exons rather than the introns, where as the reads found in the control FlipIn libraries tend to be found in the introns. Reads for Nxf1 are equally likely to be found in introns as exons. Secondly for Alyref and Chtop there are enrichments at the 3' and 5' ends of the genes respectively. The profiles for Nxf1 and FlipIn are noiser due to the lower number of reads.

