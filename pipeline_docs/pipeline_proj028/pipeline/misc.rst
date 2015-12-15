Long Expressed genes
-----------------------


.. report:: Expression.LongExpressedGenes
   :render: table
   :force:

   100 longest genes with expression over 10 RPKM



.. report:: Expression.LongCDSExpressedGenes
   :render: table
   :force:

   100 most highly expressed genes with total exonic length over 3.5kb



In the below table genes are ranked by a combination of the expression level (RPKM), the enrichment
of tags in the 3' UTR (fold enrichment over expected) and the average count of tags in the 3' UTR.

.. report:: GeneProfiles.HighlyExpressed3BiasGenes
   :render: table
   :force:

   Highly expressed genes with ChTop 3' bias
