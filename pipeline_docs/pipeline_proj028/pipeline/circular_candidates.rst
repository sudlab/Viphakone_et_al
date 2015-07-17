Candidatas for cricularisation via Alyref and Chtop
====================================================

Clip sites in the 3' and 5' of transcripts were counted, and enrichments based on utr length caluclated for each transcript. And score calculated as follows:

#. Counts and enrichments were averaged across replicates, excluding replicate 1, and any replicate of a transcript were the total clip count accross the transcirpt was zero.
#. Counts and enrichments were ranked and log transformed
#. Alyref score was calculated as the sum of the logranks for 5' utr enrichment and 5' utr clip site count.
#. Chopt score was calculated as the sum of the logranks for the 3' utr enrichment and 3' utr clip count.
#. Overall score was calculated as the sum of the Alyref score and the Chtop score.


The top twenty transcripts are shown below

.. report:: CircularCandidates.CircularCandidateTable
   :render: table

   Transcript scores for being bound by ALyref at 5' and Chtop at 3' ends

And the clip sites in the region around these genes:


.. report:: CircularCandidates.CircularCandidates
   :display: png,png,200
   :render: gallery-plot

   Plots in the region of top ten genes


