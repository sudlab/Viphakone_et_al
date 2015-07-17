Executive Summary
==================


Wellcome to the report for the project specific pipeline for proj028. This we are now moving on from the production pipeline part of the project common to all iCLIP projects, to those parts specific to this project. Hence the new report. The old report is still availible at the same address as before. So to is an updated via with all the plot updated to incorporate the new data. This is availible `Here <https://www.cgat.org/downloads/N6Cduavf7p/iCLIP_fullrun2/report/analysis.html>`_


In this report we look at the new data. Examine if the patterns we observed in the old dataset are still present, look at enrichment around the STOP codon, look at the relationship between expression and clip coverage,

New data
---------

*   The new data seems to be mostly of a high quality it sucessfully mapped and processed.
*   The new replicate, replicate 4 seems to be of a high quality, at least from an initial look
*   There is a slightly worrying trend in the data that suggests that the UMIs are failing to allow
    us to remove duplicates due to sequencing errors in the UMI itself.
*   Full details at :ref:`mapping`


Alyref and Chtop enrichements
------------------------------

*   The patterns of Alyref being enriched at the three prime end and Chtop at the five prime end are
    preserved in the new data
*   One of the control samples also has an enrichment at the 5' end, which is slightly distrubing.
*   Full details at :ref:`geneprofiles`
*   I also examined if this enrichment was at the STOP codon, as per WTAP and the methylation story.

    *  There is an enrichment at the STOP and START codons for Chtop
    *  This enrichment goes away when you normalise for UTR length or RNA expression
    *  This normalisation wasn't performed on the original methylation data, so might also have an effect there.
    *  Full details at :ref:`utrprofiles`


Correlation with expression levels
-----------------------------------

*   I looked at the corrolation of the number of clip tags on a transcript vs the expression of that transcript
*   There is a good correlation for Alyref and Chtop, a much reduced correlation for Nxf1 and FLAG
*   The reduced correlation looks to be due to a higher number of expressed transcripts with no clip tags.
*   See :ref:`rnaseq`

Calling significant clusters
-----------------------------

*   I have been working on calling significant clusters of clip tags
*   The data from this should be availible soon.

Identifying highly expressed transcripts bound at 5' end by Alfref and 3' end by Chtop
----------------------------------------------------------------------------------------

*   This list should be availible very soon.

