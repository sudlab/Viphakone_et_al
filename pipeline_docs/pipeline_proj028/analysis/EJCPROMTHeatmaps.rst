PROMPT Heatmaps for eIF4A3 and PTB.
======================================

In order to see if the apparent binding of Alyref to PROMPT sequences upstream of TSSs is specific or not, the analysis was repeated with the eIF4A3 and PTB data (and other EJC components. Full results avaiable on the 'Heatmaps' pipeline page).

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: Alyref-FLAG.union.reverse,eIF4A3-GFP.union.reverse,PTB-GFP.union.reverse

   Reverse strand reads aligned at start of genes


We can see that there is signal in both eIF4A3 and PTB. It is not as strong as for Alyref for eIF4A3, but is is stronger for PTB.

This suggests that whatever is causing the signal in Alyref is also causing signal in eIF4A3 and PTB as well. It has been suggested that this is "genomic background". I don't think this can be the case, because genomic background would not be specific to -ve strand reads and would not be consentrated at TSSs fading and fade with distance from the TSS. Thus, a more likely explaination is non-specific RNA binding, and that these factors are susceptable to non-specific binding, just as Alyref is.
