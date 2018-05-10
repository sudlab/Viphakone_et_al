Chtop and CFIm68 (CPSF6) at PolyA sites
========================================

Nico asked me to look at the relationship between CPSF6 and ChTop around PolyA sites as he thought he had observed a pattern. Specifically that Chtop was co-inciding with the CPSF6 signal at upstream PolyA sites, but not at the downstream ones where they seemed to occupy mutally exclusive spaces.

We have already seen from the exon length analysis that to an extend the profile of CPSF6 looks like that of Chtop, but this is only if you only include last exons.

Here I have created metagenes in windows 200bp upstream and 100bp downstream of Paperclip (essentially iCLIP for PABP) which identifies polyA sites. I have selected genes that have exactly 2 polyA sites so that the suppressed site can be compared to the none suppressed site. I selected genes where the upstrema site was less than half as strong as the downstream site becuase it would allow a clear distinction between sites that were being used and site that were being suppressed. Thus I selected sites where the distil polyA had at least twice the score of the proximal one.

These two filters left approximately 1000 genes to work with.

In the plot below, note that the distil site is on the left and the proximal on the right. The dashed line indicates the location of the paper clip site. Note that while the CPSF6 signal is in HEK293 cells it is PAR-CLIP rather than iCLIP.

.. report:: Paperclip.UnNormedAPAMetagenes
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=base, y=density, col=track) + facet_grid(track~pA,scale="free_y") + geom_line() + scale_x_continuous(breaks=c(0,10,20,30), labels=c(-200,-100,0,100), name ="Position relative to poly-A") + geom_vline(xintercept=20, lty=2, alpha=0.5) + theme_bw() +guides(color=FALSE)

   Per-gene normed metagene around weak 5' and strong 3' alternative polyA sites


The long and the short of it is that CPSF6 signal peaks just before the distal and proximal polyA sites. Chtop drops rapid just before the distil site, as expected because this is the end of the transcript. At the proximal polyA there is a pretty flat profile, in keeping with the lower usage of this polyA site. There is some suggestion that there might be a dip at or just upstream of the polyA site. This could be due to the footprint of other complexes.

I don't think this provides much in the way of evidence for the position or profile of Chtop being important at lower usage proixmal polyA sites.

						       
