Metagenes around paperclip sites
================================

Single polyA windows
--------------------

I got 200bp upstream and 100bp downstream around paperclip sites where a flatten gene contains only one poly-A site.

.. report:: Paperclip.SinglePolyAMetagenes
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=base, y=density, col=track) + facet_grid(track~.,scale="free_y") + geom_line() + scale_x_continuous(breaks=c(0,20,40,60), labels=c(-200,-100,0,100), name ="Position relative to poly-A") + geom_vline(xintercept=40, lty=2, alpha=0.5) + theme_bw() +guides(color=FALSE)

   Metagenes around consitutive poly-A sites.


Windows around alternate polyA sites
------------------------------------

Selected transcripts that overlapped two paperclip polyA sites and where the 3' site was the stronger of the two and created metagenes around them.

Note that the 3prime and 5prime sites are the wrong way around on the below because 3 comes before 5 in an alpha numeric sort and its not worth the effort to fix.

Dotted line is paperclip site.

.. report:: Paperclip.NormedAPAMetagenes
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=base, y=density, col=track) + facet_grid(track~pA,scale="free_y") + geom_line() + scale_x_continuous(breaks=c(0,10,20,30), labels=c(-200,-100,0,100), name ="Position relative to poly-A") + geom_vline(xintercept=20, lty=2, alpha=0.5) + theme_bw() +guides(color=FALSE)

   Per-gene normed metagene around weak 5' and strong 3' alternative polyA sites


To look at the effect of the normalization, I redid this but without the per gene normalization.


.. report:: Paperclip.UnNormedAPAMetagenes
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=base, y=density, col=track) + facet_grid(track~pA,scale="free_y") + geom_line() + scale_x_continuous(breaks=c(0,10,20,30), labels=c(-200,-100,0,100), name ="Position relative to poly-A") + geom_vline(xintercept=20, lty=2, alpha=0.5) + theme_bw() +guides(color=FALSE)

   Per-gene normed metagene around weak 5' and strong 3' alternative polyA sites


Comparing effects of ChTop knockdown and CPSF6 knockdown
=========================================================

.. report:: Paperclip.AseqVsDapars
   :render: table
   :transform: pandas
   :tf-statement: query("padj < 0.05 and PDUI_Group_diff < -0.25 and Cntrl-siCPSF6 < -0.25", engine="python")

   Genes which show proximal APA on both ChTop and CPSF6 knockdown


.. report:: Paperclip.AseqVsDaparsFishers
   :render: table
   

   Fishers test on association between significant APA on chtop KD and direction of change on CPSF6 KD


.. report:: Paperclip.AseqVsDapars
   :render: r-ggplot
   :statement: aes(x=Cntrl-siCPSF6, y=PDUI_Group_diff) + geom_point(size=0.5, alpha=0.25, pch=20) + geom_smooth(method="lm") + coord_fixed() + theme_bw(base_size=10) + xlab("Change in distil usage on CPSF6 Knockdown") + ylab("Change is distil usage on Chtop knockdown")

   Correlation of change in distil usage on Chtop knock down, compared to CPSF6 knockdown.


.. report:: Paperclip.AseqVsDapars
   :render: r-ggplot
   :statement: aes(x=Cntrl-siCPSF6, y=PDUI_Group_diff*ifelse(padj<0.05,1,NA), col=Cntrl-siCPSF6 < -0.25 & PDUI_Group_diff < -2.5) + geom_point(size=0.5, alpha=0.25, pch=20) + geom_smooth(method="lm", group=1) + coord_fixed() + theme_bw(base_size=10) + xlab("Change in distil usage on CPSF6 Knockdown") + ylab("Change is distil usage on Chtop knockdown") + scale_col_discrete(name="Selected for follow-up") + theme(legend.position="bottom")

   Correlation of change in distil usage on Chtop knock down, compared to CPSF6 knockdown, showing only genes that with significant shift to proximal usage on Chtop knockdown
   

last exon Metagenes of CPSF6 around clusters
=============================================


.. report:: Paperclip.CPSF_around_clusters
   :render: r-ggplot
   :groupby: all
   :statement: aes(base, count, colour=track) + geom_line() + facet_wrap(~track, scale="free_y") + ylim(c(0,NA)) + geom_vline(xintercept=c(0,35), lty=2, col="grey50") + theme_bw() + scale_x_continuous(labels=c(-100,-50,0,0,50,100), breaks=c(-100,-50,0,32,82,132), name="Position relative to cluster") + guides(color=FALSE)

   Metagenes of CPSF6 PAR-CLIP on last exons around clusters


Last exon metagenes of CPSF6 around clip sites
===============================================

.. report:: Paperclip.CPSF_around_sites
   :render: r-ggplot
   :groupby: all
   :statement: aes(base, count, col=track) + geom_line() + facet_wrap(~track, scale="free_y") + ylim(c(0,NA)) + theme_bw() + guides(color=FALSE) + xlim(-100,100)

   Metagenes of CPSF6 PAR-CLIP around sites of TREX iCLIP in last exons

.. report:: Paperclip.CPSF_around_apa_sites
   :render: r-ggplot
   :groupby: all
   :statement: aes(base, count, ymin=q2_5, ymax=q97_2, col=track) + geom_ribbon(col=NA, fill="grey", alpha=0.5) + geom_line(lwd=1) + facet_wrap(~track, scale="free_y") + ylim(c(0,NA)) + theme_bw() + guides(color=FALSE) + xlab("Distance from clip site") + geom_vline(xintercept=0, lty=2, col="grey50") + scale_color_manual(values=c("FlipIn-FLAG"=rgb(0.35,0.7,0.9), "Chtop-FLAG"=rgb(0,0.6,0.5), "Alyref-FLAG"=rgb(0.8,0.4,0), "Nxf1-FLAG"=rgb(0.8,0.6,0.7)))

   Metagenes of CPSF6 PAR-CLIP around sites of TREX iCLIP in last exons that show apa


