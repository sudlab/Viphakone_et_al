Expression
===========


Does expression predict the probability of being CLIP'ed?

Expression was divided into bins and the percentage of genes in each expression bin with at least one clip tag was calculated.

Expression was first calculated as the RNA-seq counts: this is a convolution of expression and length:

.. report:: Expression.ProbOfClipByExpression
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=expression+1, y=clip) + geom_point() + facet_wrap(~track) + xlab("Expression decile") + scale_x_continuous(breaks=1:10)

   Effect of RNA-seq count rank on probability of CLIP

Then expression was calculated as sailfish TPMs: this excludes length from the influance

.. report:: Expression.ProbOfClipTPM
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=expression, y=clip) + geom_point() + facet_wrap(~track)+ xlab("Expression decile") + scale_x_continuous(breaks=1:10)

   Effect of RNA-seq TPM rank on probability of CLIP.


Transcript Chuncks
--------------------

Each transcript was divided into non-overlapping chunks, such that each chunk reprsented a unqiue exon part or intron. Number of RNA seq reads or iCLIP reads in each chunk was counted, both for Choromatin associated RNA, or total RNA


.. report:: Expression.TranscriptChunks
   :render: r-ggplot
   :tracks: r(_union)
   :groupby: track
   :statement: aes(RNA+1,iCLIP+1, color = factor(Constitutive_exon), fill = factor(Constitutive_exon)) + stat_binhex(aes(alpha=log2(..count..)),color=NA) + geom_smooth() + scale_x_log10() + scale_y_log10() + theme_bw() + scale_color_discrete(labels=c("1"="Exon","0"="Intron"),name="") + scale_fill_manual( values=c("1"="#FFD500","0"="#002BFF"), labels=c("1"="Exon","0"="Intron"), name = "") + theme(legend.position="bottom", aspect.ratio=1)  + facet_wrap(~slice) + scale_alpha(range=c(0,0.7), guide=F)

   RNA-seq vs iCLIP for each transcript chunk


Looks like you actaully get more iCLIPs per RNAseq in introns than in exons  for Alyref, but this is because iCLIP is nuclear, and the total RNA is cellular, therefore under estimting the amount of intronic RNA? The opposite applies for the chromatin associated RNA.



