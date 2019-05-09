This morning Nico sent me a detailed email with some thoughts on some of the things we've been thinking about. The questions require thoughtful answers, with plots and discussion, and as such arn't really suitable for emails, so I'm going to write report pages on them, as this allows me to mix text and plots easily and edit things as things get more up to date.

This has taken me all day and I havn't been able to look at the Chtop stuff. But I know that I havn't really done much with this yet, but I will.

First up:

Detained/retained Introns
===========================


Nico asked::

    What is the proportion of transcripts having Alyref clip tags on their
    introns (or a high density of clip tags in their introns if you prefer) that
    have these introns fully excised upon Alyref knockdown in Conrad’s data?

So somethings to bear in mind when answering this question are that in RNAseq data there is no such thing as "fully-excised". What we did was to identify introns where there was a reliable increase in the amount of splicing as measured by a change in the number of reads that map to the intron to those that map to the exons. 

We should also bear in mind that not having evidence of a change is not the same thing as having evidence of no change, and thus those things that are we did not find to have significantly changed may still be changed.

Secondly, there is no clear basis to call things as high density. Our cluster calling allows for calling cluster if there is a spacial irregularity in the tags within the intron, but not if the whole intron is just dense. 

With that said, the below plot shows the relationship between the density of clip tags and the fraction of annotated retained/detained introns that have a signiicant increase in splicing efficiency. The black shows those that have any evidence of increase in efficiency , the red those that show evidence of a great than 1.5 fold increase in splicing efficiency. 

.. report:: RetainedIntrons.RetainedIntronsExpressionVsClips
   :render: r-ggplot
   :groupby: slice
   :slices: Alyref-FLAG.union
   :statement: aes(cut(log2((clip_tags/genomicData_width) + 1e-5),breaks=20)) + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & log2fold_Control_Alyref>0), col=T),fun.y=mean, geom="point") + stat_summary(mapping=aes(y=as.numeric(padj<0.1 & is.finite(padj) & log2fold_Control_Alyref>0.58), col=F),fun.y=mean, geom="point")+ facet_grid(track~., scale="free_y") +theme_bw() + scale_color_manual(values=c("FALSE"="red", "TRUE"="black"), labels=c("FALSE"="FC>1.5", "TRUE"="FC>0"), name="Fold Change") + scale_x_discrete(labels=NULL, name = "log2 (clip tags/genome size") + ylab("Fraction significant increase in splicing efficiency")

   Relationship between clip tag density and probabilty of change in splicing efficiency on Alyref KD

* the first thing to note is that indeed, not all highly bound introns have a change in splicing efficiency on Alyref knockdown. In fact the fraction is quite small. 
* Note also that the fraction is higher in nuclear vs cyto/total. One possiblity for this is that there is less intron retension in the cytoplasm to start with because NMD degrades the retained intron transcripts
*  There is a relationship however between Alyref density and the probablity of there being a change in splicing efficiency: in general the largest change is on those pieces of sequence that have an intermediate density of Alyref. One hypothesis for this is that the RNAi is only a knock down, not a knock out, and so those with the highest density of Alyref still have some Alyref after knockdown, where as those that have intermediate levels are more likely to be distrupted by the change. 
* I don't believe this speaks either way to the question of whether Alyref is loaded at transcription or simply binds to inefficiently spliced "introns".

Nico asked::

    What happens to pre-mRNAs bearing high density of Alyref intronic clip tags
    in Conrad’s data CYTO plus or minus ALYREF RNAi (if there are any in CYTO)? 
    Considering that the method they are using to extract RNA results also in
    traces of pre-mRNAs in  CYTO (I know because we are using it in the lab now),
    I suppose for this that it is the Nuclear/Cyto ratios of these transcripts 
    (spliced and unspliced) instead which should be used.

For an answer for the first part of this question, see above. The pattern in the cyto is the same as for the nuclear or total, but the absolute numbers are smaller. I'm not sure what is meant by the second half of the question. Is the question here does knock down of alyref leads to less export of unspliced transcripts to the cytoplasm? I.e. are there introns, with clip tags, that don't change levels in the nucleus on knock down, but do change levels in the cytoplasm? I.e. the splicing isn't changed, but the export is? I could have a look at this, but my feeling is that it is going to be confounded by degredation in the cytoplasm. 


Nico asked::

    Once normalised to expression levels, are all the detained introns of Sharp’s
    paper transcripts covered with a high density of Alyref clip tags?

A high density compared to what? The following plot examines the density of iCLIP tags, normalised to RNA expression on regions of sequence that are within a transcribed region, but not included in any Ensembl transcripts. I call these constituative introns. These have been broken down into those that overlap with one of Sharpes detained introns, and those that don't.

Note that this is computed using total RNA rather than the PolyA selected rna used by Conrad, so should caputre unprocesed and excised intron lariates as well as processed mature RNA.

.. report:: Expression.DetainedChunks
   :render: r-ggplot
   :tracks: Alyref_FLAG_union
   :slices: total
   :statement: aes(x=(iCLIP+1)/(RNA+1), fill = factor(Constitutive_exon)) + geom_density(data=rframe[rframe$iCLIP + rframe$RNA > 0,],alpha=0.5) + scale_x_log10() + scale_fill_discrete(labels=c("1"="Detained Intron", "0"="Not Detained"), name = "Detained?") + theme_bw() + theme(legend.position="bottom")

   Density of iCLIP sites normalised to expression for detained and non-detained introns.


As you can see, once normalised for expression, all introns, detained or not, have more or less the same density of iCLIP sites once normalised to expression level. Of course the detained introns, on average, have more CLIP tags than the rest because, on average, they are expressed to a higher level. Thus if you select on the density of Alyref iCLIP sites, you are simply indirectly selecting on the expression level. If you don't use expression data, then you will simply catch those introns that are in highly expressed transcripts. 

Indeed, if we look at the different between exons and introns we can see that the density of CLIP tags on introns is actaully slightly higher than the density of CLIP tags on exons. 

.. report:: Expression.TranscriptChunks
   :render: r-ggplot
   :tracks: Alyref_FLAG_union
   :slices: total
   :statement: aes(x=(iCLIP+1)/(RNA+1), fill = factor(Constitutive_exon)) + geom_density(alpha=0.5) + scale_x_log10() + scale_fill_discrete(labels=c("1"="Exon", "0"="Intron"), name = "") + theme_bw() + theme(legend.position="bottom") 

   Expression normalised density of iCLIP tags on introns vs exons

iCLIP for nuclear localised transcripts examines nuclear localised RNA, where as total RNAseq assays RNA level across the whole cell. Thus, this result suggests that we have iCLIP tags on intron RNA that is present in the nucleus but not resident in the cytoplasm.

Nico asked::

    Have you tried to look at whether the transcripts bearing high intronic
    density of Alyref clip tags:
    are functionally related  (GO analysis) ?
    share any common features for the splice sites surrounding these introns
    bound by Alyref (canonical/non-canonical sites, U12-dependent?)
    share features for the introns themselves ? size of the introns (although
    this doesn’t seem to influence the splicing efficiency according to Singh,
    J. & Padgett, R. A. NSMB (2009)), size of exons surrounding them (weak exon
    definition?)?
    show any bias in transcription rates which could have an influence on
    splicing efficiency and therefore on Alyref binding to introns after its
    release by 3’ end processing? 


I havn't. Mostly because there doesn't seem to be a seperate set of introns with a higher density of clip-tags than other introns, once expression has been taken into account. Instead the high numbers of tags found on introns seems to be more or less evenly distributed across introns of highly expressed genes. Of these, the only ones that would be vaguely trivial would be the GO analysis and the transcription rate analysis. This is because of a problem with how the introns are defined. Consider the example below::

    transcript 1:            |>>>>>>>>>>>>>>>>|--------------------|>>>>>>>>>>>>>>>|
    transcript 2:            |>>>>>>>>>>>>>>>>>>>>>>>>>|-----------|>>>>>>>>>>>>>>>|
    exon/intron annotations: |------E1--------|---E2---|----E3-----|------E4-------|

    constituative exons: E1, E4
    constituative introns: E3


So what do we define here as splice sites bounding the intron? Is is the end of E1 or the end of E2? How long is the intron? We could find examples that are a constituative exon followed by a consitutative intron.  These analyses could be done.

Nico suggested::

    So I think that if you could compare transcriptome-wide U2af65 deposition
    on introns and look if there is a correlation with Alyref intronic clip tags,
    that would be interesting and could give a nice explanation for why we have
    those clip tags on introns. 

Are we talking about a correlation between the number of tags on each intron or their location. If it is the former then this is doable. However, I supect the result will be that we either have a corrolation or an anti-corrolation with Alyref. Since Alyref deposition is correlated with expression level, and if U2af65 deposition is also correlated with expression level, we are going to have trouble sorting out the direction of causation. Thus we might find a connect between U2af65 and intron retension that may or may not support the previous studies, but I don't know what we'll probe about Alyref. 

If the latter, then that is harder: its difficult to know what the negative control would be. I.e. how often would be expect them to overlap if they were random. I'll have to think about it. 

Conclusions for the retained intron section
--------------------------------------------

I think the evidence is pretty strong that what determines the number of Alyref clip tags on an intron is the level of intron RNA available. I don't think that the evidence has much to say either way about whether this is due to co-transcriptional loading, or Alyref binding to slowly splicing introns because it is hanging around. In my mind the Alyref knockdown data does suggest that when you KD Alyref, you do get less of these inefficiently spliced introns in the cytoplasm. Whether this is because Alyref is actively proventing splicing, or whether this is because they are no longer exported and so hand around longer in the nucleus, leading to a more oppotunity for splicing is not clear from this data. There may be some more analyses that could be suggestive in this question, but it can probably only be decisively answered by molecular biology expreiments. 

Alyref and EJC
==================

Nico said::

    As I have told you, there is now some iCLIP data for eIF4A3 that you can use
    to match with our iCLIP data and see how our factors get deposited with
    respect to the EJC.

We have looked at the metagene profile around exon-exon boundaries. This is what we see:

.. report:: GeneProfiles.TranscriptomeExonBoundaryProfiles
   :render: r-ggplot
   :groupby: all
   :slices: union
   :statement: aes(x=position, y=density) + stat_summary(fun.y="mean", geom = "line") + facet_wrap(~track, scale="free_y") + theme_bw() + geom_vline(xintercept=c(0,-24), lty=2, lwd=0.5)

   Metagene profile centred on exon-exon junctions


Note that there is a peak around -30bp from the junction, just up stream of EJC. Of course this was also the case in Hauer et al before their correction.  These look quite like the figure from Hauer et al for eIF4A3 and SRSF3

.. image:: http://www.nature.com/ncomms/2015/150811/ncomms8921/images/ncomms8921-f2.jpg


I plan to test their correction on our data at some point. Note also the dip at zero. I wouldn't stake my life on it, but this might suggest that we are not getting binding exactly at the junction: to me this suggests that Alyref is binding the unspliced transcript.

Nico asked::

    What is the metagene profile of Alyref, Chtop, and Nxf1 over intronless genes?
    Is it like over Exons-only section? Is that information embedded in the current
    metagene profiles?

We have look at the metagene profiles from the single exon genes, but there weren't sufficient of them to get any signal above the noise. 

Read-through
================

Nico suggested::

   However, possible read-throughs induced by Alyref kn is worth checking I think.

Looking at read through and inaccurate termination is on my list of things to do. In fact once I've finished looking at how different categories of sequence effect expression normalised binding density, it is probably the next thing I'll look at on this project. But I havn't looked in much detail yet. What I have done is look at simple metagene profiles for RNAseq in knockout and wildtype samples.

.. report:: GeneProfiles.StubbsProfiles
   :render: r-ggplot
   :groupby: all
   :statement: aes(bin,area, col=condition) + geom_line() + facet_grid(fraction~.) + geom_vline(xintercept=c(1000,2000),lty=2) + theme_bw() + scale_x_continuous(breaks=c(500,1500,2500), labels=c("upstream","CDS","downstream"), name = "") + scale_y_continuous(labels=NULL, name="relative read density")

   Metagene profiles from the conrad data


As you can see there is no evidence here for a general lenghening of transcripts (as would be evidenced by an increaed signal in the downstream region). That doesn't mean something won't emerge from a closer examination. 

