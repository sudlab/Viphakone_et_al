Processing Index
================

.. report:: processing_index.ProcessingIndex
   :render: r-ggplot
   :groupby: all
   :statement: aes(x=track, y=processing_index) + geom_bar(stat="identity") + theme_bw(base_size=14) + coord_flip() + ylim(-10,10) + xlab("Sample") + ylab("Processing Index") + facet_grid(factor~.) + geom_hline(yintercept=0)

   Processing Index 
