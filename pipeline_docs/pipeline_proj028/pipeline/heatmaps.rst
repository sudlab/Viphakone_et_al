.. _pipeline-heatmaps:

Heatmaps
=========

Start aligned, forward reads
-----------------------------

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(forward)


Start aligned, reverse reads
-----------------------------

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(reverse)


End aligned, gene length sorted
----------------------------------

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(FLAG.+end_aligned.length)


End aligned, 3' UTR length sorted
------------------------------------

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: column-2
   :tracks: r(FLAG.+end_aligned.3utr)


RNA Seq depth
--------------

.. report:: Heatmaps.TrackerHeatmaps
   :render: gallery-plot
   :layout: row
   :tracks: r(star)


Normalised
-----------

.. report:: Heatmaps.NormedHeatmaps
   :render: gallery-plot
   :layout: column-3

   RNA normalised heatmaps

First Exon only
---------------

.. report:: Heatmaps.FirstExonHeatmaps
   :render: gallery-plot
   :layout: column-2

   Unnormalised heatmaps of start aligned first exons


.. report:: Heatmaps.NormedFirstExonHeatmaps
   :render: gallery-plot
   :layout: column-2

   Normalised heatmaps of start aligned first exons
