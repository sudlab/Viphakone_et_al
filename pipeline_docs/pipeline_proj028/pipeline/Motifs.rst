Motif Searching
===============


Reproducible peaks are those found in 2 or more replicates
Union tracks are peaks found in one or more replicates

Discriminative Searches against FlipIn
--------------------------------------

These analyses use Dreme to find motifs that are more common in the test than in the equivalent FlipIn sample.

Union of all clusters
+++++++++++++++++++++

.. report:: Motifs.AllAgainstFlipIn
   :render: table
   :force:
   :tracks: r(union)
   :large-html-class: sortable
   :groupby: track

   Discriminative searches of motifs against FlipIn samples in the union of all clusters


Reproducible clusters
++++++++++++++++++++++

.. report:: Motifs.AllAgainstFlipIn
   :render: table
   :force:
   :tracks: r(reproducible)
   :large-html-class: sortable
   :groupby: track

   Discriminative searches of motifs against FlipIn samples in reproducible clusters


Indevidual replicates
+++++++++++++++++++++++


.. report:: Motifs.AllAgainstFlipIn
   :render: table
   :force:
   :tracks: r(FLAG-R[0-9]+)
   :large-html-class: sortable
   :groupby: track

   Discriminative searches of motifs against FlipIn samples in reproducible cluster



Location of Motif hits
+++++++++++++++++++++++

.. report:: Motifs.DREMELocations
   :render: table
   :groupby: none
   :separate:

   Location of hits to enriched motifs


Discriminative Searches in retained introns
--------------------------------------------

Clusters within retained introns were tested against FlipIn clusters within retained introns

Union of all clusters
+++++++++++++++++++++

.. report:: Motifs.RetainedIntronDreme
   :render: table
   :force:
   :tracks: r(union)
   :large-html-class: sortable
   :groupby: track

   Discriminative searches of motifs against FlipIn samples in the union of all clusters


Reproducible clusters
++++++++++++++++++++++

.. report:: Motifs.RetainedIntronDreme
   :render: table
   :force:
   :tracks: r(reproducible)
   :large-html-class: sortable
   :groupby: track

   Discriminative searches of motifs against FlipIn samples in reproducible clusters


Indevidual replicates
+++++++++++++++++++++++


.. report:: Motifs.RetainedIntronDreme
   :render: table
   :force:
   :tracks: r(FLAG-R[0-9]+)
   :large-html-class: sortable
   :groupby: track

   Discriminative searches of motifs against FlipIn samples in reproducible clusters
