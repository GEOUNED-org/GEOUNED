Universe Fill and Transform
===========================

Components within a CAD geometry can be specified to be filled with a given universe number by assigned directly in the tree CAD by including text with _fXX_ within the component name. Where XX is the universe number to be used to fill the cell.
A transform number can also be specified during by the additional of text with _tXX_ within the component (or a subcomponent) name. Where XX is the transform number to be applied to the filled universe. It is important to note that the transform number must be defined in the MCNP input file using the appropriate transform cards.

.. image:: /images/fill_tag.png
   :alt: Example of a CAD tree with a _fXX_ and _tXX_ tags. Bodies within blanket_level_04_f3_ are filled with universe 3; the universe filling bodies within 014_t1_ is transformed by transform 1, 024_t2_ is transformed by transform 2, 034_t3_ is transformed by transform 3 and 044_t4_ is not transformed.
   :align: center
   :width: 50%

The figure above shows an example of a CAD tree with _fXX_ and _tXX_ tags. Bodies within blanket_level_04_f3_ are filled with universe 3. The universe filling bodies within 014_t1_ is transformed by transform 1, 024_t2_ is transformed by transform 2, 034_t3_ is transformed by transform 3 and 044_t4_ is not transformed.
