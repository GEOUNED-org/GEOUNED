Material Allocation
===================

The material of the components can be assigned directly in the tree CAD by including a text with _mXX_ within the component name.

.. image:: /images/mat_tag.png
   :alt: Example of a CAD tree with a _mXX_ tag.
   :align: center
   :width: 50%

The material assigned to the corresponding cell should be defined in the material definition file whose name is provided by the matFile keyword.
This file has the following format:

.. image:: /images/mat_list.png
   :alt: Example of a material list.
   :align: center
   :width: 50%

Where # is used for comments and the rest of the lines specifies the id of the material, the nominal density and the text to be included as a comment.
The nominal density value is multiplied by -1 so the criterium is inverted with respect to MCNP one (i.e. g/cm3 positive and atm/b/cm negative).
The nominal density can be changed by the inclusion of a multiplication factor using the following text in the CAD tree _mXX_dXX.XX_.
For example, in our case _m01_d2.0_ will produce a cell with material 1 and density 2.0.
