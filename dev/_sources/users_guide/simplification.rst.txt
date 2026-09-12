Simplification of Geometry
===========================

In most applications of MC codes, the starting point for the geometry is an engineering CAD model that must be simplified for radiation transport.
This preparation and simplification step is critical in the CSG approach, in order to obtain models that can be handled by conversion tools such as GEOUNED.

Ideally, the objective of the simplification process is to obtain a model that is as simple as possible while remaining representative in terms of radiation transport physics.
This means that any detail with no impact on the radiation transport should be removed. The resulting model is usually called the neutronic CAD model.

In reality, models are frequently not clean enough: there is often interference between solids, and some solids relevant for radiation transport are not defined (e.g., fluids inside pipes).
These models should be cleaned and prepared prior to the simplification process.

Moreover, in the case of a Monte Carlo Constructive Solid Geometry (MC-CSG) approach, such as the one followed by GEOUNED, there are geometry limitations in the geometry definition with respect to the B-rep approach followed by CAD tools.
These limitations are directly related to the types of surfaces that can be used.
In most codes, this means that the types of surfaces that can be used are limited to planes, second-order surfaces defined by quadratic equations, and a restricted kind of torus (those whose axis of revolution is parallel to the X, Y or Z axis).
In this sense, spline surfaces, typically used in B-rep approach, are not allowed in MC-CSG geometry and must be replaced by a set of allowed surfaces.

Starting from the engineering CAD model, the usual steps followed are:

1. Cleaning of the CAD model:

In this step, the details that are not relevant for radiation transport are removed. What is non-relevant depends on the objective of the analysis to be performed and should be decided by the analyst.
For example, in problems where radiation attenuation is relevant, the amount of material and the streaming paths should be preserved.

2. Removing non-allowed surfaces:

As mentioned above, when splines or other disallowed surfaces are present in the CAD model, they must be replaced by allowed surfaces. This process can sometimes be very tedious and time-consuming for the analyst.
Once this step is finished, a neutronic model ready to be converted should be available.
The simplification carried out up to this point can change the amount of material with respect to the engineering CAD model, so density correction factors may be required.

3. Make the model digestible for the conversion tool:

Some solids can be too complex for the conversion process. It is usually necessary, and recommended, to simplify these solids by splitting them.
These solids can be detected in the conversion tool because they take too long to convert or cause the code to crash directly.

4. Conversion and review:

Once the model has been converted, the resulting solids should be checked. In complex models, this is usually done by visual inspection and volume control.
In this sense, GEOUNED provides a way to perform volume control directly by producing a radiation source suitable for stochastic volume calculation, along with the corresponding tally.
