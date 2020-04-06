![Erosion on surface](https://raw.githubusercontent.com/a3-Unibo/a3_thesis2020/master/%40%20media/erosion.jpg)

# a3_thesis2020
  
Codes used in the Thesis sessions 2020.  
  
Please note that this material is a compendium to in-person teaching, so, several implied instructions, premises and cautions given during the dat-to-day development have not been included. - tutor: [Alessio Erioli](https://www.unibo.it/sitoweb/alessio.erioli/)  
  
  
## Installations
These are the required installations:

### Programs
. [Rhinoceros 3D v6](https://www.rhino3d.com) (includes Grasshopper)  
. [Visual Studio Community 2019](https://visualstudio.microsoft.com/vs/) (install **.NET desktop development** workload)  

### Rhino plugins
. [topmostviewport](https://www.food4rhino.com/app/topmost-viewport) - a useful floating topmost window if you have a single screen  
  
### Grasshopper Components
. [Anemone](https://www.food4rhino.com/app/anemone)  
. [Fattener](https://discourse.mcneel.com/t/skeleton-fattener-mesh-cage-morph/74766) - requires forum registration  
. [FileToScript2](https://drive.google.com/open?id=1PZIlEkYBvyzUqHkfgsY78s1Z6T_nf6wI) - if you have an older Rhino SR and the .gha is not compatible you will have to compile it from Visual Studio (instructions are included)  
. [Human](https://www.food4rhino.com/app/human)  
. [Impala](https://www.food4rhino.com/app/impala)  
. [Meshedit](https://www.food4rhino.com/app/meshedit)  
. [MeshTools](https://www.grasshopper3d.com/forum/topics/mesh-pipe) - requires forum registration  
. [Plankton](https://github.com/meshmash/Plankton) - read the "Using Plankton with Grasshopper" paragraph  
  
### Grasshopper User Objects
. [froGH 1.3 ALPHA](https://github.com/Co-de-iT/froGH/releases) -
extract the froGH folder from the zip file and place it in the Grasshopper User Objects folder: _%AppData%\Grasshopper\UserObjects_  
  
### dll Libraries
. [Noises.dll](https://drive.google.com/open?id=1UXI6UHtCaZFw4csWIghDHwlBSObBp31G) - library with Simplex Noise generation functions, it can be used to embed Noise calculations (including Curl Noise, which is based on Simplex Noise) in a custom C# script.  
Extract the .dll from the zip file and place it in the Grasshopper Libraries folder: _%AppData%\Grasshopper\Libraries_  

---

## @ utilities

Contains .gha assemblies and general purpose .gh definitions used in the workshop.

**3Dpeople_20200120.3dm** - 3D people as meshes in 3 different resolutions  
**base meshes.gh** - reference mesh models that can be use in exercises  
**bitmaps in CS.gh** - generate, load and save bitmaps from within C# components  
**Curl noise.gh** - uses Noises.dll library to compute SimplexNoise and CurlNoise functions  
**Custom material preview.gh**  
**Custom material preview.3dm**  
Use custom-generated materials in Grasshopper preview (works for Rhino rendered modes)  
**interpolate mesh data.gh** - interpolate scalar and vector data while performing a Catmull-Clark subdivision of a mesh - sometimes Millipede can be slow on big geometries. This definition allows the use of a lower-resolution mesh for faster analysis and interpolate data to use on a high-res mesh  
**M00_Millipede FEM field.gh** - simple use of Millipede Grasshopper plugin to generate a scalar and vector field of structural information over a FEM model of a mesh surface  
**M01_Millipede graphics generator.gh** - generates and bakes geometry for 3 different diagrams of Millipede generated data  
*Millipede_data.ghdata* - this file is a sample of how data is passed between M00 and M01  
**Symmetric noise - displace.gh** - Computes Noise functions in a symmetric fashion along meshes and generates a displaced mesh accordingly  
**Util_Clipping plane - Turntable base.3dm**  
**Util-01_clipping plane anim.gh**  
**Util-02_turntable.gh**  
**Util-03_view capture.gh**  
**Util-04_Anemone run, turntable and capture.gh**  
These files are helpers to generate, respectively: an animation of a moving clipping plane (for a model tomography), a turntable of one or more geometries, capture a Rhino viewport from Grasshopper, and record an animation of an Anemone loop as it executes  
  
### @ utilities/Display Modes
Contains a bunch of customized Display Modes for Rhino 6 - they can be installed in Rhino from:  
_Tools > Options > View > Display Modes > Import_  
  
### @ utilities/GH reference
Contains reference definitions for the required Grasshopper knowledge to step beyond basics.  
**EXPL_00 - DataTrees.gh** - Data Tree breakdown: structure, commands, navigation - anything you need to know to master Data Trees in Grasshopper.  
  
### @ utilities/Mesh Modeling
Rhino files and Grasshopper definitions for basic Mesh modeling (low poly to subdivision techniques)  
  
---
## codes
  
The codes are organized as follows:  
  
### GH_CSharp
The introduction to C# has been moved to a separated repo: [intro to C# in Grasshopper](https://github.com/Co-de-iT/CS-intro)  
  
### GH_Intuition
This folder contains some introductory iterative strategies developed with standard components + Anemone plug-in for a more intuitive approach
  
**01-00_iterative strategies - intuition.gh** - introduction to iterative strategies in Grasshopper  
**01-00-bis_gradient descent - intuition.gh** - simple, geometry based gradient descent algorithm in Grasshopper  
**01-01_environment and field - intuition.gh** - reading information from an environment/field  
  
### VS_Codes
This folder contains Visual Studio Projects and related .gh files for the following examples:  
  
**AgentSystemFlock** - implementation of classic Craig Reynolds flocking model  
**GradientDescent** - simple gradient descent algorithm on mesh (with mesh erosion)  
**MTSerialization** - advantages of parallelization (of a for loop) and conversion to GH_Types to speed up execution  
**Stigmergy** - example of stigmergy on Mesh, with evaporation and diffusion  
  
---
  
#### copyright
This material is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). ![cc button](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)