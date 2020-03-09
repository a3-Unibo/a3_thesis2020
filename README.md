# a3_thesis2020
  
Codes used in the Thesis sessions 2020.  
  
Please note that this material is a compendium to in-person teaching workshop days, so many implied instructions, premises and cautions given during the dat-to-day development have not been included. - tutor: [Alessio Erioli](https://www.unibo.it/sitoweb/alessio.erioli/)  
  
Required programs: [Rhinoceros 3D v6](https://www.rhino3d.com) (includes Grasshopper), [Visual Studio Community 2019](https://visualstudio.microsoft.com/vs/) (install **.NET desktop development** workload).
  
## Installations
These are the required installations:

### Rhino
. [topmostviewport](https://www.food4rhino.com/app/topmost-viewport)
  
### Components
. [Anemone](https://www.food4rhino.com/app/anemone)  
. [Fattener](https://discourse.mcneel.com/t/skeleton-fattener-mesh-cage-morph/74766) - requires forum registration  
. [Human](https://www.food4rhino.com/app/human)  
. [Impala](https://www.food4rhino.com/app/impala)  
. [Meshedit](https://www.food4rhino.com/app/meshedit)  
. [MeshTools](https://www.grasshopper3d.com/forum/topics/mesh-pipe) - requires forum registration  
. [Plankton](https://github.com/meshmash/Plankton) - read the "Using Plankton with Grasshopper" paragraph  
  
### User Objects
froGH WIP - get froGH 1.3 ALPHA from [https://github.com/Co-de-iT/froGH/releases](https://github.com/Co-de-iT/froGH/releases)
extract the froGH folder from the zip file and place it in the Grasshopper User Objects folder: %AppData%\Grasshopper\UserObjects
  
### dll Libraries
. [Noises.dll](https://drive.google.com/open?id=1UXI6UHtCaZFw4csWIghDHwlBSObBp31G) - library with Simplex Noise generation functions, it can be used to embed Noise calculations (including Curl Noise, which is based on Simplex Noise) in a custom C# script.  
Extract the .dll from the zip file and place it in the Grasshopper Libraries folder: %AppData%\Grasshopper\Libraries

---

## @ utilities

Contains .gha assemblies and .dll libraries and general purpose .gh definitions used in the workshop.

**3Dpeople_20200120** - 3D people as meshes in 3 different resolutions  
  
**bitmaps in CS** - generate, load and save bitmaps from within C# components
  
**Curl noise.gh** - uses Noises.dll library to compute SimplexNoise and CurlNoise functions  
  
**Custom material preview.gh**  
**Custom material preview.3dm**  
Use custom-generated materials in Grasshopper preview (works for Rhino rendered modes)  
  
**M00_Millipede FEM field.gh** - simple use of Millipede Grasshopper plugin to generate a scalar and vector field of structural information over a FEM model of a mesh surface  
**M01_Millipede graphics generator.gh** - generates and bakes geometry for 3 different diagrams of Millipede generated data  
*Millipede_data.ghdata* - this file is a sample of how data is passed between M00 and M01  
**interpolate mesh data.gh** - interpolate scalar and vector data while performing a Catmull-Clark subdivision of a mesh - sometimes Millipede can be slow on big geometries. This definition allows the use of a lower-resolution mesh for faster analysis and interpolate data to use on a high-res mesh  
  
**Symmetric noise - displace.gh** - Computes Noise functions in a symmetric way along meshes and generates a displaced mesh accordingly  
  
**Util_Clipping plane - Turntable base.3dm**  
**Util-01_clipping plane anim.gh**  
**Util-02_turntable.gh**  
**Util-03_view capture.gh**  
**Util-04_Anemone run, turntable and capture**  
These files are helpers to generate, respectively: an animation of a moving clipping plane (for a model tomography), a turntable of one or more geometries, capture a Rhino viewport from Grasshopper, and record an animation of an Anemone loop as it executes  
  
**base meshes.gh** - reference mesh models that can be use in exercises  
<br>

### @ utilities/Display Modes
Contains a bunch of customized Display Modes for Rhino 6 - they can be installed in Rhino from:  
_Tools > Options > View > Display Modes > Import_  


### @ utilities/Mesh Modeling
Rhino files and Grasshopper definitions for basic Mesh modeling (low poly to subdivision techniques)  

---
## codes

These folders contain the codes, organized as follows:
  
**GH_<something>** - all things Grasshopper-focused: intuition and C# introductory codes  
**VS_Codes** - all codes developed for complex strategies with Visual Studio as IDE <coming soon>  
  
### GH_CSharp
This folder contains all the Grasshopper definitions with a progressive introduction to C#.
  
**CS_00_intro.gh** - introduction to C# programming in Grasshopper  
**CS_01_data 01.gh** - data types in C# - part 1  
**CS_02_data 02.gh** - data types in C# - part 2 - loops and conditional statements  
**CS_03_functions.gh** - functions in C#  
**CS_04_classes.gh** - classes and objects in C#  
**CS_05_0_gradient descent.gh** - gradient descent example in C#  
**CS_05_1_gradient descent - erosion.gh** - gradient descent erosion example in C#  
**CS_06_delegates example.gh** - explanation of delegates, anonymous functions and lambda syntax in C#  
**CS_07_RTree point search.gh** - using RTree data structure in C# - simple example of nearest neighbours search  
