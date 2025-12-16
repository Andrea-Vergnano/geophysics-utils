Hello!

This folder collects the scripts I wrote about mesh generation for finite element processing software. I mainly use GMSH, an open source software available here https://gmsh.info/


Scripts in this folder:

1) CaveGMesher3D: This R script is meant to produce a GMSH software script (.geo), for creating the mesh of a cave or mine or tunnel with an elliptical shape, which semiaxes length is allowed to vary in space. This mesh is optimized for later use in ResIPy, an open source software for the inversion of Electrical Resistivity Tomography surveys.