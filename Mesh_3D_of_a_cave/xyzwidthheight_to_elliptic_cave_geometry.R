#This R script was written by Andrea Vergnano, 2024, for geophysical research purposes. It successfully runs in R version 4.3.2, and RStudio build 494.
#If you want to reuse it for any purpose, please cite the scientific paper of my research group about the Pugnetto cave (Italy).
#If you have any question, feel free to contact me on Github or at my email, currently andrea.vergnano@unito.it.

### WHAT IS THIS SCRIPT FOR? ###

#This R script is meant to produce a GMSH software script (.geo), for creating the mesh of a cave or mine or tunnel with an elliptical shape, which semiaxes length is allowed to vary in space.
#This mesh is optimized for later use in ResIPy, an open source software for the inversion of Electrical Resistivity Tomography surveys.

#This script takes as input a .csv file with 5 columns, named: X, Y, Z, width, height.
#Each row represents a point on the floor of the cave (which, in our use, corresponds to the positions of the electrodes in an ERT survey)
#X,Y.Z are the UTM positions of those points, whereas width and height are the dimensions of the cave in that point.

#The output ".geo" file can be opened by Gmsh open source software. To produce the mesh, in gmsh, just click on Mesh, 3D (then, optionally, Refine by splitting to have a finer mesh), Save mesh.

#The mesh (.msh file) produced in Gmsh can then be imported in ResIPy software.

### INPUT PARAMETERS ###

### !!! CHECK THEM !!! ### ### CHECK ALSO THE END OF THE SCRIPT - WHERE MESH SIZE PARAMETERS ARE GIVEN TO GMSH ###

#import csv file. put here the name of the file. No filepath is needed as long as it is in your R working directory.
XYZD <- read.csv("XYZ_height_width_electrodes_8march2024.csv")

#select the name of the output file. warning: if already present in the folder, it will be overwritten!
output_name="cave_geometry.geo"

#select the investigation depth, to create a volume that extends these meters up and down the max and min altitude of the input points:
inv_depth=30 #m

### END INPUT PARAMETERS ### ### CHECK ALSO THE END OF THE SCRIPT - WHERE MESH SIZE PARAMETERS ARE GIVEN TO GMSH ###




#write first line of the output file, in order to select OpenCASCADE kernel in gmsh.
cat(sprintf("SetFactory(\"OpenCASCADE\");\n"),file = output_name, append = FALSE)


#calculate the direction versors at each electrodes
for (j in 2:(length(XYZD$X)-1)){
  XYZD$dirX[j]=XYZD$X[j]-XYZD$X[j-1]
  XYZD$dirY[j]=XYZD$Y[j]-XYZD$Y[j-1]
  XYZD$dirZ[j]=XYZD$Z[j]-XYZD$Z[j-1]
  XYZD$angleZ[j]=atan(  XYZD$dirZ[j]/(XYZD$dirX[j]^2+XYZD$dirY[j]^2)^0.5)
  
}

##first and last position equal to second and secondfromlast.
XYZD$dirX[1]=XYZD$dirX[2]
XYZD$dirY[1]=XYZD$dirY[2]
XYZD$dirZ[1]=XYZD$dirZ[2]
XYZD$angleZ[1]=XYZD$angleZ[2]


XYZD$dirX[length(XYZD$X)]=XYZD$dirX[j]
XYZD$dirY[length(XYZD$X)]=XYZD$dirY[j]
XYZD$dirZ[length(XYZD$X)]=XYZD$dirZ[j]
XYZD$angleZ[length(XYZD$X)]=XYZD$angleZ[j]



#smooth direction versors
XYZD$dirX=smooth(XYZD$dirX)
XYZD$dirY=smooth(XYZD$dirY)
XYZD$dirZ=smooth(XYZD$dirZ)
XYZD$angleZ=smooth(XYZD$angleZ)


###For each point, create an ellipse with the axis parallel to the longitudinal direction of the cave: ###


#calculate c, the geometrical parameter useful for the subsequent for cycle, and derived from simple geometrical considerations. We basically want to create ellipses which are perpendicular to the direction of the line connecting the electrodes.
#Each ellipse have 5 important points: the center (which is of coordinates X,Y,Z+height/2), the lower point (X,Y,Z), the upper point (X,Y,Z+height), the left point (Xl,Yl,Z+height/2), and the right point(Xr,Yr, Z+height/2)
#Finding the left and right points require simple geometrical consideration. We know that Xr-X is proportional to the versor dirY calculated above (Xr-X=dirY*c), whereas Yr-Y is proportional to the versor dirX ( Yr-Y=dirX*c). 
#We also know that the square root of ((Xr-X)^2+(Yr-Y)^2) is equal to width/2, so we can solve the three equation for c, and find the formula of c written in the following:
#A drawing helps to understand.

XYZD$c=smooth(XYZD$width/(2*(XYZD$dirX^2+XYZD$dirY^2)^0.5))


#This for cycle writes into the output file the lines regarding the ellipses
 for (i in  1:length(XYZD$X)) {

     counter=i*10+1000 ###!!! check this number, if you have many points and the count goes too up and something overlaps.  
  print(i)
  ### computes the position of 5 points: the center and the top, low, right and left points of the ellipse, in a way they are perpendicular to the direction of the electrodes.
  cat(sprintf("\nPoint(%i) = {%.01f,%.01f,%.01f};  Point(%i) = {%.01f,%.01f,%.01f};  Point(%i) = {%.01f,%.01f,%.01f};  Point(%i) = {%.01f,%.01f,%.01f};  Point(%i) = {%.01f,%.01f,%.01f};",counter,   XYZD$X[i]+XYZD$dirY[i]*XYZD$c[i],XYZD$Y[i]-XYZD$dirX[i]*XYZD$c[i],XYZD$Z[i]+XYZD$height[i]/2, counter+1, XYZD$X[i],XYZD$Y[i],XYZD$Z[i]+XYZD$height[i], counter+2,   XYZD$X[i]-XYZD$dirY[i]*XYZD$c[i],XYZD$Y[i]+XYZD$dirX[i]*XYZD$c[i],XYZD$Z[i]+XYZD$height[i]/2, counter+3, XYZD$X[i],XYZD$Y[i],XYZD$Z[i], counter+4, XYZD$X[i],XYZD$Y[i],XYZD$Z[i]+XYZD$height[i]/2), file = output_name, append = TRUE)
  
  ## computes the four 90-degrees ellipse arcs 
  cat(sprintf("\nEllipse(%i) = {%i,%i,%i};  Ellipse(%i) = {%i,%i,%i};  Ellipse(%i) = {%i,%i,%i};  Ellipse(%i) = {%i,%i,%i};",counter,counter, counter+4,counter+1   , counter+1, counter+1,counter+4, counter+2,   counter+2, counter+2, counter+4,counter+3,     counter+3, counter+3,counter+4,counter  ), file = output_name, append = TRUE)
  
  ## put together the ellipse arcs into curve loops, which are then needed by the Ruled Thrusections gmsh command to create a volume.
  cat(sprintf("\nCurve Loop(%i) ={%i,%i,%i,%i};\n",i,counter,counter+1,counter+2,counter+3), file = output_name, append = TRUE)
  
  }

#The Ruled ThruSections is a gmsh command to produce a volume from a series of Curve Loops
cat(sprintf("\n\nRuled ThruSections(2) = {1:%i};",i),file = output_name, append = TRUE)




#draw the box outside the tunnel

maxX=max(XYZD$X)+inv_depth-10
minX=min(XYZD$X)-inv_depth+10
maxY=max(XYZD$Y)+inv_depth-10
minY=min(XYZD$Y)-inv_depth+10
maxZ=max(XYZD$Z)+inv_depth
minZ=min(XYZD$Z)-inv_depth


#Write the gmsh script line regarding the box, which will be a volume object in gmsh
cat(sprintf("\n\nBox(1) = {%.01f, %.01f, %.01f, %.01f, %.01f, %.01f};",minX,minY,minZ,maxX-minX,maxY-minY,maxZ-minZ),file = output_name, append = TRUE)

#Compute the boolean difference between the box and the pipe, to remove the cave volume from the box volume.
cat(sprintf("\n\nBooleanDifference{ Volume{1}; Delete; }{ Volume{2};Delete;}"),file = output_name, append = TRUE)#1 is the name of the cave volume, which was defined as volume 1.


#Create Physical volumes, useful for later processing (e.g. in ResIPy, to recognize the different volumes and possibility to assign them properties such as initial resistivity...)
cat(sprintf("\n\nPhysical Volume(\"box\", 1) = {1};\n"),file = output_name, append = TRUE)


# write on the .geo file also the part relative to the mesh size field.

### ALL OF THEM ARE MESH PARAMETER YOU CAN CHANGE ACCORDING TO YOUR PROBLEM! ###

cat(sprintf("\nField[1] = Distance;

Field[1].SurfacesList ={1:190};\n
Field[1].Sampling=20;\n

Field[2] = Threshold;\n
Field[2].InField = 1;\n
Field[2].SizeMin = 1;\n
Field[2].SizeMax = 8;\n
Field[2].DistMin = 3;\n
Field[2].DistMax = 30;\n
Mesh.Algorithm = 6; //5 handles better complex fields but mesh is worse.\n

Mesh.MeshSizeFromPoints = 0;\n
Mesh.MeshSizeFromCurvature = 0;\n
Mesh.MeshSizeExtendFromBoundary = 0;\n
Mesh.MeshSizeMax=8;\n

Field[2].StopAtDistMax = 1;\n

Background Field = 2;"),file = output_name, append = TRUE)

### END OF THE SCRIPT ###

### HAVE A NICE DAY ###