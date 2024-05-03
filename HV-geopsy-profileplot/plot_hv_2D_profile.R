### SCRIPT FOR  plotting a 2D hvsr profile taking the topography from a dem. Iinput: a list of hvfiles to be added to the profile, a dem csv xyz columns file.



  
  ### LIBRARIES Needed packages. Install them before running this script. ###
  if (!require("sp")) install.packages("sp")
  if (!require("readr")) install.packages("readr")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("tcltk")) install.packages("tcltk")
  if (!require("akima")) install.packages("akima")
  if (!require("scales")) install.packages("scales")
  if (!require("metR")) install.packages("metR")
  if (!require("kriging")) install.packages("kriging")
  
  library(sp) #required for LongLatToUtm
  library(readr)
  library(ggplot2)
  library(tcltk)
  library(akima)
  library(scales)
  library(metR)
  library(kriging)
  ### FUNCTIONS ###
  
  LongLatToUTM<-function(x,y,zone){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
  }
  

  
  ### INPUT PARAMETERS !!! CHECK THEM!!! ###
  
  recover_coords=TRUE #to recover missing coordinates in the hv files. require the correspondingly named .SAF file in the same folder of the .hv file
  recover_coords_from_file=TRUE #if you have a separate file with the coordinates of each HVSR points, with matching names
  utmzone=38 # for eventual UTM conversion performed when giving external coordinates with recover coords from file.
  Vs0=150 #constant velocity of shear waves in the top layer.
  expseth=0.3 #velocity increase factor with depth
  investigation_depth=100  
  precise_elevation=TRUE #TRUE gets higher-quality elevation profile, needed if the input dem is really coarse, but can increase drastically computation time if input DEM is high-resolution.
  distance_step=5 # x resolution in meters of the final graph
  z_step=2 #z resolution in meters
  char_to_cut=10 #Characters at the end of each file name that you want to cut from the labels in the graph. (example: if namefiles are named all like HV1_final.hv and you want to label them like HV1, this variable should be equal to 9, that is 9 characters at the end of the namefile to be cut)
  n_ite=3 #how much smoothing you want to the surface topography in the graph. Put 0 to have no smoothing. 
   ### END INPUT PARAMETERS ### 
  
  ### SEE GRAPHICS PARAMETERS AT THE END OF THE SCRIPT, BEFORE THE PLOT FUNCTION GGPLOT ###
  
  
  
  
      
  ### INPUT DATA ###
      
      Filtershv<- matrix(c("hv", "*.hv", "R code", ".s",
                           "Text", ".txt", "All files", "*"),
                         4, 2, byrow = TRUE)
      hvdata_list=tk_choose.files(filters=Filtershv)
      hvdata_list=basename(hvdata_list)
      
      DEMfile = tk_choose.files(caption="choose the DEM file which contains the x y z columns")
      
      
      if(recover_coords_from_file==TRUE){
    localization = tk_choose.files(caption="choose the file which contains the x y z information")
    xyz=read.csv(localization)
  }
  x=1
  DEM=read.table(DEMfile, quote="\"", sep=" ",comment.char="")
  
  
  ### END INPUT DATA ###
  
  ### initialize empty dataframe ###
  profile_df=data.frame(x=NULL,y=NULL,z=NULL,hv=NULL)
  
  
  
  ### FOR CYCLE THAT THAT PROCESSES EACH FILE AND STORES EVERY DATA IN A DATAFRAME ###
  
  for (data in hvdata_list){
    
    #read the -hv data
    
    print(data)
    
    #read the .hv data file
    raw_data <- read_delim(data, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
    
    
    #read the location for each data file
   
     f0_file=read_lines(data)
    
    
    #Coordinates from hv file (in UTM):
    
    UTM_index=which(grepl("Position", f0_file, fixed=TRUE))
    UTM_raw=f0_file[UTM_index]
    East=as.numeric(data.frame(a = unlist(strsplit(UTM_raw, split="[ \t]")))[3,1])
    North=as.numeric(data.frame(a = unlist(strsplit(UTM_raw, split="[ \t]")))[4,1])
    
    ### RECOVER MISSING COORDINATES ###
    
    if(recover_coords==TRUE&recover_coords_from_file==FALSE)
    {
      
      if(East==0&North==0){
        print("recovering missing coordinates from .SAF file")
        testo=readLines(sprintf("%s.SAF",substr(data,1,nchar(data)-10)),n=50) # opens the .saf file correspondent to the hv file
        # testo=testo[1:100] #to reduce the search for coordinates in the first lines, to reduce computation cost.
        Latitude_line=testo[grepl("EVT_Y", testo)]
        Latitude= as.numeric(substr(Latitude_line,8,nchar(Latitude_line)))
        Longitude_line=testo[grepl("EVT_X", testo)]
        Longitude= as.numeric(substr(Longitude_line,8,nchar(Longitude_line)))
        UTM=round(LongLatToUTM(Longitude,Latitude,utmzone),2)
        East=UTM[[2]]
        North=UTM[[3]]
        
      }
    }
    
    if(recover_coords_from_file==TRUE)
    {
      
      if(East==0&North==0){
        
        print("recovering missing coordinates from .csv file")
        
        ### extract geolocalization, based on the first 5 chars of the string, to avoid mismatch between, for example, points called L1_01B and L1_01, which have actually the same position. ###
        
        x=xyz$East_wgs84_utm38S[substr(xyz$TEM,1,5)==substr(data,1,5)]
        y=xyz$North_wgs84_utm38S[substr(xyz$TEM,1,5)==substr(data,1,5)]
        
        Longitude=max(x)
        Latitude=max(y)# to delete eventual repeated values of x y
        
        if((Latitude+Longitude)<200){# a simple recognition if we are in cartesian or geographical coordinates.
          UTM=round(LongLatToUTM(Longitude,Latitude,utmzone),2)
          East=UTM[[2]]
          North=UTM[[3]]
        }else{
          East=Longitude
          North=Latitude
        }
        
        
      }
    }
    
    ### END RECOVER MISSING COORDINATES ###
    
    

    #Calculate altitude from DEM based on minimum distance
    
    Altitude=DEM$V3[((DEM$V2-North)^2+(DEM$V1-East)^2)^0.5==min((DEM$V2-North)^2+(DEM$V1-East)^2)^0.5]
  
    #Better: calculate altitude from DEM based on distance from the nearest 4 points

    if(precise_elevation==TRUE){
      DEM$rownumber=as.numeric(row.names(DEM))
      DEM=DEM[order(((DEM$V2-North)^2+(DEM$V1-East)^2)^0.5),]
      Altitude=weighted.mean(DEM$V3[1:16],1/((DEM$V2[1:16]-North)^2+(DEM$V1[1:4]-East)^2))
      
      }
    #convert each f0 data to altitude
    
    # #iteration to solve Vs=V0(1+Z)^a
    # 
    # depth_firsttempt=(Vs0/(4*raw_data$X1))
    # Vs_firsttempt=Vs0*(1+depth_firsttempt)^expseth
    # dif=2
    # Vs_iter=Vs_firsttempt
    # while(dif>1){
    # Vs_prior=Vs_iter 
    # Vs_iter_meanup=Vs_iter
    # depth_iter=(Vs_iter/(4*raw_data$X1))
    # Vs_iter=Vs0*(1+depth_iter)^expseth
    # Vs_iter=Vs_iter*Vs0/(Vs_iter[length(Vs_iter)])#normalize to make the vs at 0 = Vs0
    # 
    # dif=sum(Vs_iter-Vs_prior)
    #     }
    #raw_data$X5=Altitude-(Vs0/(4*raw_data$X1)) #estimates the altitude of each hvsr frequency, based on the Vs/4f0 formula
    # raw_data$X5=Altitude-depth_iter #estimates the altitude of each hvsr frequency, based on the Vs/4f0 formula
    # 
    depth_Seth=((((Vs0*(1-expseth))/(4*raw_data$X1))+1)^(1/(1-expseth)))-1
    raw_data$X5=Altitude-depth_Seth
    #Add the data to a dataframe
    profile_df_temp=data.frame(name=substr(data,1,nchar(data)-char_to_cut),x=East,y=North,z=raw_data$X5,hv=raw_data$X2,top_altitude=Altitude)
    profile_df=rbind(profile_df,profile_df_temp)
    #
  }
  
  ### END FOR CYCLE THAT PROCESSES EACH FILE ###
  
  
  
  ### TRANSFORM 3D points to 2D profile ###
  
  #Create the interpolant line between the points, which will be the reference line of the 2D profile
  x0=unique(profile_df$x)
  y0=unique(profile_df$y)
  
  profile_line=lm(y0 ~ x0)#linear model of the y0 variable of x0
  intercept=profile_line$coefficients[[1]]
  angcoeff=profile_line$coefficients[[2]]
  
  #convert the explicit function y=mx+q to ax+by+c=0
  a=angcoeff
  b=-1
  c=intercept
  
  #For each point, find the closest point to the interpolant line
  
  x_on_line=((b*(b*x0-a*y0))-a*c)/(a^2+b^2)
  y_on_line=((a*(-b*x0+a*y0))-b*c)/(a^2+b^2)
  
  #order this points to assign them a distance from a point 0
  #choose as starting point the one with least x (which is the point with the lowest East coordinate)
  zero_point_index=which(x_on_line==min(x_on_line))
  zero_point_y=y_on_line[zero_point_index]
  zero_point_x=x_on_line[zero_point_index]
  distances=((x_on_line-zero_point_x)^2+(y_on_line-zero_point_y)^2)^0.5
  
  #add the column "distance" to the data frame based on the distance from the zero point on the interpolant line.
  i=0
  for (point in profile_df$x){
    i=i+1
    profile_df$distance[i]=distances[(x0==profile_df$x[i])&(y0==profile_df$y[i])]
  }

  

  #create a regular grid for plotting a contour:
  
  #find the coordinates at every distance - to calculate the altitude at every distance
  xgrid=seq(0,max(distances),by=distance_step)
  coords_on_distance_x=vector()
  coords_on_distance_y=vector()
  Altitude_on_distance=vector()
  coords_on_distance_x[1]=zero_point_x
  coords_on_distance_y[1]=zero_point_y
  Altitude_on_distance[1]=DEM$V3[((DEM$V2-zero_point_y)^2+(DEM$V1-zero_point_x)^2)^0.5==min((DEM$V2-zero_point_y)^2+(DEM$V1-zero_point_x)^2)^0.5]
  for (kl in 2:length(xgrid)){
  coords_on_distance_x[kl]=coords_on_distance_x[kl-1]+distance_step*cos(atan(angcoeff))
  #nearest neighbor interpolation of altitude of each point of the grid, based on the DEM
  coords_on_distance_y[kl]=coords_on_distance_y[kl-1]+distance_step*sin(atan(angcoeff))
  Altitude_on_distance[kl]= DEM$V3[((DEM$V2-coords_on_distance_y[kl])^2+(DEM$V1-coords_on_distance_x[kl])^2)^0.5==min((DEM$V2-coords_on_distance_y[kl])^2+(DEM$V1-coords_on_distance_x[kl])^2)^0.5]
  
  } #now all xgrid points can be related to a Altitude and coordinates.
 
 
  ##spline interpolation of altitude for each grid point. 
  if(precise_elevation==TRUE){
    
    for(gridpoint in 1:length(coords_on_distance_x)){
      
      DEM=DEM[order(((DEM$V2-coords_on_distance_y[gridpoint])^2+(DEM$V1-coords_on_distance_x[gridpoint])^2)^0.5),]
      Altitude_on_distance[gridpoint]=weighted.mean(DEM$V3[1:16],1/((DEM$V2[1:16]-coords_on_distance_y[gridpoint])^2+(DEM$V1[1:16]-coords_on_distance_x[gridpoint])^2)^0.5)
    }

    }
  
  ## SMOOTHING TOPOGRAPHY to have nicer plot, 
  if(n_ite>0){
    Altitude_on_distance=smooth(Altitude_on_distance) 
    
  for(ite in 1:n_ite){
  for (t in 2:(length(Altitude_on_distance)-1)){
    Altitude_on_distance[t]=(Altitude_on_distance[t]+Altitude_on_distance[t-1]+Altitude_on_distance[t+1])/3
  }
  }
  }
  max_hv=max(profile_df$hv)
  
  

  #interpolate without topography
  ygrid=seq(-investigation_depth,0,by=z_step)
  profile_grid=akima::interp(profile_df$distance,profile_df$z-profile_df$top_altitude,profile_df$hv,xgrid,ygrid)
   profile_grid=as.data.frame(akima::interp2xyz(profile_grid))
   
  #if interpolated without topography - put again topography on the grid, based on the distance
    profile_grid$y=profile_grid$y+Altitude_on_distance[match(profile_grid$x,xgrid)]


    # create a regular grid for a plot
    ygridn=seq(min(profile_df$top_altitude)-investigation_depth,max(profile_df$top_altitude)+z_step,by=z_step)
    profile_grid=profile_grid[is.na(profile_grid$z)==FALSE,]
    fld = interp(profile_grid$x,profile_grid$y,profile_grid$z,xgrid,ygridn,extrap = TRUE)
    fld=as.data.frame(akima::interp2xyz(fld))
for (jl in 1:length(fld$x)){
  if (fld$y[jl]>Altitude_on_distance[xgrid==fld$x[jl]]){
    fld[jl,3]=NA
  }
  if (fld$y[jl]<(Altitude_on_distance[xgrid==fld$x[jl]]-investigation_depth)){
    fld[jl,3]=NA
  }
 
}    
  ### plot the dataframe 
    
   ### INPUT GRAPHICS PARAMETERS ###
   
  # jet_colors <-  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100)
    
  jet_colors <-  colorRampPalette(c( "#00007F", "#007FFF", "aquamarine2","cornsilk", "yellow2", "red", "#7F0000"))(100)

   # max_colorscale=5 #comment next line that begins with the same variable name if you want to select manually a value of max colorscale
 max_colorscale=max_hv-1# the -1 after max_hv is to create a color scale that stops at h/v = max(H/v)-1

   breaks_var=round((seq(1, max_colorscale, length.out = 15)),1) #this length.out value should be the same of the bins in the MetR::geom_contour_fill in the ggplot function.
  label_offset=18 # increase to put labels in a higher position.
  text_size=20
  legendtext_size=14

  
### END INPUT GRAPHICS PARAMETERS ###
  
  
 ggplot2::ggplot()+
   
    metR::geom_contour_fill(aes(x=fld$x,y=fld$y,z=fld$z),bins=15)+
   metR::geom_contour2(aes(x=fld$x,y=fld$y,z=fld$z),bins=15,size=0.05)+
   
    scale_fill_stepsn(colors=jet_colors,limits=c(min(breaks_var), max(breaks_var)),breaks=breaks_var,oob=squish)+
     annotate(angle=270,size=5,geom="text",label=unique(profile_df$name),x=unique(profile_df$distance),y=unique(profile_df$top_altitude)+label_offset)+
    geom_line(color="white",aes(x=xgrid,y=Altitude_on_distance),linewidth=2.5)+
    geom_line(color="white",aes(x=xgrid,y=Altitude_on_distance-investigation_depth),linewidth=2.5)+
    geom_line(color="brown",aes(x=xgrid,y=Altitude_on_distance),linewidth=0.5)+
    theme_bw()+
    theme(text=element_text(size=text_size),legend.text=element_text(size=legendtext_size),legend.key.size = unit(1.4, 'cm'))+
    ylim(c(min(profile_df$top_altitude)-investigation_depth,max(profile_df$top_altitude)+30))+
    # coord_equal()+
    labs(x="Distance [m]",y= "Altitude [m]",fill="H/V")

  ggsave(gsub(":","-",sprintf("profile_plot%s.jpg",Sys.time())),width=12,height=6)
    
  
  
  ## eventual save of some data##
  
#   
# write.table(profile_df,"profilo_pozzi_madagascar1.txt")
# write.table(data.frame(distance=xgrid,altitude=Altitude_on_distance),"profilo_pozzi_madagascar1_topo.txt",row.names = FALSE)

### END OF THE SCRIPT ###



  
  ##old code ###
  
  
  
  #  profile_grid$z[is.na(profile_grid$z)]=0
  # 
  # #create the surface line (old way)
  # profile_df_ordered=profile_df[order(profile_df$distance),]
  # surface_topo_x=unique(profile_df_ordered$distance)
  # surface_topo_z=unique(profile_df_ordered$top_altitude,MARGIN=0)
  # surface_topo_x_resolution=2
  # 
  # 
  #                                             
  #   #prova kriging krig=kriging(profile_df$distance, profile_df$z, profile_df$hv, pixels=300)
  # xgrid=seq(0,max(distances),by=distance_step)
  # ygrid=seq(max(profile_df$top_altitude)-investigation_depth,max(profile_df$top_altitude),by=z_step)
  # profile_grid=akima::interp(profile_df$distance,profile_df$z,profile_df$hv,xgrid,ygrid)
  #  
  #  profile_grid=as.data.frame(akima::interp2xyz(profile_grid))
  #  # metR::geom_contour_fill(aes(x=krig$map$x,y=krig$map$y,z=krig$map$pred),bins=20)+
  
  # geom_point(aes(x=profile_grid$x,y=profile_grid$y,fill=profile_grid$z))+
  
  # metR::geom_contour2(aes(x=fld$x,y=fld$y,z=fld$z),bins=10,color="lightgrey")+
  # geom_contour_filled(aes(x=fld$x,y=fld$y,z=fld$z))+
  # geom_tile(aes(x=fld$x,y=fld$y,fill=fld$z))+
  # scale_fill_stepsn(colors=jet_colors,limits=c(1, max_hv),breaks=round(seq(1,max_hv,length.out=20),1),oob=squish)+#        round(exp(seq(log(1), log(max_hv), length.out = 15)),1)scale_fill_stepsn(colors=jet_colors,limits=c(1, max_hv),breaks=round(exp(seq(log(1), log(max_hv), length.out = 20)),1),oob=squish)+
  # geom_point(aes(x=profile_df$distance,y=profile_df$z,color=profile_df$hv))+ ## To just plot columns of data at each survey position
  # scale_color_gradientn(colors=jet_colors,limits=c(0, max_hv))+
  # scale_y_continuous(limits=c(max(profile_df$top_altitude)-250,max(profile_df$top_altitude)+20))+
  #
  
  
  # Altitude_on_distance= akima::interp(DEM$V1,DEM$V2,DEM$V3,coords_on_distance_x,coords_on_distance_y,extrap=TRUE,linear=FALSE)
  # Altitude_on_distance=diag(Altitude_on_distance$z)
  
  # Altitude=akima::interp(DEM$V1,DEM$V2,DEM$V3,East,North)[[3]] too slow
  