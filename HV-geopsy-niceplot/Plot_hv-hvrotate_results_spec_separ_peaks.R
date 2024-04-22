
#This script produces plots of HVSR data
#Written in January 2024 by Andrea Vergnano (andrea.vergnano@unito.it, cesare.comina@unito.it). You can use, share and modify it as long as you properly refer to the author (Creative Commons BY licence).

#Input: a series of .hv (+.log) and .grid files produced by Geopsy-hv software (with -hv and -rotate options), , and .spec files produced with the -spectrum option.



#Needed packages. Install them before running this script.
if (!require("readr")) install.packages("readr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("scales")) install.packages("scales")
if (!require("metR")) install.packages("metR")
if (!require("patchwork")) install.packages("patchwork")
if (!require("pracma")) install.packages("pracma")
if (!require("tcltk")) install.packages("tcltk")
if (!require("sp")) install.packages("sp")


library(readr)
library(ggplot2)
library(scales)
library(metR)
library(patchwork)
library(pracma)
library(tcltk)


### FUNCTIONS ###

library(sp) #required for LongLatToUtm
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}



### INPUT PARAMETERS ###
### !!!!!!!!!!!!!!!! ###
### !! CHECK THEM !! ###

f0_highpass=0 #consider as valid only values of f0 (and f1 too) > f0_highpass [Hz]. Put 0 to have no effect.
rotate90=FALSE # if you switched N and E components in geopsy.
separate_peaks=0.5 #put 0 for no effect. Value is in Hz, it is the minimum difference between peaks you expect. It is needed to not consider nearby local maxima as separate peaks.
spec=TRUE # process also .spec data which contain spectral amplitudes for each component (Z, E, N)
show_f2=TRUE # to show also f2, other than f0 and f1, in the h/v plot.
plot_georef=TRUE# produces also georeferenced rotation plots.
dimensions_georef_plot=100 # in m, dimensions of the georeferenced png rotate plot files. It actually comes 40% smaller for some reasons.
recover_coords=TRUE #to recover missing coordinates in the hv files. requires the correspondently named .SAF file in the same folder of the .hv file
recover_coords_from_file=FALSE
utmzone=33 # for eventual UTM conversion performed when giving external coordinates with recover coords from file.
cut_char=10 #characters to be cut from filename for entitling the graph. Also, one option is to go where the composite plot is created and just put the i variable, so the name of the plots is just 1,2,3,4. In any case, the software takes the files in alphabetical order, therefore, if they are named after the date-hour, they should come out in the correct order. anyanyway, there are the coordinates that prevent any mistake positioning or confusion between points.

### ramp palette for next graphs. ### 
#jet_colors <- colorRampPalette(c("bisque","darksalmon" ,"darkseagreen", "deepskyblue4","darkblue","black"))(42)
#jet_colors <- colorRampPalette(c("blue","darkcyan","lightgreen","lightyellow","pink","darkred"))(82)
#jet_colors <- colorRampPalette(c("black","green" ,"yellow", "pink","red"))(42)
# jet_colors <-  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(82)
jet_colors <-  colorRampPalette(c("black", "#00007F", "#007FFF", "aquamarine2","cornsilk", "yellow2", "red", "#7F0000"))(100)

### END INPUT PARAMETERS ###



### INPUT DATA ###

#Select all files with .hv extension in the working directory. They have to have a matching .log and a matching .grid file in the same directory
 #HVfiles=list.files(pattern=".hv$")

#Allow the user to select the hv files they are interested into. The script wants in the same folder also the .log files and .grid files.

Filtershv<- matrix(c("hv", "*.hv", "R code", ".s",
                      "Text", ".txt", "All files", "*"),
                    4, 2, byrow = TRUE)
HVfiles=tk_choose.files(filters=Filtershv)

#set the working directory where the files are
setwd(dirname(HVfiles[1]))
HVfiles=basename(HVfiles) #write the name of the files for better final output

if(recover_coords_from_file==TRUE){
localization = tk_choose.files(caption="choose the file which contains the x y z information")
xyz=read.csv(localization)
}


### END INPUT DATA ###



#Initialize list that will store information about each data file, if its data respect the 9 Sesame criteria. 

Sesame_criteria=vector("list",length(HVfiles))
names(Sesame_criteria)=HVfiles



####calculate the max peak amplitude of all the dataset to allow for using a common color scale and plot scale

#initialize a variable that will store the max peak amplitudes of each data file. This is useful to set up a certain axis-scale or color-scale limit in the plots
max_peak_amplitude=vector()

for (element in HVfiles )
  
{ 
    print(element)
  
  #read the .hv data file
  raw_data <- read_delim(element, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
  
  ############## WARNING!!! CHECK NEXT LINE, NOT FOR ALL INSTRUMENTS OR CASES!!! check f0_highpass among Input parameters above. ####
  max_peak_amplitude=append(max_peak_amplitude,max(raw_data$X2[raw_data$X1>=f0_highpass]))
}

#calculate the median peak value for scaling the next rotation images in the same way (same color scale)
#max_plot_scale=median(max_peak_amplitude[-c(which.max(max_peak_amplitude),which.min(max_peak_amplitude))])

max_plot_scale=max(max_peak_amplitude) # with this option, not the median, but the max amplitude value is used

max_plot_scale=5#to set it manually


#########################################ANALYSIS of .hv files: ############################

#this for cycle verify the SESAME criteria one .hv input file at a time, and plots their results on a Frequency-amplitude graph, saving the image in the current working directory

# the variable i is used as counter during the next for cycle.
i=0 
print("start processing .hv files...")
for (element in HVfiles )
  
  { 
i=i+1
  print(sprintf("processing .hv file %d out of %d, named %s",i,length(HVfiles),element))
  
#read the .hv data file
 raw_data <- read_delim(element, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
 
############## WARNING!!! CHECK NEXT LINE, NOT FOR ALL INSTRUMENTS OR CASES!!!
 #should not be useful anymore after introduction of f0_highpass input parameter.
# raw_data=raw_data[raw_data$X1>0.5,] #select data only greater than 0.5 Hz, given the instrument flat response only after 2 Hz
 
 #####################SESAME criteria verification: ######################
 
 #load .log file, named exactly as the .hv file (except the extension), containing infos about timewindows used in the processing performed by Geopsy (it will be useful for criterion #1)
 log_data_raw=read_lines(sprintf("%slog",substr(element,1,nchar(element)-2)))
 
 #find the part of the file related to window length
 skip_index=which(grepl("Window length", log_data_raw, fixed=TRUE) & grepl("Start time", log_data_raw, fixed=TRUE)) 
 # skip_index=which(log_data_raw=="# Start time 	 End Time 	 Window length") #alternatively to the code line above, if you are sure that all your file contains exactly that title of the timewindow section of log file
 
 #extract information about time windows from the .log file. It will be useful in Criterion #1
  log_data<- read_table(sprintf("%slog",substr(element,1,nchar(element)-2)),col_names = FALSE, skip = skip_index,comment="#",show_col_types = FALSE)
  

  #extract f0 and a0 from the .hv file. These code lines first find the .hv-file-line where f0 is reported, and then select the f0 value and the 2 values of f0+-standard deviation.
  f0_file=read_lines(element)
  f0_index=which(grepl("f0 from windows", f0_file, fixed=TRUE))
  f0_raw=f0_file[f0_index]
  f0=as.numeric(data.frame(a = unlist(strsplit(f0_raw, "\t")))[2,1])
  f0stdup=as.numeric(data.frame(a = unlist(strsplit(f0_raw, "\t")))[4,1]) #std dev of the frequency f0.
  f0stddown=as.numeric(data.frame(a = unlist(strsplit(f0_raw, "\t")))[3,1])

  #similarly to the previous lines, this code lines extract the information about a0, the peak amplitude at frequency f0.
    a0_index=which(grepl("f0 amplitude", f0_file, fixed=TRUE))
  a0_raw=f0_file[a0_index]
  a0=as.numeric(data.frame(a = unlist(strsplit(a0_raw, "\t")))[2,1])
  
  #find second peak f1, the second local maximum of the function
 local_maxima_index=which(c(0,0,diff(sign(diff(raw_data$X2))))==-2&raw_data$X1>=f0_highpass) #c(0,0,diff(...)) is because due to double diff, the first logical vector has shrinked by 2. 
 
 
 ### this difficult to interpret code is needed to eliminate the nearby peaks. It takes the maxima, and understand which maxima is too near from each other. Then, it select as wrong the least-amplitude one.
 if(any(diff(raw_data$X1[local_maxima_index])<separate_peaks))
    {
      wrong_maxima=local_maxima_index[which(diff(raw_data$X1[local_maxima_index])<1)]
      
      pip=0
      for (pippo in wrong_maxima)
      {
        pip=pip+1
        potentially_wrong=c(pippo,local_maxima_index[which(local_maxima_index==pippo)+1])
      
        if(raw_data$X2[pippo]>raw_data$X2[local_maxima_index[which(local_maxima_index==pippo)+1]])
        {
          wrong_maxima[(wrong_maxima==pippo)]=local_maxima_index[which(local_maxima_index==pippo)+1]
        }
      }
      
      wrong_maxima_index=match(wrong_maxima,local_maxima_index)
      local_maxima_index=local_maxima_index[-(wrong_maxima_index)]
      
  
 }
 
  local_maxima_df=data.frame( frequency=raw_data$X1[local_maxima_index], amplitude=raw_data$X2[local_maxima_index])
local_maxima_df=local_maxima_df[order(local_maxima_df$amplitude,decreasing=TRUE),]
f1=0
a1=0
if (length(local_maxima_df$frequency)>1){f1=local_maxima_df$frequency[2]}
if (length(local_maxima_df$frequency)>1) {a1=local_maxima_df$amplitude[2]}

### find first peak f0 similarly to f1 ### - not necessary since it is already found by geopsy and saved on the .hv file.
f0=0
a0=0
f0=local_maxima_df$frequency[1]
a0=local_maxima_df$amplitude[1]

#find third peak f0 similarly to f1
f2=0
a2=0
if (length(local_maxima_df$frequency)>2){f2=local_maxima_df$frequency[3]}
if (length(local_maxima_df$frequency)>2) {a2=local_maxima_df$amplitude[3]}

  
  ############################
  
 #1 criterion: for  every window of length twl, twl>10/f0. This stores a vector with 1 or 0, according to if the timewindows respect the criteria or not. It computes the percentage.
  
  crit1temptemp=vector()
  for (twl in log_data$X3)  #logdata$X3 is the list of the timewindows, which generally have the same length (but it is not necessary)

       {
         if (twl>10/f0)
         {
           crit1temptemp=c(crit1temptemp,1)
         }
         else 
        {
          crit1temptemp=c(crit1temptemp,0)
          
        }
          
  }
  
  #in the crit1temp variable, the percentage of respect of the 1st criterion is saved. If it is 100%, it just writes "1" on this variable, consistently with other criteria.
  crit1temp=sprintf("%d %%",sum(crit1temptemp)/length(crit1temptemp)*100)
  if(crit1temp=="100 %"){crit1temp=1}
  
  ##############################
  
  #2 criterion: f0*(sum of twl) > 200
  crit2temp=0
  
  twl_sum=sum(log_data$X3)
  if (twl_sum*f0>200){crit2temp=1}
  
  #############################################
  
  #3 criterion: dev st in the range f0/2 ...2*f0 is  less than 2
  crit3temp=0
  
  dev_st_up=raw_data$X4-raw_data$X2
  dev_st_down=raw_data$X2-raw_data$X3
  
  if ((all(dev_st_up[f0/2<raw_data$X1 & raw_data$X1<2*f0]<2))&(all(dev_st_down[f0/2<raw_data$X1 & raw_data$X1<2*f0]<2)))
  {
    crit3temp=1
  }
  
     #####4-8 criteria are about the clarity of the f0-amplitude peak
  
  #4 criterion: in the range 0.25*f0 .... f0 there is at least one amplitude < 0.5 a0
  crit4temp=0
  
  if (any(raw_data$X2[f0/4<raw_data$X1 & raw_data$X1<f0]<0.5*a0))
  {
    crit4temp=1
  }
  
  ########################################
  #5 criterion: specular to #4 criterion, range: f0 ... 4*f0
  
  crit5temp=0
  
  if (any(raw_data$X2[f0<raw_data$X1 & raw_data$X1<4*f0]<0.5*a0))
  {
    crit5temp=1
  }
  
  ##################################
  
  #6 criterion: amplitude a0 of the HV ratio is greater than 2 
  crit6temp=0
  if(a0>2){crit6temp=1}
  
  ######################################
  
  #7 criterion: in the amplitude confidence curves, there is a peak near f0, more specifically in the range of frequency values 5% near f0.
  crit7temp=0
  
  #second derivative to find maxima in confidence curves. For how it is written the following lines algorithm, it may not work if the maxima peak has a plateau (e.g. due to signal saturation)
  maxima_stddown=which(diff(sign(diff(raw_data$X3)))==-2) 
  maxima_stdup=which(diff(sign(diff(raw_data$X4)))==-2) 
  
  if((any(abs(raw_data$X1[maxima_stdup]-f0)<0.05*f0))&(any(abs(raw_data$X1[maxima_stddown]-f0)<0.05*f0)))
  {
    crit7temp=1
  }
  
  ####################################################
  
  #8 criterion: the standard deviation of f0 (not the standard deviation of a0) is within a certain range, according to a tabella in the sesame specification, here implemented as a series of if, else if statements.
  

  crit8temp=0
  if (f0>2){
  if((f0stdup-f0<0.05*f0)&(f0-f0stddown<0.05*f0))
  {
    crit8temp=1
  }
  }  else if (f0>1){
    if((f0stdup-f0<0.1*f0)&(f0-f0stddown<0.1*f0))
    {
      crit8temp=1
    }
  } else if (f0>0.5){
    if((f0stdup-f0<0.15*f0)&(f0-f0stddown<0.15*f0))
    {
      crit8temp=1
    }
  } else if (f0>0.2){
    if((f0stdup-f0<0.2*f0)&(f0-f0stddown<0.2*f0))
    {
      crit8temp=1
    }
  } else {
    if((f0stdup-f0<0.25*f0)&(f0-f0stddown<0.25*f0))
    {
      crit8temp=1
    }
  }
  ################################################################
  #9 criterion: standard deviation of a0 <1.58
  
  crit9temp=0
  
  f0_near=raw_data$X1[  which(abs(raw_data$X1 - f0) == min(abs(raw_data$X1 - f0)))   ] #find the nearest value to f0 in the frequency points (first column of .hv file), to allow for an easier subsequent if comparison.
  f0_near=f0_near[1]
  
  if (((dev_st_up[raw_data$X1==f0_near]<1.58))&((dev_st_down[raw_data$X1==f0_near]<1.58)))
  {
    crit9temp=1
  }
  
  
  #Coordinates from hv file (in UTM):
  
  UTM_index=which(grepl("Position", f0_file, fixed=TRUE))
  UTM_raw=f0_file[UTM_index]
  East=as.numeric(data.frame(a = unlist(strsplit(UTM_raw, split="[ \t]")))[3,1])
  North=as.numeric(data.frame(a = unlist(strsplit(UTM_raw, split="[ \t]")))[4,1])
  
  ### RECOVER MISSING COORDINATES ###
  
  if(recover_coords==TRUE&recover_coords_from_file==FALSE)
  {
    
    if(East==0&North==0){
      
      testo=readLines(sprintf("%s",substr(element,1,nchar(element)-3)),n=50) # opens the .saf file correspondent to the hv file
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
      
      
      ### extract geolocalization, based on the first 5 chars of the string, to avoid mismatch between, for example, points called L1_01B and L1_01, which have actually the same position. ###
      
      x=xyz$East_wgs84_utm38S[substr(xyz$TEM,1,5)==substr(element,1,5)]
      y=xyz$North_wgs84_utm38S[substr(xyz$TEM,1,5)==substr(element,1,5)]
     
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
  
  
  ### save criteria results on the Sesame criteria list, that will update for every file, to produce then a full .csv spreadsheet as report.
 
  #SAVE LATER!! Sesame_criteria[[i]]=c(crit1temp,crit2temp,crit3temp,crit4temp,crit5temp,crit6temp,crit7temp,crit8temp,crit9temp,f0,a0)
  ################################################################
  
  #plot a simple graph of amplitude versus frequency and the standard deviation.
  
p1=ggplot(raw_data)+ 
  geom_line(aes(x=X1,y=X2,linetype="Average amplitude"),color="black")+
  geom_line(aes(x=X1,y=X3,linetype="Standard deviation"),color="black")+
  geom_line(aes(x=X1,y=X4,linetype="Standard deviation"),color="black")+
  scale_x_log10(limits=c(0.5,30),breaks=c(1,2,3,4,5,10,20,30))+
#  scale_y_continuous(limits=c(0,max(raw_data$X4[raw_data$X1>0.5])))+ #this condition is because we are visualizing data only from 0.5 Hz on
  #  scale_y_continuous()+
   scale_y_continuous(limits=c(0,max_plot_scale))+
 labs(x="Frequency [Hz]",y= "H/V Ratio [-]",linetype=sprintf("f0 = %.01f Hz \nA0 = %.01f\n\nf1 = %.01f Hz \nA1 = %.01f\n\nf2 = %.01f\nA2 = %.01f",round(f0,1),round(a0,1),round(f1,1),round(a1,1),round(f2,1),round(a2,1)))+
# labs(x="Frequency [Hz]",y= "H/V Ratio [-]",title=sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-7)))+ # for saving plot alone, it requires the title
      theme_bw()
  
#   #save the image on the current working folder
# ggsave(sprintf("%s.png",substr(element,1,nchar(element)-11)), width = 9, height = 6) # for geagps .SAF files
# #ggsave(sprintf("%s.png",substr(element,1,nchar(element)-7)), width = 9, height = 6) # for minishark .txt files


#update the peak amplitude vector with the value of the current data file

  
  

### .GRID FILE ANALYSIS : directions of maximum amplitudes, polarization ###

  
  
  
print(sprintf("processing .grid file %d out of %d, named %s",i,length(HVfiles),sprintf("%sgrid",substr(element,1,nchar(element)-2))))


#read data files

raw_data_grid <- read_delim(sprintf("%sgrid",substr(element,1,nchar(element)-2)), delim = " ", escape_double = FALSE, col_names = TRUE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
raw_data_grid=as.data.frame(raw_data_grid)
colnames(raw_data_grid)=c("x","y","H/V SR")

########WARNING!!! not applicable for every case or instrument. It selects only data at frequencies greater than 0.5
raw_data_grid=raw_data_grid[raw_data_grid$x>=0.5,]


# ### !!!!!!! ROTATE 90 DEGREES THE DATA !!!!!######

if(rotate90==TRUE){
  raw_data_grid$y=raw_data_grid$y+90
  raw_data_grid$y[raw_data_grid$y>180]=raw_data_grid$y[raw_data_grid$y>180]-180
  raw_data_grid_temp0=raw_data_grid[raw_data_grid$y==180,]
  raw_data_grid_temp0$y=0
  raw_data_grid=rbind(raw_data_grid,raw_data_grid_temp0)
}
# ### END ROTATE 90 DEGREES THE DATA ###



#get the polarization direction at f0

#get a value of f0 on the gridded values nearby
f0_near_rot=  f0_near=raw_data_grid$x[  which(abs(raw_data_grid$x - f0) == min(abs(raw_data_grid$x - f0)))   ] #find the nearest value to f0 in the frequency points (first column of .hv file), to allow for an easier subsequent if comparison.
f0_near_rot=f0_near_rot[1]

#compute the maximum value of amplitude at f0, in which direction is it?
pol_direction_temp=raw_data_grid[raw_data_grid$x==f0_near_rot & raw_data_grid$y!=180,]#180 frequency is to be avoided to avoi
pol_direction=(pol_direction_temp$y[pol_direction_temp$`H/V SR`==max(pol_direction_temp$`H/V SR`)])
pol_direction=max(pol_direction) #this serves only to transform a 2-element vector in a single value

#save the polarization information on Sesame criteria output matrix

Sesame_criteria[[i]]=c(crit1temp,crit2temp,crit3temp,crit4temp,crit5temp,crit6temp,crit7temp,crit8temp,crit9temp,f0,a0,f1,a1,pol_direction,East,North)

###plot frequency - azimuth cartesian graph


### Color scale correction: putting the values higher than the max value of the color scale equal to the max value of the color scale. If don't do this, the graphs will be badly colored, due to a bug or I don't know.
#raw_data_grid$`H/V SR`[raw_data_grid$`H/V SR`>max_plot_scale]=max_plot_scale

breaks_var=round((seq(1, max_plot_scale, length.out = 15)),1) #this length.out value should be the same of the bins in the MetR::geom_contour_fill in the ggplot function.


p2=ggplot(raw_data_grid,aes(x,y))+
  metR::geom_contour_fill(aes(x=x,y=y,z=`H/V SR`))+ #interpolation is to avoid pixelated image
  labs(x="Frequency [Hz]",y= "Azimuth [degrees]",fill="H/V")+
#  labs(x="Frequency [Hz]",y= "Azimuth [degrees]",title=sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-9)))+ #with the title
 # scale_fill_continuous(limits=c(0,max_plot_scale))
#  scale_fill_stepsn(colors=jet_colors,limits=c(0,max_plot_scale),oob=squish,n.breaks=10)+ #stapped color map
  # scale_fill_gradientn(colors=jet_colors,limits=c(0,max_plot_scale),oob=squish)+ #continuous color map
  # scale_fill_gradientn(colors=colors,oob=squish)+ #without imposed scale color limit.
  scale_x_log10(limits=c(0.5,30),breaks=c(1,2,3,4,5,10,20,30))+
  # scale_fill_stepsn(colors=jet_colors,limits=c(0,max_plot_scale),breaks=round(logspace(0,log10(max_plot_scale),15),1),oob=squish)+
  scale_fill_stepsn(colors=jet_colors,limits=c(min(breaks_var), max(breaks_var)),breaks=breaks_var,oob=squish)+
  
  scale_y_continuous(breaks=c(0,45,90,135,180))+
  theme_bw()+
  theme(legend.text=element_text(size=8),legend.key.height = unit(1, 'cm'))+
  geom_vline(xintercept=c(1,2,3,4,5,10),colour="grey",linetype="dotted")+
geom_hline(yintercept=c(0,45,90,135,180),colour="grey",linetype="dotted")




#save the ggplot plot

# ggsave(sprintf("%s_Rotate.png",substr(element,1,nchar(element)-13)), width = 9, height = 6)#for .SAF geagps files 
#ggsave(sprintf("%s_Rotate.png",substr(element,1,nchar(element)-9)), width = 9, height = 6)# for minishark .txt files


###plot frequency - azimuth in a  polar heatmap using polarPlot (package openair)


# assign mirrored values to 185-355 angles.
raw_data_180_360=raw_data_grid[which((raw_data_grid$y !=0)),]
raw_data_180_360$y=180+raw_data_180_360$y
raw_data_360=rbind(raw_data_grid,raw_data_180_360)
colnames(raw_data_360)=c("Frequency[Hz]","y","H/V SR")


p3=ggplot(raw_data_360)+
  metR::geom_contour_fill(aes(x=`Frequency[Hz]`,y=y,z=`H/V SR`))+ #interpolation is to avoid pixelated image
  #  labs(x="Frequency [Hz]",y= "Azimuth [degrees]",title=sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-9)))+ #with the title
  scale_fill_stepsn(colors=jet_colors,limits=c(0,max_plot_scale),breaks=round(logspace(0,log10(max_plot_scale),15),1),oob=squish)+

  coord_polar(theta="y")+
  scale_y_continuous(breaks=c(0,90,180,270),labels=c("N","E","S","W"))+
  labs(x=" ",y=sprintf("Polarization direction = %i°",pol_direction),fill="H/V")+
  
  geom_hline(yintercept = seq(0, 360, by = 90), colour = "grey", size = 0.35,linetype="solid")+
  geom_vline(xintercept = c(1,2,5,10,20),colour = "grey", size = 0.35,linetype="dotted")+
  scale_x_log10(limits=c(0.5,30),breaks=c(0.5,2,5,10,20,30),labels=c("0","2","5","10","20","30"))+
  annotate("text",y=0,x=c(1,2,5,10,20,30),label=c("1","2","5","10","20","30 Hz"),alpha=1,color="darkorange",size=3,fontface="bold")+
  theme_bw()+
  theme(legend.text=element_text(size=8),legend.key.height = unit(1, 'cm'),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  geom_segment(aes(x=0,y=pol_direction, xend=3, yend=pol_direction),arrow=arrow(length=unit(0.3,"cm")),color="red")+ #these two lines plot the arrow in the polarization direction
  geom_segment(aes(x=0,y=pol_direction+180, xend=3, yend=pol_direction+180),arrow=arrow(length=unit(0.3,"cm")),color="red")
# annotate(geom="label", y = 0, x = c(0,2, 5, 10, 20), label = c('0 Hz', '2','5', '10', '20'),fill="white",size=2,alpha=0.6) +

#ggsave(sprintf("%s_Rotate_polar.png",substr(element,1,nchar(element)-9)), width = 9, height = 6)# for minishark .txt files

if(spec==TRUE){
  raw_data_specZ <- read_delim(sprintf("%s_Z.spec",substr(element,1,nchar(element)-3)), delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
  raw_data_specE <- read_delim(sprintf("%s_E.spec",substr(element,1,nchar(element)-3)), delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
  raw_data_specN <- read_delim(sprintf("%s_N.spec",substr(element,1,nchar(element)-3)), delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE,show_col_types = FALSE)
  
  if(rotate90==TRUE){
    raw_data_specNtemp=raw_data_specN
    raw_data_specN=raw_data_specE
    raw_data_specE=raw_data_specNtemp
  }
  
  
  p4=ggplot()+ 
    geom_line(aes(x=raw_data_specZ$X1,y=raw_data_specZ$X2,color="Z"),size=1)+
    geom_line(aes(x=raw_data_specE$X1,y=raw_data_specE$X2,color="E"),size=1)+
    geom_line(aes(x=raw_data_specN$X1,y=raw_data_specN$X2,color="N"),size=1)+
    scale_x_log10(limits=c(0.5,30),breaks=c(1,2,3,4,5,10,20,30))+
    #  scale_y_continuous(limits=c(0,max(raw_data$X4[raw_data$X1>0.5])))+ #this condition is because we are visualizing data only from 0.5 Hz on
    #  scale_y_continuous()+
    scale_y_continuous()+
    scale_color_manual(values=c("Z" = "blue", "E" = "red", "N" = "green"),name="Component")+
    labs(x="Frequency [Hz]",y= "Power spectral density [(counts)^2/Hz]",fill="Component")+
    # labs(x="Frequency [Hz]",y= "H/V Ratio [-]",title=sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-7)))+ # for saving plot alone, it requires the title
    theme_bw()
}

### MAKE A COMPOSITE PLOT ###

if(spec==TRUE)
{
  # (p1+p4)/(p2+p3)+plot_annotation(title = sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-3)))
  
  (p1+p4)/(p2+p3)+plot_annotation(title = sprintf("H/V ratio of point %i",i))
}else{
# p1/(p2+p3)+plot_annotation(title = sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-cut_char))) 
  p1/(p2+p3)+plot_annotation(title = sprintf("H/V ratio of point %i",i))
  
}
ggsave(sprintf("%s.png",substr(element,1,nchar(element)-3)), width = 9, height = 6)

### END COMPOSITE PLOT ###



###SAVE GEOREFERENCED PNG with pgw for plotting georeferenced images in qgis directly ###
if(plot_georef==TRUE){
  
  widthplot=6
  heightplot=6
  dpiplot=300
  
  #produces a plot similar to the polar coordinates one, but with transparent background
  p3+
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel background
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #avoid plotting any text or grid related to axes.
      axis.line=element_blank(),axis.text.x=element_blank(),
      axis.text.y=element_blank(),axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none") #transparent legend panel
  
  #save the plot on a png file
  ggsave(sprintf("%s_Rotate_polar.png",substr(element,1,nchar(element)-3)), width = widthplot, height = heightplot,dpi=dpiplot,bg="transparent")# for minishark .txt files
  
  #produces a tfw-like formatted files, by making 6 lines, as indicated in wikipedia: https://en.wikipedia.org/wiki/World_file
  map_units=dimensions_georef_plot #m, 
  A=map_units/(widthplot*dpiplot) #pixel size in the x-direction in map units/pixel
  D=0 #rotation about y-axis
  B=0 #rotation about x-axis
  E=-map_units/(heightplot*dpiplot)# pixel size in the y-direction in map units, almost always negative
  C=East-map_units/2 #x-coordinate of the center of the upper left pixel
  F=North+map_units/2 # y-coordinate of the center of the upper left pixel
  
  cat(sprintf("%f\n%f\n%f\n%f\n%f\n%f",A,D,B,E,C,F),file=sprintf("%s_Rotate_polar.pgw",substr(element,1,nchar(element)-3)))
}#end if
###END GEOREFERENCED PNG ###

}#end for

#save sesame criteria
Sesame_criteria=as.data.frame(Sesame_criteria)
write.csv(Sesame_criteria,"f0_A0_Sesame_criteria_report.csv",row.names=c("criterion 1","criterion 2","criterion 3","criterion 4","criterion 5","criterion 6","criterion 7","criterion 8","criterion 9","f0","A0","f1","A1","polarization direction","East","North"))






# print("Creating plots of noise polarization...")
# 
# j=0#used to count during the next for cycle
# 
# GRIDfiles=list.files(pattern=".grid$")
# 
# for (element in GRIDfiles){
#   
#   j=j+1
#   
#  
#    #polar plot and export
#   # tiff(sprintf("%s_Rotate_polar.tiff",substr(element,1,nchar(element)-13)), units="in", width=5, height=4, res=300, compression = 'lzw')#for .SAF geagps files
#  #  tiff(sprintf("%s_Rotate_polar.tiff",substr(element,1,nchar(element)-9)), units="in", width=5, height=4, res=300, compression = 'lzw')# for minishark .txt files
   
   # pp= polarPlot(mydata=raw_data_360, x="Frequency[Hz]",wd="y",pollutant="H/V SR",cols = "jet",alpha=1,limits=c(0,max_plot_scale),key.position = "bottom", key.header = sprintf("H/V ratio. Measurement number %s",substr(element,1,nchar(element)-9)), key.footer = NULL)
 #pp= polarPlot(mydata=raw_data_360, x="Frequency[Hz]",wd="y",pollutant="H/V SR",cols = "jet",alpha=1,key.position = "bottom", key.header = sprintf("H/V ratio of point %s",substr(element,1,nchar(element)-9)), key.footer = NULL)
  #all plot having the same color scale, to better comparison.
  #      dev.off()
     
        
          # }

print("You should find the output images and the report about SESAME criteria in the same folder of your data files.")
#thank you for using this script. Please contact the authors if you have any question or suggestion. 

