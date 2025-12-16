##SequencERT_Ws2D_V1.0.R##

## Release notes ##

#V 1.0 December 2025. Cleaning of this old script (first written in June 2024) and publishing on github

## End of Release notes


## introduction ##

### This script helps to create standard Wenner Schlumberger sequences optimized for a multichannel ERT resistivimeter, such as the IRIS Syscal Pro
# it also create a roll-along sequence.


### WARNING! ###

# To create a sequence that the instrument (in my case the Syscal Pro) can read as an optimized one, there is the need to create 
#"gapfiller" measurements, which are not wenner shlumeberger ones.
#They can be removed in postprocessing, after field survey.
#To check which measurements are not in wenner schlumberger configuration, just compute the mean of AB and the mean of MN, and if are not the same, they are not in WS configuration

### END WARNING ###

##end introduction



###INPUT PARAMETERS. CHECK THEM!!!###

channels=72 #number of channels or electrode of your instrument
roll_channels=24 #how many channels you want to go on when doing a roll along.
opt_ch=10 #number of measurements that the instrument can do at the same time.
min_aperture=3 #in WS the minimum aperture possible is 3
max_aperture=71 #max AB aperture, in channels. Error if greater than (channels-1)
max_meas_per_ap=5 #max measurements per aperture (not always this value will be reached anyway, for small apertures). Consider that 5 measurements will actually be 9, because 4 gapfillers will be needed, and then those 9 measurements will be performed at the same time.
  
# Set resolution of the survey. 100 is max
res_near_surface=35 #%       #this is a percentage
res_at_medium_depth=26 #%    #this is a percentage
res_at_depth=25 #%           #this is a percentage
increase_res_with_depth= 50 #% #the higher this percentage, the higher resolution at depth, but the higher number of measurements. 

time_per_measurement=2 #seconds needed for each current injection

### END INPUT PARAMETERS ###


###calculate the optimal array of apertures according to the resolution percentages in input. ###

ap_del=vector()

surf_mean_ap=max_aperture/5/2
med_mean_ap=(max_aperture/5+max_aperture/2)/2
depth_mean_ap=(max_aperture/2+max_aperture)/2
store_del=vector()
del_counter=0
for (i in min_aperture:max_aperture){
  if(del_counter>100){
    del_counter=del_counter-100
      ap_del=append(ap_del,i)
  }
  del_counter=del_counter+ ((100-res_near_surface)/(abs(surf_mean_ap-i)+1)^2+(100-res_at_medium_depth)/(abs(med_mean_ap-i)+1)^2+(100-res_at_depth)/(abs(depth_mean_ap-i)+1)^2)/(1/(abs(surf_mean_ap-i)+1)^2+1/(abs(med_mean_ap-i)+1)^2+1/(abs(depth_mean_ap-i)+1)^2)
store_del=append(store_del,((100-res_near_surface)/(abs(surf_mean_ap-i)+1)^2+(100-res_at_medium_depth)/(abs(med_mean_ap-i)+1)^2+(100-res_at_depth)/(abs(depth_mean_ap-i)+1)^2)/(1/(abs(surf_mean_ap-i)+1)^2+1/(abs(med_mean_ap-i)+1)^2+1/(abs(depth_mean_ap-i)+1)^2))
  }

plot(100-store_del)

aperture_array=c(min_aperture:max_aperture)
aperture_array=aperture_array[!aperture_array %in% ap_del]







#create a void dataframe to store the sequence in.

sequ=data.frame(A=integer(),B=integer(),M=integer(),N=integer())

##compute the WS sequence 

for (ap in aperture_array) {  #for every aperture spacing, in ascending order (min aperture=3, which is ABMN= 1 4 2 3 , max aperture =max_aperture)...
  

  for (i in seq(1,channels-ap,by=ceiling(ap/(increase_res_with_depth*max_aperture/100) ))){#for every AB of that aperture... (this ceiling is to account for larger apertures of the electrodes, and decrease resolution at depth)
    
 
    
    for (j in 1:max(min(max_meas_per_ap,ap-3),1)){ ### for every measurement you want to perform with that exact AB location...
      
      A=i
      B=i+ap
      M=i+j*ceiling((ap)/(max_meas_per_ap*4)) #this ceiling of is to account for larger apertures of the electrodes, while maintaining the optimal 5 measurements per time - which means 9 measurements includeing gapfillers.
      N=i+ap-j*ceiling((ap)/(max_meas_per_ap*4)) #max_meas_per_ap*4 was set to allow the maximum desired measurements per aperture to be reached the most often possible. (e.g. with a max_meas_per_ap of 5, which means 9 measurements including gapfillers, it will be possible to space those 5 measurements 2 instead of 1, maintaining 5 total measurements, only with an aperture greater than 20. This is the reason of max_meas_per_ap*4)
      
      if(M<N){sequ[nrow(sequ) + 1,] = list(A,B,M,N)} #to avoid creating quadripoles with M=N or M>N (this last case means that they are equal to quadripoles already computed)
    }
    
  }
  
    
  
}

### add gapfillers ###



sequ2=sequ[1,]

for (j in 2:length(sequ[,2])){
  if((sequ$A[j]==sequ$A[j-1])&(sequ$B[j]==sequ$B[j-1])&(sequ$M[j]!=sequ$N[j-1])){
      A=sequ$A[j]
      B=sequ$B[j]
      M=sequ$N[j-1]
      N=sequ$M[j]

    sequ2=rbind(sequ2,list(A,B,M,N))
  }
  sequ2=rbind(sequ2,sequ[j,])

}

sequ=sequ2
###calculate acquisition time based on number of measurements and optimization degree.

meas_count=0
counter=0

for (j in 2:length(sequ[,1])){
  if((sequ$A[j]==sequ$A[j-1])&(sequ$B[j]==sequ$B[j-1]&(counter<10))){
    counter=counter+1
  }else{
    meas_count=meas_count+1
    counter=0
  }
}

acquisition_time=meas_count*time_per_measurement

opt_degree=length(sequ[,1])/(meas_count*opt_ch)*100


#output a basic report about the created sequence

print(sprintf("Sequence created with %i quadripoles and %i current injections.", length(sequ[,1]),meas_count))
print(sprintf("Your estimated acquisition time is %.01f minutes",acquisition_time/60))
print(sprintf("Optimization percentage of this sequence is %.01f %%",opt_degree))







### In this part of the script, the roll-along sequence is created, equal to the sequence but without those measurements that would be exactly the same already done, to save time and space.
#attention has been paid that the roll along sequence is still optimized in this way.

length_sequ=length(sequ[,1])
sequ_roll_temp=sequ+roll_channels
sequ_all_temp=rbind(sequ,sequ_roll_temp)
sequ_all_temp=sequ_all_temp[!duplicated(sequ_all_temp),] #removes the duplicated measurements from the roll along sequence
sequ_roll_temp2=sequ_all_temp[(length(sequ[,1])+1):length(sequ_all_temp[,1]),]
sequ_roll=sequ_roll_temp2-roll_channels


#calculate acquisition time of roll along sequence

meas_count_roll=0
counter=0

for (j in 2:length(sequ_roll[,1])){
  if((sequ_roll$A[j]==sequ_roll$A[j-1])&(sequ_roll$B[j]==sequ_roll$B[j-1]&(counter<10))){
    counter=counter+1
  }else{
    meas_count_roll=meas_count_roll+1
    counter=0
  }
}
acquisition_time_roll=meas_count_roll*time_per_measurement

opt_degree_roll=length(sequ_roll[,1])/(meas_count_roll*opt_ch)*100


#output a basic report about the created sequence
print(sprintf("Roll along sequence created with %i quadripoles and %i current injections.", length(sequ_roll[,1]),meas_count_roll))
print(sprintf("Your estimated acquisition time is %.01f minutes",acquisition_time_roll/60))
print(sprintf("Optimization percentage of this sequence is %.01f %%",opt_degree_roll))


##export output

#first, create an output folder
ifelse(!dir.exists(file.path("output")),
       dir.create(file.path("output")),
       "Directory Exists")



#export .csv files containing the sequences. This can be imported in Resipy - Forward modeling
write.table(sequ,sprintf("output/WS%io_%i_%i_%i_%i.csv",channels,res_near_surface,res_at_medium_depth,res_at_depth,increase_res_with_depth),row.names = FALSE,quote=FALSE,sep=",")
write.table(sequ_roll,sprintf("output/WS%ioR%i_%i_%i_%i_%i.csv",channels,roll_channels,res_near_surface,res_at_medium_depth,res_at_depth,increase_res_with_depth),row.names = FALSE,quote=FALSE,sep=",")

### export .txt files containing the sequences ###
#These can be opened in Electre pro to import them into the Syscal Pro.
#In electre pro, please save them to .sqz file type, in order to save the desired Q factor, stacking options, voltage options...

#write normal sequence
ele_pos=data.frame(num=1:channels,X=0:(channels-1),Y=0,Z=0)
colnames(ele_pos)[1]="#"
write.table(ele_pos,sprintf("output/WS%io_%i_%i_%i_%i.txt",channels,res_near_surface,res_at_medium_depth,res_at_depth,increase_res_with_depth),row.names = FALSE,quote=FALSE,sep="\t")
sequ_txt=cbind(c(1:length(sequ[,1])),sequ)
colnames(sequ_txt)[1]="#"
write.table(sequ_txt,sprintf("output/WS%io_%i_%i_%i_%i.txt",channels,res_near_surface,res_at_medium_depth,res_at_depth,increase_res_with_depth),row.names = FALSE,quote=FALSE,sep="\t",append = TRUE)

#write roll along sequence
write.table(ele_pos,sprintf("output/WS%ioR%i_%i_%i_%i_%i.txt",channels,roll_channels,res_near_surface,res_at_medium_depth,res_at_depth,increase_res_with_depth),row.names = FALSE,quote=FALSE,sep="\t")
sequ_txt_roll=cbind(c(1:length(sequ_roll[,1])),sequ_roll)
colnames(sequ_txt_roll)[1]="#"
write.table(sequ_txt_roll,sprintf("output/WS%ioR%i_%i_%i_%i_%i.txt",channels,roll_channels,res_near_surface,res_at_medium_depth,res_at_depth,increase_res_with_depth),row.names = FALSE,quote=FALSE,sep="\t",append = TRUE)



###plot pseudosection - not accurate on depths ###

xc=(sequ$A+sequ$B)/2
xp=(sequ$M+sequ$N)/2
xrhoa=(xc+xp)/2
zrhoa=(abs(sequ$A-xp))/2
plot(xrhoa,-zrhoa)


### END OF THE SCRIPT. have a nice day###


