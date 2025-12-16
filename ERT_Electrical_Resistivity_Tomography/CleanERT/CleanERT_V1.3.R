###ClenERT_V1.3.R###


##Release notes ###
#1.0 April-May 2025. Basic functions like prosys III. 

#1.1 July 2025.   Added tool to create electrode indexes from their positions. They are needed e.g. for the protocol.dat format of Resipy - R2 - R3
#                 Added my email
#
#1.2 November 2025. Minor fixes
#
#1.3 December 2025 Minor fixes, updated introduction

##End of Release notes##



#introduction

#Hello! This script allows to clean raw ERT data (such as ProsysIII software) and prepare them for inversion.
#It was written by Andrea Vergnano. andrea.vergnano@posteo.eu, andrea.vergnano@unito.it.
#This script has functions stored in separate files that have to be in the same folder and working directory.

#This script allows you to:
#update electrode positions
#update electrode numbers (indexes)
#remove outlier data
#plot main variables as histograms - and consequently filter data outside the main distribution (for example, only data with recorded voltage from 0.1 to 1000 mV)
#merge close electrodes
#identify bad electrodes based on the electrodes you filtered more frequently
#plot and export pseudosection (2D or 3D)
#export to protocol.dat format to be inverted in ResIPy.
#put multiple protocol.dat files together.

#It is still semiautomatic, that is, the suggested use is to select only 1 section of the code at a time, run it, and see the results in the variable list
#You may need to change parameters, such as filter parameters, to use it accordingly to your need. Each section is commmented to be understood.


##Suggested citation for acknowledging this work:
##... We are submitting a paper that used this script to Engineering Geology journal, about Murisengo gypsum quarry.

#end of introduction


### Function list - stored in the same folder as this script ###

#1)
source("ReadProsys.R") #This function reads a .bin binary file as saved from Syscal Pro (Iris instruments) through the prosysII or III software.

### End Function list ###



##Library list ###

library(ggplot2) #to create plots
library(scales)  # makes pretty labels on the x-axis
library(patchwork) #for creating subplots
library(plotly) #to create 3D plots of pseudosection.
library(RColorBrewer)
library(dplyr)
### End libraries ###

##Read the .bin data file

filename <- file.choose()  # This returns the full file path
raw_data=ReadProsys(filename)

#Read an already processed file in the format of r2 protocol.dat
# raw_data=read.table(choose.files(),sep="\t")
# colnames(raw_data)=c("num","Xa","Xb","Xm","Xn","VsuI")
## Data processing, step by step ##

clean_data=raw_data #this clean_data variable will be updated at each step. 
#The result of each processing step is saved also in a separate variable to allow inspection
clean_data$num=c(1:length(clean_data[,1]))

##remove beginning or end measurements - store into variable clean_data_removetails ##

initial_rows_to_remove=0
final_rows_to_remove=0

clean_data=clean_data[(1+initial_rows_to_remove):(length(clean_data[,1])-final_rows_to_remove),]
clean_data_removetails=clean_data
 

##Create A, B,M, N, containing the electrode numbers (generally they are similar to Xz, Xb,Xm,Xn)

add_initial_sequence_number=0

clean_data$A=clean_data$Xa+add_initial_sequence_number
clean_data$B=clean_data$Xb+add_initial_sequence_number
clean_data$M=clean_data$Xm+add_initial_sequence_number
clean_data$N=clean_data$Xn+add_initial_sequence_number


## Adjust electrode numbering for roll-along or Streamer sequences. (when they are not similar to Xz Xb Xm, Xn)

# #Specify how many measurements are there per acquisition (how many measurements in the sequence you used)
# #you can put the total measurements if there is just one acquisition
# meas_per_acq=318 #this 318 is for a custom streamer sequence.
# 
# ##retrieve the electrode indexes instead of positions (it happens when you put the correct ele positions in the sequence instead of just 1,2,3...)
#   ele_positions=unique(as.vector(cbind(clean_data$A[1:meas_per_acq],clean_data$B[1:meas_per_acq],clean_data$M[1:meas_per_acq],clean_data$N[1:meas_per_acq])))
#   
#   ##Old for loop that explains that the following faster, vectorized commands do.
#   # for (ii in 1:length(clean_data[,1])){ #this for cycle replaces the electrode positions with their indexes.
#   #   #print(ii)
#   #   clean_data$A[ii]=which(ele_positions==clean_data$A[ii])
#   #   clean_data$B[ii]=which(ele_positions==clean_data$B[ii])
#   #   clean_data$M[ii]=which(ele_positions==clean_data$M[ii])
#   #   clean_data$N[ii]=which(ele_positions==clean_data$N[ii])
#   #   
#   # }
#   
# ##Retrieve the indexes of ABMN in ele_positions - basically converts positions into indexes. Needed if you start with ele positions instead of indexes in ABMN. see above non vectorized for loop and comments above.
#   clean_data$A= match(clean_data$A,ele_positions)
#   clean_data$B= match(clean_data$B,ele_positions)
#   clean_data$M= match(clean_data$M,ele_positions)
#   clean_data$N= match(clean_data$N,ele_positions)
#   
#   ## assign updated electrode indexes to roll-along or streamer measurements after the first one.
#   #this is a vectorized operation based on the meas_per_acqu and lenght(ele_positions) variables. The latter is basically how many electrodes do you have per acquisition
#   clean_data$A= clean_data$A+floor((which(clean_data$A==clean_data$A)-1)/meas_per_acq)*length(ele_positions)
#   clean_data$B= clean_data$B+floor((which(clean_data$B==clean_data$B)-1)/meas_per_acq)*length(ele_positions)
#   clean_data$M= clean_data$M+floor((which(clean_data$M==clean_data$M)-1)/meas_per_acq)*length(ele_positions)
#   clean_data$N= clean_data$N+floor((which(clean_data$N==clean_data$N)-1)/meas_per_acq)*length(ele_positions)
#   
# 



##Reverse eventual cables that were acquired with a reverse order:

ele_inv=c(49:72)
ele_invinv=rev(ele_inv)
#correct ap


for (ii in 1:length(clean_data$A)){

  if(any(clean_data$A[ii]==ele_inv)){
    clean_data$A[ii]=ele_invinv[ele_inv==clean_data$A[ii]]+0.001
  }
  if(any(clean_data$B[ii]==ele_inv)){
    clean_data$B[ii]=ele_invinv[ele_inv==clean_data$B[ii]]+0.001
  }
  if(any(clean_data$M[ii]==ele_inv)){
    clean_data$M[ii]=ele_invinv[ele_inv==clean_data$M[ii]]+0.001
  }
  if(any(clean_data$N[ii]==ele_inv)){
    clean_data$N[ii]=ele_invinv[ele_inv==clean_data$N[ii]]+0.001
  }

}

ele_inv=c(73:96)
ele_invinv=rev(ele_inv)
#correct ap


for (ii in 1:length(clean_data$A)){

  if(any(clean_data$A[ii]==ele_inv)){
    clean_data$A[ii]=ele_invinv[ele_inv==clean_data$A[ii]]+0.001
  }
  if(any(clean_data$B[ii]==ele_inv)){
    clean_data$B[ii]=ele_invinv[ele_inv==clean_data$B[ii]]+0.001
  }
  if(any(clean_data$M[ii]==ele_inv)){
    clean_data$M[ii]=ele_invinv[ele_inv==clean_data$M[ii]]+0.001
  }
  if(any(clean_data$N[ii]==ele_inv)){
    clean_data$N[ii]=ele_invinv[ele_inv==clean_data$N[ii]]+0.001
  }

}

clean_data$A=floor(clean_data$A)
clean_data$B=floor(clean_data$B)
clean_data$M=floor(clean_data$M)
clean_data$N=floor(clean_data$N)

# 
# 
# #
# 
ele_num_to_add=48*0 #in the CSV file, to adjust the ele numbering of this survey.
clean_data$A=clean_data$A+ele_num_to_add
clean_data$B=clean_data$B+ele_num_to_add
clean_data$M=clean_data$M+ele_num_to_add
clean_data$N=clean_data$N+ele_num_to_add
# 
# raw_data=clean_data #so now raw_data contain only the correct measurements, and have ABMN columns that are needed for the understanding of faulty electrodes.
# 

  
  
  
## Read the real position of each electrode from a four-column csv file (columns should be electrode number,x,y,z)

ele_pos=read.csv(file.choose())
# ele_pos=read.csv("Crescentino ert3d_all_nxyz_local.csv")




### WARNING ! update ele_pos$Name for streamer data ### put number of previous measurements instead of 0.
# ele_pos$Name=ele_pos$Name-(0)*length(ele_positions)


##check z if it is correct (sometimes there are errors in z taken with gps or by sampling a raster DTM)
plot(ele_pos$z)
plot(smooth(smooth(ele_pos$z,twiceit = TRUE)))

#perform a little smoothing on z
ele_pos$z=smooth(ele_pos$z)


clean_data$Xa=ele_pos[match(clean_data$A,ele_pos$Name),2]
clean_data$Xb=ele_pos[match(clean_data$B,ele_pos$Name),2]
clean_data$Xm=ele_pos[match(clean_data$M,ele_pos$Name),2]
clean_data$Xn=ele_pos[match(clean_data$N,ele_pos$Name),2]

clean_data$Ya=ele_pos[match(clean_data$A,ele_pos$Name),3]
clean_data$Yb=ele_pos[match(clean_data$B,ele_pos$Name),3]
clean_data$Ym=ele_pos[match(clean_data$M,ele_pos$Name),3]
clean_data$Yn=ele_pos[match(clean_data$N,ele_pos$Name),3]

clean_data$Za=ele_pos[match(clean_data$A,ele_pos$Name),4]
clean_data$Zb=ele_pos[match(clean_data$B,ele_pos$Name),4]
clean_data$Zm=ele_pos[match(clean_data$M,ele_pos$Name),4]
clean_data$Zn=ele_pos[match(clean_data$N,ele_pos$Name),4]


clean_data_elepos=clean_data




## Calculate the geometric factor k and, subsequently, the apparent rho (which in the raw_data is calculated based on the raw_data Xz,Xb,Xm,Xn that may have been changed by the user in the previous processing step)


### Calculate k ###

clean_data$rAM=((clean_data$Xa-clean_data$Xm)^2+ (clean_data$Ya-clean_data$Ym)^2 + (clean_data$Za-clean_data$Zm)^2)^0.5
clean_data$rBM=((clean_data$Xb-clean_data$Xm)^2+ (clean_data$Yb-clean_data$Ym)^2 + (clean_data$Zb-clean_data$Zm)^2)^0.5
clean_data$rAN=((clean_data$Xa-clean_data$Xn)^2+ (clean_data$Ya-clean_data$Yn)^2 + (clean_data$Za-clean_data$Zn)^2)^0.5
clean_data$rBN=((clean_data$Xb-clean_data$Xn)^2+ (clean_data$Yb-clean_data$Yn)^2 + (clean_data$Zb-clean_data$Zn)^2)^0.5

clean_data$k=2*pi/(1/clean_data$rAM-1/clean_data$rBM-1/clean_data$rAN+1/clean_data$rBN)
clean_data$Rho=clean_data$Vmn/clean_data$Iab*clean_data$k
# clean_data$Rho=clean_data$VsuI*clean_data$k

clean_data_k=clean_data

# clean_data=clean_data_k
### REMOVE SINGLE FAULTY ELECTRODES 
# faul_ele_remove=c(291, 306, 331 ,334)
# clean_data=clean_data[!(clean_data$A%in%faul_ele_remove | clean_data$B%in%faul_ele_remove | clean_data$M%in%faul_ele_remove | clean_data$N%in%faul_ele_remove),]


# ## Check electrodes by calculating the reciprocal error (for reciprocal measurements with current and potential electrodes switched, like 1-2-3-4; 3-4-1-2)
# 
# #First, let's order ABMN in a way A is always lower than B and M is always lower than N
# 
# clean_data$Aa=apply(data.frame(clean_data$A,clean_data$B), 1, FUN = min)# that 1 means row-wise. if 2, the function is applied column wise
# clean_data$Bb=apply(data.frame(clean_data$A,clean_data$B), 1, FUN = max)
# clean_data$Mm=apply(data.frame(clean_data$M,clean_data$N), 1, FUN = min)
# clean_data$Nn=apply(data.frame(clean_data$M,clean_data$N), 1, FUN = max)
# 
# #Now we should find the reciprocal values and set them apart in a different list, in order to calculate their rho app difference.
# 
# cd_spot_dup1=data.frame(AA=append(clean_data$Aa,clean_data$Mm),BB=append(clean_data$Bb,clean_data$Nn),MM=append(clean_data$Mm,clean_data$Aa),NN=append(clean_data$Nn,clean_data$Bb),Rho=append(clean_data$Rho,clean_data$Rho))
# cd_spot_dup2=cd_spot_dup1[duplicated(cd_spot_dup1[,1:4])| duplicated(cd_spot_dup1[,1:4],fromLast = TRUE),]
# #finding duplicates works now, but it may be enhanced to be more general. now it only finds two duplicates, but they may exist more.
# 
# # cd_spot_dup2=cd_spot_dup1[ duplicated(cd_spot_dup1[,1:4]),]
# # cd_spot_dup3=unique(cd_spot_dup2)
# # cd_spot_dup4=cd_spot_dup1[cd_spot_dup1[,1:4] %in% cd_spot_dup3[,1:4],]
# 
# #calculate error
# cd_spot_dup3=cd_spot_dup2 %>%
#   arrange(AA,BB,MM,NN, Rho) %>%
#   group_by(AA,BB,MM,NN) %>%
#   mutate(Error = abs((last(Rho) - first(Rho))/(last(Rho)+first(Rho)+0.01)))
# 
# cd_spot_dup3=cd_spot_dup3[!duplicated(cd_spot_dup3[,1:4]),]
# 
# #calculate error associated with each electrode:
# list_ele=sort(unique(append(cd_spot_dup3$AA,cd_spot_dup3$BB)))
# 
# ele_error=data.frame(ele=list_ele,error=0)
# for (ii in 1:length(cd_spot_dup3$AA)){
#   ele_error$error[cd_spot_dup3$AA[ii]]= ele_error$error[cd_spot_dup3$AA[ii]]+cd_spot_dup3$Error[ii]
#   ele_error$error[cd_spot_dup3$BB[ii]]= ele_error$error[cd_spot_dup3$BB[ii]]+cd_spot_dup3$Error[ii]
#   ele_error$error[cd_spot_dup3$MM[ii]]= ele_error$error[cd_spot_dup3$MM[ii]]+cd_spot_dup3$Error[ii]
#   ele_error$error[cd_spot_dup3$NN[ii]]= ele_error$error[cd_spot_dup3$NN[ii]]+cd_spot_dup3$Error[ii]
#   
# }
# 
# plot(ele_error$ele,ele_error$error)
# #filter bad electrodes
# bad_electrodes=c()
# clean_data=clean_data[!((clean_data$A %in% bad_electrodes)|(clean_data$B %in% bad_electrodes)|(clean_data$M %in% bad_electrodes)|(clean_data$N %in% bad_electrodes)),]
# clean_data_bad_ele=clean_data
# #end checking electrode errors by means of reciprocal values.
# 


## Plot histogram of several variables and allow the user to filter for outliers

## plot histogram of k ##

breaks=c(0,10,100,1000,10000,100000,1000000)

plotk=ggplot() + 
  geom_histogram(aes(x=abs(clean_data$k))) + 
  ggplot2::ggtitle(sprintf("%s < k < %s",round(min(abs(clean_data$k)),3),round(max(abs(clean_data$k)),3)))+
  
  scale_x_log10(
    breaks = breaks,
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  theme_bw()
## Histogram of Rho
plotRho=ggplot() + 
  geom_histogram(aes(x=(clean_data$Rho))) + 
  # scale_x_log10(
   ggplot2::ggtitle(sprintf("%s < Rho < %s",round(min(clean_data$Rho),3),round(max(clean_data$Rho),3)))+
    # labels = scales::trans_format("log10", scales::math_format(10^.x))
  # )+
  theme_bw()


## Histogram of Dev. Rho

plotDevRho=ggplot() + 
  geom_histogram(aes(x=(clean_data$DevRho)+0.1)) +
  ggplot2::ggtitle(sprintf("%s < Dev Rho < %s",round(min(clean_data$DevRho),3),round(max(clean_data$DevRho),3)))+
  
   scale_x_log10(
  
   labels = scales::trans_format("log10", scales::math_format(10^.x))
   )+
  theme_bw()

## Histogram of Iab
plotIab=ggplot() + 
  geom_histogram(aes(x=(clean_data$Iab))) + 
  ggplot2::ggtitle(sprintf("%s < Iab < %s",round(min(clean_data$Iab),3),round(max(clean_data$Iab),3)))+
  
  # scale_x_log10(
  
  # labels = scales::trans_format("log10", scales::math_format(10^.x))
  # )+
  theme_bw()

## Histogram of Vmn
plotVmn=ggplot() + 
  geom_histogram(aes(x=abs(clean_data$Vmn))) + 
  ggplot2::ggtitle(sprintf("%s < Vmn < %s",round(min(abs(clean_data$Vmn)),3),round(max(abs(clean_data$Vmn)),3)))+
  
   scale_x_log10(
  
   labels = scales::trans_format("log10", scales::math_format(10^.x))
 )+
  theme_bw()

#Make a single plot with the many histograms.
(plotRho+plotDevRho)/(plotIab+plotVmn)/plotk


##Filter the values based on the histograms

#Filter on Rho
min_rho=1
max_rho=3000
clean_data=clean_data[(clean_data$Rho)>=min_rho & (clean_data$Rho)<=max_rho,]

#Filter on Dev Rho
min_devrho=0
max_devrho=10
clean_data=clean_data[(clean_data$DevRho)>=min_devrho & (clean_data$DevRho)<=max_devrho,]

# #Filter on Iab
min_iab=1
max_iab=1000
clean_data=clean_data[(clean_data$Iab)>=min_iab & (clean_data$Iab)<=max_iab,]
# 
# #Filter on Vmn
min_vmn=0.5
max_vmn=15000 ##Syscal Pro we used has an overload protection system that does not allow to record voltages > 15V
clean_data=clean_data[(abs(clean_data$Vmn))>=min_vmn & (abs(clean_data$Vmn))<=max_vmn,]

# #Filter on k
min_k=0
max_k=10^5
clean_data=clean_data[abs(clean_data$k)>=min_k & abs(clean_data$k)<=max_k,]





##Plot pseudosection ##

## CALCULATE APPROXIMATE POSITIONS OF PSEUDOSECTION POINTS, BY AVERAGING X, Y AND A CUSTOM FORMULA FOR Z ###

clean_data$pseudoX=(clean_data$Xa+clean_data$Xb+clean_data$Xm+clean_data$Xn)/4
clean_data$pseudoY=(clean_data$Ya+clean_data$Yb+clean_data$Ym+clean_data$Yn)/4
clean_data$pseudoZ=(clean_data$Za+clean_data$Zb+clean_data$Zm+clean_data$Zn)/4

#calculate an estimated investigation depth that goes well both for Wenner-Schlumberger and Dipole Dipole and free 3D measurements
#Ztemp[i]= (  (((XA-mean(c(XN,XM)))^2  +  (YA-mean(c(YN,YM)))^2)^0.5) +  (((XB-mean(c(XN,XM)))^2  +  (YB-mean(c(YN,YM)))^2)^0.5) )/2 * 1/5 #average of spacing times a constant
clean_data$pseudoDepth= clean_data$pseudoZ-(  (((clean_data$Xa-clean_data$Xm)^2  +  (clean_data$Ya-clean_data$Ym)^2)^0.5) +  (((clean_data$Xb-clean_data$Xn)^2  +  (clean_data$Yb-clean_data$Yn)^2)^0.5) )/2 * 1/4 #average of spacing times a constant

# 
# ### ADD SOME RANDOM NOISE TO AVOID HAVING MULTIPLE SUPERIMPOSED POINTS IN THE FINAL VISUALIZATION ###
# Xtemp[i]=Xtemp[i]+runif(1,-Ztemp[i]/5,Ztemp[i]/5)
# Ytemp[i]=Ytemp[i]+runif(1,-Ztemp[i]/5,Ztemp[i]/5)
# Ztemp[i]=Ztemp[i]+runif(1,-Ztemp[i]/5,Ztemp[i]/5)




##Exterminate bad data points ##

fig_baddata <- ggplot(clean_data, aes(label = num)) + 
  geom_point(aes(x = Rho, y= pseudoDepth))+
  scale_x_log10()+
  theme_bw()+
  labs(title="Exterminate bad data points")+
  xlab("Apparent resistivity (Ohm * m)")+
  ylab("Pseudo elevation (m)")

ggplotly(fig_baddata)


#plot pseudosection. colorscales for Plotly available at https://plotly.com/r/builtin-colorscales/

fig <- plot_ly(clean_data, x = ~pseudoX, y = ~pseudoY, z = ~pseudoDepth,
               type = 'scatter3d', mode = 'markers',
marker = list(color = ~Rho, colorscale = 'Rainbow', showscale = TRUE))

# fig <- fig %>% add_markers()


fig <- fig %>% layout(scene = list(xaxis = list(title = 'East'),
                                   yaxis = list(title = 'North'),
                                   zaxis = list(title = 'Elevation')),
                      annotations = list(
                        x = 1.13,
                        y = 1.05,
                        text = 'Resistivity (Ohm * m)',
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                      ))
##add electrode positions
# fig <- fig %>%
#   add_markers(data=ele_pos_after_cleaning, name=~Name, x = ~x, y = ~y,z = ~z, marker=list(color="yellow",symbol="square"))

fig

##Export a pseudosection xyz points to visualize in other software

pseudosection=round(data.frame(clean_data$pseudoX,y=clean_data$pseudoY,z=clean_data$pseudoDepth,rho=clean_data$Rho),3)
write.csv(pseudosection,sprintf("output/pseudosection_%s.csv",basename(filename)),row.names = FALSE)


### Understand faulty electrodes based on what you filtered ###

# 
# #update raw_data and clean_data to eliminate the faulty electrodes, then repeat the count
# faulty_ele=c(205,216,217,235,237,238)
# all_ele=(unique(c(raw_data$A,raw_data$B,raw_data$M,raw_data$N)))
#  faulty_ele=all_ele
#    faulty_ele=sus_eles
# faulty_ele_combs=t(combn(faulty_ele,4))
# faulty_ele_combs=cbind(faulty_ele_combs)
# #function to speed up next for cycle
# 
# fun_countele_unfilt <- function(x, character = FALSE) {
#   sum(ele_unfilt==x)
# }
# 
# fun_countele_filt <- function(x, character = FALSE) {
#   sum(ele_filt==x)
# }
# 
# 
# ## for cycle to understand the best combination of electrodes to remove.
# variation=vector()
# for(cc in 1:length(faulty_ele_combs[,1])){
#   print(cc)
#   raw_dataF=raw_data[!((raw_data$A%in%faulty_ele_combs[cc,])|(raw_data$B%in%faulty_ele_combs[cc,])|(raw_data$M%in%faulty_ele_combs[cc,])|(raw_data$N%in%faulty_ele_combs[cc,])),]
#   clean_dataF=clean_data[!((clean_data$A%in%faulty_ele_combs[cc,])|(clean_data$B%in%faulty_ele_combs[cc,])|(clean_data$M%in%faulty_ele_combs[cc,])|(clean_data$N%in%faulty_ele_combs[cc,])),]
#   
#   
#   #Count how many of each electodes are there in the unfiltered file
#   
#   ele_unfilt=(c(raw_dataF$A,raw_dataF$B,raw_dataF$M,raw_dataF$N))
#   ele_unfilt_num=sort(unique(ele_unfilt))
#   
#   ele_unfilt_stats=data.frame(ele=ele_unfilt_num)
#   ele_unfilt_stats$num=apply(ele_unfilt_stats,1,fun_countele_unfilt)
#     
#   # ele_unfilt_stats=data.frame() #older slower version, not vectorized
#   # for (ele in ele_unfilt_num){
#   #   count=sum(ele_unfilt==ele)
#   #   ele_unfilt_stats=rbind(ele_unfilt_stats,c(ele,count))
#   # }
#   # 
#   
#   #Count how many of each electodes are there in the filtered file
#   ele_filt=(c(clean_dataF$A,clean_dataF$B,clean_dataF$M,clean_dataF$N))
#   ele_filt_num=sort(unique(ele_filt))
#   
#   ele_filt_stats=data.frame(ele=ele_filt_num)
#   ele_filt_stats$num=apply(ele_filt_stats,1,fun_countele_filt)
#   
#   # ele_filt_stats=data.frame()
#   # for (ele in ele_filt_num){
#   #   count=sum(ele_filt==ele)
#   #   ele_filt_stats=rbind(ele_filt_stats,c(ele,count))
#   # }
#   
#   # ele_filt_stats=ele_filt_stats[-15,]
#   ### calculate  the percentage filtering of each electrode
#   # variation[cc]=sum((ele_unfilt_stats[,2]-ele_filt_stats[,2])/ele_unfilt_stats[,2]*100) #alternative function to minimize: percentage
#   variation[cc]=sum((ele_unfilt_stats[,2]-ele_filt_stats[,2]))
#   
# 
# }
# 
# var_df=data.frame(n=1:length(variation),var=variation)
# var_df=var_df[order(var_df$var),]
# plot(var_df)
# sprintf("the best combination is")
# print(faulty_ele_combs[which.min(variation),])
# sprintf("the minimum variation is")
# print(min(variation))
# 
# suspected_eles_indexes=var_df$n[1:10]
# sus_eles=sort(unique(c(faulty_ele_combs[suspected_eles_indexes,])))
# 
# #plot the variation in order to spot eventual faulty electrodes (those that have the highest variation, often an outlier with respect to the other electrodes)
# #this plots shows on the y axis the variation percentage, and each electrode shows its number and how many of them were filtered.
# 
# # plot(ele_unfilt_stats[,1],variation)
# ggplot()+
#   annotate("text",ele_unfilt_stats[,1],variation,label=sprintf("%s-%s",ele_unfilt_stats[,1],ele_unfilt_stats[,2]-ele_filt_stats[,2]),size=3,color="black")+
#   theme_bw()
# sum(ele_unfilt_stats[,2]-ele_filt_stats[,2])/4 #this sum is the minimum measurements to take out
# 






# 
# raw_dataF=raw_data[!((raw_data$A%in%faulty_ele)|(raw_data$B%in%faulty_ele)|(raw_data$M%in%faulty_ele)|(raw_data$N%in%faulty_ele)),]
# clean_dataF=clean_data[!((clean_data$A%in%faulty_ele)|(clean_data$B%in%faulty_ele)|(clean_data$M%in%faulty_ele)|(clean_data$N%in%faulty_ele)),]
# 
# 
# #Count how many of each electodes are there in the unfiltered file
# 
# ele_unfilt=(c(raw_dataF$A,raw_dataF$B,raw_dataF$M,raw_dataF$N))
# ele_unfilt_num=sort(unique(ele_unfilt))
# 
# ele_unfilt_stats=data.frame()
# for (ele in ele_unfilt_num){
#   count=sum(ele_unfilt==ele)
#   ele_unfilt_stats=rbind(ele_unfilt_stats,c(ele,count))
# }
# 
# #Count how many of each electodes are there in the filtered file
# ele_filt=(c(clean_dataF$A,clean_dataF$B,clean_dataF$M,clean_dataF$N))
# ele_filt_num=sort(unique(ele_filt))
# 
# ele_filt_stats=data.frame()
# for (ele in ele_filt_num){
#   count=sum(ele_filt==ele)
#   ele_filt_stats=rbind(ele_filt_stats,c(ele,count))
# }
# 
# # ele_filt_stats=ele_filt_stats[-15,]
# ### calculate  the percentage filtering of each electrode
# variation=(ele_unfilt_stats[,2]-ele_filt_stats[,2])/ele_unfilt_stats[,2]*100
# 
# #plot the variation in order to spot eventual faulty electrodes (those that have the highest variation, often an outlier with respect to the other electrodes)
# #this plots shows on the y axis the variation percentage, and each electrode shows its number and how many of them were filtered.
# 
# # plot(ele_unfilt_stats[,1],variation)
# ggplot()+
#   annotate("text",ele_unfilt_stats[,1],variation,label=sprintf("%s-%s",ele_unfilt_stats[,1],ele_unfilt_stats[,2]-ele_filt_stats[,2]),size=3,color="black")+
#   theme_bw()
# sum(ele_unfilt_stats[,2]-ele_filt_stats[,2])/4 #this sum is the minimum measurements to take out

### End understand faulty electrodes ###




##here you can save temporarily a dataframe of clean_data, to append later
# clean_data_all=data.frame()

# clean_data_3D3=clean_data
# clean_data_3D1$Nn<-NULL to remove a column from a dataframe
# clean_data_all=rbind(clean_data_3D1,clean_data_3D2,clean_data_3D3)
# clean_data=clean_data_all




### MERGE CLOSE  ELECTRODES ###
# 
# clean_data_backup=clean_data
# clean_data=clean_data_backup
# ele_pos_after_cleaning_backup=ele_pos_after_cleaning
# ele_pos_after_cleaning=ele_pos_after_cleaning_backup

#This option is to write a protocol.dat file in which the close electrodes are indicated by a single index. This is to avoid difficulties in mesh creation because there are too-close electrodes.
#For example, one could group all the electrodes within a 0.2 m range.

merge_range = 0.2 # meters

#update the topography file: it may not contain all the points, because some electrodes could have been totally rejected in the previous filtering process.
ele_pos_after_cleaning=ele_pos[ele_pos$Name %in% as.vector(rbind(clean_data$A,clean_data$B,clean_data$M,clean_data$N)) ,]


### preliminary - do also later -  Rename all electrode from 1 to n , avoiding any hole in the data - make ele_pos_after_filtering$name be equal to 1:length(ele_pos_after_filtering$name ###

clean_data$A= match(clean_data$A,ele_pos_after_cleaning$Name)
clean_data$B= match(clean_data$B,ele_pos_after_cleaning$Name)
clean_data$M= match(clean_data$M,ele_pos_after_cleaning$Name)
clean_data$N= match(clean_data$N,ele_pos_after_cleaning$Name)

ele_pos_after_cleaning$Name=1:length(ele_pos_after_cleaning$Name)



#search for close electrodes and make a list of close electrode numbers.
list_close_ele=list()
countlist=0

for (iii in 1:length(ele_pos_after_cleaning$Name)){
  if(ele_pos_after_cleaning$Name[iii] %in% unlist(list_close_ele)==0){ #only if they are not already in the list. unlist is needed to transform the list in a simple vector
    countlist=countlist+1
    temp_dist=((ele_pos_after_cleaning$x-ele_pos_after_cleaning$x[iii])^2+(ele_pos_after_cleaning$y-ele_pos_after_cleaning$y[iii])^2)^0.5
    list_close_ele[[countlist]]=which(temp_dist<merge_range)[(which(temp_dist<merge_range)%in% unlist(list_close_ele))==0]
    
    if(length(list_close_ele[[countlist]])>1){
      for(jjj in list_close_ele[[countlist]]){ #one should create a self-nested function for this.
        temp_dist=((ele_pos_after_cleaning$x-ele_pos_after_cleaning$x[jjj])^2+(ele_pos_after_cleaning$y-ele_pos_after_cleaning$y[jjj])^2)^0.5
        list_close_ele[[countlist]]=append( list_close_ele[[countlist]],which(temp_dist<merge_range)[(which(temp_dist<merge_range)%in% unlist(list_close_ele))==0])
        
        
        if(length(list_close_ele[[countlist]])>1){
          for(kkk in list_close_ele[[countlist]]){ #one should create a self-nested function for this.
            temp_dist=((ele_pos_after_cleaning$x-ele_pos_after_cleaning$x[kkk])^2+(ele_pos_after_cleaning$y-ele_pos_after_cleaning$y[kkk])^2)^0.5
            list_close_ele[[countlist]]=append( list_close_ele[[countlist]],which(temp_dist<merge_range)[(which(temp_dist<merge_range)%in% unlist(list_close_ele))==0])
            
            
            if(length(list_close_ele[[countlist]])>1){
              for(www in list_close_ele[[countlist]]){ #one should create a self-nested function for this.
                temp_dist=((ele_pos_after_cleaning$x-ele_pos_after_cleaning$x[www])^2+(ele_pos_after_cleaning$y-ele_pos_after_cleaning$y[www])^2)^0.5
                list_close_ele[[countlist]]=append( list_close_ele[[countlist]],which(temp_dist<merge_range)[(which(temp_dist<merge_range)%in% unlist(list_close_ele))==0])
                
              }
            }
            
            
          }
        }
        
        
      }
    }
    
  }
}


# Now, we have a list of all the INDEXES of the electrodes to be merged. To know the name of the electrode,
#we should check the ele_pos_after_cleaning$Name at that index. 

#We now have to rename all the electrodes in each element of the list as the first of the n elements, and also update the 
#ele_pos_after_cleaning variable, using e.g. the average position of these points.

for(element in list_close_ele){
if(length(unlist(element))>1){#if there is at least 2 elements in the list entry - it means there are close electrodes to be merged
  clean_data$A[clean_data$A %in% unlist(element)]=unlist(element)[1]
  clean_data$B[clean_data$B %in% unlist(element)]=unlist(element)[1]
  clean_data$M[clean_data$M %in% unlist(element)]=unlist(element)[1]
  clean_data$N[clean_data$N %in% unlist(element)]=unlist(element)[1]
  
  ele_pos_after_cleaning$x[unlist(element)]=mean(ele_pos_after_cleaning$x[unlist(element)])
  ele_pos_after_cleaning$y[unlist(element)]=mean(ele_pos_after_cleaning$y[unlist(element)])
  ele_pos_after_cleaning$z[unlist(element)]=mean(ele_pos_after_cleaning$z[unlist(element)])
 
   ele_pos_after_cleaning$Name[unlist(element)]=unlist(element)[1]
  
}
  
}

length(unique(as.vector(rbind(clean_data$A,clean_data$B,clean_data$M,clean_data$N))))==length(list_close_ele)

#now we have to eliminate from the electrode position list all the duplicate electrodes.
# ele_pos_after_cleaning=ele_pos_after_cleaning[(!duplicated(ele_pos_after_cleaning$x) & !duplicated(ele_pos_after_cleaning$y)),] #not working properli probably due to numerical approximations.
ele_pos_after_cleaning=ele_pos_after_cleaning[(!duplicated(ele_pos_after_cleaning$Name)),]

### END MERGE CLOSE ELECTRODES ###


### Rename all electrode from 1 to n , avoiding any hole in the data - make ele_pos_after_filtering$name be equal to 1:length(ele_pos_after_filtering$name ###
#needed even if you did already before the merging of electrodes
clean_data$A= match(clean_data$A,ele_pos_after_cleaning$Name)
clean_data$B= match(clean_data$B,ele_pos_after_cleaning$Name)
clean_data$M= match(clean_data$M,ele_pos_after_cleaning$Name)
clean_data$N= match(clean_data$N,ele_pos_after_cleaning$Name)

ele_pos_after_cleaning$Name=1:length(ele_pos_after_cleaning$Name)



#export topography file (it may have smoothed z compared to the input file)
write.csv(ele_pos_after_cleaning,sprintf("output/ele_pos_%s.csv",basename(filename)),row.names = FALSE)

###   WARNING!!!  ### 

# The exported topography may not contain all the points, because some electrode could have been totally rejected in the process, or merged.
#This should not a problem in importing into resipy, since the electrode numbering between protocol.dat file and topo file should match.




##Export into protocol.dat format (ResIPy or R2 or R3)

write.table(data.frame(1:length(clean_data$A),clean_data$A,clean_data$B,clean_data$M,clean_data$N,clean_data$Vmn/clean_data$Iab),sprintf("output/protocol_%s_.dat",basename(filename)),row.names = FALSE,sep="\t",col.names = FALSE)
# write.table(data.frame(1:length(clean_data$A),clean_data$A,clean_data$B,clean_data$M,clean_data$N,clean_data$VsuI),sprintf("protocol_%s_.dat",basename(filename)),row.names = FALSE,sep="\t",col.names = FALSE)

# 
# ## concatenate multiple protocol.dat##
conc_files=choose.files()
output_protocol=data.frame()
for (fil in conc_files){
  output_protocol=rbind(output_protocol,read.csv(fil,sep='\t',col.names = c("num","a","b","m","n","dVsuI")))
}
  output_protocol$num=1:length(output_protocol$num)

length(sort(unique(c(output_protocol$a,output_protocol$b,output_protocol$m,output_protocol$n)))) #to check that every electrode is still there#
#if not, consider to add fake measurements in order to make resipy positioning work.

  write.table(output_protocol,sprintf("%output/s_all.dat",basename(filename)),row.names = FALSE,sep="\t",col.names = FALSE)
# #   
# Export into Res2dinv .dat format