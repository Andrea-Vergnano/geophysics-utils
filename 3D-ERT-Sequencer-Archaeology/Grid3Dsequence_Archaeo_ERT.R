
#Hello! This script creates an ERT sequence for a 3D survey setup with electrodes at the surface located in a m x n grid.
#This sequence is particularly optimized if you have a georesistivimeter such as IRIS Syscal Pro that measures 10 potential differences simultaneously,
#speeding up the survey.
#The goal of an archaeological 3D survey is to acquire as many data as possible near the surface, until 1-2 meters depth, for archaeological purposes. Therefore,
#we generally want to focus on very shallow investigation depths, with a very dense data distibution near the surface.

library(readr)
#First, we select the survey geometry.

### INPUT PARAMETERS ###
m=12 #electrodes in a row
n=12 #electrodes in a column
ele_dist=2 #distance between electodes in the grid
num_pot=20 #number of potential measurements per each current dipole (a multiple of 10 is suggested to obtain maximum optimization in Syscal Pro georesistivimeter or similar)
meandering=TRUE #if TRUE, the script will order the electrode positions as a snake. If you instead prefer to make the survey grid with each cable's lowest electrode number always start from the same side of the grid, put FALSE.
random_ele_pot=FALSE #experimental: just pick random potential electrodes assigned to each current dipole


##Since the focus of ERT surveys in archaeology is often to scan the first meters of subsoil,the program focuses on potential electrodes nearby and around the current electrodes
##First, the program sets the AB spacings, which can be set in "ws_aperture_list_hor" and "ws_aperture_list_ver" variables. # Optionally, you can add AB with a diagonal spacing (AB adjacent in diagonal)
## Second, according to those spacings, it will create a certain number (=num_pot) of potential electrodes to be measured during that AB current injection


#Please, set to TRUE the following variables if you want to add them to the sequence.
diagonal=FALSE
ws=TRUE

##Please, set the following variables to generate the AB current electrodes with different spacings between them. For a simple case, use just ws_aper... = c(1)
#If you want AB to have a large aperture set for example =c(10). You can set a free combination of AB apertures, for example =c(1,5,8), and they can be different in horizontal and vertical directions.

# ws_aperture should be less than m or n! - it represents how much spacing between A and B current electrodes for 3D configuration. Practically, it is equal to B-A in the case the electrode grid has a spacing of 1m
ws_aperture_list_hor=c(1,2,3) #here we define how many different AB apertures to be selected in horizontal direction
ws_aperture_list_ver=c(1,2,3)#here we define how many different AB apertures to be selected in vertical direction

### END INPUT PARAMETERS ###


#Calculation of the electrode grid based on the input parameters.
ele_pos=data.frame(num=1:(m*n),x=floor((0:(m*n-1))/m),y=0:(m*n-1)-m*floor((0:(m*n-1))/m)) #calculates the electrodes position in the grid
ele_pos$x=ele_pos$x*ele_dist #multiply the electrode positions for the spacing between electrodes.
ele_pos$y=ele_pos$y*ele_dist #multiply the electrode positions for the spacing between electrodes.



#creation of the possible AB (current electrodes) combinations

print("creating the current electrodes")

#initialize the variables

dd_diag_right=data.frame()
dd_diag_left=data.frame()

ws_ap_hor=data.frame()
ws_ap_ver=data.frame()


if(diagonal==TRUE){
  
  dd_diag_right=data.frame(A=1:(m*(n-1)-1),B=(m+2):(m*n))
  dd_diag_left=data.frame(A=(m+1):(m*n-1),B=2:(m*(n-1)))
}

ws_all=data.frame()
if(ws==TRUE) {

  for (ws_aperture in ws_aperture_list_hor){
 ws_ap_hor=data.frame(A=1:(m*n-ws_aperture),B=(ws_aperture+1):(m*n)) #wenner shlumberger current electrodes in x direction, only spacings = ele_dist*3 are considered (the minimum wenner-schlumberger aperture)
 ws_ap_hor=ws_ap_hor[(ws_ap_hor$A%%m %in% c(1:(m-ws_aperture))),]#remove some AB that we do not want because they are far away between each other
 ws_all=rbind(ws_all,ws_ap_hor)
  }
  
  for (ws_aperture in ws_aperture_list_ver){
    
 ws_ap_ver=data.frame(A=1:(m*n),B=1:(m*n)+m*(ws_aperture)) #wenner shlumberger current electrodes in y direction, only spacings = ele_dist*3 are considered (the minimum wenner-schlumberger aperture)
 ws_ap_ver=ws_ap_ver[ws_ap_ver$B<=m*n,] #remove some AB that we do not want because they are far away between each other
 ws_all=rbind(ws_all,ws_ap_ver)
  }
}

 
#This following line puts together all the current AB dipoles that we created above and that we want to calculate the corresponding potential MN electrodes of.
#If you set to FALSE some ot these AB dipoles in the input parameters, they will not be added anyway to the "ap" variable because they will be empty variables. Instead. if you modified the script and
# you add other kinds of AB dipoles, do not forget to add them here below too.
 

ap=rbind(dd_diag_right,dd_diag_left,ws_all) #new code just uses ws_all for everything.  It still does not support diagonal AB apertures.




#initialize the data frame that will contain the final sequence
abmn=data.frame(A=0,B=0,M=0,N=0)


#The following for loop is the creation of the MN potential electrodes associated with each current electrode.

print("creating the potential electrodes")

if(random_ele_pot==FALSE){ #this is the normal case, if random is FALSE
for (i in 1:length(ap[,1])){ # For each AB current dipole:
  
  #determine the x y position of the center of the dipole
  xA=ele_pos$x[ele_pos$num==ap[i,1]]
  xB=ele_pos$x[ele_pos$num==ap[i,2]]
  yA=ele_pos$y[ele_pos$num==ap[i,1]]
  yB=ele_pos$y[ele_pos$num==ap[i,2]]
  xAB=(xA+xB)/2
  yAB=(yA+yB)/2
  
  #determine the distance matrix between the current dipole and the potential potential electrodes
  
  dists=((ele_pos$x-xAB)^2+(ele_pos$y-yAB)^2)^0.5
  dists[(ele_pos$num==ap[i,1]) | (ele_pos$num==ap[i,2])]=9999*ele_dist # set a high value of distance to the current electrodes in order to not consider them in the following calculations of the potential electrodes
  

  
  ele_pot=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes associated to each current dipole
  M=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes M associated to each current dipole
  N=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes N associated to each current dipole
  
  
  ele_pot[1]= which(dists==min(dists))[1] #The first chosen potential electrodes is one of the electrodes nearest to the center of the current dipole, but not the electrodes of the current dipole itself.
  
  #find all the other potential electrodes and store them temporarily in M and N variable
  for (j in 2:(num_pot+1))  {
   
    #calculate another distance matrix, which is the distance of each electrode from the previous potential electrode
    xJ=ele_pos$x[ele_pos$num==ele_pot[j-1]]
    yJ=ele_pos$y[ele_pos$num==ele_pot[j-1]]
    distj=((ele_pos$x-xJ)^2+(ele_pos$y-yJ)^2)^0.5
    distj[which(ele_pos$num %in% ele_pot[1:j-1])]=9999*ele_dist #set a high number for the electrodes already used (which are stored in ele_pot), so they will not be chosen again
    
    #select the electrode nearest to the current dipole among the nearests to the previous potential electrode, and not already used.
    #this algorithm basically will select all the potential electrodes around the first potential electrodes, going farther and farther from the center of the current dipole (xAB,yAB)
    ele_pot[j]= ele_pos$num[ele_pos$num==which(dists+distj*2==min(dists+distj*2))[1]][1]
   
    #create the M and N vectors to be used with the current AB dipole (which is ap[i,])
      M[j-1] =ele_pot[j-1]
      N[j-1]=ele_pot[j]
    
  }
 
   #create a temporary abmn dataframe containing the quadripoles related to the current AB current electrode
  abmntemp=data.frame(A=ap[i,1],B=ap[i,2],M=M,N=N)
  
  #update the general ABMN dataframe. At the end of the for loop, which end is the next line, this will be the final electrode sequence.
   abmn=rbind(abmn,abmntemp) 
}

abmn=abmn[-1,] #remove the first line that was just 0 because I put it because I always need to find some tricks since I cannot code so well.

}else{ #if random_ele_pot==TRUE
  
  for (i in 1:length(ap[,1])){ # For each AB current dipole:
    
    possible_mn=c(1:(m*n))[! c(1:(m*n)) %in% c(ap$A[i],ap$B[i])]
    #generate five random integers between 1 and 20 (sample without replacement)
      
    mnrand=sample(possible_mn, num_pot+1, replace=FALSE) #if we use replace=FALSE then we do not allow the same integer to be generated more than once.
    
    
    mrand=mnrand[1:num_pot]
    nrand=mnrand[2:(num_pot+1)]
    abmntemp=data.frame(A=ap$A[i],B=ap$B[i],M=mrand,N=nrand)
    abmn=rbind(abmn,abmntemp) 
  }
  abmn=abmn[-1,] #remove the first line that was just 0 because I put it because I always need to find some tricks since I cannot code so well.
  
}

###reorder the electrode positions and the AB dipoles in order to make it meandering, like a snake. 
#Basically, every second row has to be inverted both in electrode positions and in the abmn
#It may take some time.


if(meandering==TRUE){
  print("reordeing the sequence to make the electrode positions meandering")
  
  abmn=as.matrix(abmn)
  for (k in 1:n){
    
    if (k%%2==0){ #this k%%2 selects every second row to be reversed
      
      #this creates the variables that store the electrodes index to be reversed.
      ele_inv=(m*(k-1)+1):(m*k)
      ele_invinv=rev(ele_inv)

      for (ii in 1:length(abmn[,1])){
        
        if(ii%%100==0){
          print(sprintf("processing line %s of %s, quadripole %s",k,n,ii))
        }
        if(any(abmn[ii,1]==ele_inv)){
          abmn[ii,1]=ele_invinv[ele_inv==abmn[ii,1]]+0.1 #This +0.001 and then the floor function is a trick that I found to easily replace one value by another. Maybe it can be done more efficiently.
        }
        if(any(abmn[ii,2]==ele_inv)){
          abmn[ii,2]=ele_invinv[ele_inv==abmn[ii,2]]+0.1
        }
        if(any(abmn[ii,3]==ele_inv)){
          abmn[ii,3]=ele_invinv[ele_inv==abmn[ii,3]]+0.1
        }
        if(any(abmn[ii,4]==ele_inv)){
          abmn[ii,4]=ele_invinv[ele_inv==abmn[ii,4]]+0.1
        }
        
      }
      abmn=floor(abmn)
      
      #Similarly, this for cycle corrects the electrode positions
      for (ii in 1:length(ele_pos$num)){
        
        if(any(ele_pos$num[ii]==ele_inv)){
          ele_pos$num[ii]=ele_invinv[ele_inv==ele_pos$num[ii]]+0.001
          
          
          
        }
        ele_pos$num=floor(ele_pos$num)
        
      }
      
    }
    
  }
  
}

ele_pos=ele_pos[order(ele_pos$num),]



#save the sequence in ResIPy and Electre formats

ele_pos$z=0

### export .csv files containing the sequences. This can be imported in Resipy - Forward modeling ###
write.table(abmn,sprintf("Seq_grid_%ix%i.csv",m,n),row.names = FALSE,quote=FALSE,sep=",")
write.table(ele_pos,sprintf("Ele_pos_grid%ix%i.csv",m,n),row.names = FALSE,quote=FALSE,sep=",")

### export .txt files containing the sequences ###
#These can be opened in Electre pro to import them into the Syscal Pro.
#In electre pro, please save them to .sqz file type, in order to save the desired Q factor, stacking options, voltage options...

colnames(ele_pos)=c("#","X","Y","Z")
write.table(ele_pos,sprintf("Grid_%ix%i.txt",m,n),row.names = FALSE,quote=FALSE,sep="\t")
sequ_txt=cbind(c(1:length(abmn[,1])),abmn)
colnames(sequ_txt)[1]="#"
write.table(sequ_txt,sprintf("Grid_%ix%i.txt",m,n),row.names = FALSE,quote=FALSE,sep="\t",append=TRUE)





###old code ###
# # 
# # # check plot the distance array
#  library(ggplot2)
# 
# ggplot()+
#   geom_point(aes(x=ele_pos$Y,y=ele_pos$X,color=(dists)),size=10)+
#   #scale_color_distiller(palette=1)+
#   annotate("text",x=ele_pos$Y,y=ele_pos$X,label=sprintf("1-%s",ele_pos$`#`),size=4,color="darkred")+
#   annotate("text",x=ele_pos$Y+1,y=ele_pos$X+1,label=sprintf("2-%s",ele_pos$`#`),size=4, color="darkgreen")+
#   
#   xlab("X distance (m)")+
#   ylab("Y distance (m)")+
#   theme_bw()+
#   theme(axis.title.x = element_text(color="black", size=24),axis.title.y = element_text(color="black", size=24))+
#  theme(  axis.text = element_text(size=20))
# 
# 
# 
# 






