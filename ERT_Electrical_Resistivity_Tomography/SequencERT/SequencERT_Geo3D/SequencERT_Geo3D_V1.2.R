## Grid3Dsequence_GeoArchaeo_ERT V1.2 ###

## Release notes ##

#V 1.0 June 2025: Creation of this script for a geologic/mining investigation with Electrical Resistivity Tomography (ERT)

#V 1.1 October 2025:
                    #now the probability function that selects MN electrodes varies with AB distance. Now you can have more skewed MN electrodes distribution when AB are close to each other, while maintaining a more uniform distribution whan AB are distant from each other
                    #now you can define a range around the center of the AB dipole in which all the potential electrodes will be surely selected by the probability function
                    #fixed sampling probability by reducing 0.001 to 0.00000001
                    #added some histogram plots to increase knowledge on the sequence created.
                    #I did this update to investigate an embankment with a 3D configuration.
#V 1.2 December 2025 #output directory updated

## End of Release notes


### Introduction

#Hello! This script creates an ERT sequence for a 3D survey setup with electrodes at the surface located in a m x n grid.
#This sequence is particularly optimized if you have a georesistivimeter such as IRIS Syscal Pro that measures 10 potential differences simultaneously,
#speeding up the survey.
#This script was adapted from another script of mine, focused on archaeological 3D surveys.
#( The goal of an archaeological 3D survey is to acquire as many data as possible near the surface, until 1-2 meters depth, for archaeological purposes. Therefore,
#we generally want to focus on very shallow investigation depths, with a very dense data distribution near the surface. )

#Instead, here we focus on geologic surveys, so we want to go deep with the investigation.
#This script is therefore optimized for 3D rectangular-shape grids (e.g. 3 cables of 48 electrodes each, or 6 x 24), and it is meant to select
#some meaningful AB spacings (see the variables ws_aperture_list...) and, for each AB spacing, a series of and MN potential electrodes that 
#are more spaced the more distant from the center of the current dipole AB. This allows to have a high investigation depth while maintaining the geometrical factor k quite low.

library(readr)

#First, we select the survey geometry.

### INPUT PARAMETERS ###
m=48 #electrodes in a row (or cable...if you have your rows exactly of the length of the cable)
n=3 #electrodes in a column = number of rows (or cables)
ele_dist_x=8 #distance between rows
ele_dist_y=1 #distance between electrodes in a row (or cable electrode spacing)
num_pot=20 #number of potential measurements per each current dipole (a multiple of 10 is suggested to obtain maximum optimization in Syscal Pro georesistivimeter or similar)
meandering=FALSE #if TRUE, the script will order the electrode positions as a snake. If you instead prefer to make the survey grid with each cable's lowest electrode number always start from the same side of the grid, put FALSE.
random_ele_pot=FALSE #experimental: just pick random potential electrodes assigned to each current dipole
delete_inline_measurements=FALSE #if you perform already line by line measurements with traditional 2D sequences, and you wanna save space and time.

##First, the program sets the AB spacings, which can be set in "ws_aperture_list_hor" and "ws_aperture_list_ver" variables. # Optionally, you can add AB with a diagonal spacing (diag_aperture_list)
## Second, according to those spacings, it will create a certain number (=num_pot) of potential electrodes to be measured during that AB current injection


#Please, set to TRUE the following variables if you want to add them to the sequence. Diagonal means AB current electrodes that are spaced diagonally within the
#electrode grid. ws stands for Wenner-Schlumberger but in this script contains also dipole dipole electrodes, basically it is all those AB dipoles which are not spaced diagonally, but only along a row or column.
diagonal=TRUE
ws=TRUE

##Please, set the following variables to generate the AB current electrodes with different spacings between them. For a simple case, use just ws_aper... = c(1)
#If you want AB to have a large aperture set for example =c(10). You can set a free combination of AB apertures, for example =c(1,5,8), and they can be different in horizontal and vertical directions.

# ws_aperture should be less than m or n! - it represents how much spacing between A and B current electrodes for 3D configuration. Practically, it is equal to B-A in the case the electrode grid has a spacing of 1m
ws_aperture_list_hor=c(1,3,45,46,47) #here we define how many different AB apertures to be selected in horizontal direction
ws_aperture_list_ver=c(1)#here we define how many different AB apertures to be selected in vertical direction

#diagonals should be a list of c(firstnumber,secondnumber), like list(c(1,2),c(2,10), c(3,20). First number should be less than n, second number less than m
#these numbers means, for example c(2,10), that the diagonal AB is made by a certain series of A electrodes, and the B electrodes are in positions, relative to A, of 2*ele_dist_x on the x direction, and of 10*ele_dist_y on the y direction
diag_aperture_list=list(c(1,8),c(2,16),c(2,47))

## input parameters for the choice of the potential electrodes.

#As stated at the beginning of the script, for each current electrode AB we need to select a meaningful series (possibly a multiple of 10, in the case of optimized Syscal Pro georesistivimeter)
#of MN potential electrodes. The idea, to minimize k but allow deep investigation at the same time, is to select MN more spaced and more spaced, as they are more distant from the center of the AB current dipole.

#At first,a random weighted distribution (with sample(x,num,prob) function of R) was tried, with the weights inversely proportional to the distance from the center of AB
#However it was too random and it was probable that the most distant electrodes were underrepresented.
#Therefore, a custom algorithm to build a "less-random" spatial distribution of electrodes was created.
#This algorithm basically divides the candidates to be potential electrodes (all electrodes but the current electrodes) in several groups (num_pot +1). Each group has more elements if it is more distant
#from the center of the AB dipole. The choice of these groups is linket to the parameters mult_fact and elev_fact. For each group, the algorithm selects just one potential electrodes, randomly, but with a weight. 
#This weight is higher for the more distant electrodes within the group (this weight is linked to the parameter elev_fact_inv). 

#in V1.1 this was converted to a range, based on the distance between A and B. 

#if you want to choose different parameters, we suggest you to try them in the commented trial section below, meant for testing and plotting the 
# quasi-random potential electodes grid. The default values of 1.5, 1.8, 11 are a good balance and worked well for a 48x3 grid, where the 3 lines were spaced 5 times the electrode spacing.

# mult_fact=1.5 # this parameter, if lower, it will skew more the electrode distribution
# elev_fact=1.8 # this parameter, if higher, it will skew more the electrode distribution
# elev_fact_inv=11 # this parameter, if higher, it will choose more consistently the most distant electrodes - the potential electrodes distribution will be "less random".

mult_fact_range=c(1.1,2) # this parameter, if lower, it will skew more the electrode distribution
elev_fact_range=c(2,1.1) # this parameter, if higher, it will skew more the electrode distribution
ele_dist_range=c(min(ele_dist_x,ele_dist_y),(((m-1)*ele_dist_y)^2+((n-1)*ele_dist_x)^2)^0.5) # a vector with two elements: the distance of the closest couple of electrodes and the distance of the most distant couple of electrodes.
elev_fact_inv=2 # this parameter, if higher, it will choose more consistently the most distant electrodes - the potential electrodes distribution will be "less random".
range_sure_ele_range=c(4,0.1) #m #in this range of distance from the center of the current dipole, all the electrodes will be picked. This range is a range, for larger aperture it shrinks, because for large apertures we are not interested in selecting all these electrodes near the AB dipole center


### END INPUT PARAMETERS ###


#Calculation of the electrode grid based on the input parameters.
ele_pos=data.frame(num=1:(m*n),x=floor((0:(m*n-1))/m),y=0:(m*n-1)-m*floor((0:(m*n-1))/m)) #calculates the electrodes position in the grid
ele_pos$x=ele_pos$x*ele_dist_x #multiply the electrode positions for the spacing between electrodes.
ele_pos$y=ele_pos$y*ele_dist_y #multiply the electrode positions for the spacing between electrodes.

# #example plot of the electrode grid
# plot(ele_pos$x,ele_pos$y)
# 


### trial quasi-random sampling on elecrodes - see later where it is really implemented in the code. use  ###
# #this lines for try the parameters of sampling and plot them visually.
# 
# y0=500
# x0=0
#  prob_samp=((ele_pos$y-y0)^2+(ele_pos$x-x0)^2)^0.5
# # prob_samp=abs((ele_pos$y-y0))#for basing distance only on y
# 
# prob_samp=max(prob_samp)*mult_fact-(prob_samp)
# prob_samp=prob_samp^elev_fact
# prob_samp=prob_samp/sum(prob_samp)*num_pot
# prob_inv=1/prob_samp^elev_fact_inv
# 
# s_i=1
# s_i_start=1
# ele_sampled=vector()
# pr_count=0
# while (s_i<=m*n){
# 
#   pr_count=pr_count+prob_samp[s_i]
# 
#   if(pr_count>1 | s_i%%m==0 ){
# 
#     if(s_i==s_i_start){
#       ele_sampled=append(ele_sampled,s_i)
#     }else{
#     # ele_sampled=append(ele_sampled,sample(c(s_i_start:s_i),1,prob=prob_samp[s_i_start:s_i]))
#     # ele_sampled=append(ele_sampled,sample(c(s_i_start:s_i),1))
#       ele_sampled=append(ele_sampled,sample(c(s_i_start:s_i),1,prob=prob_inv[s_i_start:s_i]))
#     }
#     s_i_start=s_i+1
#     pr_count=pr_count-1
# 
#   }
# 
#   s_i=s_i+1
# }
# plot(ele_pos$x[ele_pos$num %in% ele_sampled],ele_pos$y[ele_pos$num %in% ele_sampled])

### End of: trial quasi-random sampling on elecrodes - see later where it is really implemented in the code.###




#creation of the possible AB (current electrodes) combinations

print("creating the current electrodes")

#initialize the variables

dd_diag_right=data.frame()
dd_diag_left=data.frame()
dd_diag_all=data.frame()

ws_ap_hor=data.frame()
ws_ap_ver=data.frame()

ws_all=data.frame()

# this for cycle creates the various AB combinations which are "diagonal". It is a basic algorithm that performs basic logic and geometrical calculations.
if(diagonal==TRUE){
  
  # dd_diag_right=data.frame(A=1:(m*(n-1)-1),B=(m+2):(m*n))
  # dd_diag_left=data.frame(A=(m+1):(m*n-1),B=2:(m*(n-1)))
  for(diag_aperture in diag_aperture_list){
    dd_diag_right=data.frame(A=1:(m*(n-diag_aperture[1])-diag_aperture[2]),B=(m*diag_aperture[1]+1+diag_aperture[2]):(m*n))
    dd_diag_right=dd_diag_right[(abs(floor((dd_diag_right$B-1)/m)-floor((dd_diag_right$A-1)/m))==diag_aperture[1]),]#this deletion is needed because the algorithm I use to create them creates too many of them
    
    dd_diag_left=data.frame(A=(m*diag_aperture[1]+1):(m*n-diag_aperture[2]),B=(1+diag_aperture[2]):(m*(n-diag_aperture[1])))
    dd_diag_left=dd_diag_left[(abs(floor((dd_diag_left$B-1)/m)-floor((dd_diag_left$A-1)/m))==diag_aperture[1]),]#this deletion is needed because the algorithm I use to create them creates too many of them
    
    
     dd_diag_all=rbind(dd_diag_all,dd_diag_left,dd_diag_right)
  }
}

# this for cycle creates the various AB combinations which are not "diagonal". It is a basic algorithm that performs basic logic and geometrical calculations.

if(ws==TRUE) {

  for (ws_aperture in ws_aperture_list_hor){
 ws_ap_hor=data.frame(A=1:(m*n-ws_aperture),B=(ws_aperture+1):(m*n)) #wenner shlumberger current electrodes in x direction, only spacings = ele_dist*3 are really WS (the minimum wenner-schlumberger aperture), if lower they are actually dipole-dipole electrodes.
 ws_ap_hor=ws_ap_hor[(ws_ap_hor$A%%m %in% c(1:(m-ws_aperture))),]#remove some AB that we do not want because they are far away between each other
 ws_all=rbind(ws_all,ws_ap_hor)
  }
  
  for (ws_aperture in ws_aperture_list_ver){
    
 ws_ap_ver=data.frame(A=1:(m*n),B=1:(m*n)+m*(ws_aperture)) #wenner shlumberger current electrodes in y direction, only spacings = ele_dist*3 are really WS (the minimum wenner-schlumberger aperture), if lower they are actually dipole-dipole electrodes.)
 ws_ap_ver=ws_ap_ver[ws_ap_ver$B<=m*n,] #remove some AB that we do not want because they are far away between each other
 ws_all=rbind(ws_all,ws_ap_ver)
  }
}

 
#This following line puts together all the current AB dipoles that we created above and that we want to calculate the corresponding potential MN electrodes of.
#If you set to FALSE some of these AB dipoles in the input parameters, they will not be added anyway to the "ap" variable because they will be empty variables. Instead. if you modified the script and
# you add other kinds of AB dipoles, do not forget to add them here below too.
 

ap=rbind(ws_all,dd_diag_all) #new code just uses ws_all and dd_diag_all




#initialize the data frame that will contain the final sequence
abmn=data.frame(A=0,B=0,M=0,N=0)


#The following for loop is the creation of the MN potential electrodes associated with each current electrode.

print("creating the potential electrodes")

if(random_ele_pot==FALSE){ #this is the normal case, if random is FALSE
for (i in 1:length(ap[,1])){ # For each AB current dipole:
  
  # i=1
  # ap[i,1]=43
  # ap[i,2]=48
  # ele_pos$num=ele_pos$`#`
  # ele_pos$y=ele_pos$Y
  # 
  # ele_pos$x=ele_pos$X
  
  #determine the x y position of the center of the dipole
  xA=ele_pos$x[ele_pos$num==ap[i,1]]
  xB=ele_pos$x[ele_pos$num==ap[i,2]]
  yA=ele_pos$y[ele_pos$num==ap[i,1]]
  yB=ele_pos$y[ele_pos$num==ap[i,2]]
  xAB=(xA+xB)/2
  yAB=(yA+yB)/2
  
    #determine the distance between A and B
  ele_dist_AB=((xA-xB)^2+(yA-yB)^2)^0.5+1
  
  ###select a series of num_pot potential electrodes distributed quasi-randomly (more uniform than random) - weighted on distance
  #we want to select the potential electrodes spread all over the domain, but more closest to the center of the current dipole, and order them cable by cable simply.
  
  
  y0=yAB
  x0=xAB

  #remove the current electrodes from the potential potential electrode list
  ele_pos_temp=ele_pos[(ele_pos$num!=ap[i,1] & ele_pos$num!=ap[i,2]),]
  
  ##Calculate a probability to pick the electrode as potential electrode for this AB current dipole 
  #(less probable if more distant from the center of the current dipole)
  
  
  #First, calculate mult_fact and elev_fact from their possible ranges, based on distance of the couple of AB current electrodes. The more distant they are, the more uniform distribution (less skewed): we want to obtain to reduce the geometrical factor k.
  #it is a linear model depending on the distance range "ele_dist_range". We want to calculate elev_fact and mult_fact as functions of ele_dist_range
  #calculate slope and q (slope and intercept of the linear model)
  
  slope= (elev_fact_range[2]-elev_fact_range[1])/(ele_dist_range[2]-ele_dist_range[1])
  q=elev_fact_range[1]-(ele_dist_range[1]*slope)
  
  elev_fact=slope*ele_dist_AB+q
  
  #same for mult_fact
  slope= (mult_fact_range[2]-mult_fact_range[1])/(ele_dist_range[2]-ele_dist_range[1])
  q=mult_fact_range[1]-(ele_dist_range[1]*slope)
  
  mult_fact=slope*ele_dist_AB+q
  
  
  #same for range_sure_ele
  slope= (range_sure_ele_range[2]-range_sure_ele_range[1])/(ele_dist_range[2]-ele_dist_range[1])
  q=range_sure_ele_range[1]-(ele_dist_range[1]*slope)
  
  range_sure_ele=slope*ele_dist_AB+q
  
  
  ##target sum of prob_samp is num_pot+1+0.001, and no prob_samp > 1 obviously.
  
  
  prob_samp=((ele_pos_temp$y-y0)^2+(ele_pos_temp$x-x0)^2)^0.5-min(ele_dist_range)
  
  
  # prob_samp=abs((ele_pos_temp$y-y0))#for basing distance only on y
  prob_samp_sure=which(prob_samp<=range_sure_ele)
  prob_samp_notsure=which(prob_samp>range_sure_ele)
  prob_samp_target_sum=num_pot+1.0000001
  prob_samp_target_sum_corr=prob_samp_target_sum-length(which(prob_samp<range_sure_ele))*1
  
  prob_samp=max(prob_samp)*mult_fact-(prob_samp) # adjust this probability with the mult_fact, which can be selected at the beginning of the script.
  prob_samp=prob_samp^elev_fact  # adjust this probability with the elev_fact, which can be selected at the beginning of the script.
  

  prob_samp[prob_samp_sure]=1
  prob_samp[prob_samp_notsure]=prob_samp[prob_samp_notsure]/sum(prob_samp[prob_samp_notsure])*prob_samp_target_sum_corr
   # prob_samp=prob_samp/sum(prob_samp)*(num_pot+1+0.0001)#this normalization will ensure that the potential electrodes will be exactly num_pot, as desired. +0.0001 avoids errors due to numerical approximations.
  
  
  # plot(prob_samp)
  prob_inv=1/prob_samp^elev_fact_inv # this prob_inv is an inverse of the prob. It is meant later, to make most distant electrodes be more represented in the final distribution. It can be selected at the beginning of the script.
  
  #initialize variables for choosing the quasi-random distance-based electrode distribution
  s_i_start=1
  ele_sampled=vector()
  pr_count=0
  
  #calculate the quasi-random distance-based electrode distribution.
  #Basically, to create a random but homogeneous distribution, this code
  #divides the electrodes in groups, which contain more electrodes if those groups are more distant from the AB center
  #therefore, we will obtain a distribution that have many electrodes near AB center, and less far away.
  #However, we want that the most distant ones are well represented in this distribution, to increase investigation depth.
  #Therefore, among each electrode group, we select most probably the more distant one (this is the meaning of prob_inv variable)
  
  for (s_i in 1:length(ele_pos_temp$num)){
    # s_i=s_i+1
    
    pr_count=pr_count+prob_samp[s_i]
    
    #each electrode has a probability. when summing adjacent electrode probabilities, at some point we reach 1. 
    #when we reach 1, we select randomly one of those electrodes and move to the next electrode.
    
    if(pr_count>=1 | ele_pos_temp$num[s_i]%%m==0 ){ #this ele_pos_temp$num[s_i]%%m==0 ensures that when the counter arrives at the end of a line, it will select one of the last electrodes of the line.
      
      if(s_i==s_i_start){ #if there is just one electrode in the group, select it
        ele_sampled=append(ele_sampled,ele_pos_temp$num[s_i])
      }else{ #if there are more than 1 electrode in the group, select it randomly, but more probably the more distant from AB (see prob_inv)
      
        ele_sampled=append(ele_sampled,ele_pos_temp$num[sample(c(s_i_start:s_i),1,prob=prob_inv[s_i_start:s_i])])
      }
      s_i_start=s_i+1
      pr_count=pr_count-1 #reset the probability after having added an electrode to the ele_sampled final list.
      
    }
    
  }
  
  M=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes M associated to each current dipole
  N=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes N associated to each current dipole

  for (j in 2:(length(ele_sampled)))  {

    #create the M and N vectors to be used with the current AB dipole (which is ap[i,])
      M[j-1] =ele_sampled[j-1]
      N[j-1]=ele_sampled[j]
  
  }
  # print(length(ele_sampled))
  length(abmn$A)
  
  ### This commented part is from the Archaeo script, meant to produce all those potential electrodes closest to the center of the current electrode.
  # #determine the distance matrix between the current dipole and the potential potential electrodes
  # 
  # dists=((ele_pos$x-xAB)^2+(ele_pos$y-yAB)^2)^0.5
  # dists[(ele_pos$num==ap[i,1]) | (ele_pos$num==ap[i,2])]=9999*ele_dist_x # set a high value of distance to the current electrodes in order to not consider them in the following calculations of the potential electrodes
  # 
  # 
  # 
  # ele_pot=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes associated to each current dipole
  # M=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes M associated to each current dipole
  # N=1:num_pot*0 #initialize the temporary storage variable for the potential electrodes N associated to each current dipole
  # 
  # 
  # ele_pot[1]= which(dists==min(dists))[1] #The first chosen potential electrodes is one of the electrodes nearest to the center of the current dipole, but not the electrodes of the current dipole itself.
  # 
  # #find all the other potential electrodes and store them temporarily in M and N variable
  # for (j in 2:(num_pot+1))  {
  #  
  #   #calculate another distance matrix, which is the distance of each electrode from the previous potential electrode
  #   xJ=ele_pos$x[ele_pos$num==ele_pot[j-1]]
  #   yJ=ele_pos$y[ele_pos$num==ele_pot[j-1]]
  #   distj=((ele_pos$x-xJ)^2+(ele_pos$y-yJ)^2)^0.5
  #   distj[which(ele_pos$num %in% ele_pot[1:j-1])]=9999*ele_dist #set a high number for the electrodes already used (which are stored in ele_pot), so they will not be chosen again
  #   
  #   #select the electrode nearest to the current dipole among the nearests to the previous potential electrode, and not already used.
  #   #this algorithm basically will select all the potential electrodes around the first potential electrodes, going farther and farther from the center of the current dipole (xAB,yAB)
  #   ele_pot[j]= ele_pos$num[ele_pos$num==which(dists+distj*2==min(dists+distj*2))[1]][1]
  #  
  #   #create the M and N vectors to be used with the current AB dipole (which is ap[i,])
  #     M[j-1] =ele_pot[j-1]
  #     N[j-1]=ele_pot[j]
  #   
  # }
  ### End of: This commented part is from the Archaeo script, meant to produce all those potential electrodes closest to the center of the current electrode.
  
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


if(delete_inline_measurements==TRUE){
  abmn=abmn[((floor((abmn$A-1)/m)==floor((abmn$B-1)/m)) & (floor((abmn$B-1)/m) ==floor((abmn$M-1)/m)) & (floor((abmn$M-1)/m)==floor((abmn$N-1)/m)))==FALSE,]
}

### save the sequence in ResIPy and Electre formats ###


#first, create an output folder
ifelse(!dir.exists(file.path("output")),
       dir.create(file.path("output")),
       "Directory Exists")



ele_pos$z=0

### export .csv files containing the sequences. This can be imported in Resipy - Forward modeling ###
write.table(abmn,sprintf("output/Seq_grid_%ix%i.csv",m,n),row.names = FALSE,quote=FALSE,sep=",")
write.table(ele_pos,sprintf("output/Ele_pos_grid%ix%i.csv",m,n),row.names = FALSE,quote=FALSE,sep=",")


sequ_txt=cbind(c(1:length(abmn[,1])),abmn)
colnames(sequ_txt)[1]="#"
### Calculate k for check ###
# sequ_txt$A=143
# sequ_txt$B=144
# sequ_txt$M=141
# sequ_txt$N=142
A=sequ_txt$A
B=sequ_txt$B
M=sequ_txt$M
N=sequ_txt$M
AA=ele_pos$x[sequ_txt$A]
MM=ele_pos$x[sequ_txt$M]

rAM=((ele_pos$x[A]-ele_pos$x[M])^2+ (ele_pos$y[A]-ele_pos$y[M])^2 + (ele_pos$z[A]-ele_pos$z[M])^2)^0.5

rAM=((ele_pos$x[sequ_txt$A]-ele_pos$x[sequ_txt$M])^2+ (ele_pos$y[sequ_txt$A]-ele_pos$y[sequ_txt$M])^2 + (ele_pos$z[sequ_txt$A]-ele_pos$z[sequ_txt$M])^2)^0.5
rBM=((ele_pos$x[sequ_txt$B]-ele_pos$x[sequ_txt$M])^2+ (ele_pos$y[sequ_txt$B]-ele_pos$y[sequ_txt$M])^2 + (ele_pos$z[sequ_txt$B]-ele_pos$z[sequ_txt$M])^2)^0.5
rAN=((ele_pos$x[sequ_txt$A]-ele_pos$x[sequ_txt$N])^2+ (ele_pos$y[sequ_txt$A]-ele_pos$y[sequ_txt$N])^2 + (ele_pos$z[sequ_txt$A]-ele_pos$z[sequ_txt$N])^2)^0.5
rBN=((ele_pos$x[sequ_txt$B]-ele_pos$x[sequ_txt$N])^2+ (ele_pos$y[sequ_txt$B]-ele_pos$y[sequ_txt$N])^2 + (ele_pos$z[sequ_txt$B]-ele_pos$z[sequ_txt$N])^2)^0.5
  
sequ_txt$k=2*pi/(1/rAM-1/rBM-1/rAN+1/rBN)


## plot histogram of k ##

library(ggplot2)
library(scales)  # makes pretty labels on the x-axis

breaks=c(0,10,100,1000,10000,100000,1000000)

ggplot() + 
  geom_histogram(aes(x=abs(sequ_txt$k))) + 
  scale_x_log10(
    breaks = breaks,
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  theme_bw()
sum(abs(sequ_txt$k[abs(sequ_txt$k)<10^10]))


#plot histogram of electrode distribution

hist(rbind(sequ_txt$A,sequ_txt$B))
hist(rbind(sequ_txt$N,sequ_txt$M))


### export .txt files containing the sequences ###
#These can be opened in Electre pro to import them into the Syscal Pro.
#In electre pro, please save them to .sqz file type, in order to save the desired Q factor, stacking options, voltage options...

colnames(ele_pos)=c("#","X","Y","Z")
write.table(ele_pos,sprintf("output/Grid_%ix%i.txt",m,n),row.names = FALSE,quote=FALSE,sep="\t")

write.table(sequ_txt[,1:5],sprintf("output/Grid_%ix%i.txt",m,n),row.names = FALSE,quote=FALSE,sep="\t",append=TRUE)


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






