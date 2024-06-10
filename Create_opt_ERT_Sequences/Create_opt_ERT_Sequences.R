
### This script helps to create standard dipole dipole sequences optimized for a multichannel ERT resistivimeter, such as the IRIS Syscal Pro
#it also create a roll-along sequence.
#each sequence (normal and roll along) also contains many repeated measurements to calculate the error model. No need to do another sequence with inverted AB and MN.

###INPUT PARAMETERS. CHECK THEM!!!###

channels=72 #number of channels or electrode of your instrument
roll_channels=24 #how many channels you want to go on when doing a roll along.
opt_ch=10 #number of measurements that the instrument can do at the same time in a dipole dipole sequence.
min_opt_ch=7 # put lower than opt_ch to add some quadripoles not completely optimized , but still performing at least min_opt_ch measurements at a time
max_a=10 #a is the spacing between the two current electrodes, and also between the two potential electrodes, in a standard dipole dipole configuration
time_per_measurement=2 #seconds needed for each current injection

### END INPUT PARAMETERS ###

#create a void dataframe to store the sequence in.

sequ=data.frame(A=integer(),B=integer(),M=integer(),N=integer())

##compute the standard dipole dipole sequence 
meas_count=0

for (a in 1:max_a) {
  
  if( (channels-min_opt_ch*a-a*2)>0 ) {
  
  for (i in 1:(channels-min_opt_ch*a-a*2)){
    
    meas_count=meas_count+1
    for (j in 1:opt_ch){
      
      A=i
      B=i+a
      M=i+a+j*a
      N=i+a+j*a+a
      sequ[nrow(sequ) + 1,] = list(A,B,M,N)
    }
    
  }
  
  }  
  
}
#remove measurements with one or more electrodes exceeding the max electrode number (=channels) (it was easier to create them and then remove them)
sequ <- sequ[!(sequ$A>channels | sequ$B>channels | sequ$M>channels | sequ$N>channels), ]



###add reverse sequence for repeated measurements
for (a in 1:max_a) {
  
  if( (channels-min_opt_ch*a-a*2)>0 ) {
    
    for (i in 1:(channels-min_opt_ch*a-a*2)){
      
      meas_count=meas_count+1
      for (j in 1:opt_ch){
        
        A=channels+1-(i)
        B=channels+1-(i+a)
        M=channels+1-(i+a+j*a)
        N=channels+1-(i+a+j*a+a)
        sequ[nrow(sequ) + 1,] = list(A,B,M,N)
      }
      
    }
    
  }  
  
}
#remove measurements with one or more electrodes whose number is lower than 1

sequ <- sequ[!(sequ$A<1 | sequ$B<1 | sequ$M<1 | sequ$N<1), ]



###calculate acquisition time based on number of measurements and optimization degree.

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
for (j in 2:length(sequ_roll[,1])){
  if((sequ_roll$A[j]==sequ_roll$A[j-1])&(sequ_roll$B[j]==sequ_roll$B[j-1])){
    
  }else{meas_count_roll=meas_count_roll+1}
}
acquisition_time_roll=meas_count_roll*time_per_measurement

opt_degree_roll=length(sequ_roll[,1])/(meas_count_roll*opt_ch)*100


#output a basic report about the created sequence
print(sprintf("Roll along sequence created with %i quadripoles and %i current injections.", length(sequ_roll[,1]),meas_count_roll))
print(sprintf("Your estimated acquisition time is %.01f minutes",acquisition_time_roll/60))
print(sprintf("Optimization percentage of this sequence is %.01f %%",opt_degree_roll))


#export .csv files containing the sequences. This can be imported in Resipy - Forward modeling
write.table(sequ,sprintf("DD%i.csv",channels),row.names = FALSE,quote=FALSE,sep=",")
write.table(sequ_roll,sprintf("DD%i_roll.csv",channels),row.names = FALSE,quote=FALSE,sep=",")

### export .txt files containing the sequences ###
#These can be opened in Electre pro to import them into the Syscal Pro.
#In electre pro, please save them to .sqz file type, in order to save the desired Q factor, stacking options, voltage options...

#write normal sequence
ele_pos=data.frame(num=1:channels,X=0:(channels-1),Y=0,Z=0)
colnames(ele_pos)[1]="#"
write.table(ele_pos,sprintf("DD%i.txt",channels),row.names = FALSE,quote=FALSE,sep="\t")
sequ_txt=cbind(c(1:length(sequ[,1])),sequ)
colnames(sequ_txt)[1]="#"
write.table(sequ_txt,sprintf("DD%i.txt",channels),row.names = FALSE,quote=FALSE,sep="\t",append = TRUE)

#write roll along sequence
write.table(ele_pos,sprintf("DD%i_roll.txt",channels),row.names = FALSE,quote=FALSE,sep="\t")
sequ_txt_roll=cbind(c(1:length(sequ_roll[,1])),sequ_roll)
colnames(sequ_txt_roll)[1]="#"
write.table(sequ_txt_roll,sprintf("DD%i_roll.txt",channels),row.names = FALSE,quote=FALSE,sep="\t",append = TRUE)



### END OF THE SCRIPT. have a nice day###