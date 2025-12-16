##ReordERT_V1.0.R##

## Release notes ##

#1.0 December 2025 Cleaned it up, commented it better, updated output directory, drafted the suggested citation, uploading it to github

## End of Release notes




### INTRODUCTION ###

### This script was first written in December 2024 by Andrea Vergnano. You can use it for any purpose by acknowledging my contribution.
# For any comments, you can mail me at andrea.vergnano@unito.it, or andrea.vergnano@posteo.eu

#Hello! This script takes any quadripole sequence for Electrical Resistivity Tomography as input, and optimizes it for:
# 1) for being able to remove the first cable as soon as possible, to begin dismantling the ERT line while the instrument is still acquiring
# 2) for not having any potential electrode just after being used as current electrode, to avoid polarization effects that induce measurements errors.


#Citation suggested: 
# Vergnano et al., 2025. Full-3D Electrical Resistivity Tomography (ERT) for mining activity planning and sinkhole assessment: a case study in the Murisengo gypsum quarry, 2025. Submitted to Engineering Geology.
# This citation is temporary, and hopefully this manuscript will be published in the next month
#

### Some further comments (in Italian)
# In pratica, prima riordina la sequenza in base all'elettrodo col numero più basso che trova in ogni gruppo di quadripoli acquisiti simultaneamente. 
# La sequenza ordinata per numero di elettrodi, ma non ancora ottimizzata per elettrodi di potenziale dopo quelli di corrente, la chiamo s. Il buffer, per ora vuoto, lo chiamo b, e la sequenza finale, per ora vuota, la chiamo so.
# Il codice ora mette tutti i gruppi di quadripoli che riesce dalla matrice b alla matrice so, se tutti gli MN rispettano la condizione di aspettare n_min iniezioni di corrente dopo essere stati usati come AB. Se non trova niente di utile in b, oppure b è vuoto (come nel caso iniziale), allora pesca da s per mettere un nuovo gruppo di misure dentro so, sempre che rispetti la condizione, altrimenti lo mette nel b. Dopo aver messo un nuovo gruppo di misure da s a so, ripete da capo a controllare tutti i gruppi di misure nel buffer b se rispettino la condizione e possano essere spostate a so. Ripete tutto questo fino a finire le misure. 
# Ne rimangono spesso alcune alla fine nel buffer che non riescono più a rispettare la condizione (non tante, generalmente meno di 10 iniezioni di corrente), che quindi appendo alla fine di so senza più controllare ... (forse si potrebbero spostare manualmente all'inizio, visto che son tutte quelle misure tipo 68 69 70 71 che rimangono al fondo, così di sicuro non interferiscono con le prime...)
# Un ciclo for prova diversi n_min (numero di iniezioni di corrente che vogliamo lasciar passare perché un elettrodo usato come A o B sia usato come M o N) e plotta un grafico da cui puoi capire quale valore di n_min ti piace di più, in base a quanto è alto e in base a quante misure restano al fondo. ripeti il codice con quell'n_min e ti salva il txt di electre.
# Conviene togliere il primo cavo non quando le misure arrivano all'elettrodo 25, ma conviene aspettare che ci sia l'elettrodo 26 o 27, visto che qualche 24 potrebbe essere rimasto nel buffer per un pochetto (non molto, dal buffer pesco sempre a partire dai primi, i più "vecchi", che essendo la sequenza s ordinata per numero di elettrodo, i primi del buffer sono anche quelli con il numero di elettrodo più basso).

### END INTRODUCTION ###




trial=1

### optional part to uncomment if you want a better optimization. Read below: ###
### UNCOMMENT AND RUN THIS PART THE FIRST TIME YOU WANT TO OPTIMIZE A CERTAIN SEQUENCE ###
### ALSO, UNCOMMENT THE RELATIVE PART AT THE END OF THE SCRIPT, BEFORE THE SAVE DATA ###
### ALSO, COMMENT n_min IN INPUT PARAMETERS ###
### This part is needed to understand the best n_min. 
### n_min is the minimum number of current injections you want to wait for using again a potential electrode (M or N) after it has been used as a current electrode (A or B).
### Run it (uncomment the following lines) the first time you want to optimize a sequence, then comment n_min in INPUT PARAMETERS. You will see a plot at the end, representing the
# 
# len_b_trial=vector()
# for (trial in 1:10){ #a trial and error method to estimate the best number of nmin: the one high enough that also leave the least number of measurements in the buffer b
#   n_min=trial
# print(sprintf("trial with n_min %s",trial))

### END trial to understand best n_min REMEMBER TO COMMENT ALSO n_min in the INPUT PARAMETERS and UNCOMMENT TTHE RELATIVE PART AT THE END OF THE SCRIPT BEFORE SAVE DATA

### INPUT PARAMETERS ###

ch_max=20 #max number of simultaneous measurements. It depends on the instrument. For example, the Syscal Pro by IRIS instruments can perform up to 10 measurements simultaneously, provided that they have the same current electrodes (A, B) and concatenated potential electrodes (M, N so that M[i]=N[i-1] )
n_min=3 #minimum number of current injections to wait for using again a potential electrode after it has been used as a current electrode.
n_min_ini=n_min #this copy variable is for knowing always what was n_min, since it is updated (lowered) near the end of the sequence.
input_sequence_file=(file.choose())
first_cable_length=24 #number of electrodes in the first cable you want to dismantle before the end of the survey

file_type="ELECTRE" #choose "ELECTRE" or "RESIPY". Electre is the electre pro txt, resipy is a simple csv with A, B,M,N columns.
howmany_ele=72 #only useful if "ELECTRE" is chosen above.
### END INPUT PARAMETERS ###


#Read data#
library(readr)
#

#the following is to import electre style txts.

if(file_type=="ELECTRE") {
  sequ <- read_delim(input_sequence_file, delim = "\t",escape_double = FALSE, trim_ws = TRUE, skip = howmany_ele+1)
  
}else{
  sequ=read.csv(input_sequence_file) #this is to import resipy-style csvs.
}

input_sequence_file_base=basename(input_sequence_file)

len=length(as.matrix(sequ[,1])) # stores the length of the input data for convenience.


### 1st step ###

### order the sequence in a way you can remove the first cable the soonest possible ###
### this part of the script was literally translated from my fortran code, so the syntax could be a little different from the rest of the code.

A=sequ$A
B=sequ$B
M=sequ$M
N=sequ$N

# Calculate the first index of each optimized group of measurements
groupcount <- 1
group <- c(1)

 c=1
 
for (i in 2:length(A)) {
 
  if (A[i] != A[i-1] | B[i] != B[i-1]| c>=ch_max | M[i] != N[i-1] ) {
    groupcount <- groupcount + 1
    group <- c(group, i)
    c=1
  } else{c=c+1}
}

group <- c(group, length(A) + 1)

# Calculate the minimum electrode index for each group of measurements
groupmin <- numeric(groupcount)
for (i in 1:groupcount) {
  groupminA <- min(A[group[i]:(group[i+1]-1)])
  groupminB <- min(B[group[i]:(group[i+1]-1)])
  groupminM <- min(M[group[i]:(group[i+1]-1)])
  groupminN <- min(N[group[i]:(group[i+1]-1)])
  groupmin[i] <- min(groupminA, groupminB, groupminM, groupminN)
}

# Rewrite the A, B, M, N arrays writing first those that have the least minimum electrode index
indmin <- 1
indmax <- 0
Aord <- numeric(length(A))
Bord <- numeric(length(B))
Mord <- numeric(length(M))
Nord <- numeric(length(N))

for (i in 1:len) {
  for (j in 1:groupcount) {
    if (groupmin[j] == i) {
      indmax <- indmin + group[j+1] - group[j] - 1
      Aord[indmin:indmax] <- A[group[j]:(group[j+1]-1)]
      Bord[indmin:indmax] <- B[group[j]:(group[j+1]-1)]
      Mord[indmin:indmax] <- M[group[j]:(group[j+1]-1)]
      Nord[indmin:indmax] <- N[group[j]:(group[j+1]-1)]
      indmin <- indmax + 1
    }
  }
}

### END  1st step ###
#update the old sequ matrix with the new ordered matrix, and store it in the variable called "s".

s=as.matrix(data.frame(A=Aord,B=Bord,M=Mord,N=Nord)) # store in s the intermediate sequence already ordered by electrode number but still not optimized for potential and current electrodes positions
so=matrix(nrow=0,ncol=4) #initialize output matrix containing the final sequence
b=matrix(nrow=0,ncol=4) #initialize a buffer matrix containing the part of the sequence that cannot be put at that time into the sequence, due to potential electrodes after current electrodes
so_g=matrix(nrow=0,ncol=3)#initialize a side variable of so that contains information about how many measurements are in each group
b_g=matrix(nrow=0,ncol=3)#initialize a side variable of b that contains information about how many measurements are in each group

### 2nd step ###
### reorder the sequence by allowing no potential electrode just after current electrode. ### 
### The number of current injections to wait before using a potential electrode after being used as current electrodes is n_min ###

i=1
ii=1

while (i<=len){ # This is the main routine that reorders the sequence. It copies quadrupole groups from s to so, but if finds a quadripole group in which even a single electrode was used as current electrode in the previous n_min quadripole groups, it stores it temporarily in the buffer b, which is scanned after every new entry from s into so
  print(sprintf("processing quadripole %s",i))
  ### Reduce n_min for late stages of the sequence. This allows to better empty the buffer b from the last remaining quadripole groups.
 # if(i>9/10*len){n_min=ceiling(n_min_ini*1/2)}
#  if(i>99/100*len){n_min=ceiling(n_min_ini*1/5)}
 # if(i>len-2*ch_max){n_min=1}
  
  
  flag1=1 #this flag indicates if at least one measurement group was taken from the buffer
  
  while(length(b[,1])>0 & flag1==1){ # Search in the buffer if there is some quadrupole groups to put into so, only if respect the condition of no potential electrodes just after being used as current electrodes
    
    flag1=0
    flag=1
    ii=1
    b_gcount=0
    
    while (ii<=length(b[,1]) & flag==1){#this while search from the first to the last measurement in the buffer and when found one, it puts it into so from b
      jj=1
      # print(ii)
      while(ii+jj<=length(b[,1]) && (b[ii,1]==b[ii+jj,1] & b[ii,2]==b[ii+jj,2] &  b[ii+jj,3]==b[ii+jj-1,4] & jj<ch_max)) { ## This while is checking accurately how many measurements after i belong to the same group ## && is short-circuited: if the first argument is false, it does not go on in evaluating the rest
        jj=jj+1
      }
      b_gcount=b_gcount+1
      kk=min(n_min,length(so_g[,1]))
      
      if( any(  b[ii:(ii+jj-1),3:4] %in% so[ so_g[(length(so_g[,1])-kk+1):length(so_g[,1]) ,3],1:2]) ){ #this if checks if the potential electrodes are present in previous groups of so, until n_max groups before.
        #do nothing and advance to the next group of buffer measurements
        ii=ii+jj# go to the next measurement (skipping all the measurement of the same group, since they were already put into so or b)
      }else{        #if there is a quadrupole group that is OK, move it from b into so
        so=rbind(so,b[ii:(ii+jj-1),])
        so_gtemp=rbind(so_g,c(length(so_g[,1])+1,jj ,0) )         #update so_groups so_g
        
        so_gtemp[length(so_gtemp[,1]),3]=sum(so_gtemp[,2])
        so_g=so_gtemp
        
        flag=0 # good, we found a measurement in the buffer to put in so, so we exit from scanning now, we update the buffer and then with flag1=1 we rescan again
        #update buffer b
        btemp=b[-(ii:(ii+jj-1)),]
        b=matrix(nrow=0,ncol=4)
        b=rbind(b,btemp)
        #update buffer groups b_g
        b_gt=matrix(nrow=0,ncol=3)
        b_gt=rbind(b_gt,b_g[-(b_g[,1]==b_gcount),])
        b_gt[,1]=c(1:length(b_gt[,1]))
        b_gt[,3]=cumsum(b_gt[,2])
        b_g=b_gt
        flag1=1 # since we found a measurement we scan again the buffer from beginning if we find another one
        #if no measurements were found, flag1 remain to 0, we exit from searching the buffer and we continue to take another measurement from s
        
        
      }
    }
    
    
    
    
  } #end while(length(b[,1])>0 & flag1=1) #end look in the buffer for something
  
  #begin taking data from the main s sequence. This part of the routine is activated only when there is no OK quadrupole groups left to take from the buffer b.
  j=1
  # print(i)
   while(i+j<=len && (s[i,1]==s[i+j,1] & s[i,2]==s[i+j,2] &  s[i+j,3]==s[i+j-1,4] & j<ch_max)) { ## This while is checking accurately how many measurements after i belong to the same group ## && is short-circuited: if the first argument is false, it does not go on in evaluating the rest
     j=j+1
   }
    #j now tells how many measurements are in this group of measurements which will be performed simultaneously by the instrument
    #i:i+j-1 are the indexes of the s measurements in this group
  
  if(length(so[,1])>0){ # if instead so is empty, just put the first quadrupole group from s to so
    k=min(n_min,length(so_g[,1]))
 
    # Check if in the previous groups in the output matrix so is already there that electrode used as current.
    
    if( any(  s[i:(i+j-1),3:4] %in% so[ so_g[(length(so_g[,1])-k+1):length(so_g[,1]) ,3],1:2]) ){ #this if checks if the potential electrodes are present in previous groups of so, until n_max groups before.
      #include in buffer in this case, and not put it into so.
      b=rbind(b,s[i:(i+j-1),])
      #update buffer info b_g
    
      b_gtemp=rbind(b_g,c(length(b_g[,1])+1,j ,0) )
      b_gtemp[length(b_gtemp[,1]),3]=sum(b_gtemp[,2])
      b_g=b_gtemp
      
      
    }else {
        #copy from s into so
     so=rbind(so,s[i:(i+j-1),])
      
        #update so_g, the variable that stores information about so.

     so_gtemp=rbind(so_g,c(length(b_g[,1])+1,j ,0) )
     so_gtemp[length(so_gtemp[,1]),3]=sum(so_gtemp[,2])
          so_g=so_gtemp     
        
      }
  
  }else{#if instead so is empty, just put the first quadrupole group from s to so
    so=rbind(so,s[i:(i+j-1),])
    so_gtemp=rbind(so_g,c(length(so_g[,1])+1,j ,0) )
    so_gtemp[length(so_gtemp[,1]),3]=sum(so_gtemp[,2])
    so_g=so_gtemp
  }
  
  i=i+j # go to the next measurement (skipping all the measurement of the same group, since they were already put into so or b)
}

# now so contains the final sequence, and b may contain some leftovers, which should be very few , and they are appended at the end of the sequence without any further scanning.

##recalculate b_g to be sure. in one case it was miscalculated before... not clear why
c=1
group <- c(1)
groupcount=1
if(length(b[,1])>1){
  for (i in 2:length(b[,1]))
    if (b[i,1] != b[i-1,1] | b[i,2] != b[i-1,2]| c>=ch_max | b[i,3] != b[i-1,4] ) {
      groupcount <- groupcount + 1
      group <- c(group, i)
      c=1
    } else{c=c+1}
  b_g[,1]=c(1:groupcount)
  b_g[,2]=c(diff(group),c)
  b_g[,3]=cumsum(b_g[,2])}else{
    b_g[,1]=1
    b_g[,2]=1
    b_g[,3]=1
  }

### append buffer leftover from the end of the sequence going reverse (before and after are reversed with respect to the common sense)###

if(length(b_g[,1])>0){
  
for (bgbg in 1:length(b_g[,1])){
  print(bgbg)
  #identify quadrupole group to be appended to the main sequence
  temp_g=matrix((b[(b_g[bgbg,3]-b_g[bgbg,2]+1):(b_g[bgbg,3]),]),ncol=4)
  
  #check where to put it into the main sequence, scanning from the beginning
  found=0
  temp_i=length(so_g[,1])-1 #this is the temporary index of where the new group from b could be inserted
  while (found==0){
    
    check1=TRUE
    check2=TRUE
    #check groups before temp_i ( we want to check a number of groups equal to n_min_ini)
    temp_so_g_before=min(length(so_g[,1]),temp_i+n_min_ini)
    temp_so_ind_bef=c((so_g[temp_i,3]+1):(so_g[temp_so_g_before,3])) #indices of so to be checked
    if(any(temp_g[,1:2] %in% so[temp_so_ind_bef,3:4])){check1=FALSE}
    
    ##check groups after temp_i ( we want to check a number of groups equal to n_min_ini) 
    temp_so_g_after=max(1,temp_i-n_min_ini)
    temp_so_ind_aft=c((so_g[temp_so_g_after,3]-so_g[temp_so_g_after,2]+1):(so_g[temp_i,3])) #indices of so to be checked
    if(any(temp_g[,3:4] %in% so[temp_so_ind_aft,1:2])){check1=FALSE}
    
    
    if(check1==FALSE || check2==FALSE)    {
      temp_i=temp_i-1
      if(temp_i==0){ # do not allow temp_i to go below 1, which would be an impossible index- In this case, reduce n_min_ini to allow faster convergence.
        temp_i=length(so_g[,1])-1
        n_min_ini=floor(n_min_ini/2)
      }
    # print("not found")
    }
    
    if(check1==TRUE && check2==TRUE)      {      #if both checks are ok, so their value is TRUE, we found the position and we can put the group from b into so, and recalculate so_g
      print("found")
      so_n=rbind(so[1:max(temp_so_ind_aft),],temp_g,so[min(temp_so_ind_bef):length(so[,1]),])
      so=so_n
      
      #update so_g
      so_gtemp=rbind(so_g[1:(temp_i),],c(temp_i,length(temp_g[,1]),0),so_g[((temp_i+1):length(so_g[,1])),] )
      so_gtemp[,3]=cumsum(so_gtemp[,2])
      so_g=so_gtemp
      found=1
      
    }  
    
    } #end while found==FALSE
} #end for bgbg in 1:length....

}#END if(length(b_g[,1])>0)...

### end appending buffer leftovers ###

### CALCULATE AFTER WHICH MEASUREMENT YOU CAN TAKE OUT THE FIRST CABLE ###
sequ_txt=cbind(c(1:len),so)
first_cable_out=max(sequ_txt[(sequ_txt[,2]<=first_cable_length |sequ_txt[,5]<=first_cable_length |sequ_txt[,3]<=first_cable_length | sequ_txt[,4]<=first_cable_length ),1])

second_cable_out=max(sequ_txt[(sequ_txt[,2]<=first_cable_length*2 |sequ_txt[,5]<=first_cable_length*2 |sequ_txt[,3]<=first_cable_length*2 | sequ_txt[,4]<=first_cable_length*2 ),1])


####UNCOMMENT IF TRIAL TO UNDERSTAND BEST n_min ###

# len_b_trial[trial]=length(b[,1])
# }#end trial for
# plot(c(1:trial),first_cable_out,type="h")

#### END UNCOMMENT IF TRIAL TO UNDERSTAND BEST n_min ###





### SAVE DATA ###

#first, create an output folder
ifelse(!dir.exists(file.path("output")),
       dir.create(file.path("output")),
       "Directory Exists")


#Save the sequence in .txt format which is readable by ELectre Pro software (Iris Instruments), which can be used to put this sequence into the georesistivimeter Syscal Pro

if(file_type=="ELECTRE"){
  ele_pos_ELECTRE=readLines(input_sequence_file,n=howmany_ele+1)
  writeLines(ele_pos_ELECTRE,sprintf("output/o-%s-1stcable-off-at-measure%s_2nd-at-%s.txt",input_sequence_file_base,first_cable_out,second_cable_out))
} else{
  ele_pos=data.frame(num=c(1:max(s[,])),x=c(1:max(s[,])),y=0,z=0)
  colnames(ele_pos)=c("#","X","Y","Z")
  write.table(ele_pos,sprintf("output/o-%s-1stcable-off-at-measure%s_2nd-at-%s.txt",input_sequence_file_base,first_cable_out,second_cable_out),row.names = FALSE,quote=FALSE,sep="\t")
  
}

sequ_txt=cbind(c(1:len),so)
colnames(sequ_txt)=c("#","A","B","M","N")
write.table(sequ_txt,sprintf("output/o-%s-1stcable-off-at-measure%s_2nd-at-%s.txt",input_sequence_file_base,first_cable_out,second_cable_out),row.names = FALSE,quote=FALSE,sep="\t",append=TRUE)

 # Save the sequence in .csv format for direct import by ResIPy software when performing a forward simulation.
write.table(so,sprintf("output/%s%s.csv",input_sequence_file_base,format(Sys.time(), "%d-%h-%y_%H.%M")),sep=",",row.names = FALSE,col.names = c("A","B","M","N"),quote=FALSE)

### END SAVE DATA ###



###
### END SCRIPT. BYE! ###



##old code##

# s=as.data.frame(s)
# for (g in 1:length(s[,1])){
#   s[g,5]=sprintf("%s %s %s %s",s[g,1],s[g,2],s[g,3],s[g,4])
#   
# }