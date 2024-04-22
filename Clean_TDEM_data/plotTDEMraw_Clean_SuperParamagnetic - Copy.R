#uncomment following line and last line to make this script a function
#plotTDEMraw_Clean_SuperParamagnetic=function(){  



# This script reads a raw data file from TEM-FAST instrument, a TDEM (Time-domain electromagnetic method) instrument.
# This raw data file contains some metadata lines, and then 4 columns, containing Channel, Time (in microseconds), Signal (in V/A), and Error (in V/A)
# This script runs in R version 4.3.2 and Rstudio 2023.09.1 Build 494.
# If you want to cite it, please cite the publication of mine (Andrea Vergnano) about TDEM and HVSR measurements in Jangany, Madagascar, 2024.

# What does this script do:
#0) Asks an input file to be selected interactively
#1) Allows the user to graphycally select which data to keep, to clean noisy parts of the signal, such as those due to superparamagnetic effects.
#2) Calculates apparent resistivity (Rhoapp), and transforms it to effective resistivity and depth, according to Meju (1998). A plot of two Meju's equations for Rho_eff is produced. A plot of Obukhov conditions is produced.
#3) Produces an output .csv file containing also the apparent and effective resistivity, and nicely formatted in columns to be imported in a custom MATLAB script to invert the TDEM data in 1D.


### LOADING REQUIRED PACKAGES ###
if (!require("readr")) install.packages("readr")
if (!require("tcltk")) install.packages("tcltk")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("shiny")) install.packages("shiny")
if (!require("clipr")) install.packages("clipr")
if (!require("data.table")) install.packages("data.table")
if (!require("viridis")) install.packages("viridis")
if (!require("paletteer")) install.packages("paletteer")
if (!require("rmatio")) install.packages("rmatio")


library(readr)
library(tcltk)
library(ggplot2)
library(shiny)
library(clipr)
library(data.table)
library(viridis)
library(paletteer)
library(rmatio)

### END LOADING REQUIRED PACKAGES ###



### DEFINITION OF THE GRAPHICAL INTERFACE TO SELECT DATA TO KEEP - BY MEANS OF A SHINYAPP ###


####Shinibrush -rectangular selection of data . Plots also e(t)*t vs t to recognize when it's flat (Superparamagnetic effect)#####
ui <- fluidPage(
  plotOutput("plot", brush = "plot_brush"),
  tableOutput("data")
)
server <- function(input, output, session) {
  output$plot <- renderPlot({ggplot2::ggplot(data=rawdata)+
     # geom_point(mapping=aes(x=Time,y=`E/I[V/A]`),color="lightblue")+  
     # geom_point(mapping=aes(x=Time,y=`E/I[V/A]`*Time),color="darkgreen")+
     # geom_point(mapping=aes(x=Time,y=(`rho_app`)),color="red")+
      
      geom_point(mapping=aes(x=Time,y=`E/I[V/A]`,color="Signal"))+  
      geom_point(mapping=aes(x=Time,y=`E/I[V/A]`*Time,color="Signal * Time"))+
      geom_point(mapping=aes(x=Time,y=(`rho_app`),color="Rho app"))+
      scale_x_log10()+
      scale_y_log10()+
      geom_errorbar(mapping=aes(x=Time,ymin=`E/I[V/A]`- `Err[V/A]`*10,ymax=`E/I[V/A]`+ `Err[V/A]`*10, color="Errorbar multiplied by 10"))+
      #theme_bw()
    #  scale_color_viridis(discrete=TRUE, option="viridis")
      scale_color_paletteer_d("fishualize::Thalassoma_pavo")+
      labs(caption="Draw a rectangle around the Signal data points to keep. \n You can check the pickings in the table that will appear below. \n If you missed some points, just adjust or redraw the rectangle. \n Tips: \n -Avoid points where the signal rises with time \n -Avoid points where Signal*time is constant in time \n -Avoid points where errorbar is too high \n -Check also apparent resistivity (rho_app) data points if they seem ok \n -When finished, close the window and the script will continue")
  }, res = 96)
  
  output$data <- renderTable({
    brushedPoints(rawdata, input$plot_brush)
   # write_clip(brushedPoints(rawdata, input$plot_brush)) to copy to the clipboard (when I did not know the global assignment <<- below)
pickings <<- (brushedPoints(rawdata, input$plot_brush)) #assigns the pickings to the global variable "pickings".# all columns of the input dataframe are saved.
  })
  
  #a little 3 lines to allow the code to run after the shiniapp(ui,server) command.
  session$onSessionEnded(function() {
    stopApp()
  })
}

### END OF SHINYAPP DEFINITION ###


###STEP 1 import the raw file #####

element=file.choose()
 
  rawdata_metadata <- (read_delim(element, delim = "\t", 
                        escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE, skip = 3))[1:4,]
 
   rawdata <- read_delim(element, delim = "\t",
                         escape_double = FALSE, trim_ws = TRUE,
                         skip =7)
   rawdata$Time=rawdata$Time/10^6
   
   
   #read some metadata (loop dimensions, current)
   TxS=as.numeric(rawdata_metadata$X2[2])
   RxA=as.numeric(rawdata_metadata$X2[2])*as.numeric(rawdata_metadata$X4[2])
   if(TxS^2 !=RxA)print("WARNING - Check loop dimensions on raw data file")
  
   #read value of I which is inside other characters in the same cell.
    I=rawdata_metadata$X6[1]
    I <- unlist(regmatches(I, gregexpr('\\(?[0-9,.]+', I)))
    I <- as.numeric(gsub('\\(', '-', gsub(',', '', I)))
   
    #check if electrical current value seems OK.
    
if(is.numeric(I)==FALSE || I>3 || I<0.5) print("WARNING - CHECK current value in raw data file")
   
    
    # Calculates apparent resistivity (rho app), according to Obukhov,1986, cited by Spichak, 2007.
   
  mu=4*pi*10^(-7) #magnetic permeability of free space in Henry/meter
  k2=(pi^0.5)/20 # it should be an instrument constant (this is good for sure for for TEM-FAST, since it comes from its manual written mostly by Spichak)
  TxR=TxS/(pi^0.5)
  signal=rawdata$`E/I[V/A]`
  times=rawdata$Time
  
  rho_app=(k2*mu^(5/2)*TxR^4/(times^(5/2)*(signal/I)))^(2/3) #according to Obukhov,1986, cited by Spichak, 2007.
  
  rawdata$rho_app=rho_app #update the raw data dataframe with values of rho_app
 
  
   
  #check if Obukhov's condition is respected Ob>>1
  
  Ob=times/(mu*TxR^2/rho_app)
  ggplot()+
    geom_point(mapping=aes(x=times,y=Ob,color="Obukhov condition "))+
    geom_point(mapping=aes(x=times,y=rho_app,color="Rho apparent"))+
    
    scale_x_log10()+
    ylim(0,300)+
    labs(caption="Obukhov condition should be >> 1 to trust Rho apparent")
  
  
  
  
  #Perform rhoapp-rhodepth transform according to Meju (equation 5, which is not the best but it is simpler)
  k=2.3
  alpha=0.15
  
  rho_eff_eq5=k*rho_app*exp(-1+alpha)
  depth_eff_eq5=(2*times*rho_app/mu)^(0.5)/k
#   ggplot()+
#     geom_point(mapping=aes(x=rho_eff_eq5,y=depth_eff_eq5,color=(times)))+
# scale_x_log10(limits=c(30,500))+
#         scale_y_reverse(limits=c(300,0))+
#           labs(caption="rho vs depth after equation 5 of Meju 1998")
#   
  #Perform rhoapp-rhodepth transform according to Meju (equation 1). Maybe...put it after selecting the good points since it needs some smoothing.
  
  T=3.9*times
  rho_app_smooth=rho_app#smooth(rho_app)
  depth_eff_eq1=(rho_app_smooth*T/(2*pi*mu))^(0.5)
  d_logrho=c(0,diff(log(rho_app)))
  d_logT=c(0,diff(log(T)))
  rho_eff_eq1=rho_app*(1+d_logrho/d_logT)
  
  # ggplot()+
  #   geom_point(mapping=aes(x=rho_eff_eq1,y=depth_eff_eq1,color=(times)))+
  #   scale_x_log10(limits=c(30,500))+
  #   scale_y_reverse(limits=c(300,0))+
  #   
  #   labs(caption="rho vs depth after equation 1 of Meju 1998")
  # 
  # 
  # #plot comparing the two equations:
  # ggplot()+
  #   geom_point(mapping=aes(x=rho_eff_eq1,y=depth_eff_eq1,color="equation 1"))+
  #   scale_x_continuous(limits=c(-300,5000))+
  #   scale_y_reverse(limits=c(300,0))+
  #   geom_point(mapping=aes(x=rho_eff_eq5,y=depth_eff_eq5,color="equation 5"))+
  #   
  #   labs(caption="rho vs depth after equation 1 and 5 of Meju 1998")
  # 
  #averaging equation 1 and 5 of Meju, 1998.
  
  rho_eff_eq1_5=rho_eff_eq5*rho_eff_eq1/((rho_eff_eq5+rho_eff_eq1)/2)
  depth_eff_eq1_5=depth_eff_eq1/2+depth_eff_eq5/2
  
  ggplot()+
    geom_point(mapping=aes(x=rho_eff_eq1_5,y=depth_eff_eq1_5,color="equation 1 and 5 averaged"))+
    scale_x_log10(limits=c(30,500))+
    scale_y_reverse(limits=c(300,0))+
    geom_point(mapping=aes(x=rho_eff_eq1,y=depth_eff_eq1,color="equation 1"))+
    geom_point(mapping=aes(x=rho_eff_eq5,y=depth_eff_eq5,color="equation 5"))+
    
    labs(caption="rho vs depth after equation 1 and 5 of Meju 1998")
  # 
  #this part of the script was meant to allow the shiniapp to work into a function, it had problems of timing
  #i=0
  # 
  # while (i <4) {
  #   i=i+1
  #   print(shinyApp(ui, server))
  #   Sys.sleep(5)
  #   a=1
  #   }
  # shinyApp(ui, server)
  # shinyApp(ui, server)
  # shinyApp(ui, server)
  
  
  
###STEP 2 Plot and select the square of data to keep, then close the windows and go to STEP 3 ##
  
print(shinyApp(ui, server))
   



   
###STEP 3 export the data to a csv. 
#THIS STEP MIGHT HAVE TO BE RUN SEPARATELY, IF THE CODE STOPS AFTER SHINYAPP. Check if the output is there, it should be already in the same folder of the raw data!
   
#rawdata_selected=read.table(text = read_clip(), header = TRUE, sep = "\t")
rawdata_selected=pickings
rawdata_selected$RxA=RxA
rawdata_selected$TxS=TxS
rawdata_selected$TxR=TxR
rawdata_selected$I=I
rawdata_selected$depth_eff_eq1=depth_eff_eq1[rawdata_selected$Channel] #it saves only those depth_eff and rho_eff correspondent to the data selected by manual picking
rawdata_selected$rho_eff_eq1=rho_eff_eq1[rawdata_selected$Channel]
rawdata_selected$depth_eff_eq5=depth_eff_eq5[rawdata_selected$Channel]
rawdata_selected$rho_eff_eq5=rho_eff_eq5[rawdata_selected$Channel]

print("You can find the output in the same folder of the raw data")

####____save a simple CSV with the rawdata_selected matrix of the tdem profile for further interpretation, geolocalization, etc.
fwrite(rawdata_selected,sprintf("%s_.csv",element),row.names=FALSE,scipen=50)

#####_____save a .csv file for TDEM1D interpretation. it is transposed since the 1D matlab code wants an input matrix aligned horizontally

#final_data_for_TDEM1D_input=t(rawdata_selected)

#fwrite(final_data_for_TDEM1D_input,sprintf("%s_.csv",element),col.names=FALSE,scipen=50) #scipen50 è necessario nelle opzioni di fwrite per avere dei numeri corretti in notazione non scientifica nell'export.




### Try to directly save a .mat file
####____save a .mat file. it loses some data when reimporting to matlab... not perfectly compatible

# mat_example$data$cparams$times=list(rawdata_selected$Time) #this list command is very important to make the save work.
# mat_example$data$tdem$signal=list(rawdata_selected$`E/I[V/A]`)
# mat_example$data$tdem$error=list(rawdata_selected$`Err[V/A]`)
# mat_example$data$config$TxS=TxS
# mat_example$data$config$RxA=RxA
# mat_example$data$config$TxR=TxR
# mat_example$data$config$I=I
# 
# write.mat(mat_example,sprintf("%s.mat",element))




#}


#trial and errors

# Sys.sleep(3)
#rm(list = ls())
# final_string=sprintf("%.08f",final_data_for_TDEM1D_input)
# write.csv(final_data_for_TDEM1D_input,sprintf("%s.csv",element),col.names=FALSE)
# write.table((final_data_for_TDEM1D_input), sep="\t", file=sprintf("%s.txt",element), row.names=FALSE)
#write.table(final_data_for_TDEM1D_input,sprintf("%s.csv",element),col.names=FALSE,scipen=50)
#row.names=c("channel","Time(micros)","E/I(V/A)","Err(V/A)")




#uncomment next line and first line to make this script a function
#} 
