#This script is used to plot figures for the software Bayescan in R.
#Copyright (C) 2010  Matthieu Foll

#This script assume assume you are in a directory that contains:
#A folder where you have the file "plot.r".
#A folder where you have the output of running Bayescan with ex: "output_fst.txt".


#Load packages
library(ggplot2)

#Use the function already available in Bayescan. First load the function in R.
source("plot_R.r")

##Run the function with the file ex: "output_fst.txt"
results1<-plot_bayescan("../results/data_wp_fst.txt", FDR=0.001) ##significance threshold

#Check results
summary(results1)


##Other way to visualize the results:
#Load your file ex: "output_fst.txt"
juniperus<-read.csv("../results/data_wp_fst.txt", 
                       sep="", header= TRUE,
                       col.names = c("prob", "log10(PO)", "qval", "alpha", "fst"))

#Create the column "bayescan number"
bayescan_number<-c(1:2904)

#add column to our data frame
juniperus1<-cbind(juniperus, bayescan_number)

#Explore with a different type of plot
ggplot(juniperus1, aes(x=bayescan_number, y=fst)) +
  geom_point(size=2, alpha=0.5,
             aes(colour = qval<= 0.001 )) + # color by significance threshold
  scale_color_manual(values = c("black", "red")) #choose the colors

