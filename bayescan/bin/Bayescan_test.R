#This script is used to plot figures from the outputs files of Bayescan software, 
#using the function plot_R Copyright (C) 2010  Matthieu Foll

#This script assume assume you are in a directory that contains:
#A folder where you have the file "plot.r".
#A folder where you have the output of running Bayescan with ex: "output_fst.txt".


#Load packages
library(ggplot2)

#Use the function already available in Bayescan. First load the function in R.
source("plot_R.r")

##Run the function with the file ex: "output_fst.txt"
pdf("results1.pdf")
results1<-plot_bayescan("../Bayescan_data/data_wp_fst.txt", FDR=0.001) ##significance threshold
dev.off()

#Check results
summary(results1)

#Check the outlier list (locus id)
head(results1$outliers)

##Another ways to visualize the results:
#Load your file ex: "output_fst.txt"
juniperus<-read.csv("../Bayescan_data/data_wp_fst.txt", 
                       sep="", header= TRUE,
                       col.names = c("prob", "log10(PO)", "qval", "alpha", "fst"))

#Create the column "bayescan number"
Bayescan_number<-c(1:2904)

#add column to our data frame
juniperus1<-cbind(juniperus, Bayescan_number)

##Explore with differents types of plots
#Option1
results2<-ggplot(juniperus1, aes(x=Bayescan_number, y=fst))+
  geom_point(size=2, alpha=0.5,
             aes(colour = qval<= 0.001 )) + # color by significance threshold
  scale_color_manual(values = c("black", "red")) #choose the colors
ggsave("results2.pdf")

#Option2
results3<-ggplot(juniperus1, aes(x=alpha, y=fst)) +
  geom_point(size=2, alpha=0.5,
             aes(colour = qval<= 0.001 )) + # color by significance threshold
  scale_color_manual(values = c("black", "red")) #choose the colors
ggsave("results3.pdf")
