####Script to run LFMM####
##This script assumes you are in the working directory where you have the files to run LFMM

#Install the packages "lfmm" and "LEA"
library(devtools)
devtools::install_github("bcm-uga/lfmm")
devtools::install_github("bcm-uga/LEA")

library(lfmm)
library(LEA)
library(raster)
library(rgdal)

#Convert our file with ext .ped to format lfmm 
#The file with ext .ped must not have headers
output = ped2lfmm("batch_1.plink_wp.ped", force = TRUE)

##Climate data##

#Choose the directory where we have the climate data layers (only those without correlated values) with ext .asc
setwd(choose.dir())###directory of the climate data layers
clim.list <- list.files(".",pattern = "*.asc$",full.names = T) 
clim.layer<- stack(clim.list) #stacks the layers into a  object

#Plotting
pdf("clim.layer.pdf")
plot(clim.layer)
dev.off()

#Back to our working directory
setwd(choose.dir())
sample.coord <-read.table("datos_coord.txt", header=T, stringsAsFactors=F)#open coordinates file
sample.coord

#Define the spatial projection system that the points are in (usually WGS84)
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  
sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('decimalLongitude','decimalLatitude')], proj4string=CRS(crs.wgs), data=sample.coord)
#Extract the climate data for each point
clim.points <- extract(clim.layer, sample.coord.sp)


#Combine the sample coordinates with the climate data points
clim.points <- cbind(sample.coord, clim.points)  
write.table(clim.points, "clim.points", sep="\t", quote=F, row.names=F)  
clim.points

### Test the colinearity with the varianze inflation factor, if VIF>10 implies multicolinearity, we must remove those variables
vif(clim.points[,4:22])
no_corr <- vifstep(clim.points[,4:22], th=10) # th=treshold
no_corr
##We stay only with the bios 17, 18, 3, 9#

#Save climate data without sample names, latitude, longitude, or column names for LFMM
clim.env <- clim.points[, c(12,13,16,22)]
colnames(clim.env) <- NULL
clim.env
write.table(clim.env, "clim.env", sep="\t", quote=F, row.names=F)


##Assesing population structure##
#Run snmf to estimate K, considering K from 1-11 which is the number of our populations:
project = NULL
project = snmf("batch_1.plink_wp.lfmm", K = 1:11, entropy = TRUE, repetitions = 10, project = "new")
pdf("sNMF1.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()
#The "best" value of K has the lowest cross-entropy value (y-axis), in this case k=3

#Plotting the ancestry proportions#
best = which.min(cross.entropy(project, K = 3))
pdf("sNMF.barchart.pdf")
barchart(project, K = 3, run = best, border = NA, space = 0, col = c("red", "blue", "green"), xlab = "Individuals", ylab = "Ancestry proportions") -> bp
axis(1, at = 1:length(bp$order), labels = bp$order, las=1, cex.axis = .3)
dev.off()



####Run LFMM#####
#WE RUN THE LFMM with k=the "best" number of clusters that we inferred above, the number of repetitions of the model to run, the number of processors to use, and the number of iterations and burn-in.
project = NULL
project = lfmm("batch_1.plink_wp.lfmm", "clim.env", K = 3, repetitions = 5, CPU = 16, iterations = 1000, burnin = 500, project = "new")

##Combine the data from the five repetitions and compute new calibrated P-values.
##Bioclim 17: Precipitation of driest quarter
z.pdry = z.scores(project, K = 3, d = 1)
z.pdry <- apply(z.pdry, 1, median)

##We need to calculate ?? (the "genomic inflation factor")
lambda.pdry = median(z.pdry^2)/qchisq(0.5, df = 1)

##Calculate the "adjusted" P-values
p.pdry.adj = pchisq(z.pdry^2/lambda.pdry, df = 1, lower = FALSE)

##Repeat with the others climates variables
##Bioclim 18: Precipitation of warmest quarter
z.pwarm = z.scores(project, K = 3, d = 2)
z.pwarm <- apply(z.pwarm, 1, median)
lambda.pwarm = median(z.pwarm^2)/qchisq(0.5, df = 1)
p.pwarm.adj = pchisq(z.pwarm^2/lambda.pwarm, df = 1, lower = FALSE)

##Bioclim 3: Isothermality
z.iso = z.scores(project, K = 3, d = 3)
z.iso <- apply(z.iso, 1, median)
lambda.iso = median(z.iso^2)/qchisq(0.5, df = 1)
p.iso.adj = pchisq(z.iso^2/lambda.iso, df = 1, lower = FALSE)

##Bioclim 9: Mean temperature of driest quearter
z.mtdry = z.scores(project, K = 3, d = 4)
z.mtdry <- apply(z.mtdry, 1, median)
lambda.mtdry = median(z.mtdry^2)/qchisq(0.5, df = 1)
p.mtdry.adj = pchisq(z.mtdry^2/lambda.mtdry, df = 1, lower = FALSE)

#Make histograms of the P-values
pdf("LFMM_P_Histograms.pdf")
par(mfrow = c(4,1))
hist(p.pdry.adj, col = "blue", main = "Pdry", xlab='')
hist(p.pwarm.adj, col = "blue", main = "Pwarm", xlab='')
hist(p.iso.adj, col = "blue", main = "Iso", xlab='')
hist(p.mtdry.adj, col = "blue", main = "Mtdry", xlab=expression(italic(P)))
dev.off()
##The "best" K and proper calibration value will lead to histograms of P-values that are flat, 
#except perhaps with an elevated frequency of very low P-values, representing the outliers

#If the model is behaving well, we can continue

###Adjustment of q values###
install_github("jdstorey/qvalue")
library(qvalue)

#We need to correct for multiple testing, for this we adjust the P-values to Q-values
q.pdry<-qvalue(p.pdry.adj)$qvalues
q.pwarm<-qvalue(p.pwarm.adj)$qvalues
q.iso<-qvalue(p.iso.adj)$qvalues
q.mtdry<-qvalue(p.mtdry.adj)$qvalues


#Manhattan plots to visually summarize large numbers of association tests#
pdf("LFMM_Manhattan.pdf")
par(mfrow = c(4,1))
plot(-log10(q.pdry), pch = 19, col = "blue", cex = .7, xlab = '', ylim=c(0,12))
plot(-log10(q.warm), pch = 19, col = "blue", cex = .7, xlab = '',ylim=c(0,12))
plot(-log10(q.iso), pch = 19, col = "blue", cex = .7, xlab = '',ylim=c(0,12))
plot(-log10(q.mtdry), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)", ylim=c(0,12))
dev.off()
## -log10(0.01) = 2, so values greater than 2 representing low Q-values less than 0.01.

##To know how many of the significant ones (Q < 0.05) are also of large effect we can look 
##at the set of significant SNPs that also have very high or very low z-scores.
sum(q.pdry<0.01 & abs(z.pdry)>2)
sum(q.pwarm<0.01 & abs(z.pwarm)>2)
sum(q.iso<0.01 & abs(z.iso)>2)
sum(q.mtdry<0.01 & abs(z.mtdry)>2)

#We can combine all the z and Q-values into a table and save for use in other software.
lfmm.results <- cbind(z.pdry, q.pdry, z.pwarm, q.pwarm, z.iso, q.iso, z.mtdry, q.mtdry)
head(lfmm.results)
write.table(lfmm.results, "lfmm.results", sep="\t", quote=F, row.names=F)
