#########################################################################
#########################################################################
####                                                                 ####
####                                                                 ####
####                                                                 ####
####                           LAND-SE                               ####
####             LANDSLIDE SUSCEPTIBILITY EVALUATION                 ####
####                           IRPI CNR                              ####
####                    MAURO ROSSI - IRPI CNR                       ####
####                   v1r0b30 - 18 January 2016                     ####
####                                                                 ####
#### Copyright (C) 2016 Mauro Rossi                                  ####
####                                                                 ####
#### This program is free software; you can redistribute it and/or   ####
#### modify it under the terms of the GNU General Public License     ####
#### as published by the Free Software Foundation; either version 2  ####
#### of the License, or (at your option) any later version. ###      ####
####                                                                 ####
#### This program is distributed in the hope that it will be useful, ####
#### but WITHOUT ANY WARRANTY; without even the implied warranty of  ####
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the ####
#### GNU General Public License for more details.                    ####
####                                                                 ####
####     Istituto di Ricerca per la Protezione Idrogeologica         ####
####              Consiglio Nazionale delle Ricerche                 ####
####                    Gruppo di Geomorfologia                      ####
####                  Via della Madonna Alta, 126                    ####
####                    06128 Perugia (Italia)                       ####
####                       +39 075 5014421                           ####
####                       +39 075 5014420                           ####
####                   mauro.rossi@irpi.cnr.it                       ####
####                  geomorfologia@irpi.cnr.it                      ####
####                                                                 ####
####           This script was prepared using R 3.0.1                ####
####         The script requires the following R packages:           ####
####                       1: MASS                                   ####
####                       2: vcd                                    ####
####                       3: verification                           ####
####                       4: perturb                                ####
####                       5: glm                                    ####
####                       6: nnet                                   ####
####                       7: raster                                 ####
####                       8: rasterVis                              ####
####                       9: RColorBrewer                           ####
####                                                                 ####
####     INPUTS: 1) data.txt file (tab delimited)                    ####
####                1st column -> identification value               ####
####                2nd column -> grouping variable                  ####
####                Other columns -> explanatory variables           ####
####                                                                 ####
####             2) validation.txt file (tab delimited)              ####
####                1st column -> identification value               ####
####                2nd column -> validation grouping variable       ####
####                Other columns -> explanatory variables           ####
####                                                                 ####
####             3) training.shp (polygons or points)                ####
####                column -> identification value                   ####
####                column -> mapping unit area                      ####
####                          (only for polygons)                    ####
####                column -> landslide area in the mapping unit     ####
####                          (only for polygons)                    ####
####                                                                 ####  
####             4) validation.shp (polygons or points)              ####
####                column -> identification value                   ####
####                column -> mapping unit area                      ####
####                          (only for polygons)                    ####
####                column -> landslide area in the mapping unit     ####
####                          (only for polygons)                    ####
####                                                                 ####                                                             
####     3) and 4) are not mandatory. They are useful to produce     ####
####     output maps (shp and pdf) and success rate and              ####
####     prediction rate curve.                                      ####
####     The configuration_spatial_data.txt file contains their      ####
####     configuration parameter.                                    ####                                                           
####                                                                 ####                                                          
####                4th column -> analysis parameter                 ####
####                              For QDA can be:                    ####
####                              SEL (eliminate dummy variables)    ####
####                              DUM (maintain dummy variables      ####
####                              trasformed in numeric              ####
####                              introducing a random variation     ####
####                              between -0.1 and +0.1 and 0 are    ####
####                              avoided sampling quantities        ####
####                              close to 0)                        ####
####                              For NNM can be:                    ####
####                              NOR (default weigth, half          ####
####                              variables in the hidden layer)     ####
####                              OPT (auto optimize the neural      ####
####                              structure:  slower)                ####
####                                                                 ####
####                                                                 ####
#########################################################################
#########################################################################


#C:\PROGRA~1\R\R-3.0.3\bin\R.exe --no-save --args -cd /media/sf_disco_dati/R/SusceptibilityAnalysis/testing -wd /media/sf_disco_dati/R/SusceptibilityAnalysis/testing < rainfall_events_commented.R > susceptibilty.log
#pars<-c("-cd","/media/sf_disco_dati/R/SusceptibilityAnalysis/testing","-wd","/media/sf_disco_dati/R/SusceptibilityAnalysis/testing")
pars <-commandArgs(trailingOnly=TRUE)

# rm(list=(ls()))
# graphics.off()
#memory.limit(size=120000)
time_start_calculation<-Sys.time()

if (length(table(pars == "-wd"))==2)
  {
  wd_selected<-pars[which(pars=="-wd")+1]
  } else
  {
  #wd_selected<-"/polygon_Random_validation_entirearea/" # linux
  wd_selected<-"C:/polygon_Random_validation_entirearea/" # windows
  }
setwd(wd_selected)

if (length(table(pars == "-cd"))==2)
  {
  cd_selected<-pars[which(pars=="-cd")+1]
  } else
  {
  #cd_selected<-"/polygon_Random_validation_entirearea/" # linux
  cd_selected<-"C:/polygon_Random_validation_entirearea/" # windows
  }



#--------------------------- PARAMETER DEFINITION ---------------------------#
load_rdata<-FALSE
name_inventory<-"inventory"
rdata_file<-paste("datatable_inventory.RData",sep="")

enable_probability_optimal_binary_classification<-FALSE
enable_probability_optimal_classification<-FALSE # This require enable_probability_optimal_binary_classification<-TRUE
type_probability_optimal_classification<-"proportional" # or "fixed": "proportional" will partition the unexplained errors (1-TPR and 1-TNR) proportionally but with the some original number of classed, "fixed" will used a fixed subdivision of the TPR and FPR range, and depending from the optima binary value should results in different number of classes


enable_screen_plotting<-FALSE
enable_multicollinearity_test<-FALSE
enable_rocplot_confidenceinterval<-FALSE
threshold_series<-seq(0,1,0.01)
#threshold_series<-NULL #Default procedure
bootstrap_constant_correction<-TRUE #To enable when explanatory variables contain mostly 0 value (blocking susceptibility model convergence)
enable_detailed_data_export<-TRUE

#--------------------------- SEED SETTING ---------------------------#
seed.value<-NULL
#seed.value<-1 # Uncomment this line if you want to fix a seed. If this is the case multiple run of the script will give always the same result.
if (is.numeric(seed.value)){seed.value<-seed.value} else {seed.value<-round(runif(1,min=1,max=10000))}


#---------------------  READ CONFIGURATION FILE ---------------------#
configuration.table<-read.table(paste(cd_selected,"LAND-SE_configuration.txt",sep=""),header = TRUE,dec=".", sep="\t")
model.type<-as.character(configuration.table$MODEL)
model.run.matrix<-as.character(configuration.table$RUN)
bootstrap.sample.values<-as.numeric(configuration.table$BOOTSTRAP_SAMPLES_ROC_CURVE)
analysis.parameter.matrix<-as.character(configuration.table$ANALYSIS_PARAMETER)
bootstrap.model.variability<-as.character(configuration.table$BOOTSTRAP_MODEL_VARIABILITY_RUN)
bootstrap.sample.model<-as.numeric(configuration.table$BOOTSTRAP_SAMPLES_MODEL_VARIABILITY)

#--------------- READ CONFIGURATION SPATIAL DATA FILE ---------------#
configuration.spatial.data.table<-read.table(paste(cd_selected,"LAND-SE_configuration_spatial_data.txt",sep=""),header = TRUE,dec=".", sep="\t",stringsAsFactors=FALSE)
spatial.data.type<-as.character(configuration.spatial.data.table$TYPE)
spatial.data.presence<-as.character(configuration.spatial.data.table$PRESENCE)
spatial.data.id<-as.character(configuration.spatial.data.table$ID_FIELD)
spatial.data.epsg<-as.numeric(configuration.spatial.data.table$EPSG)

#------------------------- READ OF THE DATA -------------------------#
if (load_rdata==FALSE)
  {
  training.table<-read.table("training.txt",header = TRUE,dec=".", sep="\t")
  } else
  {
  load(file=rdata_file)
  }
      
#names(training.table)
#dim(training.table)
training.table<-na.omit(training.table)
#dim(training.table)

explanatory.variables<-training.table[3:dim(training.table)[2]]
#str(explanatory.variables)
#range(explanatory.variables)
#dim(explanatory.variables)

data.variables<-training.table[,2:dim(training.table)[2]]
#data.variables<-training.table[,2:58]
#dim(data.variables)

grouping.variable<-as.factor(training.table[,2])
#str(grouping.variable)
identification.value<-training.table[,1]
#str(identification.value)

## Histogram of posterior grouping variable
breaks.histogram.values<-c(0,0.2,0.45,0.55,0.8,1)

color_ramp_fun<-colorRamp(c(rgb(38,115,0,max=255),rgb(255,255,0,max=255),rgb(255,0,0,max=255)))
#color_ramp_fun(breaks.histogram.values)

color_ramp_palette_fun<-colorRampPalette(c(rgb(38,115,0,max=255),rgb(255,255,0,max=255),rgb(255,0,0,max=255)))
#color_ramp_palette_fun(length(breaks.histogram.values)-1)



if (enable_screen_plotting==TRUE)
{
dev.new()
hist(training.table[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of grouping variable", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
}
pdf(file = "GroupingVariable_Histogram.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
hist(training.table[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of grouping variable", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
hist(training.table[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of grouping variable", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))

dev.off()

#-------------------  READ OF THE VALIDATION DATA -------------------#
if (load_rdata==FALSE)
  {
  validation.table<-read.table("validation.txt",header = TRUE,dec=".", sep="\t")
  }
      
#dim(validation.table)

validation.table<-na.omit(validation.table)
#dim(validation.table)

validation.explanatory.variables<-validation.table[3:dim(validation.table)[2]]
#str(validation.explanatory.variables)
#range(validation.explanatory.variables)
#dim(validation.explanatory.variables)

validation.variables<-validation.table[,2:dim(validation.table)[2]]
#dim(validation.variables)

validation.grouping.variable<-as.factor(validation.table[,2])
#str(validation.grouping.variable)
validation.identification.value<-validation.table[,1]
#str(validation.identification.value)

## Histogram of validation grouping variable

if (enable_screen_plotting==TRUE)
{
dev.new()
hist(validation.table[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of validation grouping variable", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
}
pdf(file = "GroupingVariable_Histogram_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
hist(validation.table[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of validation grouping variable", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
dev.off()


############ Selezioe soooinsieme dati
#############
#indexselection<-sample.int(n=length(grouping.variable), size = 10000)
#identification.value<-identification.value[indexselection]
#training.table<-training.table[indexselection,]
#grouping.variable<-grouping.variable[indexselection]
#explanatory.variables<-explanatory.variables[indexselection,]
#print(table(colSums(explanatory.variables)==0))
#colSums(explanatory.variables)
#names(explanatory.variables)
################
################




#---------------------- READ SPATIAL DATA FILE ----------------------#
if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POLYGONS")
  {
  library(rgdal)
  shape_training<-readOGR(dsn="training.shp",layer="training")
  index_col_shape_training<-match(configuration.spatial.data.table[c(3,5,6)],colnames(shape_training@data))
  index_col_shape_training<-index_col_shape_training[which(is.finite(index_col_shape_training))]
  shape_training@data<-shape_training@data[,index_col_shape_training]
  #shape_training@data<-shape_training@data[order(shape_training@data[,index_col_shape_training[1]]),index_col_shape_training]

  
  
  shape_validation<-readOGR(dsn="validation.shp",layer="validation")
  index_col_shape_validation<-match(configuration.spatial.data.table[c(3,5,6)],names(shape_validation@data))
  index_col_shape_validation<-index_col_shape_validation[which(is.finite(index_col_shape_validation))]
  shape_validation@data<-shape_validation@data[,index_col_shape_validation]
  #shape_validation@data<-shape_validation@data[order(shape_validation@data[,index_col_shape_training[1]]),index_col_shape_validation]

  
  shape_merge_field<-names(shape_training@data)[1]
  
  #spplot(shape_validation)
    
  #shape_validation@data
  #slot(shape_validation, "data") # equivalent to the previous
  #slot(shape_validation, "polygons")
  
  # Matching DATA ID and Shapefile ID and adding variable columuns
  #if(identical(as.character(identification.value),as.character(shape_training@data[,shape_merge_field])) & identical(as.character(validation.identification.value),as.character(shape_validation@data[,shape_merge_field])))
  if(identical(sort(as.character(identification.value)),sort(as.character(shape_training@data[,shape_merge_field]))) & identical(sort(as.character(validation.identification.value)),sort(as.character(shape_validation@data[,shape_merge_field]))))
    {
	shape_training@data <- merge(x=shape_training@data,y=training.table,by.x=shape_merge_field, by.y=names(training.table)[1], all.x=T, sort=F)	
	shape_validation@data <- merge(x=shape_validation@data,y=validation.table,by.x=shape_merge_field, by.y=names(validation.table)[1], all.x=T, sort=F)	
	#shape_training@data<-cbind(shape_training@data, training.table[,-c(1,2)])
    #shape_validation@data<-cbind(shape_validation@data, validation.table[,-c(1,2)])
	} else
    {
    print("ERROR - Data ID and Shapefile ID doesn't match")
    }
  # In case of mismatching  due to order of IDs, do not change the ID in shape data but please refere to http://mapserver.org/utilities/sortshp.html
  
  
  #is.projected(shapefile.subdivision)
  #proj4string(shapefile.subdivision)
  #str(shapefile.subdivision)
  #summary(shapefile.subdivision)
  #coordinates(shapefile.subdivision)
  
  #### spplot parameter training
  # Determination of position of map items
  longitude.map.extent.training<-bbox(shape_training)[1,2]-bbox(shape_training)[1,1]
  latitude.map.extent.training<-bbox(shape_training)[2,2]-bbox(shape_training)[2,1]
  arrow.scale.value.training<-round(longitude.map.extent.training,-log10(longitude.map.extent.training))/10
  scale.scale.value.training<-round(longitude.map.extent.training,-log10(longitude.map.extent.training))/10*2
  arrow.lon.training<-bbox(shape_training)[1,2]-arrow.scale.value.training
  arrow.lat.training<-bbox(shape_training)[2,2]-arrow.scale.value.training*1.5
  scale.lon.training<-bbox(shape_training)[1,2]-scale.scale.value.training*1.5
  scale.lat.training<-bbox(shape_training)[2,1]+arrow.scale.value.training/2
  text.scale1.lon.training<-bbox(shape_training)[1,2]-scale.scale.value.training*1.5
  text.scale1.lat.training<-bbox(shape_training)[2,1]+arrow.scale.value.training
  text.scale2.lon.training<-text.scale1.lon.training+scale.scale.value.training
  text.scale2.lat.training<-text.scale1.lat.training
  
  # Definition of map items
  arrow.training <- list("SpatialPolygonsRescale", layout.north.arrow(),  offset = c(arrow.lon.training,arrow.lat.training), scale = arrow.scale.value.training, which = 1)
  scale.training <- list("SpatialPolygonsRescale", layout.scale.bar(),  offset = c(scale.lon.training,scale.lat.training), scale = scale.scale.value.training, fill=c("transparent","black"), which = 1)
  text.scale1.training <- list("sp.text", c(text.scale1.lon.training,text.scale1.lat.training), "0", which = 1, cex=0.7)
  text.scale2.training <- list("sp.text", c(text.scale2.lon.training,text.scale2.lat.training), paste(scale.scale.value.training," m",sep=""), which = 1, cex=0.7)
  
  #### spplot parameter validation
  longitude.map.extent.validation<-bbox(shape_validation)[1,2]-bbox(shape_validation)[1,1]
  latitude.map.extent.validation<-bbox(shape_validation)[2,2]-bbox(shape_validation)[2,1]
  arrow.scale.value.validation<-round(longitude.map.extent.validation,-log10(longitude.map.extent.validation))/10
  scale.scale.value.validation<-round(longitude.map.extent.validation,-log10(longitude.map.extent.validation))/10*2
  arrow.lon.validation<-bbox(shape_validation)[1,2]-arrow.scale.value.validation
  arrow.lat.validation<-bbox(shape_validation)[2,2]-arrow.scale.value.validation*1.5
  scale.lon.validation<-bbox(shape_validation)[1,2]-scale.scale.value.validation*1.5
  scale.lat.validation<-bbox(shape_validation)[2,1]+arrow.scale.value.validation/2
  text.scale1.lon.validation<-bbox(shape_validation)[1,2]-scale.scale.value.validation*1.5
  text.scale1.lat.validation<-bbox(shape_validation)[2,1]+arrow.scale.value.validation
  text.scale2.lon.validation<-text.scale1.lon.validation+scale.scale.value.validation
  text.scale2.lat.validation<-text.scale1.lat.validation
  
  # Definition of map items
  arrow.validation <- list("SpatialPolygonsRescale", layout.north.arrow(),  offset = c(arrow.lon.validation,arrow.lat.validation), scale = arrow.scale.value.validation, which = 1)
  scale.validation <- list("SpatialPolygonsRescale", layout.scale.bar(),  offset = c(scale.lon.validation,scale.lat.validation), scale = scale.scale.value.validation, fill=c("transparent","black"), which = 1)
  text.scale1.validation <- list("sp.text", c(text.scale1.lon.validation,text.scale1.lat.validation), "0", which = 1, cex=0.7)
  text.scale2.validation <- list("sp.text", c(text.scale2.lon.validation,text.scale2.lat.validation), paste(scale.scale.value.validation," m",sep=""), which = 1, cex=0.7)
  
  # Plot and export of maps
  breaks.map.susceptibility<-c(0,0.2,0.45,0.55,0.8,1.0001)
  breaks.map.uncertainty<-c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)
  breaks.map.matching.code<-c(1,2,3,4,5)
  color.vector.susceptibility<-color_ramp_palette_fun(length(breaks.histogram.values)-1)
  require(RColorBrewer)
  color.vector.uncertainty<-rev(brewer.pal(length(breaks.map.uncertainty)-1,"YlOrRd"))
  color.vector.matching<-c(rgb(38,115,0,max=255),rgb(233,255,190,max=255),rgb(215,215,158,max=255),rgb(115,115,0,max=255))
  }

if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POINTS")
	{
	library(rgdal)
	shape_training<-readOGR(dsn="training.shp",layer="training")
	shape_validation<-readOGR(dsn="validation.shp",layer="validation")
	shape_merge_field<-(configuration.spatial.data.table$ID_FIELD)
	#spplot(shape_validation)
	 
	#shape_validation@data
	#slot(shape_validation, "data") # equivalent to the previous
	#slot(shape_validation, "polygons")
		  
	## Matching DATA ID and Shapefile ID and adding variable columuns
	#if(identical(as.character(identification.value),as.character(shape_training@data[,shape_merge_field])) & identical(as.character(validation.identification.value),as.character(shape_validation@data[,shape_merge_field])))
	#{
	#shape_training@data <- merge(x=shape_training@data,y=training.table,by.x=shape_merge_field, by.y=names(training.table)[1], all.x=T, sort=F)	
	#shape_validation@data <- merge(x=shape_validation@data,y=validation.table,by.x=shape_merge_field, by.y=names(validation.table)[1], all.x=T, sort=F)	
	##shape_training@data<-cbind(shape_training@data, training.table[,-c(1,2)])
	##shape_validation@data<-cbind(shape_validation@data, validation.table[,-c(1,2)])
	#} else
	#{
	#print("ERROR - Data ID and Shapefile ID doesn't match")
	#}
		  
		  
	#is.projected(shapefile.subdivision)
	#proj4string(shapefile.subdivision)
	#str(shapefile.subdivision)
	#summary(shapefile.subdivision)
	#coordinates(shapefile.subdivision)
	 
	#### spplot parameter training
	# Determination of position of map items
	longitude.map.extent.training<-bbox(shape_training)[1,2]-bbox(shape_training)[1,1]
	latitude.map.extent.training<-bbox(shape_training)[2,2]-bbox(shape_training)[2,1]
	arrow.scale.value.training<-round(longitude.map.extent.training,-log10(longitude.map.extent.training))/10
	scale.scale.value.training<-round(longitude.map.extent.training,-log10(longitude.map.extent.training))/10*2
	arrow.lon.training<-bbox(shape_training)[1,2]-arrow.scale.value.training
	arrow.lat.training<-bbox(shape_training)[2,2]-arrow.scale.value.training*1.5
	scale.lon.training<-bbox(shape_training)[1,2]-scale.scale.value.training*1.5
	scale.lat.training<-bbox(shape_training)[2,1]+arrow.scale.value.training/2
	text.scale1.lon.training<-bbox(shape_training)[1,2]-scale.scale.value.training*1.5
	text.scale1.lat.training<-bbox(shape_training)[2,1]+arrow.scale.value.training
	text.scale2.lon.training<-text.scale1.lon.training+scale.scale.value.training
	text.scale2.lat.training<-text.scale1.lat.training
	
	# Definition of map items
	arrow.training <- list("SpatialPolygonsRescale", layout.north.arrow(),  offset = c(arrow.lon.training,arrow.lat.training), scale = arrow.scale.value.training, which = 1)
	scale.training <- list("SpatialPolygonsRescale", layout.scale.bar(),  offset = c(scale.lon.training,scale.lat.training), scale = scale.scale.value.training, fill=c("transparent","black"), which = 1)
	text.scale1.training <- list("sp.text", c(text.scale1.lon.training,text.scale1.lat.training), "0", which = 1, cex=0.7)
	text.scale2.training <- list("sp.text", c(text.scale2.lon.training,text.scale2.lat.training), paste(scale.scale.value.training," m",sep=""), which = 1, cex=0.7)
	
	#### spplot parameter validation
	longitude.map.extent.validation<-bbox(shape_validation)[1,2]-bbox(shape_validation)[1,1]
	latitude.map.extent.validation<-bbox(shape_validation)[2,2]-bbox(shape_validation)[2,1]
	arrow.scale.value.validation<-round(longitude.map.extent.validation,-log10(longitude.map.extent.validation))/10
	scale.scale.value.validation<-round(longitude.map.extent.validation,-log10(longitude.map.extent.validation))/10*2
	arrow.lon.validation<-bbox(shape_validation)[1,2]-arrow.scale.value.validation
	arrow.lat.validation<-bbox(shape_validation)[2,2]-arrow.scale.value.validation*1.5
	scale.lon.validation<-bbox(shape_validation)[1,2]-scale.scale.value.validation*1.5
	scale.lat.validation<-bbox(shape_validation)[2,1]+arrow.scale.value.validation/2
	text.scale1.lon.validation<-bbox(shape_validation)[1,2]-scale.scale.value.validation*1.5
	text.scale1.lat.validation<-bbox(shape_validation)[2,1]+arrow.scale.value.validation
	text.scale2.lon.validation<-text.scale1.lon.validation+scale.scale.value.validation
	text.scale2.lat.validation<-text.scale1.lat.validation
	
	# Definition of map items
	arrow.validation <- list("SpatialPolygonsRescale", layout.north.arrow(),  offset = c(arrow.lon.validation,arrow.lat.validation), scale = arrow.scale.value.validation, which = 1)
	scale.validation <- list("SpatialPolygonsRescale", layout.scale.bar(),  offset = c(scale.lon.validation,scale.lat.validation), scale = scale.scale.value.validation, fill=c("transparent","black"), which = 1)
	text.scale1.validation <- list("sp.text", c(text.scale1.lon.validation,text.scale1.lat.validation), "0", which = 1, cex=0.7)
	text.scale2.validation <- list("sp.text", c(text.scale2.lon.validation,text.scale2.lat.validation), paste(scale.scale.value.validation," m",sep=""), which = 1, cex=0.7)
	  
	# Plot and export of maps
	breaks.map.susceptibility<-c(0,0.2,0.45,0.55,0.8,1.0001)
	breaks.map.uncertainty<-c(0,0.001,0.005,0.01,0.05,0.1,0.5,1)
	breaks.map.matching.code<-c(1,2,3,4,5)
	color.vector.susceptibility<-color_ramp_palette_fun(length(breaks.histogram.values)-1)
	require(RColorBrewer)
	color.vector.uncertainty<-rev(brewer.pal(length(breaks.map.uncertainty)-1,"YlOrRd"))
	color.vector.matching<-c(rgb(38,115,0,max=255),rgb(233,255,190,max=255),rgb(215,215,158,max=255),rgb(115,115,0,max=255))
	}
  
  
#---------------------  COLLINEARITY EVALUATION ---------------------#

# Colldiag is an implementation of the regression collinearity diagnostic procedures found in
# Belsley, Kuh, and Welsch (1980). These procedures examine the ?conditioning? of the matrix
# of independent variables.

# Colldiag computes the condition indexes of the matrix. If the largest condition index
# (the condition number) is large (Belsley et al suggest 30 or higher), then there may be
# collinearity problems. All large condition indexes may be worth investigating.

# Colldiag also provides further information that may help to identify the source of these problems,
# the variance decomposition proportions associated with each condition index. If a large condition
# index is associated two or more variables with large variance decomposition proportions, these
# variables may be causing collinearity problems. Belsley et al suggest that a large proportion is 50
# percent or more.

if(enable_multicollinearity_test==TRUE)
	{
	#load collinearity package (perturb)
	library(perturb)
	#colnames(explanatory.variables)
	collinearity.test<-colldiag(explanatory.variables)
	#collinearity.test$condindx 
	#collinearity.test$pi 
	range(collinearity.test$condindx)
	
	if(range(collinearity.test$condindx)[2] >= 30) {      
		collinearity.value<-"Some explanatory variables are collinear"
		} else {      
		collinearity.value<-"Explanatory variables are not collinear"
		}
	
	print(collinearity.test,fuzz=.5)
	collinearity.evaluation.matrix<-print(collinearity.test,fuzz=.5)
	write.table("COLLINEARITY ANALYSIS RESULT",file="result_Collinearity_Analysis.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("EXPLANATION",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("This analysis was performed with Colldiag an implementation of the regression collinearity diagnostic procedures found in Belsley, Kuh,
	and Welsch (1980). These procedures examine the ?conditioning? of the matrix of independent variables. The procedure computes the condition
	indexes of the matrix. If the largest condition index (the condition number) is large (Belsley et al suggest 30 or higher), then there may be
	collinearity problems. All large condition indexes may be worth investigating. The procedure also provides further information that may help to
	identify the source of these problems, the variance decomposition proportions associated with each condition index. If a large condition
	index (> 30) is associated with two or more variables with large variance decomposition proportions, these variables may be causing collinearity problems.
	Belsley et al suggest that a large proportion is 50 percent or more.",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t",
	row.names=FALSE, col.names=FALSE)
	write.table("",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("RESULTS",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table(paste("Largest condition index (the condition number) =",range(collinearity.test$condindx)[2]),file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table(collinearity.value,file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table("Matrix of the variance decomposition proportions associated with each condition index (1st column)",file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table(rbind(colnames(collinearity.evaluation.matrix),collinearity.evaluation.matrix),file="result_Collinearity_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
	}



#----------------------- BEST MODEL SELECTION -----------------------#

#    ### BEST MODEL SELECTION USING AKAIKE INFORMATION CRITERION (AIC)
#    library(MASS)
#    
#    ### LINEAR MODEL
#    #Definition of lm General Model
#    model.lm<-lm(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables) 
#    #Definition of lm Upper Model
#    upper.model.lm<-as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+")))
#    #Definition of lm Lower Model
#    lower.model.lm<-as.formula(paste(names(data.variables)[1],"~",names(data.variables[2])))
#    
#    variable.selection.AIC.lm<-stepAIC(model.lm, scope=list(upper.model.lm,lower.model.lm), direction="both")
#    variable.selection.AIC.lm$call
#    variable.selection.AIC.lm$anova
#    
#    
#    ### GENERALIZED LINEAR MODEL
#    #Definition of glm General Model
#    model.glm<-glm(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables) 
#    
#    #Definition of glm Upper Model
#    upper.model.glm<-as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+")))
#    #Definition of glm Lower Model
#    lower.model.glm<-as.formula(paste(names(data.variables)[1],"~",names(data.variables[2])))
#    
#    variable.selection.AIC.glm<-stepAIC(model.glm, scope=list(upper.model.glm,lower.model.glm), direction="both")
#    variable.selection.AIC.glm$call
#    variable.selection.AIC.glm$anova 
#    
#    
#    
#    ### BEST MODEL SELECTION USING USCHI'S CLASSIFICATION PERFORMANCE MEASURES (UCPM) 
#    library(klaR)
#    variable.selection.ucpm.lda<-stepclass(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables,method="lda",direction="backward",improvement=0.01, criterion="CFvec")
#    variable.selection.ucpm.lda$call
#    variable.selection.ucpm.lda$model
#    variable.selection.ucpm.lda$process
#    names(variable.selection.ucpm.lda)
#    variable.selection.ucpm.lda$formula
#    
#    names(data.variables[,2:dim(data.variables)[2]])
#    
#    ### BEST MODEL SELECTION USING USCHI'S CLASSIFICATION PERFORMANCE MEASURES (UCPM) 
#    library(klaR)
#    variable.selection.wilks<-greedy.wilks(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables,niveau=0.1)
#    names(variable.selection.wilks)
#    
#    ### LINEAR MODEL
#    #Definition of lm General Model
#    model.lm<-lm(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables) 


#load package (MASS)
library(MASS)

#------------------- LINEAR DISCRIMINANT ANALISYS -------------------#

##### Linear Discriminant Analisys Run
if(model.run.matrix[1] == "YES")
  {

  #if(class(result.lda[1]) == "NULL") # Other Error selection criteria. In this case the IF istruction must be put at the end of the script 
  if (class(try(lda(explanatory.variables, grouping.variable, tol=0.001, method="moment")))=="try-error")  
    { 
    lda(explanatory.variables, grouping.variable, tol=0.001, method="moment")
    write.table("Linear Discriminant Analysis was not completed",file="Error_LDA_Analysis.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="Error_LDA_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("Error LOG",file="Error_LDA_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind("Message",rev(1:length(as.vector(.Traceback)))," ->",as.vector(.Traceback)),file="Error_LDA_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    }
    
  ##### Linear Discriminant Analisys using data.frame
  result.lda<-NULL
  result.lda<-lda(explanatory.variables, grouping.variable, tol=0.001, method="moment") 

  # Result Predicted
  predict.result.lda<-predict(result.lda)
  str(predict.result.lda)
  
  ### predict.result.lda$class is obtained rounding the posterior probability associated to 1 (predict.result.lda$posterior[,2]) 
  length(predict.result.lda$class[predict.result.lda$class==0])
  length(training.table[,2][training.table[,2]==0])
  
  cross.classification.lda<-table(grouping.variable,predict.result.lda$class,dnn=c("Observed","Predicted"))
  rownames(cross.classification.lda)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.lda)<-list("No Landslide","Landslide") # Predicted    
  str(cross.classification.lda)
  
  # Assignation of a matching code between observed and predicted values
  result.lda.matching.code<-paste(grouping.variable,as.numeric(levels(predict.result.lda$class))[predict.result.lda$class],sep="")
  result.lda.matching.code<-gsub("00","1",result.lda.matching.code)
  result.lda.matching.code<-gsub("01","2",result.lda.matching.code)
  result.lda.matching.code<-gsub("10","3",result.lda.matching.code)
  result.lda.matching.code<-gsub("11","4",result.lda.matching.code)
  result.lda.matching.code<-as.numeric(result.lda.matching.code)


  ##### Linear Discriminant Analisys using formula
  #read of the data
  #training.table<-read.table("data_provafrax_Export_Staffora_filtered.txt",header = TRUE,dec=".", sep="\t")
  #str(training.table)
  #training.table.names<-(colnames(training.table))
  #length(training.table$FRAXD2)
  
  #load package (MASS)
  #library(MASS)
  
  ##### Linear Discriminant Analisys
  #result.lda<-lda(FRAXD2 ~ ORDER + MAGN + LINK_LEN + AREAT_C + SLO_ARE + R + ELV_M + ELV_STD + SLO_ANG + SLO_ANG2 + ANG_STD + LNK_ANG + SLO_LEN + LEN_STD + ANGLE1 + ANGLE2 + ANGLE3 + CONC + CONV + COV_COC + COC_COV + RET + IRR + CC + AG_VA_PA + ALB_ZEB + ALLUVIO + AR_BIS + AR_R_M_P + AR_SCA + AT_PA_CA + CA_AN_CA + CA_PEN + DETRITO + MR_AN_LO + MR_B_R_C + MR_BOSM + MR_P_R_B + SACONG + INSIDE + ALV + BD + BMD + BR + INC + PRA + RIM + SEM + URB + VIG + REG + FRA + TRA + CAO + TR1 + TR2 + TR3 + FRA_OLD, training.table, tol=0.001, method="moment") 
  #str(result.lda)
  
  #cross.classification.lda<-table(training.table$FRAXD2[-result.lda$na.action],predict.result.lda$class,dnn=c("Observed","Predicted"))
  #rownames(cross.classification.lda)<-list("Landslides","No Landslides") # Observed
  #colnames(cross_classification_lda)<-list("Landslides","No Landslides") # Predicted    
  #str(cross.classification.lda)

  
  
  #Elaboration of Coefficient of association for contingency table 
  #load package (vcd)  
  library(vcd)
  
  #help(package=vcd)         
  contingency.table.lda<-table2d_summary(cross.classification.lda)
  test.table.lda<-assocstats(cross.classification.lda)
  #co_table(cross.classification.lda, margin=1)
  #mar_table(cross.classification.lda) 
  #structable(cross.classification.lda)

  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  fourfold(round(cross.classification.lda/sum(cross.classification.lda)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  #fourfold(cross.classification.lda, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    
  #dev.new()
  #spine(cross.classification.lda)
  #dev.new()       
  #mosaic(cross.classification.lda)
  #dev.new()       
  #pairs(cross.classification.lda)
  #dev.new()
  #sieve(cross.classification.lda)
  #dev.new()
  #strucplot(cross.classification.lda)
  #dev.new()       
  #cotabplot(cross.classification.lda)
  #dev.new()       
  #doubledecker(cross.classification.lda)
  #dev.new()       
  #grid_barplot(cross.classification.lda)
  }
  
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #A ROC curve plots the false alarm rate against the hit rate
  #for a probablistic forecast for a range of thresholds. 
  
  #load package (verification)  
  library(verification)
  
  #verify function
  #Based on the type of inputs, this function calculates a range of verification statistics and skill scores.
  #Additionally, it creates a verify class object that can be further analyzed.
  
  ##### ROC PLOT OBS - BINARY PREDICTION
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(as.numeric(training.table$FRAXD2[-result.lda$na.action]),as.numeric(predict.result.lda$class),main = "ROC PLOT: LINEAR DISCRIMINANT ANALYSIS MODEL - BINARY PREDICTED", binormal = TRUE, plot = "both")

  ##### ROC PLOT OBS - POSTERIOR PROBABILITY ASSOCIATED TO 1                                                                                 
  ## 1st method
  #dev.new()
  #roc.plot(training.table[,2],predict.result.lda$posterior[,2],main = "ROC PLOT: LINEAR DISCRIMINANT ANALYSIS MODEL", binormal = TRUE, plot = "both")
  #}
  
  # 2nd method using verify function
  verification.results.lda<-verify(training.table[,2],predict.result.lda$posterior[,2],frcst.type="prob", obs.type="binary")
  summary(verification.results.lda)
  #str(verification.results.lda)
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(verification.results.lda, main = "ROC PLOT: LINEAR DISCRIMINANT ANALYSIS MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #}
  area.under.roc.curve.lda<-roc.area(training.table[,2],predict.result.lda$posterior[,2])

  ## showing confidence intervals.  MAY BE SLOW
  
  if (cross.classification.lda[1,2]==0 | cross.classification.lda[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  roc.plot(verification.results.lda, main = "ROC PLOT: LINEAR DISCRIMINANT ANALYSIS MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[1] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.lda$A,2),";  Sample size = ",area.under.roc.curve.lda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[1], sep=""), side=3, col="red", cex=0.8)
  ## Histogram of posterior probability
  dev.new()
  hist(predict.result.lda$posterior[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Linear Disciminant Analysis susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  }
  pdf(file = "result_LDA_Histogram.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  hist(predict.result.lda$posterior[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Linear Disciminant Analysis susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  dev.off()

  # EXPORT OF PLOT FOR LDA MODEL
  pdf(file = "result_LDA_FourfoldPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.lda/sum(cross.classification.lda)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  #fourfold(cross.classification.lda, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  dev.off()
  
  #pdf(file = "result_LDA_ROCPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.results.lda, main = "ROC PLOT: LINEAR DISCRIMINANT ANALYSIS MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #dev.off()
  
  pdf(file = "result_LDA_ROCPlot_bootstrap.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.results.lda, main = "ROC PLOT: LINEAR DISCRIMINANT ANALYSIS MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[1] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.lda$A,2),";  Sample size = ",area.under.roc.curve.lda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[1], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  ## BOOTSTRAP PROCEDURE FOR THE ESTIMATION OF MODEL PREDICTION VARIABILITY
  if(bootstrap.model.variability[1] == "YES")
    {
    bootstrap.sample.model.lda<-bootstrap.sample.model[1]

    matrix.bootstrap.model.lda<-matrix(data=NA, nrow=dim(training.table)[1], ncol=(bootstrap.sample.model.lda*3)+1)
    colnames(matrix.bootstrap.model.lda)<-rep("na",(bootstrap.sample.model.lda*3)+1)
    matrix.bootstrap.model.lda[,1]<-identification.value
    colnames(matrix.bootstrap.model.lda)[1]<-"ID"
    name.sel.run<-paste(rep("ID_Selection_Run",bootstrap.sample.model.lda),1:bootstrap.sample.model.lda,sep="_")
    colnames(matrix.bootstrap.model.lda)[seq(2,(bootstrap.sample.model.lda*3)-1,3)]<-name.sel.run
    name.prob.run<-paste(rep("Probability_Run",bootstrap.sample.model.lda),1:bootstrap.sample.model.lda,sep="_")
    colnames(matrix.bootstrap.model.lda)[seq(3,(bootstrap.sample.model.lda*3),3)]<-name.prob.run
    name.pred.run<-paste(rep("Prediction_Run",bootstrap.sample.model.lda),1:bootstrap.sample.model.lda,sep="_")
    colnames(matrix.bootstrap.model.lda)[seq(4,(bootstrap.sample.model.lda*3)+1,3)]<-name.pred.run

    selection.index<-NULL
    library(MASS)
    #Bootstrap procedure
    for (count.boot in 1:bootstrap.sample.model.lda)
        {
        selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
        matrix.bootstrap.model.lda[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
        explanatory.variables.bootstrap.model.lda<-training.table[selection.index,3:dim(training.table)[2]]
        grouping.variable.bootstrap.model.lda<-as.factor(training.table[selection.index,2])
        #result.bootstrap.model.lda<-lda(explanatory.variables.bootstrap.model.lda, grouping.variable.bootstrap.model.lda, tol=0.001, method="moment")
        while(inherits(try(result.bootstrap.model.lda<-lda(explanatory.variables.bootstrap.model.lda, grouping.variable.bootstrap.model.lda, tol=0.001, method="moment"),silent=TRUE),what="try-error"))
			{
			print(paste("Count boot: ",count.boot," - Boostrap while resampling",sep=""))
			selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
			matrix.bootstrap.model.lda[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
			explanatory.variables.bootstrap.model.lda<-training.table[selection.index,3:dim(training.table)[2]]
			if(bootstrap_constant_correction==TRUE)
				{
				print("Performing bootstrap 0 value correction")
				indexbootstrapcosntant<-which(explanatory.variables.bootstrap.model.lda==0,arr.ind=TRUE)
				explanatory.variables.bootstrap.model.lda[indexbootstrapcosntant]<-runif(length(indexbootstrapcosntant)/2, min = 0.00001, max = 0.01)
				}
			grouping.variable.bootstrap.model.lda<-as.factor(training.table[selection.index,2])
			}
		#matrix.bootstrap.model.lda[as.numeric(names(table(selection.index))),(count.boot*3)]<-predict(result.bootstrap.model.lda,newdata=explanatory.variables[as.numeric(names(table(selection.index))),])$posterior[,2]
        matrix.bootstrap.model.lda[,(count.boot*3)+1]<-predict(result.bootstrap.model.lda,newdata=explanatory.variables)$posterior[,2]
		matrix.bootstrap.model.lda[as.numeric(names(table(selection.index))),(count.boot*3)]<-matrix.bootstrap.model.lda[as.numeric(names(table(selection.index))),(count.boot*3)+1]
		}
    # Export of bootstrap sample
    write.table(matrix.bootstrap.model.lda,file="result_LDA_BootstrapSamples.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)

    ID.bootstrap.model.lda.count<-numeric(length=dim(training.table)[1])
    #Probability (selected values)
    bootstrap.model.lda.probability.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.probability.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.probability.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.probability.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.probability.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.probability.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)

    #Prediction (all values)
    bootstrap.model.lda.prediction.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.prediction.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.prediction.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.prediction.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.prediction.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.lda.prediction.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
#    for (count.row.variability in 1:dim(training.table)[1])
#        {
#        # Statistics on boostrapped probability
#        ID.bootstrap.model.lda.count[count.row.variability]<-length(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(2,(bootstrap.sample.model.lda*3)-1,3)]))
#        bootstrap.model.lda.probability.mean[count.row.variability]<-mean(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(3,(bootstrap.sample.model.lda*3),3)]))
#        bootstrap.model.lda.probability.sd[count.row.variability]<-sd(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(3,(bootstrap.sample.model.lda*3),3)]))
#        bootstrap.model.lda.probability.min[count.row.variability]<-min(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(3,(bootstrap.sample.model.lda*3),3)]))
#        bootstrap.model.lda.probability.max[count.row.variability]<-max(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(3,(bootstrap.sample.model.lda*3),3)]))
#        bootstrap.model.lda.probability.sderror[count.row.variability]<-bootstrap.model.lda.probability.sd[count.row.variability]/ID.bootstrap.model.lda.count[count.row.variability]
#        bootstrap.model.lda.probability.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(3,(bootstrap.sample.model.lda*3),3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
#        # Statistics on boostrapped prediction
#        bootstrap.model.lda.prediction.mean[count.row.variability]<-mean(matrix.bootstrap.model.lda[count.row.variability,seq(4,(bootstrap.sample.model.lda*3)+1,3)])
#        bootstrap.model.lda.prediction.sd[count.row.variability]<-sd(matrix.bootstrap.model.lda[count.row.variability,seq(4,(bootstrap.sample.model.lda*3)+1,3)])
#        bootstrap.model.lda.prediction.min[count.row.variability]<-min(matrix.bootstrap.model.lda[count.row.variability,seq(4,(bootstrap.sample.model.lda*3)+1,3)])
#        bootstrap.model.lda.prediction.max[count.row.variability]<-max(matrix.bootstrap.model.lda[count.row.variability,seq(4,(bootstrap.sample.model.lda*3)+1,3)])
#        bootstrap.model.lda.prediction.sderror[count.row.variability]<-bootstrap.model.lda.prediction.sd[count.row.variability]/bootstrap.sample.model.lda
#        bootstrap.model.lda.prediction.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.lda[count.row.variability,seq(4,(bootstrap.sample.model.lda*3)+1,3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
#        }


		fun_length<-function(x) {length_var<-length(x[which(is.finite(x))]); return(length_var)}
		fun_quantile<-function(x) {quantile_var<-t(quantile(x[which(is.finite(x))],probs=c(0,0.05,0.25,0.5,0.75,0.95,1))); return(quantile_var)}
		ID.bootstrap.model.lda.count<-apply(matrix.bootstrap.model.lda[,grep("ID_Selection",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=fun_length)
		bootstrap.model.lda.probability.mean<-apply(matrix.bootstrap.model.lda[,grep("Probability",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=mean,na.rm = TRUE)
		bootstrap.model.lda.probability.sd<-apply(matrix.bootstrap.model.lda[,grep("Probability",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=sd,na.rm = TRUE)
		bootstrap.model.lda.probability.min<-apply(matrix.bootstrap.model.lda[,grep("Probability",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=min,na.rm = TRUE)
		bootstrap.model.lda.probability.max<-apply(matrix.bootstrap.model.lda[,grep("Probability",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=max,na.rm = TRUE)
		bootstrap.model.lda.probability.sderror<-bootstrap.model.lda.probability.sd/bootstrap.sample.model.lda
		bootstrap.model.lda.probability.quantiles<-apply(matrix.bootstrap.model.lda[,grep("Probability",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=fun_quantile)
		bootstrap.model.lda.prediction.mean<-apply(matrix.bootstrap.model.lda[,grep("Prediction",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=mean,na.rm = TRUE)
		bootstrap.model.lda.prediction.sd<-apply(matrix.bootstrap.model.lda[,grep("Prediction",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=sd,na.rm = TRUE)
		bootstrap.model.lda.prediction.min<-apply(matrix.bootstrap.model.lda[,grep("Prediction",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=min,na.rm = TRUE)
		bootstrap.model.lda.prediction.max<-apply(matrix.bootstrap.model.lda[,grep("Prediction",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=max,na.rm = TRUE)
		bootstrap.model.lda.prediction.sderror<-bootstrap.model.lda.prediction.sd/bootstrap.sample.model.lda
		bootstrap.model.lda.prediction.quantiles<-apply(matrix.bootstrap.model.lda[,grep("Prediction",colnames(matrix.bootstrap.model.lda))],MARGIN=1,FUN=fun_quantile)
		

		
		
		
    # Export of bootstrap sample statistics
    write.table(cbind("ID","LDA_NumberSelectedSamples","LDA_Probability_Mean","LDA_Probability_Sd","LDA_Probability_Min","LDA_Probability_Max","LDA_Probability_Sderror","LDA_Probability_Quantiles_0","LDA_Probability_Quantiles_0.05","LDA_Probability_Quantiles_0.25","LDA_Probability_Quantiles_0.5","LDA_Probability_Quantiles_0.75","LDA_Probability_Quantiles_0.95","LDA_Probability_Quantiles_1","LDA_Prediction_Mean","LDA_Prediction_Sd","LDA_Prediction_Min","LDA_Prediction_Max","LDA_Prediction_Sderror","LDA_Prediction_Quantiles_0","LDA_Prediction_Quantiles_0.05","LDA_Prediction_Quantiles_0.25","LDA_Prediction_Quantiles_0.5","LDA_Prediction_Quantiles_0.75","LDA_Prediction_Quantiles_0.95","LDA_Prediction_Quantiles_1"),file="result_LDA_BootstrapStatistics.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind(identification.value,ID.bootstrap.model.lda.count,bootstrap.model.lda.probability.mean,bootstrap.model.lda.probability.sd,bootstrap.model.lda.probability.min,bootstrap.model.lda.probability.max,bootstrap.model.lda.probability.sderror,t(bootstrap.model.lda.probability.quantiles),bootstrap.model.lda.prediction.mean,bootstrap.model.lda.prediction.sd,bootstrap.model.lda.prediction.min,bootstrap.model.lda.prediction.max,bootstrap.model.lda.prediction.sderror,t(bootstrap.model.lda.prediction.quantiles)),file="result_LDA_BootstrapStatistics.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    
    if (enable_screen_plotting==TRUE)
    {
    dev.new()
    plot(bootstrap.model.lda.probability.mean,bootstrap.model.lda.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="LDA BOOTSTRAP: Mean Probability vs Mean Prediction")
    abline(a=0,b=1,col="red",lty=1,lwd=1)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lda,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    }
    
    pdf(file = "result_LDA_BootstrapMeansComparison.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(bootstrap.model.lda.probability.mean,bootstrap.model.lda.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="LDA BOOTSTRAP: Mean Probability vs Mean Prediction")
    abline(a=0,b=1,col="red",lty=1,lwd=1)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lda,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    dev.off()

    #if (enable_screen_plotting==TRUE)
    #{
    #dev.new()
    #double.sd.histogram.variability<-hist(bootstrap.model.lda.probability.sd*2,breaks=seq(0,1,0.05),labels=TRUE)
    #plot(double.sd.histogram.variability$counts, seq(0,0.95,0.05), type="S",ylim=c(0,1),labels=TRUE)
    #}
    
    # BOOTSTRAPPED PROBABILITY - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.probability.lda<-cbind(bootstrap.model.lda.probability.mean,2*bootstrap.model.lda.probability.sd)
    parabola.probability.lda<-na.omit(parabola.probability.lda[order(parabola.probability.lda[,1]),])
    colnames(parabola.probability.lda)<-c("abscissa","ordinate")

    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.probability.lda <- nls(parabola.probability.lda[,"ordinate"] ~ coeff.a*(parabola.probability.lda[,"abscissa"]^2) + (-1)*coeff.a*parabola.probability.lda[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.probability.lda<-predict(fit.parabola.probability.lda)
    #coef(fit.parabola.probability.lda)
    
    if (enable_screen_plotting==TRUE)
    {
    dev.new()
    plot(parabola.probability.lda[,"abscissa"],parabola.probability.lda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="LDA Model Probability Variability (Bootstrap)")
    lines(parabola.probability.lda[,"abscissa"],value.parabola.probability.lda,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.lda),3),coeff.b= -round(coef(fit.parabola.probability.lda),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_LDA_BootstrapProbabilityVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.probability.lda[,"abscissa"],parabola.probability.lda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="LDA Model Probability Variability (Bootstrap)")
    lines(parabola.probability.lda[,"abscissa"],value.parabola.probability.lda,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.lda),3),coeff.b= -round(coef(fit.parabola.probability.lda),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
    
    # BOOTSTRAPPED PREDICTION - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.prediction.lda<-cbind(bootstrap.model.lda.prediction.mean,2*bootstrap.model.lda.prediction.sd)
    parabola.prediction.lda<-parabola.prediction.lda[order(parabola.prediction.lda[,1]),]
    colnames(parabola.prediction.lda)<-c("abscissa","ordinate")

    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.prediction.lda <- nls(parabola.prediction.lda[,"ordinate"] ~ coeff.a*(parabola.prediction.lda[,"abscissa"]^2) + (-1)*coeff.a*parabola.prediction.lda[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.prediction.lda<-predict(fit.parabola.prediction.lda)
    #coef(fit.parabola.prediction.lda)
    
    if (enable_screen_plotting==TRUE)
    {
    dev.new()
    plot(parabola.prediction.lda[,"abscissa"],parabola.prediction.lda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="LDA Model Prediction Variability (Bootstrap)")
    lines(parabola.prediction.lda[,"abscissa"],value.parabola.prediction.lda,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.lda),3),coeff.b= -round(coef(fit.parabola.prediction.lda),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_LDA_BootstrapPredictionVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.prediction.lda[,"abscissa"],parabola.prediction.lda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="LDA Model Prediction Variability (Bootstrap)")
    lines(parabola.prediction.lda[,"abscissa"],value.parabola.prediction.lda,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.lda),3),coeff.b= -round(coef(fit.parabola.prediction.lda),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
  }

  ## Sensitivity, Specificity, Cohens kappa plot
  dev.new()
  roc.plot.lda.series<-roc.plot(verification.results.lda,binormal=bin,plot=FALSE,show.thres=plot_thres,thresholds=threshold_series)
  dev.off()
  #str(roc.plot.lda.series)
  #roc.plot.lda.series$plot.data
  #str(roc.plot.lda.series$plot.data)

  ###########################################
  # min(abs(TPR - (1-FPR)))
  if(enable_probability_optimal_binary_classification==FALSE)
    {
    dev.new()
    roc.plot.lda.series<-roc.plot(verification.results.lda,binormal=bin,plot=FALSE,show.thres=FALSE,thresholds=threshold_series)
    dev.off()
    } else
    {
    dev.new()
    roc.plot.lda.series<-roc.plot(verification.results.lda,binormal=bin,plot=FALSE,show.thres=FALSE)
    dev.off()  
    }
  
  if(enable_probability_optimal_binary_classification==TRUE)
    {
    lda.probability.classification.optimal<-data.frame(prob_thres=roc.plot.lda.series$plot.data[,1,1],tpr=roc.plot.lda.series$plot.data[,2,1],fpr=roc.plot.lda.series$plot.data[,3,1],tnr=(1-roc.plot.lda.series$plot.data[,3,1]),diff_abs_tpr_tnr=abs(roc.plot.lda.series$plot.data[,2,1]-(1-roc.plot.lda.series$plot.data[,3,1])),optimal_sel=NA,breaks_sel=NA)
    index.lda.filter<-which(lda.probability.classification.optimal$prob_thres>0 & lda.probability.classification.optimal$prob_thres<1) # removing strnge thresh values
    lda.probability.classification.optimal<-rbind(c(0,1,1,0,1,NA,NA),lda.probability.classification.optimal[index.lda.filter,],c(1,0,0,1,1,NA,NA))
    lda.optimal.index<-which(lda.probability.classification.optimal$diff_abs_tpr_tnr==min(lda.probability.classification.optimal$diff_abs_tpr_tnr))
    lda.probability.classification.optimal$optimal_sel[lda.optimal.index]<-TRUE
    lda.probability.optimal.binary.threshold<-lda.probability.classification.optimal$prob_thres[lda.optimal.index]
    
    ### Generating the optimal fourfold plot
    cross.classification.lda.optimal<-table(grouping.variable,predict.result.lda$posterior[,2]>lda.probability.optimal.binary.threshold,dnn=c("Observed","Predicted"))
    rownames(cross.classification.lda.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.lda.optimal)<-list("No Landslide","Landslide") # Predicted    
    str(cross.classification.lda.optimal)
    # Assignation of a matching code between observed and predicted values
    result.lda.matching.code.optimal<-paste(grouping.variable,as.numeric(predict.result.lda$posterior[,2]>lda.probability.optimal.binary.threshold),sep="")
    result.lda.matching.code.optimal<-gsub("00","1",result.lda.matching.code.optimal)
    result.lda.matching.code.optimal<-gsub("01","2",result.lda.matching.code.optimal)
    result.lda.matching.code.optimal<-gsub("10","3",result.lda.matching.code.optimal)
    result.lda.matching.code.optimal<-gsub("11","4",result.lda.matching.code.optimal)
    result.lda.matching.code.optimal<-as.numeric(result.lda.matching.code.optimal)
    
    if (enable_screen_plotting==TRUE)
      {
      dev.new()
      fourfold(round(cross.classification.lda.optimal/sum(cross.classification.lda.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
      }
    # EXPORT OF PLOT FOR LDA MODEL
    pdf(file = "result_LDA_FourfoldPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.lda.optimal/sum(cross.classification.lda.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.lda.optimal, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    dev.off()
    
    ### Optimal susceptibility classes identification 
    if(enable_probability_optimal_classification==TRUE)
      {
      lda.unexplained.errors<-round(1-(lda.probability.classification.optimal$tpr[lda.optimal.index]+lda.probability.classification.optimal$tpr[lda.optimal.index])/2,2)
  
      if(type_probability_optimal_classification=="proportional")
        {
        lda.unexplained.errors.partition<-1-(lda.unexplained.errors/((length(breaks.histogram.values)/2)))*(1:((length(breaks.histogram.values)/2)))
        
        for(count_part in 1:(length(lda.unexplained.errors.partition)-1))
          {
          #count_part<-1
          unexplained.errors.partition.sel<-lda.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((lda.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          lda.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((lda.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          lda.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
          }
        lda.breaks.histogram.values.optimal<-c(0,lda.probability.classification.optimal$prob_thres[which(lda.probability.classification.optimal$breaks_sel==TRUE)],1)
        #lda.probability.classification.optimal[which(lda.probability.classification.optimal$breaks_sel==TRUE),]
        }
      
      if(type_probability_optimal_classification=="fixed")
        {
        step.lda.unexplained.fixed<-0.1
        if(lda.unexplained.errors<=0.1) step.lda.unexplained.fixed<-0.05
        if(lda.unexplained.errors<=0.05) step.lda.unexplained.fixed<-0.025
        lda.unexplained.errors.partition<-seq(step.lda.unexplained.fixed,1-step.lda.unexplained.fixed,step.lda.unexplained.fixed)[seq(step.lda.unexplained.fixed,1-step.lda.unexplained.fixed,step.lda.unexplained.fixed)>(1-lda.unexplained.errors)]
  
        for(count_part in 1:(length(lda.unexplained.errors.partition)-1))
          {
          #count_part<-1
          unexplained.errors.partition.sel<-lda.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((lda.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          lda.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((lda.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          lda.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
          }
        lda.breaks.histogram.values.optimal<-c(0,lda.probability.classification.optimal$prob_thres[which(lda.probability.classification.optimal$breaks_sel==TRUE)],1)
        #lda.probability.classification.optimal[which(lda.probability.classification.optimal$breaks_sel==TRUE),]
        }
      
      if (enable_screen_plotting==TRUE)
        {
        dev.new()
        hist(predict.result.lda$posterior[,2], breaks=lda.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal LDA susceptibility", col=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1))
        plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("LDA OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
        mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
        mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
        for (count in 1:(length(lda.breaks.histogram.values.optimal)-1))
          {
          #count=1
          polygon(c(lda.breaks.histogram.values.optimal[count:(count+1)],rev(lda.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1)[count])  
          }
        polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
        lines(lda.probability.classification.optimal$prob_thres,lda.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
        lines(lda.probability.classification.optimal$prob_thres,lda.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
        index_points_plot<-c(1,which(lda.probability.classification.optimal$breaks_sel==TRUE),dim(lda.probability.classification.optimal)[1])
        points(lda.probability.classification.optimal[index_points_plot[1:floor(length(lda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
        text(lda.probability.classification.optimal[index_points_plot[1:floor(length(lda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(lda.probability.classification.optimal[index_points_plot[1:floor(length(lda.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(lda.breaks.histogram.values.optimal)/2)-1)))
        points(lda.probability.classification.optimal[index_points_plot[ceiling(length(lda.breaks.histogram.values.optimal)/2):length(lda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
        text(lda.probability.classification.optimal[index_points_plot[ceiling(length(lda.breaks.histogram.values.optimal)/2):length(lda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(lda.probability.classification.optimal[index_points_plot[ceiling(length(lda.breaks.histogram.values.optimal)/2):length(lda.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(lda.breaks.histogram.values.optimal)-ceiling(length(lda.breaks.histogram.values.optimal)/2)),3))
        }
  
      pdf(file = "result_LDA_Histogram_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      hist(predict.result.lda$posterior[,2], breaks=lda.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal LDA susceptibility", col=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1))
      dev.off()
  
      pdf(file = "result_LDA_ModelEvaluationPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("LDA OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
      mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
      mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
      for (count in 1:(length(lda.breaks.histogram.values.optimal)-1))
        {
        #count=1
        polygon(c(lda.breaks.histogram.values.optimal[count:(count+1)],rev(lda.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1)[count])  
        }
      polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
      lines(lda.probability.classification.optimal$prob_thres,lda.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
      lines(lda.probability.classification.optimal$prob_thres,lda.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
      index_points_plot<-c(1,which(lda.probability.classification.optimal$breaks_sel==TRUE),dim(lda.probability.classification.optimal)[1])
      points(lda.probability.classification.optimal[index_points_plot[1:floor(length(lda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(lda.probability.classification.optimal[index_points_plot[1:floor(length(lda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(lda.probability.classification.optimal[index_points_plot[1:floor(length(lda.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(lda.breaks.histogram.values.optimal)/2)-1)))
      points(lda.probability.classification.optimal[index_points_plot[1+ceiling(length(lda.breaks.histogram.values.optimal)/2):length(lda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
      text(lda.probability.classification.optimal[index_points_plot[1+ceiling(length(lda.breaks.histogram.values.optimal)/2):length(lda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(lda.probability.classification.optimal[index_points_plot[1+ceiling(length(lda.breaks.histogram.values.optimal)/2):length(lda.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(lda.breaks.histogram.values.optimal)-1-ceiling(length(lda.breaks.histogram.values.optimal)/2)),3))
      points(lda.probability.classification.optimal[lda.optimal.index,c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(lda.probability.classification.optimal[lda.optimal.index,c("prob_thres","tpr")],labels=with(round(lda.probability.classification.optimal[lda.optimal.index,c("prob_thres","tpr")],3),paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(4))
      dev.off()
    }
    }
    
  ###########################################

  contingency.table.matrix.lda<-matrix(nrow=dim(roc.plot.lda.series$plot.data)[1],ncol=8)
  colnames(contingency.table.matrix.lda)<-c("Threshold","TP","TN","FP","FN","TPR","FPR","COHEN_KAPPA")
  contingency.table.matrix.lda[,1]<-roc.plot.lda.series$plot.data[,1,1]
  contingency.table.matrix.lda[,6]<-roc.plot.lda.series$plot.data[,2,1]
  contingency.table.matrix.lda[,7]<-roc.plot.lda.series$plot.data[,3,1]
  values.observed<-training.table[,2]
  values.predicted<-predict.result.lda$posterior[,2]
  for (count.threshold.series in 1:dim(roc.plot.lda.series$plot.data)[1])
      {
      value.threshold<-contingency.table.matrix.lda[count.threshold.series,1]
      values.probability.reclassified<-NULL
      values.probability.reclassified<-as.numeric(values.predicted>value.threshold) 
      #sum(values.probability.reclassified-round(values.predicted)) # Check sum: It has to be 0 if threshold is equal to 1
      series.pasted<-paste(values.observed,values.probability.reclassified,sep="")
      series.pasted<-gsub("00","1",series.pasted)
      series.pasted<-gsub("01","2",series.pasted)
      series.pasted<-gsub("10","3",series.pasted)
      series.pasted<-gsub("11","4",series.pasted)
      series.pasted<-as.numeric(series.pasted)
      TP<-as.numeric(sum(series.pasted>=4)) # True Positive
      FN<-as.numeric(sum(series.pasted>=3 & series.pasted<4)) # False Negative
      FP<-as.numeric(sum(series.pasted>=2 & series.pasted<3)) # False Positive
      TN<-as.numeric(sum(series.pasted>=1 & series.pasted<2)) # True Negative              
      #TPR<-TP/(TP+FN) # Hit Rate or True Positive Rate or Sensitivity - Assigned before the for cicle using rocplot data
      #FPR<-FP/(FP+TN) # False Alarm Rate or False Positive Rate or 1-Specificity
      # Cohen's Kappa = (agreement-chance)/(1-chance)  where agreement=(TP+TN)/(TP+TN+FP+FN) and chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
      agreement=(TP+TN)/(TP+TN+FP+FN)
      chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
      cohen.kappa.value<-(agreement-chance)/(1-chance)
      #Other
      #library(vcd)
      #cohen.kappa.value<-Kappa(cross.classification.table)
      contingency.table.matrix.lda[count.threshold.series,2]<-TP
      contingency.table.matrix.lda[count.threshold.series,3]<-TN
      contingency.table.matrix.lda[count.threshold.series,4]<-FP
      contingency.table.matrix.lda[count.threshold.series,5]<-FN
      contingency.table.matrix.lda[count.threshold.series,8]<-cohen.kappa.value
      print(value.threshold)
      }

  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(roc.plot.lda.series$plot.data[,1,1],roc.plot.lda.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="LDA MODEL EVALUATION PLOT")
  lines(roc.plot.lda.series$plot.data[,1,1],1-roc.plot.lda.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.lda.series$plot.data[,1,1], contingency.table.matrix.lda[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  }
  pdf(file = "result_LDA_ModelEvaluationPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(roc.plot.lda.series$plot.data[,1,1],roc.plot.lda.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="LDA MODEL EVALUATION PLOT")
  lines(roc.plot.lda.series$plot.data[,1,1],1-roc.plot.lda.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.lda.series$plot.data[,1,1], contingency.table.matrix.lda[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  dev.off()

  ## VALIDATION OF LDA MODEL (Matching LDA posterior probability results and validation grouping variable)
  
  # Result Predicted
  predict.result.lda.validation<-predict(result.lda,newdata=validation.explanatory.variables)
  #str(predict.result.lda.validation)
  #summary(predict.result.lda.validation)
  cross.classification.validation.lda<-table(validation.grouping.variable,predict.result.lda.validation$class,dnn=c("Observed","Predicted"))
  rownames(cross.classification.validation.lda)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.validation.lda)<-list("No Landslide","Landslide") # Predicted
  str(cross.classification.validation.lda)
  #cross.classification.lda.validation<-table(validation.grouping.variable,predict.result.lda$class,dnn=c("Observed","Predicted"))
  
 
  #Elaboration of Coefficient of association for contingency table
  #load package (vcd)
  library(vcd)

  #help(package=vcd)
  contingency.table.validation.lda<-table2d_summary(cross.classification.validation.lda)
  test.table.validation.lda<-assocstats(cross.classification.validation.lda)

  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  fourfold(round(cross.classification.validation.lda/sum(cross.classification.validation.lda)*100,2), std="margin", main="VALIDATION LDA MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  #fourfold(cross.classification.validation.lda, std="margin", main="VALIDATION LDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  }
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #load package (verification)
  library(verification)

  # 2nd method using verify function
  verification.validation.lda<-verify(validation.table[,2],predict.result.lda.validation$posterior[,2], frcst.type="prob", obs.type="binary")
  #summary(verification.validation.lda)

  # showing confidence intervals.  MAY BE SLOW
  area.under.roc.curve.validation.lda<-roc.area(validation.table[,2],predict.result.lda.validation$posterior[,2])
  
  if (cross.classification.validation.lda[1,2]==0 | cross.classification.validation.lda[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  roc.plot(verification.validation.lda, main = "ROC PLOT: VALIDATION LDA MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[1] , alpha = 0.05, extra=TRUE, legend=TRUE, show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.validation.lda$A,2),";  Sample size = ",area.under.roc.curve.validation.lda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[1], sep=""), side=3, col="red", cex=0.8)
  }

  # EXPORT OF PLOT FOR VALIDATION OF LDA MODEL

  pdf(file = "result_LDA_FourfoldPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.validation.lda/sum(cross.classification.validation.lda)*100,2), std="margin", main="VALIDATION LDA MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  #fourfold(cross.classification.validation.lda, std="margin", main="VALIDATION LDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  dev.off()

  #pdf(file = "result_LDA_ROCPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.validation.lda, main = "ROC PLOT: VALIDATION LDA MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.validation.lda<-roc.area(training.table[,2],predict.result.lda.validation$posterior[,2])
  #dev.off()

  pdf(file = "result_LDA_ROCPlot_bootstrap_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.validation.lda, main = "ROC PLOT: VALIDATION LDA MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[1] , alpha = 0.05, extra=TRUE, legend=TRUE, show.thres=plot_thres,thresholds=NULL)
  mtext(paste("ROC area = ",round(area.under.roc.curve.validation.lda$A,2),";  Sample size = ",area.under.roc.curve.validation.lda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[1], sep=""), side=3, col="red", cex=0.8)
  dev.off()

  # Assignation of a matching code between observed and predicted values calculated using the validation dataset
  validation.lda.matching.code<-paste(validation.grouping.variable,as.numeric(levels(predict.result.lda.validation$class))[predict.result.lda.validation$class],sep="")
  validation.lda.matching.code<-gsub("00","1",validation.lda.matching.code)
  validation.lda.matching.code<-gsub("01","2",validation.lda.matching.code)
  validation.lda.matching.code<-gsub("10","3",validation.lda.matching.code)
  validation.lda.matching.code<-gsub("11","4",validation.lda.matching.code)
  validation.lda.matching.code<-as.numeric(validation.lda.matching.code)

  ##########################################
  if(enable_probability_optimal_binary_classification==TRUE)
    {
    cross.classification.validation.lda.optimal<-table(validation.grouping.variable,as.numeric(predict.result.lda.validation$posterior[,2]>lda.probability.optimal.binary.threshold),dnn=c("Observed","Predicted"))
    rownames(cross.classification.validation.lda.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.validation.lda.optimal)<-list("No Landslide","Landslide") # Predicted
    
    #Different plots for contingency table
    if (enable_screen_plotting==TRUE)
      {
      dev.new()
      fourfold(round(cross.classification.validation.lda.optimal/sum(cross.classification.validation.lda.optimal)*100,2), std="margin", main="VALIDATION LDA MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
      #fourfold(cross.classification.validation.lda, std="margin", main="VALIDATION LDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
      }
    if (cross.classification.validation.lda.optimal[1,2]==0 | cross.classification.validation.lda.optimal[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
    pdf(file = "result_LDA_FourfoldPlot_Validation_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.validation.lda.optimal/sum(cross.classification.validation.lda.optimal)*100,2), std="margin", main="VALIDATION LDA MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.lda, std="margin", main="VALIDATION LDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    dev.off()
    
    
    # Assignation of a optimal matching code between observed and predicted values calculated using the validation dataset
    validation.lda.matching.code.optimal<-paste(validation.grouping.variable,as.numeric(predict.result.lda.validation$posterior[,2]>lda.probability.optimal.binary.threshold),sep="")
    validation.lda.matching.code.optimal<-gsub("00","1",validation.lda.matching.code.optimal)
    validation.lda.matching.code.optimal<-gsub("01","2",validation.lda.matching.code.optimal)
    validation.lda.matching.code.optimal<-gsub("10","3",validation.lda.matching.code.optimal)
    validation.lda.matching.code.optimal<-gsub("11","4",validation.lda.matching.code.optimal)
    validation.lda.matching.code.optimal<-as.numeric(validation.lda.matching.code.optimal)
    }
  #########################################
  

  # EXPORT OF LDA MODEL RESULTS
  write.table("RESULTS OF LINEAR DISCRIMINANT ANALYSIS",file="result_LDA.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("LDA MODEL OUTPUTS",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Prior Probabilities",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(names(result.lda$prior),result.lda$prior),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Total number",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind("N",result.lda$N),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Counts",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(names(result.lda$counts),result.lda$counts),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Means",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(t(rbind(c("",colnames(result.lda$means)),cbind(rownames(result.lda$means),result.lda$means))),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Discriminant function coefficients",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  #write.table(cbind(rownames(result.lda$scaling),result.lda$scaling),file="result_LDA.tsv", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  scaling.order<-result.lda$scaling[order(result.lda$scaling),]
  scaling.matrix.ordered<-cbind(names(scaling.order),scaling.order)
  write.table(scaling.matrix.ordered,file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE MODEL RESULT",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.lda$table[,1,],contingency.table.lda$table[,2,],contingency.table.lda$table[,3,])),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE VALIDATION",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.validation.lda$table[,1,],contingency.table.validation.lda$table[,2,],contingency.table.validation.lda$table[,3,])),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("MATCHING CODE DEFINITION",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(c("","OBSERVED NO LANDSLIDES: 0","OBSERVED LANDSLIDES: 1"), c("PREDICTED NO LANDSLIDES: 0","00 -> Code 1","10 -> Code 3"), c("PREDICTED LANDSLIDES: 1","01 -> Code 2","11 -> Code 4")),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  ########
  if(enable_probability_optimal_binary_classification==FALSE) 
    {
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","LDA MODEL POSTERIOR PROBABILITY","LDA MODEL CLASSIFICATION","LDA MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict.result.lda$posterior[,2],as.numeric(levels(predict.result.lda$class))[predict.result.lda$class],result.lda.matching.code)),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","LDA VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.lda.validation$posterior[,2],as.numeric(levels(predict.result.lda.validation$class))[predict.result.lda.validation$class],validation.lda.matching.code)),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    } else
    {
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    if(enable_probability_optimal_classification==TRUE) 
      {
      write.table(paste("OPTIMAL SUSCEPTIBILITY PARTITION -> Method: ",type_probability_optimal_classification,sep=""),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(lda.probability.classification.optimal[index_points_plot,c("prob_thres","tnr","tpr")]),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
      } else
      {
      write.table(paste("OPTIMAL SUSCEPTIBILITY BINARY PARTITION",sep=""),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(lda.probability.classification.optimal[lda.optimal.index,c("prob_thres","tnr","tpr")]),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
      }
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","LDA MODEL POSTERIOR PROBABILITY","LDA MODEL CLASSIFICATION","LDA MODEL RESULT MATCHING CODE","LDA OPTIMAL MODEL CLASSIFICATION","LDA OPTIMAL MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict.result.lda$posterior[,2],as.numeric(levels(predict.result.lda$class))[predict.result.lda$class],result.lda.matching.code,as.numeric(predict.result.lda$posterior[,2]>lda.probability.optimal.binary.threshold),result.lda.matching.code.optimal)),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","LDA VALIDATION MATCHING CODE","OPTIMAL VALIDATION CLASSIFICATION","OPTIMAL LDA VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.lda.validation$posterior[,2],as.numeric(levels(predict.result.lda.validation$class))[predict.result.lda.validation$class],validation.lda.matching.code,as.numeric(predict.result.lda.validation$posterior[,2]>lda.probability.optimal.binary.threshold),validation.lda.matching.code.optimal)),file="result_LDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    }
  ########

  # PLOT AND EXPORT OF LDA MODEL MAPS
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POLYGONS")
	    {
      ################################################
      ### Prima di cancellare testare con polygoni
      ##############################################
	    #shape_training_lda<-shape_training
	    #result_training_lda_shape<-cbind(identification.value,training.table[,2],predict.result.lda$posterior[,2],as.numeric(levels(predict.result.lda$class))[predict.result.lda$class],result.lda.matching.code)
	    #colnames(result_training_lda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH")
		  #shape_training_lda@data <- merge(x=shape_training_lda@data,y=result_training_lda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
		  ##writeOGR(shape_training_lda,dsn="result_LDA_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
	    #writeOGR(shape_training_lda,dsn="result_LDA_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
	            
	    #shape_validation_lda<-shape_validation
	    #result_validation_lda_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.lda.validation$posterior[,2],as.numeric(levels(predict.result.lda.validation$class))[predict.result.lda.validation$class],validation.lda.matching.code)
	    #colnames(result_validation_lda_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
	  	#shape_validation_lda@data <- merge(x=shape_validation_lda@data,y=result_validation_lda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
	    ##writeOGR(shape_validation_lda,dsn="result_LDA_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
	  	#writeOGR(shape_validation_lda,dsn="result_LDA_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
	    ##################################################################
      ### DA qui Nuovo
      ##################################################################   
      #### aggiungere a poligoni colonna incertezza ed esportazione pdf incertezza
	    shape_training_lda<-shape_training
	    result_training_lda_shape<-cbind(identification.value,training.table[,2],predict.result.lda$posterior[,2],as.numeric(levels(predict.result.lda$class))[predict.result.lda$class],result.lda.matching.code,ID.bootstrap.model.lda.count,bootstrap.model.lda.probability.mean,bootstrap.model.lda.probability.sd,bootstrap.model.lda.probability.min,bootstrap.model.lda.probability.max,bootstrap.model.lda.probability.sderror,t(bootstrap.model.lda.probability.quantiles),bootstrap.model.lda.prediction.mean,bootstrap.model.lda.prediction.sd,bootstrap.model.lda.prediction.min,bootstrap.model.lda.prediction.max,bootstrap.model.lda.prediction.sderror,t(bootstrap.model.lda.prediction.quantiles))
	    colnames(result_training_lda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LDA_SAMP","LDA_PMean","LDA_PSd","LDA_PMin","LDA_PMax","LDA_PSder","LDA_PQ0","LDA_PQ_005","LDA_PQ_025","LDA_PQ05","LDA_PQ_075","LDA_PQ095","LDA_PQ1","LDA_PrMean","LDA_PrSd","LDA_PrMin","LDA_PrMax","LDA_PrSder","LDA_PrQ0","LDA_PrQ005","LDA_PrQ025","LDA_PrQ05","LDA_PrQ075","LDA_PrQ095","LDA_PrQ1")
	    ###########
	    if(enable_probability_optimal_binary_classification==TRUE) 
	      {
	      result_training_lda_shape<-cbind(result_training_lda_shape,as.numeric(predict.result.lda$posterior[,2]>lda.probability.optimal.binary.threshold),result.lda.matching.code.optimal)
	      colnames(result_training_lda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LDA_SAMP","LDA_PMean","LDA_PSd","LDA_PMin","LDA_PMax","LDA_PSder","LDA_PQ0","LDA_PQ_005","LDA_PQ_025","LDA_PQ05","LDA_PQ_075","LDA_PQ095","LDA_PQ1","LDA_PrMean","LDA_PrSd","LDA_PrMin","LDA_PrMax","LDA_PrSder","LDA_PrQ0","LDA_PrQ005","LDA_PrQ025","LDA_PrQ05","LDA_PrQ075","LDA_PrQ095","LDA_PrQ1","OPT_CLASS","OPT_MATCH")
	      }
	    #############
	    
	    shape_training_lda@data <- merge(x=shape_training_lda@data,y=result_training_lda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
	    #writeOGR(shape_training_lda,dsn="result_LDA_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
	    writeOGR(shape_training_lda,dsn="result_LDA_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
	    # WARNING: The validation does't have the unvertainty estimation: probably this can be daone using the parabolic error function 
	    shape_validation_lda<-shape_validation
	    result_validation_lda_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.lda.validation$posterior[,2],as.numeric(levels(predict.result.lda.validation$class))[predict.result.lda.validation$class],validation.lda.matching.code)
	    colnames(result_validation_lda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
	    ###########
	    if(enable_probability_optimal_binary_classification==TRUE) 
	      {
	      result_validation_lda_shape<-cbind(result_validation_lda_shape,as.numeric(predict.result.lda.validation$posterior[,2]>lda.probability.optimal.binary.threshold),validation.lda.matching.code.optimal)
	      colnames(result_validation_lda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
	      }
	    ############
	    shape_validation_lda@data <- merge(x=shape_validation_lda@data,y=result_validation_lda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
	    shape_validation_lda@data <- cbind(shape_validation_lda@data,PROB_SDMOD=(coefficients(fit.parabola.probability.lda)*(shape_validation_lda@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.lda)*shape_validation_lda@data$VAL_PROB))
	    #writeOGR(shape_validation_lda,dsn="result_LDA_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
	    writeOGR(shape_validation_lda,dsn="result_LDA_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
		
	    # Plot and export of maps
	    # LDA Susceptibility
	    #dev.new()
	    pdf(file = "result_LDA_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
	    print(spplot(obj=shape_training_lda, zcol=c("MOD_PROB"), names.attr=c("LDA MODEL PROBABILITY"), main="LDA MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
	    dev.off()
      
  
	    # LDA Model Matching Code
	    #dev.new()
	    pdf(file = "result_LDA_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
	    print(spplot(obj=shape_training_lda, zcol=c("MOD_MATCH"), names.attr=c("LDA MODEL MATCHING CODE"), main="LDA MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
	    dev.off()
	    
	    # LDA Model Uncertainty
	    #dev.new()
	    pdf(file = "result_LDA_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
	    print(spplot(obj=shape_training_lda, zcol=c("LDA_PrSd"), names.attr=c("LDA MODEL UNCERTAINTY"), main="LDA MODEL UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
	    dev.off()

	    # LDA Validation Susceptibility
	    #dev.new()
	    pdf(file = "result_LDA_Validation_Susceptibility_Map.pdf")
	    print(spplot(obj=shape_validation_lda, zcol=c("VAL_PROB"), names.attr=c("LDA VALIDATION PROBABILITY"), main="LDA VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
	    dev.off()
	    
	    # LDA Validation Matching Code
	    #dev.new()
	    pdf(file = "result_LDA_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
	    print(spplot(obj=shape_validation_lda, zcol=c("VAL_MATCH"), names.attr=c("LDA VALIDATION MATCHING CODE"), main="LDA VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
	    dev.off()
	    
	    # LDA Validation Uncertainty
	    #dev.new()
	    pdf(file = "result_LDA_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
	    print(spplot(obj=shape_validation_lda, zcol=c("PROB_SDMOD"), names.attr=c("LDA VALIDATION UNCERTAINTY"), main="LDA VALIDATION UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
	    dev.off()

      
		  ########### TBT
		  if(enable_probability_optimal_binary_classification==TRUE) 
  		  {
		    # LDA Model Matching Code
		    #dev.new()
		    pdf(file = "result_LDA_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		    print(spplot(obj=shape_training_lda, zcol=c("OPT_MATCH"), names.attr=c("LDA MODEL MATCHING CODE"), main="LDA MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
		    dev.off()
        
		    # LDA Model Matching Code
		    #dev.new()
		    pdf(file = "result_LDA_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		    print(spplot(obj=shape_validation_lda, zcol=c("OPT_MATCH"), names.attr=c("LDA VALIDATION MATCHING CODE"), main="LDA VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
		    dev.off()
		    		    
		    if(enable_probability_optimal_classification==TRUE)
		      {
		      pdf(file = "result_LDA_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		      print(spplot(obj=shape_training_lda, zcol=c("MOD_PROB"), names.attr=c("LDA MODEL PROBABILITY"), main="LDA MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=(lda.breaks.histogram.values.optimal)+c(rep(0,(length(lda.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1)))
		      dev.off()
          
		      pdf(file = "result_LDA_Validation_Susceptibility_Map_Optimal.pdf")
		      print(spplot(obj=shape_validation_lda, zcol=c("VAL_PROB"), names.attr=c("LDA VALIDATION PROBABILITY"), main="LDA VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=(lda.breaks.histogram.values.optimal)+c(rep(0,(length(lda.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1)))
		      dev.off()
          }
		    }
		  ###########
      
	    
	 
	    ### Success & prediction rate curve
	    susc_values_lda<-c(1,0.8,0.55,0.45,0.2,0)
	    
	    # ordering data for susceptibility
	    ind_col_succ_pred_rate<-which(colnames(shape_training_lda@data) %in% c(configuration.spatial.data.table[c(5,6)],"MOD_PROB"))
	    shape_training_lda@data<-shape_training_lda@data[order(shape_training_lda@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
	    shape_training_lda@data[,1]<-cumsum(shape_training_lda@data[,1])/sum(shape_training_lda@data[,1])*100
	    shape_training_lda@data[,2]<-cumsum(shape_training_lda@data[,2])/sum(shape_training_lda@data[,2])*100
	    
	    ind_pro_fun_mod<-approxfun(shape_training_lda@data[,3],shape_training_lda@data[,1])
	    area_susc_values_lda_mod<-c(0,ind_pro_fun_mod(susc_values_lda[-c(1,length(susc_values_lda))]),100)
	    
	    #dev.new()
	    pdf(file = "result_LDA_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
	    plot(0,0,col="transparent",main="LDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
	    for (count in 1:(length(susc_values_lda)-1))
	      {
	      #count=1
	      polygon(c(area_susc_values_lda_mod[count:(count+1)],rev(area_susc_values_lda_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
	      }
	    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
	    lines(shape_training_lda@data[,1],shape_training_lda@data[,2],col="black")
	    dev.off()
      
      
		  ########### TBT
		  if(enable_probability_optimal_classification==TRUE) 
		  {
		    area_susc_values_lda_mod_optimal<-c(0,ind_pro_fun_mod(rev(lda.breaks.histogram.values.optimal)[-c(1,length(lda.breaks.histogram.values.optimal))]),100)
		    pdf(file = "result_LDA_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		    plot(0,0,col="transparent",main="LDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
		    for (count in 1:(length(rev(lda.breaks.histogram.values.optimal))-1))
		    {
		      #count=1
		      polygon(c(area_susc_values_lda_mod_optimal[count:(count+1)],rev(area_susc_values_lda_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1))[count])  
		    }
		    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
		    lines(shape_training_lda@data[,1],shape_training_lda@data[,2],col="black")
		    dev.off()
		  }
		  ############
	    
	    
	    ind_col_succ_pred_rate<-which(colnames(shape_validation_lda@data) %in% c(configuration.spatial.data.table[c(5,6)],"VAL_PROB"))
	    shape_validation_lda@data<-shape_validation_lda@data[order(shape_validation_lda@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
	    shape_validation_lda@data[,1]<-cumsum(shape_validation_lda@data[,1])/sum(shape_validation_lda@data[,1])*100
	    shape_validation_lda@data[,2]<-cumsum(shape_validation_lda@data[,2])/sum(shape_validation_lda@data[,2])*100
	    
	    ind_pro_fun_val<-approxfun(shape_validation_lda@data[,3],shape_validation_lda@data[,1])
	    area_susc_values_lda_val<-c(0,ind_pro_fun_val(susc_values_lda[-c(1,length(susc_values_lda))]),100)
	    
	    #dev.new()
	    pdf(file = "result_LDA_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
	    plot(0,0,col="transparent",main="LDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
	    for (count in 1:(length(susc_values_lda)-1))
	      {
	      #count=1
	      polygon(c(area_susc_values_lda_val[count:(count+1)],rev(area_susc_values_lda_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
	      }
	    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
	    lines(shape_validation_lda@data[,1],shape_validation_lda@data[,2],col="black")
	    dev.off()
            
		  ########### TBT
		  if(enable_probability_optimal_classification==TRUE) 
		    {
		    area_susc_values_lda_val_optimal<-c(0,ind_pro_fun_val(rev(lda.breaks.histogram.values.optimal)[-c(1,length(lda.breaks.histogram.values.optimal))]),100)
		    pdf(file = "result_LDA_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		    plot(0,0,col="transparent",main="LDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
		    for (count in 1:(length(rev(lda.breaks.histogram.values.optimal))-1))
		      {
		      #count=1
		      polygon(c(area_susc_values_lda_val_optimal[count:(count+1)],rev(area_susc_values_lda_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1))[count])  
		      }
		    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
		    lines(shape_validation_lda@data[,1],shape_validation_lda@data[,2],col="black")
		    dev.off()
		    }
		  ############
		  }
	

	if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POINTS")
		{
		shape_training_lda<-shape_training
		result_training_lda_shape<-cbind(identification.value,training.table[,2],predict.result.lda$posterior[,2],as.numeric(levels(predict.result.lda$class))[predict.result.lda$class],result.lda.matching.code,ID.bootstrap.model.lda.count,bootstrap.model.lda.probability.mean,bootstrap.model.lda.probability.sd,bootstrap.model.lda.probability.min,bootstrap.model.lda.probability.max,bootstrap.model.lda.probability.sderror,t(bootstrap.model.lda.probability.quantiles),bootstrap.model.lda.prediction.mean,bootstrap.model.lda.prediction.sd,bootstrap.model.lda.prediction.min,bootstrap.model.lda.prediction.max,bootstrap.model.lda.prediction.sderror,t(bootstrap.model.lda.prediction.quantiles))
		colnames(result_training_lda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LDA_SAMP","LDA_PMean","LDA_PSd","LDA_PMin","LDA_PMax","LDA_PSder","LDA_PQ0","LDA_PQ_005","LDA_PQ_025","LDA_PQ05","LDA_PQ_075","LDA_PQ095","LDA_PQ1","LDA_PrMean","LDA_PrSd","LDA_PrMin","LDA_PrMax","LDA_PrSder","LDA_PrQ0","LDA_PrQ005","LDA_PrQ025","LDA_PrQ05","LDA_PrQ075","LDA_PrQ095","LDA_PrQ1")
		###########
		if(enable_probability_optimal_binary_classification==TRUE) 
		  {
		  result_training_lda_shape<-cbind(result_training_lda_shape,as.numeric(predict.result.lda$posterior[,2]>lda.probability.optimal.binary.threshold),result.lda.matching.code.optimal)
		  colnames(result_training_lda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LDA_SAMP","LDA_PMean","LDA_PSd","LDA_PMin","LDA_PMax","LDA_PSder","LDA_PQ0","LDA_PQ_005","LDA_PQ_025","LDA_PQ05","LDA_PQ_075","LDA_PQ095","LDA_PQ1","LDA_PrMean","LDA_PrSd","LDA_PrMin","LDA_PrMax","LDA_PrSder","LDA_PrQ0","LDA_PrQ005","LDA_PrQ025","LDA_PrQ05","LDA_PrQ075","LDA_PrQ095","LDA_PrQ1","OPT_CLASS","OPT_MATCH")
		  }
		#############
		
		shape_training_lda@data <- merge(x=shape_training_lda@data,y=result_training_lda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
		#writeOGR(shape_training_lda,dsn="result_LDA_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
		writeOGR(shape_training_lda,dsn="result_LDA_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
		
		
		# WARNING: The validation does't have the unvertainty estimation: probably this can be daone using the parabolic error function 
		shape_validation_lda<-shape_validation
		result_validation_lda_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.lda.validation$posterior[,2],as.numeric(levels(predict.result.lda.validation$class))[predict.result.lda.validation$class],validation.lda.matching.code)
		colnames(result_validation_lda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
		###########
		if(enable_probability_optimal_binary_classification==TRUE) 
		  {
		  result_validation_lda_shape<-cbind(result_validation_lda_shape,as.numeric(predict.result.lda.validation$posterior[,2]>lda.probability.optimal.binary.threshold),validation.lda.matching.code.optimal)
		  colnames(result_validation_lda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
		  }
		############
		
		shape_validation_lda@data <- merge(x=shape_validation_lda@data,y=result_validation_lda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
		shape_validation_lda@data <- cbind(shape_validation_lda@data,PROB_SDMOD=(coefficients(fit.parabola.probability.lda)*(shape_validation_lda@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.lda)*shape_validation_lda@data$VAL_PROB))
		#writeOGR(shape_validation_lda,dsn="result_LDA_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
		writeOGR(shape_validation_lda,dsn="result_LDA_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
		
		require(raster)
		
		# Plot and export of maps
		# LDA Susceptibility
		#dev.new()
		pdf(file = "result_LDA_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
		layer_gridded<-shape_training_lda
		layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_PROB"])
		gridded(layer_gridded)<-TRUE
		layer_gridded_raster<-raster(layer_gridded)
		#res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
		#zoom(layer_gridded_raster)		
		dev.off()
		
		writeRaster(layer_gridded_raster, filename="result_LDA_Model_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
		###########
		if(enable_probability_optimal_classification==TRUE) 
		  {
		  pdf(file = "result_LDA_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
  	  print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1),breaks=round(lda.breaks.histogram.values.optimal,3)))
		  dev.off()
		  }
		############
		
		
		
		# LDA Model Matching Code
		#dev.new()
		pdf(file = "result_LDA_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
		layer_gridded<-shape_training_lda
		layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_MATCH"])
		gridded(layer_gridded)<-TRUE
		layer_gridded_raster<-raster(layer_gridded)
		#res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
		print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
		legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
		dev.off()

    writeRaster(layer_gridded_raster, filename="result_LDA_Model_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
		
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
      {
      #dev.new()
      pdf(file = "result_LDA_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_lda
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
      print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_LDA_Model_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_LDA_Model_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_lda
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_LDA_Model_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
      }
    ############
    
  	# LDA Model uncertainity
		#dev.new()
		pdf(file = "result_LDA_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
		layer_gridded<-shape_training_lda
		layer_gridded@data<-as.data.frame(layer_gridded@data[,"LDA_PrSd"])
		gridded(layer_gridded)<-TRUE
		layer_gridded_raster<-raster(layer_gridded)
		#res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
		dev.off()
    
		writeRaster(layer_gridded_raster, filename="result_LDA_Model_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
		

    # LDA Validation Susceptibility
		#dev.new()
		pdf(file = "result_LDA_Validation_Susceptibility_Map.pdf")
		layer_gridded<-shape_validation_lda
		layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_PROB"])
		gridded(layer_gridded)<-TRUE
		layer_gridded_raster<-raster(layer_gridded)
		#res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
		#zoom(layer_gridded_raster)		
		dev.off()
		
		writeRaster(layer_gridded_raster, filename="result_LDA_Validation_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
		    
		###########
		if(enable_probability_optimal_classification==TRUE) 
		  {
		  pdf(file = "result_LDA_Validation_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		  print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1),breaks=round(lda.breaks.histogram.values.optimal,3)))
		  dev.off()
	  	}
		############
		
		# LDA Validation Matching Code
		#dev.new()
		pdf(file = "result_LDA_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
		layer_gridded<-shape_validation_lda
		layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_MATCH"])
		gridded(layer_gridded)<-TRUE
		layer_gridded_raster<-raster(layer_gridded)
		#res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
		legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    #zoom(layer_gridded_raster)  	
		dev.off()
    
		writeRaster(layer_gridded_raster, filename="result_LDA_Validation_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)

    ########
		if(enable_probability_optimal_binary_classification==TRUE) 
		{
		  #dev.new()
		  pdf(file = "result_LDA_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		  layer_gridded<-shape_validation_lda
		  layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
		  gridded(layer_gridded)<-TRUE
		  layer_gridded_raster<-raster(layer_gridded)
		  #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		  res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		  print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
		  legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
		  #zoom(layer_gridded_raster)  	
		  dev.off()
		  writeRaster(layer_gridded_raster, filename="result_LDA_Validation_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

		  #dev.new()
		  pdf(file = "result_LDA_Validation_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		  layer_gridded<-shape_validation_lda
		  layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
		  gridded(layer_gridded)<-TRUE
		  layer_gridded_raster<-raster(layer_gridded)
		  #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		  res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		  print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
		  legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
		  dev.off()
		  
		  writeRaster(layer_gridded_raster, filename="result_LDA_Validation_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
		  
      
		  }
		#######
    
		# LDA Model uncertainity
		#dev.new()
		pdf(file = "result_LDA_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
		layer_gridded<-shape_validation_lda
		layer_gridded@data<-as.data.frame(layer_gridded@data[,"PROB_SDMOD"])
		gridded(layer_gridded)<-TRUE
		layer_gridded_raster<-raster(layer_gridded)
		#res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
		res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
		print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
		#zoom(layer_gridded_raster)		
		dev.off()
    
		writeRaster(layer_gridded_raster, filename="result_LDA_Validation_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
	
		### Success & prediction rate curve
		susc_values_lda<-c(1,0.8,0.55,0.45,0.2,0)
		
		# ordering data for susceptibility
		ind_col_succ_pred_rate<-which(colnames(shape_training_lda@data) %in% c("MOD_GROUP","MOD_PROB"))
		shape_training_lda@data<-cbind(PIXEL_AREA=rep(configuration.spatial.data.table[c(8)]^2,dim(shape_training_lda@data)[1]),shape_training_lda@data[order(shape_training_lda@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
		shape_training_lda@data$MOD_GROUP<-shape_training_lda@data$MOD_GROUP*configuration.spatial.data.table[c(8)]^2
		shape_training_lda@data[,1]<-cumsum(shape_training_lda@data[,1])/sum(shape_training_lda@data[,1])*100
		shape_training_lda@data[,2]<-cumsum(shape_training_lda@data[,2])/sum(shape_training_lda@data[,2])*100
		
		ind_pro_fun_mod<-approxfun(shape_training_lda@data[,3],shape_training_lda@data[,1])
		area_susc_values_lda_mod<-c(0,ind_pro_fun_mod(susc_values_lda[-c(1,length(susc_values_lda))]),100)
		
		#dev.new()
		pdf(file = "result_LDA_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		plot(0,0,col="transparent",main="LDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
		for (count in 1:(length(susc_values_lda)-1))
			{
			#count=1
			polygon(c(area_susc_values_lda_mod[count:(count+1)],rev(area_susc_values_lda_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
			}
		polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
		lines(shape_training_lda@data[,1],shape_training_lda@data[,2],col="black")
		dev.off()
		
		###########
		if(enable_probability_optimal_classification==TRUE) 
  		{
		  area_susc_values_lda_mod_optimal<-c(0,ind_pro_fun_mod(rev(lda.breaks.histogram.values.optimal)[-c(1,length(lda.breaks.histogram.values.optimal))]),100)
		  pdf(file = "result_LDA_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
		  plot(0,0,col="transparent",main="LDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
		  for (count in 1:(length(rev(lda.breaks.histogram.values.optimal))-1))
		    {
		    #count=1
		    polygon(c(area_susc_values_lda_mod_optimal[count:(count+1)],rev(area_susc_values_lda_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1))[count])  
		    }
		  polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
		  lines(shape_training_lda@data[,1],shape_training_lda@data[,2],col="black")
		  dev.off()
	  	}
		############
		
		
		ind_col_succ_pred_rate<-which(colnames(shape_validation_lda@data) %in% c("VAL_GROUP","VAL_PROB"))
		shape_validation_lda@data<-cbind(PIXEL_AREA=rep(as.numeric(configuration.spatial.data.table[c(8)])^2,dim(shape_validation_lda@data)[1]),shape_validation_lda@data[order(shape_validation_lda@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
		shape_validation_lda@data$VAL_GROUP<-shape_validation_lda@data$VAL_GROUP*as.numeric(configuration.spatial.data.table[c(8)])^2
		shape_validation_lda@data[,1]<-cumsum(shape_validation_lda@data[,1])/sum(shape_validation_lda@data[,1])*100
		shape_validation_lda@data[,2]<-cumsum(shape_validation_lda@data[,2])/sum(shape_validation_lda@data[,2])*100
	
		ind_pro_fun_val<-approxfun(shape_validation_lda@data[,3],shape_validation_lda@data[,1])
		area_susc_values_lda_val<-c(0,ind_pro_fun_val(susc_values_lda[-c(1,length(susc_values_lda))]),100)
		
		#dev.new()
		pdf(file = "result_LDA_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		plot(0,0,col="transparent",main="LDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
		for (count in 1:(length(susc_values_lda)-1))
			{
			#count=1
			polygon(c(area_susc_values_lda_val[count:(count+1)],rev(area_susc_values_lda_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
			}
		polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
		lines(shape_validation_lda@data[,1],shape_validation_lda@data[,2],col="black")
		dev.off()
		
		###########
		
		if(enable_probability_optimal_classification==TRUE) 
	  	{
		  area_susc_values_lda_val_optimal<-c(0,ind_pro_fun_val(rev(lda.breaks.histogram.values.optimal)[-c(1,length(lda.breaks.histogram.values.optimal))]),100)
		  pdf(file = "result_LDA_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		  plot(0,0,col="transparent",main="LDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
		  for (count in 1:(length(rev(lda.breaks.histogram.values.optimal))-1))
		    {
		    #count=1
		    polygon(c(area_susc_values_lda_val_optimal[count:(count+1)],rev(area_susc_values_lda_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lda.breaks.histogram.values.optimal)-1))[count])  
		    }
		  polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
		  lines(shape_validation_lda@data[,1],shape_validation_lda@data[,2],col="black")
		  dev.off()
		  }
		############
		}
		
	}


#----------------- QUADRATIC DISCRIMINANT ANALISYS ------------------#

if(model.run.matrix[2] == "YES")
  {
  # Changing of Dummy Explanatory Variables in Numeric variable
  if (analysis.parameter.matrix[2] == "DUM")
  {
    print("The Quadratic Discriminant Analsysis (QDA) will be performed using dummy variables, but a random variation in these variables will be introduced")
    for (count.variables in 1:dim(explanatory.variables)[2])
    {
      #print(range(explanatory.variables[,count.variables]))
      if (min(explanatory.variables[,count.variables])==0 & max(explanatory.variables[,count.variables])==1)  
      { 
        set.seed(seed.value)
        explanatory.variables[,count.variables]<-explanatory.variables[,count.variables]+runif(dim(explanatory.variables)[1],-0.1,0.1)  
        validation.explanatory.variables[,count.variables]<-explanatory.variables[,count.variables]
      }   
    }
    print("Performing 0 value correction")
    indexbootstrapcosntant_training<-which(explanatory.variables==0,arr.ind=TRUE)
    explanatory.variables[indexbootstrapcosntant_training]<-runif(length(indexbootstrapcosntant_training)/2, min = 0.00001, max = 0.01)
    indexbootstrapcosntant_validation<-which(validation.explanatory.variables==0,arr.ind=TRUE)
    validation.explanatory.variables[indexbootstrapcosntant_validation]<-runif(length(indexbootstrapcosntant_validation)/2, min = 0.00001, max = 0.01)
  }
  
  
  if (analysis.parameter.matrix[2] == "SEL")
  {
    print("The Quadratic Discriminant Analsysis (QDA) will be performed excluding dummy variables")
    index.variables.dummy<-NULL
    for (count.variables in 1:dim(explanatory.variables)[2])
    {
      print(range(explanatory.variables[,count.variables]))
      if (min(explanatory.variables[,count.variables])==0 & max(explanatory.variables[,count.variables])==1)  
      { 
        index.variables.dummy<-c(index.variables.dummy,count.variables)
      }      
    }
    if(length(index.variables.dummy)>0) {explanatory.variables<-explanatory.variables[,-index.variables.dummy]}
    if(length(index.variables.dummy)>0) {validation.explanatory.variables<-validation.explanatory.variables[,-index.variables.dummy]}
    #str(explanatory.variables)
  }
  
  
  if (class(try(qda(explanatory.variables, grouping.variable, method="moment")))=="try-error")  
  { 
    #qda(explanatory.variables, grouping.variable, method="moment")
    write.table("Quadratic Discriminant Analysis was not completed",file="Error_QDA_Analysis.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="Error_QDA_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("Error LOG",file="Error_QDA_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind("Message",rev(1:length(as.vector(.Traceback)))," ->",as.vector(.Traceback)),file="Error_QDA_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } 
  
  ##### Quadratic Discriminant Analisys using data.frame
  result.qda<-NULL
  result.qda<-qda(explanatory.variables, grouping.variable, method="mle")
  
  # Result Predicted
  predict.result.qda<-predict(result.qda)
  str(predict.result.qda)
  
  cross.classification.qda<-table(grouping.variable,predict.result.qda$class,dnn=c("Observed","Predicted"))
  rownames(cross.classification.qda)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.qda)<-list("No Landslide","Landslide") # Predicted    
  str(cross.classification.qda)
  
  # Assignation of a matching code between observed and predicted values
  result.qda.matching.code<-paste(grouping.variable,as.numeric(levels(predict.result.qda$class))[predict.result.qda$class],sep="")
  result.qda.matching.code<-gsub("00","1",result.qda.matching.code)
  result.qda.matching.code<-gsub("01","2",result.qda.matching.code)
  result.qda.matching.code<-gsub("10","3",result.qda.matching.code)
  result.qda.matching.code<-gsub("11","4",result.qda.matching.code)
  result.qda.matching.code<-as.numeric(result.qda.matching.code)
  
  #Elaboration of Coefficient of association for contingency table 
  #load package (vcd)  
  library(vcd)
  
  #help(package=vcd)         
  contingency.table.qda<-table2d_summary(cross.classification.qda)
  test.table.qda<-assocstats(cross.classification.qda)
  
  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.qda/sum(cross.classification.qda)*100,2), std="margin", main="QUADRATIC DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.qda, std="margin", main="QUADRATIC DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  }
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #A ROC curve plots the false alarm rate against the hit rate
  #for a probablistic forecast for a range of thresholds. 
  
  #load package (verification)  
  library(verification)
  
  #verify function
  #Based on the type of inputs, this function calculates a range of verification statistics and skill scores.
  #Additionally, it creates a verify class object that can be further analyzed.
  
  
  ##### ROC PLOT OBS - POSTERIOR PROBABILITY ASSOCIATED TO 1                                                                                 
  ## 1st method
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(training.table[,2],predict.result.qda$posterior[,2],main = "ROC PLOT: QUADRATIC DISCRIMINANT ANALYSIS MODEL", binormal = TRUE, plot = "both")
  #}
  
  # 2nd method using verify function
  verification.results.qda<-verify(training.table[,2],predict.result.qda$posterior[,2], frcst.type="prob", obs.type="binary")
  #summary(verification.results.qda)
  #str(verification.results.qda)
  
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(verification.results.qda, main = "ROC PLOT: QUADRATIC DISCRIMINANT ANALYSIS MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #}
  area.under.roc.curve.qda<-roc.area(training.table[,2],predict.result.qda$posterior[,2])
  
  ## showing confidence intervals.  MAY BE SLOW
  
  if (cross.classification.qda[1,2]==0 | cross.classification.qda[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.results.qda, main = "ROC PLOT: QUADRATIC DISCRIMINANT ANALYSIS MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[2] , alpha = 0.05, extra=TRUE, legend=TRUE, show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.qda$A,2),";  Sample size = ",area.under.roc.curve.qda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[2], sep=""), side=3, col="red", cex=0.8)
    ## Histogram of posterior probability
    dev.new()
    hist(predict.result.qda$posterior[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Quadratic Disciminant Analysis susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  }
  pdf(file = "result_QDA_Histogram.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  hist(predict.result.qda$posterior[,2], breaks=breaks.histogram.values,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Quadratic Disciminant Analysis susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  dev.off()
  
  
  # EXPORT OF PLOT FOR QDA MODEL
  pdf(file = "result_QDA_FourfoldPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.qda/sum(cross.classification.qda)*100,2), std="margin", main="QUADRATIC DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  #fourfold(cross.classification.qda, std="margin", main="QUADRATIC DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  dev.off()
  
  #pdf(file = "result_QDA_ROCPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.results.qda, main = "ROC PLOT: QUADRATIC DISCRIMINANT ANALYSIS MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #dev.off()
  
  pdf(file = "result_QDA_ROCPlot_bootstrap.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.results.qda, main = "ROC PLOT: QUADRATIC DISCRIMINANT ANALYSIS MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[2] , alpha = 0.05, extra=TRUE, legend=TRUE, show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.qda$A,2),";  Sample size = ",area.under.roc.curve.qda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[2], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  ## BOOTSTRAP PROCEDURE FOR THE ESTIMATION OF MODEL PREDICTION VARIABILITY
  if(bootstrap.model.variability[2] == "YES")
  {
    bootstrap.sample.model.qda<-bootstrap.sample.model[2]
    
    matrix.bootstrap.model.qda<-matrix(data=NA, nrow=dim(explanatory.variables)[1], ncol=(bootstrap.sample.model.qda*3)+1)
    colnames(matrix.bootstrap.model.qda)<-rep("na",(bootstrap.sample.model.qda*3)+1)
    matrix.bootstrap.model.qda[,1]<-identification.value
    colnames(matrix.bootstrap.model.qda)[1]<-"ID"
    name.sel.run<-paste(rep("ID_Selection_Run",bootstrap.sample.model.qda),1:bootstrap.sample.model.qda,sep="_")
    colnames(matrix.bootstrap.model.qda)[seq(2,(bootstrap.sample.model.qda*3)-1,3)]<-name.sel.run
    name.prob.run<-paste(rep("Probability_Run",bootstrap.sample.model.qda),1:bootstrap.sample.model.qda,sep="_")
    colnames(matrix.bootstrap.model.qda)[seq(3,(bootstrap.sample.model.qda*3),3)]<-name.prob.run
    name.pred.run<-paste(rep("Prediction_Run",bootstrap.sample.model.qda),1:bootstrap.sample.model.qda,sep="_")
    colnames(matrix.bootstrap.model.qda)[seq(4,(bootstrap.sample.model.qda*3)+1,3)]<-name.pred.run
    
    selection.index<-NULL
    library(MASS)
    #Bootstrap procedure
    for (count.boot in 1:bootstrap.sample.model.qda)
    {
      selection.index<-sample(1:dim(explanatory.variables)[1], replace=TRUE, prob=NULL)
      matrix.bootstrap.model.qda[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
      explanatory.variables.bootstrap.model.qda<-explanatory.variables[selection.index,]
      grouping.variable.bootstrap.model.qda<-as.factor(training.table[selection.index,2])
      #result.bootstrap.model.qda<-qda(explanatory.variables.bootstrap.model.qda, grouping.variable.bootstrap.model.qda, tol=0.001, method="moment")
      while(inherits(try(result.bootstrap.model.qda<-qda(explanatory.variables.bootstrap.model.qda, grouping.variable.bootstrap.model.qda, tol=0.001, method="moment"),silent=TRUE),what="try-error"))
      {
        print(paste("Count boot: ",count.boot," - Boostrap while resampling",sep=""))
        selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
        matrix.bootstrap.model.qda[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
        explanatory.variables.bootstrap.model.qda<-training.table[selection.index,3:dim(training.table)[2]]
        if(bootstrap_constant_correction==TRUE)
        {
          print("Performing bootstrap 0 value correction")
          indexbootstrapcosntant<-which(explanatory.variables.bootstrap.model.qda==0,arr.ind=TRUE)
          explanatory.variables.bootstrap.model.qda[indexbootstrapcosntant]<-runif(length(indexbootstrapcosntant)/2, min = 0.00001, max = 0.01)
        }
        grouping.variable.bootstrap.model.qda<-as.factor(training.table[selection.index,2])
      }
      matrix.bootstrap.model.qda[,(count.boot*3)+1]<-predict(result.bootstrap.model.qda,newdata=explanatory.variables)$posterior[,2]
      matrix.bootstrap.model.qda[as.numeric(names(table(selection.index))),(count.boot*3)]<-matrix.bootstrap.model.qda[as.numeric(names(table(selection.index))),(count.boot*3)+1]  	
    }
    # Export of bootstrap sample
    write.table(matrix.bootstrap.model.qda,file="result_QDA_BootstrapSamples.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    
    ID.bootstrap.model.qda.count<-numeric(length=dim(training.table)[1])
    #Probability (selected values)
    bootstrap.model.qda.probability.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.probability.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.probability.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.probability.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.probability.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.probability.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #Prediction (all values)
    bootstrap.model.qda.prediction.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.prediction.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.prediction.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.prediction.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.prediction.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.qda.prediction.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #    for (count.row.variability in 1:dim(training.table)[1])
    #        {
    #        # Statistics on boostrapped probability
    #        ID.bootstrap.model.qda.count[count.row.variability]<-length(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(2,(bootstrap.sample.model.qda*3)-1,3)]))
    #        bootstrap.model.qda.probability.mean[count.row.variability]<-mean(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(3,(bootstrap.sample.model.qda*3),3)]))
    #        bootstrap.model.qda.probability.sd[count.row.variability]<-sd(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(3,(bootstrap.sample.model.qda*3),3)]))
    #        bootstrap.model.qda.probability.min[count.row.variability]<-min(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(3,(bootstrap.sample.model.qda*3),3)]))
    #        bootstrap.model.qda.probability.max[count.row.variability]<-max(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(3,(bootstrap.sample.model.qda*3),3)]))
    #        bootstrap.model.qda.probability.sderror[count.row.variability]<-bootstrap.model.qda.probability.sd[count.row.variability]/ID.bootstrap.model.qda.count[count.row.variability]
    #        bootstrap.model.qda.probability.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(3,(bootstrap.sample.model.qda*3),3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        # Statistics on boostrapped prediction
    #        bootstrap.model.qda.prediction.mean[count.row.variability]<-mean(matrix.bootstrap.model.qda[count.row.variability,seq(4,(bootstrap.sample.model.qda*3)+1,3)])
    #        bootstrap.model.qda.prediction.sd[count.row.variability]<-sd(matrix.bootstrap.model.qda[count.row.variability,seq(4,(bootstrap.sample.model.qda*3)+1,3)])
    #        bootstrap.model.qda.prediction.min[count.row.variability]<-min(matrix.bootstrap.model.qda[count.row.variability,seq(4,(bootstrap.sample.model.qda*3)+1,3)])
    #        bootstrap.model.qda.prediction.max[count.row.variability]<-max(matrix.bootstrap.model.qda[count.row.variability,seq(4,(bootstrap.sample.model.qda*3)+1,3)])
    #        bootstrap.model.qda.prediction.sderror[count.row.variability]<-bootstrap.model.qda.prediction.sd[count.row.variability]/bootstrap.sample.model.qda
    #        bootstrap.model.qda.prediction.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.qda[count.row.variability,seq(4,(bootstrap.sample.model.qda*3)+1,3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        }
    
    
    fun_length<-function(x) {length_var<-length(x[which(is.finite(x))]); return(length_var)}
    fun_quantile<-function(x) {quantile_var<-t(quantile(x[which(is.finite(x))],probs=c(0,0.05,0.25,0.5,0.75,0.95,1))); return(quantile_var)}
    ID.bootstrap.model.qda.count<-apply(matrix.bootstrap.model.qda[,grep("ID_Selection",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=fun_length)
    bootstrap.model.qda.probability.mean<-apply(matrix.bootstrap.model.qda[,grep("Probability",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.qda.probability.sd<-apply(matrix.bootstrap.model.qda[,grep("Probability",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.qda.probability.min<-apply(matrix.bootstrap.model.qda[,grep("Probability",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.qda.probability.max<-apply(matrix.bootstrap.model.qda[,grep("Probability",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.qda.probability.sderror<-bootstrap.model.qda.probability.sd/bootstrap.sample.model.qda
    bootstrap.model.qda.probability.quantiles<-apply(matrix.bootstrap.model.qda[,grep("Probability",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=fun_quantile)
    bootstrap.model.qda.prediction.mean<-apply(matrix.bootstrap.model.qda[,grep("Prediction",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.qda.prediction.sd<-apply(matrix.bootstrap.model.qda[,grep("Prediction",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.qda.prediction.min<-apply(matrix.bootstrap.model.qda[,grep("Prediction",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.qda.prediction.max<-apply(matrix.bootstrap.model.qda[,grep("Prediction",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.qda.prediction.sderror<-bootstrap.model.qda.prediction.sd/bootstrap.sample.model.qda
    bootstrap.model.qda.prediction.quantiles<-apply(matrix.bootstrap.model.qda[,grep("Prediction",colnames(matrix.bootstrap.model.qda))],MARGIN=1,FUN=fun_quantile)
    
    
    # Export of bootstrap sample statistics
    write.table(cbind("ID","QDA_NumberSelectedSamples","QDA_Probability_Mean","QDA_Probability_Sd","QDA_Probability_Min","QDA_Probability_Max","QDA_Probability_Sderror","QDA_Probability_Quantiles_0","QDA_Probability_Quantiles_0.05","QDA_Probability_Quantiles_0.25","QDA_Probability_Quantiles_0.5","QDA_Probability_Quantiles_0.75","QDA_Probability_Quantiles_0.95","QDA_Probability_Quantiles_1","QDA_Prediction_Mean","QDA_Prediction_Sd","QDA_Prediction_Min","QDA_Prediction_Max","QDA_Prediction_Sderror","QDA_Prediction_Quantiles_0","QDA_Prediction_Quantiles_0.05","QDA_Prediction_Quantiles_0.25","QDA_Prediction_Quantiles_0.5","QDA_Prediction_Quantiles_0.75","QDA_Prediction_Quantiles_0.95","QDA_Prediction_Quantiles_1"),file="result_QDA_BootstrapStatistics.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind(identification.value,ID.bootstrap.model.qda.count,bootstrap.model.qda.probability.mean,bootstrap.model.qda.probability.sd,bootstrap.model.qda.probability.min,bootstrap.model.qda.probability.max,bootstrap.model.qda.probability.sderror,t(bootstrap.model.qda.probability.quantiles),bootstrap.model.qda.prediction.mean,bootstrap.model.qda.prediction.sd,bootstrap.model.qda.prediction.min,bootstrap.model.qda.prediction.max,bootstrap.model.qda.prediction.sderror,t(bootstrap.model.qda.prediction.quantiles)),file="result_QDA_BootstrapStatistics.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    
    if (enable_screen_plotting==TRUE)
    {
      #dev.new()
      #double.sd.histogram.variability<-hist(bootstrap.model.qda.probability.sd*2,breaks=seq(0,1,0.05),labels=TRUE)
      #plot(double.sd.histogram.variability$counts, seq(0,0.95,0.05), type="S",ylim=c(0,1), labels=TRUE)
      dev.new()
      plot(bootstrap.model.qda.probability.mean,bootstrap.model.qda.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="QDA BOOTSTRAP: Mean Probability vs Mean Prediction")
      abline(a=0,b=1,col="red",lty=1,lwd=1)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.qda,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    }
    
    pdf(file = "result_QDA_BootstrapMeansComparison.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(bootstrap.model.qda.probability.mean,bootstrap.model.qda.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="QDA BOOTSTRAP: Mean Probability vs Mean Prediction")
    abline(a=0,b=1,col="red",lty=1,lwd=1)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.qda,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    dev.off()
    
    
    # BOOTSTRAPPED PROBABILITY - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.probability.qda<-cbind(bootstrap.model.qda.probability.mean,2*bootstrap.model.qda.probability.sd)
    parabola.probability.qda<-na.omit(parabola.probability.qda[order(parabola.probability.qda[,1]),])
    colnames(parabola.probability.qda)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.probability.qda <- nls(parabola.probability.qda[,"ordinate"] ~ coeff.a*(parabola.probability.qda[,"abscissa"]^2) + (-1)*coeff.a*parabola.probability.qda[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.probability.qda<-predict(fit.parabola.probability.qda)
    #coef(fit.parabola.probability.qda)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.probability.qda[,"abscissa"],parabola.probability.qda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="QDA Model Probability Variability (Bootstrap)")
      lines(parabola.probability.qda[,"abscissa"],value.parabola.probability.qda,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.qda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.qda),3),coeff.b= -round(coef(fit.parabola.probability.qda),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_QDA_BootstrapProbabilityVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.probability.qda[,"abscissa"],parabola.probability.qda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="QDA Model Probability Variability (Bootstrap)")
    lines(parabola.probability.qda[,"abscissa"],value.parabola.probability.qda,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.qda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.qda),3),coeff.b= -round(coef(fit.parabola.probability.qda),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
    
    # BOOTSTRAPPED PREDICTION - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.prediction.qda<-cbind(bootstrap.model.qda.prediction.mean,2*bootstrap.model.qda.prediction.sd)
    parabola.prediction.qda<-parabola.prediction.qda[order(parabola.prediction.qda[,1]),]
    colnames(parabola.prediction.qda)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.prediction.qda <- nls(parabola.prediction.qda[,"ordinate"] ~ coeff.a*(parabola.prediction.qda[,"abscissa"]^2) + (-1)*coeff.a*parabola.prediction.qda[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.prediction.qda<-predict(fit.parabola.prediction.qda)
    #coef(fit.parabola.prediction.qda)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.prediction.qda[,"abscissa"],parabola.prediction.qda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="QDA Model Prediction Variability (Bootstrap)")
      lines(parabola.prediction.qda[,"abscissa"],value.parabola.prediction.qda,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.qda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.qda),3),coeff.b= -round(coef(fit.parabola.prediction.qda),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_QDA_BootstrapPredictionVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.prediction.qda[,"abscissa"],parabola.prediction.qda[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="QDA Model Prediction Variability (Bootstrap)")
    lines(parabola.prediction.qda[,"abscissa"],value.parabola.prediction.qda,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.qda,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.qda),3),coeff.b= -round(coef(fit.parabola.prediction.qda),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
  }
  
  ## Sensitivity, Specificity, Cohens kappa plot
  dev.new()
  roc.plot.qda.series<-roc.plot(verification.results.qda,binormal=bin,plot=FALSE,show.thres=plot_thres,thresholds=threshold_series)
  dev.off()
  #str(roc.plot.qda.series)
  #roc.plot.qda.series$plot.data
  #str(roc.plot.qda.series$plot.data)
  
  ###########################################
  # min(abs(TPR - (1-FPR)))
  if(enable_probability_optimal_binary_classification==FALSE)
  {
    dev.new()
    roc.plot.qda.series<-roc.plot(verification.results.qda,binormal=bin,plot=FALSE,show.thres=FALSE,thresholds=threshold_series)
    dev.off()
  } else
  {
    dev.new()
    roc.plot.qda.series<-roc.plot(verification.results.qda,binormal=bin,plot=FALSE,show.thres=FALSE)
    dev.off()  
  }
  
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    qda.probability.classification.optimal<-data.frame(prob_thres=roc.plot.qda.series$plot.data[,1,1],tpr=roc.plot.qda.series$plot.data[,2,1],fpr=roc.plot.qda.series$plot.data[,3,1],tnr=(1-roc.plot.qda.series$plot.data[,3,1]),diff_abs_tpr_tnr=abs(roc.plot.qda.series$plot.data[,2,1]-(1-roc.plot.qda.series$plot.data[,3,1])),optimal_sel=NA,breaks_sel=NA)
    index.qda.filter<-which(qda.probability.classification.optimal$prob_thres>0 & qda.probability.classification.optimal$prob_thres<1) # removing strnge thresh values
    qda.probability.classification.optimal<-rbind(c(0,1,1,0,1,NA,NA),qda.probability.classification.optimal[index.qda.filter,],c(1,0,0,1,1,NA,NA))
    qda.optimal.index<-which(qda.probability.classification.optimal$diff_abs_tpr_tnr==min(qda.probability.classification.optimal$diff_abs_tpr_tnr))
    qda.probability.classification.optimal$optimal_sel[qda.optimal.index]<-TRUE
    qda.probability.optimal.binary.threshold<-qda.probability.classification.optimal$prob_thres[qda.optimal.index]
    
    ### Generating the optimal fourfold plot
    cross.classification.qda.optimal<-table(grouping.variable,predict.result.qda$posterior[,2]>qda.probability.optimal.binary.threshold,dnn=c("Observed","Predicted"))
    rownames(cross.classification.qda.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.qda.optimal)<-list("No Landslide","Landslide") # Predicted    
    str(cross.classification.qda.optimal)
    # Assignation of a matching code between observed and predicted values
    result.qda.matching.code.optimal<-paste(grouping.variable,as.numeric(predict.result.qda$posterior[,2]>qda.probability.optimal.binary.threshold),sep="")
    result.qda.matching.code.optimal<-gsub("00","1",result.qda.matching.code.optimal)
    result.qda.matching.code.optimal<-gsub("01","2",result.qda.matching.code.optimal)
    result.qda.matching.code.optimal<-gsub("10","3",result.qda.matching.code.optimal)
    result.qda.matching.code.optimal<-gsub("11","4",result.qda.matching.code.optimal)
    result.qda.matching.code.optimal<-as.numeric(result.qda.matching.code.optimal)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.qda.optimal/sum(cross.classification.qda.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    }
    # EXPORT OF PLOT FOR QDA MODEL
    pdf(file = "result_QDA_FourfoldPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.qda.optimal/sum(cross.classification.qda.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.qda.optimal, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    dev.off()
    
    ### Optimal susceptibility classes identification 
    if(enable_probability_optimal_classification==TRUE)
    {
      qda.unexplained.errors<-round(1-(qda.probability.classification.optimal$tpr[qda.optimal.index]+qda.probability.classification.optimal$tpr[qda.optimal.index])/2,2)
      
      if(type_probability_optimal_classification=="proportional")
      {
        qda.unexplained.errors.partition<-1-(qda.unexplained.errors/((length(breaks.histogram.values)/2)))*(1:((length(breaks.histogram.values)/2)))
        
        for(count_part in 1:(length(qda.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-qda.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((qda.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          qda.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((qda.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          qda.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        qda.breaks.histogram.values.optimal<-c(0,qda.probability.classification.optimal$prob_thres[which(qda.probability.classification.optimal$breaks_sel==TRUE)],1)
        #qda.probability.classification.optimal[which(qda.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if(type_probability_optimal_classification=="fixed")
        {
        step.qda.unexplained.fixed<-0.1
        if(qda.unexplained.errors<=0.1) step.qda.unexplained.fixed<-0.05
        if(qda.unexplained.errors<=0.05) step.qda.unexplained.fixed<-0.025
        qda.unexplained.errors.partition<-seq(step.qda.unexplained.fixed,1-step.qda.unexplained.fixed,step.qda.unexplained.fixed)[seq(step.qda.unexplained.fixed,1-step.qda.unexplained.fixed,step.qda.unexplained.fixed)>(1-qda.unexplained.errors)]
                
        for(count_part in 1:(length(qda.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-qda.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((qda.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          qda.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((qda.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          qda.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        qda.breaks.histogram.values.optimal<-c(0,qda.probability.classification.optimal$prob_thres[which(qda.probability.classification.optimal$breaks_sel==TRUE)],1)
        #qda.probability.classification.optimal[which(qda.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if (enable_screen_plotting==TRUE)
      {
        dev.new()
        hist(predict.result.qda$posterior[,2], breaks=qda.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal QDA susceptibility", col=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1))
        plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("QDA OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
        mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
        mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
        for (count in 1:(length(qda.breaks.histogram.values.optimal)-1))
        {
          #count=1
          polygon(c(qda.breaks.histogram.values.optimal[count:(count+1)],rev(qda.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1)[count])  
        }
        polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
        lines(qda.probability.classification.optimal$prob_thres,qda.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
        lines(qda.probability.classification.optimal$prob_thres,qda.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
        index_points_plot<-c(1,which(qda.probability.classification.optimal$breaks_sel==TRUE),dim(qda.probability.classification.optimal)[1])
        points(qda.probability.classification.optimal[index_points_plot[1:floor(length(qda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
        text(qda.probability.classification.optimal[index_points_plot[1:floor(length(qda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(qda.probability.classification.optimal[index_points_plot[1:floor(length(qda.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(qda.breaks.histogram.values.optimal)/2)-1)))
        points(qda.probability.classification.optimal[index_points_plot[ceiling(length(qda.breaks.histogram.values.optimal)/2):length(qda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
        text(qda.probability.classification.optimal[index_points_plot[ceiling(length(qda.breaks.histogram.values.optimal)/2):length(qda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(qda.probability.classification.optimal[index_points_plot[ceiling(length(qda.breaks.histogram.values.optimal)/2):length(qda.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(qda.breaks.histogram.values.optimal)-ceiling(length(qda.breaks.histogram.values.optimal)/2)),3))
      }
      
      pdf(file = "result_QDA_Histogram_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      hist(predict.result.qda$posterior[,2], breaks=qda.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal QDA susceptibility", col=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1))
      dev.off()
      
      pdf(file = "result_QDA_ModelEvaluationPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("QDA OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
      mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
      mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
      for (count in 1:(length(qda.breaks.histogram.values.optimal)-1))
      {
        #count=1
        polygon(c(qda.breaks.histogram.values.optimal[count:(count+1)],rev(qda.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1)[count])  
      }
      polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
      lines(qda.probability.classification.optimal$prob_thres,qda.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
      lines(qda.probability.classification.optimal$prob_thres,qda.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
      index_points_plot<-c(1,which(qda.probability.classification.optimal$breaks_sel==TRUE),dim(qda.probability.classification.optimal)[1])
      points(qda.probability.classification.optimal[index_points_plot[1:floor(length(qda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(qda.probability.classification.optimal[index_points_plot[1:floor(length(qda.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(qda.probability.classification.optimal[index_points_plot[1:floor(length(qda.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(qda.breaks.histogram.values.optimal)/2)-1)))
      points(qda.probability.classification.optimal[index_points_plot[1+ceiling(length(qda.breaks.histogram.values.optimal)/2):length(qda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
      text(qda.probability.classification.optimal[index_points_plot[1+ceiling(length(qda.breaks.histogram.values.optimal)/2):length(qda.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(qda.probability.classification.optimal[index_points_plot[1+ceiling(length(qda.breaks.histogram.values.optimal)/2):length(qda.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(qda.breaks.histogram.values.optimal)-1-ceiling(length(qda.breaks.histogram.values.optimal)/2)),3))
      points(qda.probability.classification.optimal[qda.optimal.index,c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(qda.probability.classification.optimal[qda.optimal.index,c("prob_thres","tpr")],labels=with(round(qda.probability.classification.optimal[qda.optimal.index,c("prob_thres","tpr")],3),paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(4))
      dev.off()
    }
  }
  
  ###########################################
  
  contingency.table.matrix.qda<-matrix(nrow=dim(roc.plot.qda.series$plot.data)[1],ncol=8)
  colnames(contingency.table.matrix.qda)<-c("Threshold","TP","TN","FP","FN","TPR","FPR","COHEN_KAPPA")
  contingency.table.matrix.qda[,1]<-roc.plot.qda.series$plot.data[,1,1]
  contingency.table.matrix.qda[,6]<-roc.plot.qda.series$plot.data[,2,1]
  contingency.table.matrix.qda[,7]<-roc.plot.qda.series$plot.data[,3,1]
  values.observed<-training.table[,2]
  values.predicted<-predict.result.qda$posterior[,2]
  for (count.threshold.series in 1:dim(roc.plot.qda.series$plot.data)[1])
  {
    value.threshold<-contingency.table.matrix.qda[count.threshold.series,1]
    values.probability.reclassified<-NULL
    values.probability.reclassified<-as.numeric(values.predicted>value.threshold) 
    #sum(values.probability.reclassified-round(values.predicted)) # Check sum: It has to be 0 if threshold is equal to 1
    series.pasted<-paste(values.observed,values.probability.reclassified,sep="")
    series.pasted<-gsub("00","1",series.pasted)
    series.pasted<-gsub("01","2",series.pasted)
    series.pasted<-gsub("10","3",series.pasted)
    series.pasted<-gsub("11","4",series.pasted)
    series.pasted<-as.numeric(series.pasted)
    TP<-as.numeric(sum(series.pasted>=4)) # True Positive
    FN<-as.numeric(sum(series.pasted>=3 & series.pasted<4)) # False Negative
    FP<-as.numeric(sum(series.pasted>=2 & series.pasted<3)) # False Positive
    TN<-as.numeric(sum(series.pasted>=1 & series.pasted<2)) # True Negative              
    #TPR<-TP/(TP+FN) # Hit Rate or True Positive Rate or Sensitivity - Assigned before the for cicle using rocplot data
    #FPR<-FP/(FP+TN) # False Alarm Rate or False Positive Rate or 1-Specificity
    # Cohen's Kappa = (agreement-chance)/(1-chance)  where agreement=(TP+TN)/(TP+TN+FP+FN) and chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    agreement=(TP+TN)/(TP+TN+FP+FN)
    chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    cohen.kappa.value<-(agreement-chance)/(1-chance)
    #Other
    #library(vcd)
    #cohen.kappa.value<-Kappa(cross.classification.table)
    contingency.table.matrix.qda[count.threshold.series,2]<-TP
    contingency.table.matrix.qda[count.threshold.series,3]<-TN
    contingency.table.matrix.qda[count.threshold.series,4]<-FP
    contingency.table.matrix.qda[count.threshold.series,5]<-FN
    contingency.table.matrix.qda[count.threshold.series,8]<-cohen.kappa.value
  }
  
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(roc.plot.qda.series$plot.data[,1,1],roc.plot.qda.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="QDA MODEL EVALUATION PLOT")
  lines(roc.plot.qda.series$plot.data[,1,1],1-roc.plot.qda.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.qda.series$plot.data[,1,1], contingency.table.matrix.qda[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  }
  pdf(file = "result_QDA_ModelEvaluationPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(roc.plot.qda.series$plot.data[,1,1],roc.plot.qda.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="QDA MODEL EVALUATION PLOT")
  lines(roc.plot.qda.series$plot.data[,1,1],1-roc.plot.qda.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.qda.series$plot.data[,1,1], contingency.table.matrix.qda[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  dev.off()

  
  
  
  
  ## VALIDATION OF QDA MODEL (Matching QDA posterior probability results and validation grouping variable)
  
  # Result Predicted
  predict.result.qda.validation<-predict(result.qda,newdata=validation.explanatory.variables)
  #str(predict.result.qda.validation)
  #summary(predict.result.qda.validation)
  cross.classification.validation.qda<-table(validation.grouping.variable,predict.result.qda.validation$class,dnn=c("Observed","Predicted"))
  rownames(cross.classification.validation.qda)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.validation.qda)<-list("No Landslide","Landslide") # Predicted
  #str(cross.classification.validation.qda)
  
  
  #Elaboration of Coefficient of association for contingency table
  #load package (vcd)
  library(vcd)
  
  #help(package=vcd)
  contingency.table.validation.qda<-table2d_summary(cross.classification.validation.qda)
  test.table.validation.qda<-assocstats(cross.classification.validation.qda)
  
  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.validation.qda/sum(cross.classification.validation.qda)*100,2), std="margin", main="VALIDATION QDA MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.qda, std="margin", main="VALIDATION QDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #load package (verification)
  library(verification)
  
  # 2nd method using verify function
  verification.validation.qda<-verify(validation.table[,2],predict.result.qda.validation$posterior[,2], frcst.type="prob", obs.type="binary")
  #summary(verification.validation.qda)
  
  # showing confidence intervals.  MAY BE SLOW
  area.under.roc.curve.validation.qda<-roc.area(validation.table[,2],predict.result.qda.validation$posterior[,2])
  
  if (cross.classification.validation.qda[1,2]==0 | cross.classification.validation.qda[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.validation.qda, main = "ROC PLOT: VALIDATION QDA MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[1] , alpha = 0.05, extra=TRUE, legend=TRUE, show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.validation.qda$A,2),";  Sample size = ",area.under.roc.curve.validation.qda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[2], sep=""), side=3, col="red", cex=0.8)
  }
  
  # EXPORT OF PLOT FOR VALIDATION OF QDA MODEL
  
  pdf(file = "result_QDA_FourfoldPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.validation.qda/sum(cross.classification.validation.qda)*100,2), std="margin", main="VALIDATION QDA MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  #fourfold(cross.classification.validation.qda, std="margin", main="VALIDATION QDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  dev.off()
  
  #pdf(file = "result_QDA_ROCPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.validation.qda, main = "ROC PLOT: VALIDATION QDA MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.validation.qda<-roc.area(validation.table[,2],predict.result.qda.validation$posterior[,2])
  #dev.off()
  
  pdf(file = "result_QDA_ROCPlot_bootstrap_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.validation.qda, main = "ROC PLOT: VALIDATION QDA MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[1] , alpha = 0.05, extra=TRUE, legend=TRUE, show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.validation.qda$A,2),";  Sample size = ",area.under.roc.curve.validation.qda$n.total,";  Bootstrap samples = ",bootstrap.sample.values[2], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  # Assignation of a matching code between observed and predicted values calculated using the validation dataset
  validation.qda.matching.code<-paste(validation.grouping.variable,as.numeric(levels(predict.result.qda.validation$class))[predict.result.qda.validation$class],sep="")
  validation.qda.matching.code<-gsub("00","1",validation.qda.matching.code)
  validation.qda.matching.code<-gsub("01","2",validation.qda.matching.code)
  validation.qda.matching.code<-gsub("10","3",validation.qda.matching.code)
  validation.qda.matching.code<-gsub("11","4",validation.qda.matching.code)
  validation.qda.matching.code<-as.numeric(validation.qda.matching.code)
  
  ##########################################
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    cross.classification.validation.qda.optimal<-table(validation.grouping.variable,as.numeric(predict.result.qda.validation$posterior[,2]>qda.probability.optimal.binary.threshold),dnn=c("Observed","Predicted"))
    rownames(cross.classification.validation.qda.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.validation.qda.optimal)<-list("No Landslide","Landslide") # Predicted
    
    #Different plots for contingency table
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.validation.qda.optimal/sum(cross.classification.validation.qda.optimal)*100,2), std="margin", main="VALIDATION QDA MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
      #fourfold(cross.classification.validation.qda, std="margin", main="VALIDATION QDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    }
    if (cross.classification.validation.qda.optimal[1,2]==0 | cross.classification.validation.qda.optimal[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
    pdf(file = "result_QDA_FourfoldPlot_Validation_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.validation.qda.optimal/sum(cross.classification.validation.qda.optimal)*100,2), std="margin", main="VALIDATION QDA MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.qda, std="margin", main="VALIDATION QDA MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    dev.off()
    
    
    # Assignation of a optimal matching code between observed and predicted values calculated using the validation dataset
    validation.qda.matching.code.optimal<-paste(validation.grouping.variable,as.numeric(predict.result.qda.validation$posterior[,2]>qda.probability.optimal.binary.threshold),sep="")
    validation.qda.matching.code.optimal<-gsub("00","1",validation.qda.matching.code.optimal)
    validation.qda.matching.code.optimal<-gsub("01","2",validation.qda.matching.code.optimal)
    validation.qda.matching.code.optimal<-gsub("10","3",validation.qda.matching.code.optimal)
    validation.qda.matching.code.optimal<-gsub("11","4",validation.qda.matching.code.optimal)
    validation.qda.matching.code.optimal<-as.numeric(validation.qda.matching.code.optimal)
  }
  #########################################
  
  
  # EXPORT OF QDA MODEL RESULTS
  write.table("RESULTS OF QUADRATIC DISCRIMINANT ANALYSIS",file="result_QDA.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("QDA MODEL OUTPUTS",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Prior Probabilities",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(names(result.qda$prior),result.qda$prior),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Total number",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind("N",result.qda$N),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Counts",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(names(result.qda$counts),result.qda$counts),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Means",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(t(rbind(c("",colnames(result.qda$means)),cbind(rownames(result.qda$means),result.qda$means))),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Discriminant function coefficients",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  #Scaling coefficients
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(paste("Coefficients Group",levels(grouping.variable)[1]),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("",colnames(result.qda$scaling[,,1])),cbind(rownames(result.qda$scaling[,,1]),result.qda$scaling[,,1])),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(paste("Coefficients Group",levels(grouping.variable)[2]),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("",colnames(result.qda$scaling[,,2])),cbind(rownames(result.qda$scaling[,,2]),result.qda$scaling[,,2])),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE MODEL RESULT",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.qda$table[,1,],contingency.table.qda$table[,2,],contingency.table.qda$table[,3,])),file="result_QDA.txt", append=TRUE, quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE VALIDATION",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.validation.qda$table[,1,],contingency.table.validation.qda$table[,2,],contingency.table.validation.qda$table[,3,])),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("MATCHING CODE DEFINITION",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(c("","OBSERVED NO LANDSLIDES: 0","OBSERVED LANDSLIDES: 1"), c("PREDICTED NO LANDSLIDES: 0","00 -> Code 1","10 -> Code 3"), c("PREDICTED LANDSLIDES: 1","01 -> Code 2","11 -> Code 4")),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  ########
  if(enable_probability_optimal_binary_classification==FALSE) 
  {
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","QDA MODEL POSTERIOR PROBABILITY","QDA MODEL CLASSIFICATION","QDA MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict.result.qda$posterior[,2],as.numeric(levels(predict.result.qda$class))[predict.result.qda$class],result.qda.matching.code)),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","QDA VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.qda.validation$posterior[,2],as.numeric(levels(predict.result.qda.validation$class))[predict.result.qda.validation$class],validation.qda.matching.code)),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } else
  {
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    if(enable_probability_optimal_classification==TRUE) 
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY PARTITION -> Method: ",type_probability_optimal_classification,sep=""),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(qda.probability.classification.optimal[index_points_plot,c("prob_thres","tnr","tpr")]),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    } else
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY BINARY PARTITION",sep=""),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(qda.probability.classification.optimal[qda.optimal.index,c("prob_thres","tnr","tpr")]),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    }
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","QDA MODEL POSTERIOR PROBABILITY","QDA MODEL CLASSIFICATION","QDA MODEL RESULT MATCHING CODE","QDA OPTIMAL MODEL CLASSIFICATION","QDA OPTIMAL MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict.result.qda$posterior[,2],as.numeric(levels(predict.result.qda$class))[predict.result.qda$class],result.qda.matching.code,as.numeric(predict.result.qda$posterior[,2]>qda.probability.optimal.binary.threshold),result.qda.matching.code.optimal)),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","QDA VALIDATION MATCHING CODE","OPTIMAL VALIDATION CLASSIFICATION","OPTIMAL QDA VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.qda.validation$posterior[,2],as.numeric(levels(predict.result.qda.validation$class))[predict.result.qda.validation$class],validation.qda.matching.code,as.numeric(predict.result.qda.validation$posterior[,2]>qda.probability.optimal.binary.threshold),validation.qda.matching.code.optimal)),file="result_QDA.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  }
  ########
  
  # PLOT AND EXPORT OF QDA MODEL MAPS
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POLYGONS")
  {
    ################################################
    ### Prima di cancellare testare con polygoni
    ##############################################
    #shape_training_qda<-shape_training
    #result_training_qda_shape<-cbind(identification.value,training.table[,2],predict.result.qda$posterior[,2],as.numeric(levels(predict.result.qda$class))[predict.result.qda$class],result.qda.matching.code)
    #colnames(result_training_qda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH")
    #shape_training_qda@data <- merge(x=shape_training_qda@data,y=result_training_qda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_qda,dsn="result_QDA_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_training_qda,dsn="result_QDA_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    #shape_validation_qda<-shape_validation
    #result_validation_qda_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.qda.validation$posterior[,2],as.numeric(levels(predict.result.qda.validation$class))[predict.result.qda.validation$class],validation.qda.matching.code)
    #colnames(result_validation_qda_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    #shape_validation_qda@data <- merge(x=shape_validation_qda@data,y=result_validation_qda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    ##writeOGR(shape_validation_qda,dsn="result_QDA_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_validation_qda,dsn="result_QDA_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ### DA qui Nuovo
    ##################################################################   
    #### aggiungere a poligoni colonna incertezza ed esportazione pdf incertezza
    shape_training_qda<-shape_training
    result_training_qda_shape<-cbind(identification.value,training.table[,2],predict.result.qda$posterior[,2],as.numeric(levels(predict.result.qda$class))[predict.result.qda$class],result.qda.matching.code,ID.bootstrap.model.qda.count,bootstrap.model.qda.probability.mean,bootstrap.model.qda.probability.sd,bootstrap.model.qda.probability.min,bootstrap.model.qda.probability.max,bootstrap.model.qda.probability.sderror,t(bootstrap.model.qda.probability.quantiles),bootstrap.model.qda.prediction.mean,bootstrap.model.qda.prediction.sd,bootstrap.model.qda.prediction.min,bootstrap.model.qda.prediction.max,bootstrap.model.qda.prediction.sderror,t(bootstrap.model.qda.prediction.quantiles))
    colnames(result_training_qda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","QDA_SAMP","QDA_PMean","QDA_PSd","QDA_PMin","QDA_PMax","QDA_PSder","QDA_PQ0","QDA_PQ_005","QDA_PQ_025","QDA_PQ05","QDA_PQ_075","QDA_PQ095","QDA_PQ1","QDA_PrMean","QDA_PrSd","QDA_PrMin","QDA_PrMax","QDA_PrSder","QDA_PrQ0","QDA_PrQ005","QDA_PrQ025","QDA_PrQ05","QDA_PrQ075","QDA_PrQ095","QDA_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_qda_shape<-cbind(result_training_qda_shape,as.numeric(predict.result.qda$posterior[,2]>qda.probability.optimal.binary.threshold),result.qda.matching.code.optimal)
      colnames(result_training_qda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","QDA_SAMP","QDA_PMean","QDA_PSd","QDA_PMin","QDA_PMax","QDA_PSder","QDA_PQ0","QDA_PQ_005","QDA_PQ_025","QDA_PQ05","QDA_PQ_075","QDA_PQ095","QDA_PQ1","QDA_PrMean","QDA_PrSd","QDA_PrMin","QDA_PrMax","QDA_PrSder","QDA_PrQ0","QDA_PrQ005","QDA_PrQ025","QDA_PrQ05","QDA_PrQ075","QDA_PrQ095","QDA_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############

    shape_training_qda@data <- merge(x=shape_training_qda@data,y=result_training_qda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_qda,dsn="result_QDA_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_qda,dsn="result_QDA_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    # WARNING: The validation does't have the unvertainty estimation: probably this can be daone using the parabolic error function 
    shape_validation_qda<-shape_validation
    result_validation_qda_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.qda.validation$posterior[,2],as.numeric(levels(predict.result.qda.validation$class))[predict.result.qda.validation$class],validation.qda.matching.code)
    colnames(result_validation_qda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_qda_shape<-cbind(result_validation_qda_shape,as.numeric(predict.result.qda.validation$posterior[,2]>qda.probability.optimal.binary.threshold),validation.qda.matching.code.optimal)
      colnames(result_validation_qda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    shape_validation_qda@data <- merge(x=shape_validation_qda@data,y=result_validation_qda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_qda@data <- cbind(shape_validation_qda@data,PROB_SDMOD=(coefficients(fit.parabola.probability.qda)*(shape_validation_qda@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.qda)*shape_validation_qda@data$VAL_PROB))
    #writeOGR(shape_validation_qda,dsn="result_QDA_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_qda,dsn="result_QDA_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ##################################################################
    
    # Plot and export of maps
    # QDA Susceptibility
    #dev.new()
    pdf(file = "result_QDA_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_qda, zcol=c("MOD_PROB"), names.attr=c("QDA MODEL PROBABILITY"), main="QDA MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # QDA Model Matching Code
    #dev.new()
    pdf(file = "result_QDA_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_qda, zcol=c("MOD_MATCH"), names.attr=c("QDA MODEL MATCHING CODE"), main="QDA MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
    
    # QDA Model Uncertainty
    #dev.new()
    pdf(file = "result_QDA_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_qda, zcol=c("QDA_PrSd"), names.attr=c("QDA MODEL UNCERTAINTY"), main="QDA MODEL UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()
    
    # QDA Validation Susceptibility
    #dev.new()
    pdf(file = "result_QDA_Validation_Susceptibility_Map.pdf")
    print(spplot(obj=shape_validation_qda, zcol=c("VAL_PROB"), names.attr=c("QDA VALIDATION PROBABILITY"), main="QDA VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # QDA Validation Matching Code
    #dev.new()
    pdf(file = "result_QDA_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_qda, zcol=c("VAL_MATCH"), names.attr=c("QDA VALIDATION MATCHING CODE"), main="QDA VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
    
    # QDA Validation Uncertainty
    #dev.new()
    pdf(file = "result_QDA_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_qda, zcol=c("PROB_SDMOD"), names.attr=c("QDA VALIDATION UNCERTAINTY"), main="QDA VALIDATION UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()

    ########### TBT
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      # QDA Model Matching Code
      #dev.new()
      pdf(file = "result_QDA_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_training_qda, zcol=c("OPT_MATCH"), names.attr=c("QDA MODEL MATCHING CODE"), main="QDA MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      # QDA Model Matching Code
      #dev.new()
      pdf(file = "result_QDA_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_validation_qda, zcol=c("OPT_MATCH"), names.attr=c("QDA VALIDATION MATCHING CODE"), main="QDA VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      if(enable_probability_optimal_classification==TRUE)
      {
        pdf(file = "result_QDA_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
        print(spplot(obj=shape_training_qda, zcol=c("MOD_PROB"), names.attr=c("QDA MODEL PROBABILITY"), main="QDA MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=(qda.breaks.histogram.values.optimal)+c(rep(0,(length(qda.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1)))
        dev.off()
        
        pdf(file = "result_QDA_Validation_Susceptibility_Map_Optimal.pdf")
        print(spplot(obj=shape_validation_qda, zcol=c("VAL_PROB"), names.attr=c("QDA VALIDATION PROBABILITY"), main="QDA VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=(qda.breaks.histogram.values.optimal)+c(rep(0,(length(qda.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1)))
        dev.off()
      }
    }
    ###########
    
    
    ### Success & prediction rate curve
    susc_values_qda<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_qda@data) %in% c(configuration.spatial.data.table[c(5,6)],"MOD_PROB"))
    shape_training_qda@data<-shape_training_qda@data[order(shape_training_qda@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_training_qda@data[,1]<-cumsum(shape_training_qda@data[,1])/sum(shape_training_qda@data[,1])*100
    shape_training_qda@data[,2]<-cumsum(shape_training_qda@data[,2])/sum(shape_training_qda@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_qda@data[,3],shape_training_qda@data[,1])
    area_susc_values_qda_mod<-c(0,ind_pro_fun_mod(susc_values_qda[-c(1,length(susc_values_qda))]),100)
    
    #dev.new()
    pdf(file = "result_QDA_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="QDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_qda)-1))
    {
      #count=1
      polygon(c(area_susc_values_qda_mod[count:(count+1)],rev(area_susc_values_qda_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_qda@data[,1],shape_training_qda@data[,2],col="black")
    dev.off()
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_qda_mod_optimal<-c(0,ind_pro_fun_mod(rev(qda.breaks.histogram.values.optimal)[-c(1,length(qda.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_QDA_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="QDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(qda.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_qda_mod_optimal[count:(count+1)],rev(area_susc_values_qda_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_qda@data[,1],shape_training_qda@data[,2],col="black")
      dev.off()
    }
    ############
    
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_qda@data) %in% c(configuration.spatial.data.table[c(5,6)],"VAL_PROB"))
    shape_validation_qda@data<-shape_validation_qda@data[order(shape_validation_qda@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_validation_qda@data[,1]<-cumsum(shape_validation_qda@data[,1])/sum(shape_validation_qda@data[,1])*100
    shape_validation_qda@data[,2]<-cumsum(shape_validation_qda@data[,2])/sum(shape_validation_qda@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_qda@data[,3],shape_validation_qda@data[,1])
    area_susc_values_qda_val<-c(0,ind_pro_fun_val(susc_values_qda[-c(1,length(susc_values_qda))]),100)
    
    #dev.new()
    pdf(file = "result_QDA_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="QDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_qda)-1))
    {
      #count=1
      polygon(c(area_susc_values_qda_val[count:(count+1)],rev(area_susc_values_qda_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_qda@data[,1],shape_validation_qda@data[,2],col="black")
    dev.off()
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_qda_val_optimal<-c(0,ind_pro_fun_val(rev(qda.breaks.histogram.values.optimal)[-c(1,length(qda.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_QDA_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="QDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(qda.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_qda_val_optimal[count:(count+1)],rev(area_susc_values_qda_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_qda@data[,1],shape_validation_qda@data[,2],col="black")
      dev.off()
    }
    ############
    
  }
  
  
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POINTS")
  {
    shape_training_qda<-shape_training
    result_training_qda_shape<-cbind(identification.value,training.table[,2],predict.result.qda$posterior[,2],as.numeric(levels(predict.result.qda$class))[predict.result.qda$class],result.qda.matching.code,ID.bootstrap.model.qda.count,bootstrap.model.qda.probability.mean,bootstrap.model.qda.probability.sd,bootstrap.model.qda.probability.min,bootstrap.model.qda.probability.max,bootstrap.model.qda.probability.sderror,t(bootstrap.model.qda.probability.quantiles),bootstrap.model.qda.prediction.mean,bootstrap.model.qda.prediction.sd,bootstrap.model.qda.prediction.min,bootstrap.model.qda.prediction.max,bootstrap.model.qda.prediction.sderror,t(bootstrap.model.qda.prediction.quantiles))
    colnames(result_training_qda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","QDA_SAMP","QDA_PMean","QDA_PSd","QDA_PMin","QDA_PMax","QDA_PSder","QDA_PQ0","QDA_PQ_005","QDA_PQ_025","QDA_PQ05","QDA_PQ_075","QDA_PQ095","QDA_PQ1","QDA_PrMean","QDA_PrSd","QDA_PrMin","QDA_PrMax","QDA_PrSder","QDA_PrQ0","QDA_PrQ005","QDA_PrQ025","QDA_PrQ05","QDA_PrQ075","QDA_PrQ095","QDA_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_qda_shape<-cbind(result_training_qda_shape,as.numeric(predict.result.qda$posterior[,2]>qda.probability.optimal.binary.threshold),result.qda.matching.code.optimal)
      colnames(result_training_qda_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","QDA_SAMP","QDA_PMean","QDA_PSd","QDA_PMin","QDA_PMax","QDA_PSder","QDA_PQ0","QDA_PQ_005","QDA_PQ_025","QDA_PQ05","QDA_PQ_075","QDA_PQ095","QDA_PQ1","QDA_PrMean","QDA_PrSd","QDA_PrMin","QDA_PrMax","QDA_PrSder","QDA_PrQ0","QDA_PrQ005","QDA_PrQ025","QDA_PrQ05","QDA_PrQ075","QDA_PrQ095","QDA_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    
    shape_training_qda@data <- merge(x=shape_training_qda@data,y=result_training_qda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_qda,dsn="result_QDA_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_qda,dsn="result_QDA_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    
    # WARNING: The validation does't have the unvertainty estimation: probabilty this can be daone using the parabolic error function 
    shape_validation_qda<-shape_validation
    result_validation_qda_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.qda.validation$posterior[,2],as.numeric(levels(predict.result.qda.validation$class))[predict.result.qda.validation$class],validation.qda.matching.code)
    colnames(result_validation_qda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_qda_shape<-cbind(result_validation_qda_shape,as.numeric(predict.result.qda.validation$posterior[,2]>qda.probability.optimal.binary.threshold),validation.qda.matching.code.optimal)
      colnames(result_validation_qda_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    
    shape_validation_qda@data <- merge(x=shape_validation_qda@data,y=result_validation_qda_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_qda@data <- cbind(shape_validation_qda@data,PROB_SDMOD=(coefficients(fit.parabola.probability.qda)*(shape_validation_qda@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.qda)*shape_validation_qda@data$VAL_PROB))
    #writeOGR(shape_validation_qda,dsn="result_QDA_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_qda,dsn="result_QDA_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    require(raster)
    
    # Plot and export of maps
    # QDA Susceptibility
    #dev.new()
    pdf(file = "result_QDA_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_qda
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_QDA_Model_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_QDA_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1),breaks=round(qda.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    # QDA Model Matching Code
    #dev.new()
    pdf(file = "result_QDA_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_qda
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
    print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    
    writeRaster(layer_gridded_raster, filename="result_QDA_Model_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_QDA_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_qda
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
      print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_QDA_Model_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_QDA_Model_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_qda
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_QDA_Model_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
    }
    ############
    
    # QDA Model uncertainity
    #dev.new()
    pdf(file = "result_QDA_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_qda
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"QDA_PrSd"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_QDA_Model_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    
    
    # QDA Validation Susceptibility
    #dev.new()
    pdf(file = "result_QDA_Validation_Susceptibility_Map.pdf")
    layer_gridded<-shape_validation_qda
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_QDA_Validation_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_QDA_Validation_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1),breaks=round(qda.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    
    
    # QDA Validation Matching Code
    #dev.new()
    pdf(file = "result_QDA_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_qda
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    #zoom(layer_gridded_raster)  	
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_QDA_Validation_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    
    ########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_QDA_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_qda
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      #zoom(layer_gridded_raster)  	
      dev.off()
      writeRaster(layer_gridded_raster, filename="result_QDA_Validation_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_QDA_Validation_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_qda
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_QDA_Validation_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
      
    }
    #######
    
    
    # QDA Model uncertainity
    #dev.new()
    pdf(file = "result_QDA_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_qda
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"PROB_SDMOD"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_QDA_Validation_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    
    ### Success & prediction rate curve
    susc_values_qda<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_qda@data) %in% c("MOD_GROUP","MOD_PROB"))
    shape_training_qda@data<-cbind(PIXEL_AREA=rep(configuration.spatial.data.table[c(8)]^2,dim(shape_training_qda@data)[1]),shape_training_qda@data[order(shape_training_qda@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_training_qda@data$MOD_GROUP<-shape_training_qda@data$MOD_GROUP*configuration.spatial.data.table[c(8)]^2
    shape_training_qda@data[,1]<-cumsum(shape_training_qda@data[,1])/sum(shape_training_qda@data[,1])*100
    shape_training_qda@data[,2]<-cumsum(shape_training_qda@data[,2])/sum(shape_training_qda@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_qda@data[,3],shape_training_qda@data[,1])
    area_susc_values_qda_mod<-c(0,ind_pro_fun_mod(susc_values_qda[-c(1,length(susc_values_qda))]),100)
    
    #dev.new()
    pdf(file = "result_QDA_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="QDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_qda)-1))
    {
      #count=1
      polygon(c(area_susc_values_qda_mod[count:(count+1)],rev(area_susc_values_qda_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_qda@data[,1],shape_training_qda@data[,2],col="black")
    dev.off()
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_qda_mod_optimal<-c(0,ind_pro_fun_mod(rev(qda.breaks.histogram.values.optimal)[-c(1,length(qda.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_QDA_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="QDA SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(qda.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_qda_mod_optimal[count:(count+1)],rev(area_susc_values_qda_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_qda@data[,1],shape_training_qda@data[,2],col="black")
      dev.off()
    }
    ############
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_qda@data) %in% c("VAL_GROUP","VAL_PROB"))
    shape_validation_qda@data<-cbind(PIXEL_AREA=rep(as.numeric(configuration.spatial.data.table[c(8)])^2,dim(shape_validation_qda@data)[1]),shape_validation_qda@data[order(shape_validation_qda@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_validation_qda@data$VAL_GROUP<-shape_validation_qda@data$VAL_GROUP*as.numeric(configuration.spatial.data.table[c(8)])^2
    shape_validation_qda@data[,1]<-cumsum(shape_validation_qda@data[,1])/sum(shape_validation_qda@data[,1])*100
    shape_validation_qda@data[,2]<-cumsum(shape_validation_qda@data[,2])/sum(shape_validation_qda@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_qda@data[,3],shape_validation_qda@data[,1])
    area_susc_values_qda_val<-c(0,ind_pro_fun_val(susc_values_qda[-c(1,length(susc_values_qda))]),100)
    
    #dev.new()
    pdf(file = "result_QDA_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="QDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_qda)-1))
    {
      #count=1
      polygon(c(area_susc_values_qda_val[count:(count+1)],rev(area_susc_values_qda_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_qda@data[,1],shape_validation_qda@data[,2],col="black")
    dev.off()
    ###########
    
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_qda_val_optimal<-c(0,ind_pro_fun_val(rev(qda.breaks.histogram.values.optimal)[-c(1,length(qda.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_QDA_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="QDA PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(qda.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_qda_val_optimal[count:(count+1)],rev(area_susc_values_qda_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(qda.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_qda@data[,1],shape_validation_qda@data[,2],col="black")
      dev.off()
    }
    ############
    
  }
  
  
  explanatory.variables<-training.table[,3:dim(training.table)[2]] # Restore to original values of explanatory variables
  validation.explanatory.variables<-validation.table[3:dim(validation.table)[2]] # Restore to original values of validation explanatory variables
  
}


#-------------------- LOGISTIC REGRESSION MODEL ---------------------#

if(model.run.matrix[3] == "YES")
  {
  #library(Zelig)
  #if (class(try(zelig(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables, model="logit",cite=FALSE)))=="try-error")
  if (class(try(glm(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables,family=binomial())))=="try-error")
  { 
    #zelig(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables, model="logit")
    write.table("Analysis based on Logistic Regression Model was not completed",file="Error_LRM_Analysis.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="Error_LRM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("Error LOG",file="Error_LRM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind("Message",rev(1:length(as.vector(.Traceback)))," ->",as.vector(.Traceback)),file="Error_LRM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } 
  
  result.lrm<-NULL
  #result.lrm<-zelig(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables, model="logit",cite=FALSE)
  result.lrm<-glm(as.formula(paste(names(data.variables)[1],"~",paste(names(data.variables[,2:dim(data.variables)[2]]),collapse= "+"))), data=data.variables,family=binomial())
  #str(result.lrm)
  #names(result.lrm)
  
  #for predicted value (posterior probablity calculated with model) result.lrm$result$fitted.values was considered
  
  # Result Predicted  
  #test.explanatory.variables<-setx(result.lrm, fn=NULL, cond=FALSE)
  
  #problem with sim: while zelig() and sex() work well, sim() doesn't work. Changing data.variables (escludind last two variables) sim() works well. Why?    
  #predict.result.lrm<-sim(result.lrm, test.explanatory.variables, num=c(2,2), prev = NULL, bootstrap=FALSE, bootfn=NULL)
  #names(predict.result.lrm)
  
  #predict.result.lrm$result$qi$ev[]
  #predict.result.lrm$result$qi$pr[]
  
  #plot(predict.result.lrm)
  
  
  ### WARNING: From Zelig for R-2.15.1 result.lrm$fitted.values is substitute with result.lrm$result$fitted.values
  # result.lrm$y with result.lrm$result$y
  # result.lrm$coefficients with result.lrm$result$coefficients
  # see str(result.lrm$result) and #str(result.lrm)
  # To run the script with R afeter version R-2.15.1 substitute result.lrm$ with result.lrm$result
  
  #cross.classification.lrm<-table(as.numeric(result.lrm$result$y),round(result.lrm$result$fitted.values),dnn=c("Observed","Predicted"))
  cross.classification.lrm<-table(grouping.variable,round(predict(result.lrm, type="response")),dnn=c("Observed","Predicted"))
  rownames(cross.classification.lrm)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.lrm)<-list("No Landslide","Landslide") # Predicted    
  str(cross.classification.lrm)
  
  # Assignation of a matching code between observed and predicted values
  #result.lrm.matching.code<-paste(grouping.variable,round(result.lrm$result$fitted.values),sep="")
  result.lrm.matching.code<-paste(grouping.variable,round(predict(result.lrm, type="response")),sep="")
  result.lrm.matching.code<-gsub("00","1",result.lrm.matching.code)
  result.lrm.matching.code<-gsub("01","2",result.lrm.matching.code)
  result.lrm.matching.code<-gsub("10","3",result.lrm.matching.code)
  result.lrm.matching.code<-gsub("11","4",result.lrm.matching.code)
  result.lrm.matching.code<-as.numeric(result.lrm.matching.code)
  
  #Elaboration of Coefficient of association for contingency table 
  #load package (vcd)  
  library(vcd)
  
  #help(package=vcd)         
  contingency.table.lrm<-table2d_summary(cross.classification.lrm)
  test.table.lrm<-assocstats(cross.classification.lrm)
  
  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.lrm/sum(cross.classification.lrm)*100,2), std="margin",  main="LOGISTIC REGRESSION MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.lrm, std="margin",  main="LOGISTIC REGRESSION MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #A ROC curve plots the false alarm rate against the hit rate
  #for a probablistic forecast for a range of thresholds. 
  
  #load package (verification)  
  library(verification)
  
  #verify function
  #Based on the type of inputs, this function calculates a range of verification statistics and skill scores.
  #Additionally, it creates a verify class object that can be further analyzed.
  
  
  ##### ROC PLOT OBS - POSTERIOR PROBABILITY ASSOCIATED TO 1                                                                                 
  
  # Method using verify function
  #verification.results.lrm<-verify(result.lrm$result$y,result.lrm$result$fitted.values, frcst.type="prob", obs.type="binary")
  verification.results.lrm<-verify(training.table[,2],predict(result.lrm, type="response"), frcst.type="prob", obs.type="binary")
  
  
  
  #str(verification.results.lrm)
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(verification.results.lrm, main = "ROC PLOT: LOGISTIC REGRESSION MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.lrm<-roc.area(result.lrm$result$y,result.lrm$result$fitted.values)
  #}
  area.under.roc.curve.lrm<-roc.area(training.table[,2],predict(result.lrm, type="response"))
  
  ## showing confidence intervals.  MAY BE SLOW
  
  if (cross.classification.lrm[1,2]==0 | cross.classification.lrm[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.results.lrm, main = "ROC PLOT: LOGISTIC REGRESSION MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[3] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.lrm$A,2),";  Sample size = ",area.under.roc.curve.lrm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[3], sep=""), side=3, col="red", cex=0.8)
    ## Histogram of posterior probability
    dev.new()                            
    #hist(result.lrm$result$fitted.values, breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
    hist(predict(result.lrm, type="response"), breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  }
  pdf(file = "result_LRM_Histogram.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #hist(result.lrm$result$fitted.values, breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  hist(predict(result.lrm, type="response"), breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  dev.off() 
  
  # EXPORT OF PLOT FOR LRM MODEL
  
  pdf(file = "result_LRM_FourfoldPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.lrm/sum(cross.classification.lrm)*100,2), std="margin",  main="LOGISTIC REGRESSION MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  #fourfold(cross.classification.lrm, std="margin",  main="LOGISTIC REGRESSION MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  dev.off()
  
  #pdf(file = "result_LRM_ROCPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.results.lrm, main = "ROC PLOT: LOGISTIC REGRESSION MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #dev.off()
  
  pdf(file = "result_LRM_ROCPlot_bootstrap.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.results.lrm, main = "ROC PLOT: LOGISTIC REGRESSION MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[3] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.lrm$A,2),";  Sample size = ",area.under.roc.curve.lrm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[3], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  ## BOOTSTRAP PROCEDURE FOR THE ESTIMATION OF MODEL PREDICTION VARIABILITY
  if(bootstrap.model.variability[3] == "YES")
  {
    bootstrap.sample.model.lrm<-bootstrap.sample.model[3]
    
    matrix.bootstrap.model.lrm<-matrix(data=NA, nrow=dim(training.table)[1], ncol=(bootstrap.sample.model.lrm*3)+1)
    colnames(matrix.bootstrap.model.lrm)<-rep("na",(bootstrap.sample.model.lrm*3)+1)
    matrix.bootstrap.model.lrm[,1]<-identification.value
    colnames(matrix.bootstrap.model.lrm)[1]<-"ID"
    name.sel.run<-paste(rep("ID_Selection_Run",bootstrap.sample.model.lrm),1:bootstrap.sample.model.lrm,sep="_")
    colnames(matrix.bootstrap.model.lrm)[seq(2,(bootstrap.sample.model.lrm*3)-1,3)]<-name.sel.run
    name.prob.run<-paste(rep("Probability_Run",bootstrap.sample.model.lrm),1:bootstrap.sample.model.lrm,sep="_")
    colnames(matrix.bootstrap.model.lrm)[seq(3,(bootstrap.sample.model.lrm*3),3)]<-name.prob.run
    name.pred.run<-paste(rep("Prediction_Run",bootstrap.sample.model.lrm),1:bootstrap.sample.model.lrm,sep="_")
    colnames(matrix.bootstrap.model.lrm)[seq(4,(bootstrap.sample.model.lrm*3)+1,3)]<-name.pred.run
    
    selection.index<-NULL
    #library(Zelig)
    #Bootstrap procedure
    for (count.boot in 1:bootstrap.sample.model.lrm)
    {
      selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
      matrix.bootstrap.model.lrm[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
      data.variables.bootstrap.model.lrm<-training.table[selection.index,2:dim(training.table)[2]]
      explanatory.variables.bootstrap.model.lrm<-training.table[selection.index,3:dim(training.table)[2]]
      grouping.variable.bootstrap.model.lrm<-as.factor(training.table[selection.index,2])
      #result.bootstrap.model.lrm<-zelig(as.formula(paste(names(data.variables.bootstrap.model.lrm)[1],"~",paste(names(data.variables.bootstrap.model.lrm[,2:dim(data.variables.bootstrap.model.lrm)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.lrm, model="logit",cite=FALSE)
      #result.bootstrap.model.lrm<-glm(as.formula(paste(names(data.variables.bootstrap.model.lrm)[1],"~",paste(names(data.variables.bootstrap.model.lrm[,2:dim(data.variables.bootstrap.model.lrm)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.lrm,family=binomial())
      while(inherits(try(result.bootstrap.model.lrm<-glm(as.formula(paste(names(data.variables.bootstrap.model.lrm)[1],"~",paste(names(data.variables.bootstrap.model.lrm[,2:dim(data.variables.bootstrap.model.lrm)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.lrm,family=binomial()),silent=TRUE),what="try-error"))
      {
        print(paste("Count boot: ",count.boot," - Boostrap while resampling",sep=""))
        selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
        matrix.bootstrap.model.lrm[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
        explanatory.variables.bootstrap.model.lrm<-training.table[selection.index,3:dim(training.table)[2]]
        if(bootstrap_constant_correction==TRUE)
        {
          print("Performing bootstrap 0 value correction")
          indexbootstrapcosntant<-which(explanatory.variables.bootstrap.model.lrm==0,arr.ind=TRUE)
          explanatory.variables.bootstrap.model.lrm[indexbootstrapcosntant]<-runif(length(indexbootstrapcosntant)/2, min = 0.00001, max = 0.01)
        }
        grouping.variable.bootstrap.model.lrm<-as.factor(training.table[selection.index,2])
      }
      excluded.variables.bootstrap.model.lrm<-which(match(result.bootstrap.model.lrm$coefficients,NA)==1)
      if (length(excluded.variables.bootstrap.model.lrm) != 0)
      {
        data.variables.bootstrap.model.lrm.selected<-data.variables.bootstrap.model.lrm[,-excluded.variables.bootstrap.model.lrm]
        #setx.data.prediction<-training.table[,2:dim(training.table)[2]][,-excluded.variables.bootstrap.model.lrm]
      } else
      {
        data.variables.bootstrap.model.lrm.selected<-data.variables.bootstrap.model.lrm
        #setx.data.prediction<-training.table[,2:dim(training.table)[2]]
      }
      #result.bootstrap.model.lrm.selected<-zelig(as.formula(paste(names(data.variables.bootstrap.model.lrm.selected)[1],"~",paste(names(data.variables.bootstrap.model.lrm.selected[,2:dim(data.variables.bootstrap.model.lrm.selected)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.lrm.selected, model="logit",cite=FALSE)
      result.bootstrap.model.lrm.selected<-glm(as.formula(paste(names(data.variables.bootstrap.model.lrm.selected)[1],"~",paste(names(data.variables.bootstrap.model.lrm.selected[,2:dim(data.variables.bootstrap.model.lrm.selected)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.lrm.selected,family=binomial())
      #x.result.bootstrap.model.lrm.selected.prediction<-setx(result.bootstrap.model.lrm.selected,data=setx.data.prediction,fn=NULL)
      #matrix.bootstrap.model.lrm[,(count.boot*3)+1]<-sim(result.bootstrap.model.lrm.selected,x=x.result.bootstrap.model.lrm.selected.prediction)$result$fitted.values
      #matrix.bootstrap.model.lrm[as.numeric(names(table(selection.index))),(count.boot*3)]<-matrix.bootstrap.model.lrm[as.numeric(names(table(selection.index))),(count.boot*3)+1]
      matrix.bootstrap.model.lrm[,(count.boot*3)+1]<-predict(result.bootstrap.model.lrm.selected,newdata=explanatory.variables,type="response")
      matrix.bootstrap.model.lrm[as.numeric(names(table(selection.index))),(count.boot*3)]<-matrix.bootstrap.model.lrm[as.numeric(names(table(selection.index))),(count.boot*3)+1]
    }
    
    # Export of bootstrap sample
    write.table(matrix.bootstrap.model.lrm,file="result_LRM_BootstrapSamples.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    
    ID.bootstrap.model.lrm.count<-numeric(length=dim(training.table)[1])
    #Probability (selected values)
    bootstrap.model.lrm.probability.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.probability.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.probability.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.probability.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.probability.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.probability.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #Prediction (all values)
    bootstrap.model.lrm.prediction.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.prediction.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.prediction.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.prediction.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.prediction.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.lrm.prediction.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #    for (count.row.variability in 1:dim(training.table)[1])
    #        {
    #        # Statistics on boostrapped probability
    #        ID.bootstrap.model.lrm.count[count.row.variability]<-length(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(2,(bootstrap.sample.model.lrm*3)-1,3)]))
    #        bootstrap.model.lrm.probability.mean[count.row.variability]<-mean(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(3,(bootstrap.sample.model.lrm*3),3)]))
    #        bootstrap.model.lrm.probability.sd[count.row.variability]<-sd(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(3,(bootstrap.sample.model.lrm*3),3)]))
    #        bootstrap.model.lrm.probability.min[count.row.variability]<-min(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(3,(bootstrap.sample.model.lrm*3),3)]))
    #        bootstrap.model.lrm.probability.max[count.row.variability]<-max(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(3,(bootstrap.sample.model.lrm*3),3)]))
    #        bootstrap.model.lrm.probability.sderror[count.row.variability]<-bootstrap.model.lrm.probability.sd[count.row.variability]/ID.bootstrap.model.lrm.count[count.row.variability]
    #  	series.matrix.bootstrap.probability<-matrix.bootstrap.model.lrm[count.row.variability,seq(4,(bootstrap.sample.model.lrm*3)+1,3)]
    #		bootstrap.model.lrm.probability.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(3,(bootstrap.sample.model.lrm*3),3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        # Statistics on boostrapped prediction
    #        bootstrap.model.lrm.prediction.mean[count.row.variability]<-mean(matrix.bootstrap.model.lrm[count.row.variability,seq(4,(bootstrap.sample.model.lrm*3)+1,3)])
    #        bootstrap.model.lrm.prediction.sd[count.row.variability]<-sd(matrix.bootstrap.model.lrm[count.row.variability,seq(4,(bootstrap.sample.model.lrm*3)+1,3)])
    #        bootstrap.model.lrm.prediction.min[count.row.variability]<-min(matrix.bootstrap.model.lrm[count.row.variability,seq(4,(bootstrap.sample.model.lrm*3)+1,3)])
    #        bootstrap.model.lrm.prediction.max[count.row.variability]<-max(matrix.bootstrap.model.lrm[count.row.variability,seq(4,(bootstrap.sample.model.lrm*3)+1,3)])
    #        bootstrap.model.lrm.prediction.sderror[count.row.variability]<-bootstrap.model.lrm.prediction.sd[count.row.variability]/bootstrap.sample.model.lrm
    #		bootstrap.model.lrm.prediction.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.lrm[count.row.variability,seq(4,(bootstrap.sample.model.lrm*3)+1,3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #		}
    
    fun_length<-function(x) {length_var<-length(x[which(is.finite(x))]); return(length_var)}
    fun_quantile<-function(x) {quantile_var<-t(quantile(x[which(is.finite(x))],probs=c(0,0.05,0.25,0.5,0.75,0.95,1))); return(quantile_var)}
    ID.bootstrap.model.lrm.count<-apply(matrix.bootstrap.model.lrm[,grep("ID_Selection",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=fun_length)
    bootstrap.model.lrm.probability.mean<-apply(matrix.bootstrap.model.lrm[,grep("Probability",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.lrm.probability.sd<-apply(matrix.bootstrap.model.lrm[,grep("Probability",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.lrm.probability.min<-apply(matrix.bootstrap.model.lrm[,grep("Probability",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.lrm.probability.max<-apply(matrix.bootstrap.model.lrm[,grep("Probability",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.lrm.probability.sderror<-bootstrap.model.lrm.probability.sd/bootstrap.sample.model.lrm
    bootstrap.model.lrm.probability.quantiles<-apply(matrix.bootstrap.model.lrm[,grep("Probability",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=fun_quantile)
    bootstrap.model.lrm.prediction.mean<-apply(matrix.bootstrap.model.lrm[,grep("Prediction",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.lrm.prediction.sd<-apply(matrix.bootstrap.model.lrm[,grep("Prediction",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.lrm.prediction.min<-apply(matrix.bootstrap.model.lrm[,grep("Prediction",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.lrm.prediction.max<-apply(matrix.bootstrap.model.lrm[,grep("Prediction",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.lrm.prediction.sderror<-bootstrap.model.lrm.prediction.sd/bootstrap.sample.model.lrm
    bootstrap.model.lrm.prediction.quantiles<-apply(matrix.bootstrap.model.lrm[,grep("Prediction",colnames(matrix.bootstrap.model.lrm))],MARGIN=1,FUN=fun_quantile)
    
    
    
    # Export of bootstrap sample statistics
    write.table(cbind("ID","LRM_NumberSelectedSamples","LRM_Probability_Mean","LRM_Probability_Sd","LRM_Probability_Min","LRM_Probability_Max","LRM_Probability_Sderror","LRM_Probability_Quantiles_0","LRM_Probability_Quantiles_0.05","LRM_Probability_Quantiles_0.25","LRM_Probability_Quantiles_0.5","LRM_Probability_Quantiles_0.75","LRM_Probability_Quantiles_0.95","LRM_Probability_Quantiles_1","LRM_Prediction_Mean","LRM_Prediction_Sd","LRM_Prediction_Min","LRM_Prediction_Max","LRM_Prediction_Sderror","LRM_Prediction_Quantiles_0","LRM_Prediction_Quantiles_0.05","LRM_Prediction_Quantiles_0.25","LRM_Prediction_Quantiles_0.5","LRM_Prediction_Quantiles_0.75","LRM_Prediction_Quantiles_0.95","LRM_Prediction_Quantiles_1"),file="result_LRM_BootstrapStatistics.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind(identification.value,ID.bootstrap.model.lrm.count,bootstrap.model.lrm.probability.mean,bootstrap.model.lrm.probability.sd,bootstrap.model.lrm.probability.min,bootstrap.model.lrm.probability.max,bootstrap.model.lrm.probability.sderror,t(bootstrap.model.lrm.probability.quantiles),bootstrap.model.lrm.prediction.mean,bootstrap.model.lrm.prediction.sd,bootstrap.model.lrm.prediction.min,bootstrap.model.lrm.prediction.max,bootstrap.model.lrm.prediction.sderror,t(bootstrap.model.lrm.prediction.quantiles)),file="result_LRM_BootstrapStatistics.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    
    if (enable_screen_plotting==TRUE)
    {
      #dev.new()
      #double.sd.histogram.variability<-hist(bootstrap.model.lrm.probability.sd*2,breaks=seq(0,1,0.05),labels=TRUE)
      #plot(double.sd.histogram.variability$counts, seq(0,0.95,0.05), type="S",ylim=c(0,1), labels=TRUE)
      dev.new()
      plot(bootstrap.model.lrm.probability.mean,bootstrap.model.lrm.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="LRM BOOTSTRAP: Mean Probability vs Mean Prediction")
      abline(a=0,b=1,col="red",lty=1,lwd=1)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lrm,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    }
    
    pdf(file = "result_LRM_BootstrapMeansComparison.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(bootstrap.model.lrm.probability.mean,bootstrap.model.lrm.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="LRM BOOTSTRAP: Mean Probability vs Mean Prediction")
    abline(a=0,b=1,col="red",lty=1,lwd=1)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lrm,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    dev.off()
    
    
    # BOOTSTRAPPED PROBABILITY - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.probability.lrm<-cbind(bootstrap.model.lrm.probability.mean,2*bootstrap.model.lrm.probability.sd)
    parabola.probability.lrm<-na.omit(parabola.probability.lrm[order(parabola.probability.lrm[,1]),])
    colnames(parabola.probability.lrm)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.probability.lrm <- nls(parabola.probability.lrm[,"ordinate"] ~ coeff.a*(parabola.probability.lrm[,"abscissa"]^2) + (-1)*coeff.a*parabola.probability.lrm[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.probability.lrm<-predict(fit.parabola.probability.lrm)
    #coef(fit.parabola.probability.lrm)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.probability.lrm[,"abscissa"],parabola.probability.lrm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="LRM Model Probability Variability (Bootstrap)")
      lines(parabola.probability.lrm[,"abscissa"],value.parabola.probability.lrm,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lrm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.lrm),3),coeff.b= -round(coef(fit.parabola.probability.lrm),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_LRM_BootstrapProbabilityVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.probability.lrm[,"abscissa"],parabola.probability.lrm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="LRM Model Probability Variability (Bootstrap)")
    lines(parabola.probability.lrm[,"abscissa"],value.parabola.probability.lrm,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lrm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.lrm),3),coeff.b= -round(coef(fit.parabola.probability.lrm),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
    
    # BOOTSTRAPPED PREDICTION - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.prediction.lrm<-cbind(bootstrap.model.lrm.prediction.mean,2*bootstrap.model.lrm.prediction.sd)
    parabola.prediction.lrm<-parabola.prediction.lrm[order(parabola.prediction.lrm[,1]),]
    colnames(parabola.prediction.lrm)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.prediction.lrm <- nls(parabola.prediction.lrm[,"ordinate"] ~ coeff.a*(parabola.prediction.lrm[,"abscissa"]^2) + (-1)*coeff.a*parabola.prediction.lrm[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.prediction.lrm<-predict(fit.parabola.prediction.lrm)
    #coef(fit.parabola.prediction.lrm)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.prediction.lrm[,"abscissa"],parabola.prediction.lrm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="LRM Model Prediction Variability (Bootstrap)")
      lines(parabola.prediction.lrm[,"abscissa"],value.parabola.prediction.lrm,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lrm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.lrm),3),coeff.b= -round(coef(fit.parabola.prediction.lrm),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_LRM_BootstrapPredictionVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.prediction.lrm[,"abscissa"],parabola.prediction.lrm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="LRM Model Prediction Variability (Bootstrap)")
    lines(parabola.prediction.lrm[,"abscissa"],value.parabola.prediction.lrm,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.lrm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.lrm),3),coeff.b= -round(coef(fit.parabola.prediction.lrm),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
  }
  
  ## Sensitivity, Specificity, Cohens kappa plot
  dev.new()
  roc.plot.lrm.series<-roc.plot(verification.results.lrm,binormal=bin,plot=FALSE,show.thres=plot_thres,thresholds=threshold_series)
  dev.off()
  #str(roc.plot.lrm.series)
  #roc.plot.lrm.series$plot.data
  #str(roc.plot.lrm.series$plot.data)
  
  
  ###########################################
  # min(abs(TPR - (1-FPR)))
  if(enable_probability_optimal_binary_classification==FALSE)
  {
    dev.new()
    roc.plot.lrm.series<-roc.plot(verification.results.lrm,binormal=bin,plot=FALSE,show.thres=FALSE,thresholds=threshold_series)
    dev.off()
  } else
  {
    dev.new()
    #roc.plot.lrm.series<-roc.plot(verification.results.lrm,binormal=bin,plot=FALSE,show.thres=FALSE)
    # the following to fasten the calculation but less accurate
    roc.plot.lrm.series<-roc.plot(verification.results.lrm,binormal=bin,plot=FALSE,show.thres=FALSE,thresholds=threshold_series)
    dev.off()  
  }
  
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    lrm.probability.classification.optimal<-data.frame(prob_thres=roc.plot.lrm.series$plot.data[,1,1],tpr=roc.plot.lrm.series$plot.data[,2,1],fpr=roc.plot.lrm.series$plot.data[,3,1],tnr=(1-roc.plot.lrm.series$plot.data[,3,1]),diff_abs_tpr_tnr=abs(roc.plot.lrm.series$plot.data[,2,1]-(1-roc.plot.lrm.series$plot.data[,3,1])),optimal_sel=NA,breaks_sel=NA)
    index.lrm.filter<-which(lrm.probability.classification.optimal$prob_thres>0 & lrm.probability.classification.optimal$prob_thres<1) # removing strnge thresh values
    lrm.probability.classification.optimal<-rbind(c(0,1,1,0,1,NA,NA),lrm.probability.classification.optimal[index.lrm.filter,],c(1,0,0,1,1,NA,NA))
    lrm.optimal.index<-which(lrm.probability.classification.optimal$diff_abs_tpr_tnr==min(lrm.probability.classification.optimal$diff_abs_tpr_tnr))
    lrm.probability.classification.optimal$optimal_sel[lrm.optimal.index]<-TRUE
    lrm.probability.optimal.binary.threshold<-lrm.probability.classification.optimal$prob_thres[lrm.optimal.index]
    ### Generating the optimal fourfold plot
    cross.classification.lrm.optimal<-table(grouping.variable,predict(result.lrm, type="response")>lrm.probability.optimal.binary.threshold,dnn=c("Observed","Predicted"))
    rownames(cross.classification.lrm.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.lrm.optimal)<-list("No Landslide","Landslide") # Predicted    
    str(cross.classification.lrm.optimal)
    # Assignation of a matching code between observed and predicted values
    result.lrm.matching.code.optimal<-paste(grouping.variable,as.numeric(predict(result.lrm, type="response")>lrm.probability.optimal.binary.threshold),sep="")
    result.lrm.matching.code.optimal<-gsub("00","1",result.lrm.matching.code.optimal)
    result.lrm.matching.code.optimal<-gsub("01","2",result.lrm.matching.code.optimal)
    result.lrm.matching.code.optimal<-gsub("10","3",result.lrm.matching.code.optimal)
    result.lrm.matching.code.optimal<-gsub("11","4",result.lrm.matching.code.optimal)
    result.lrm.matching.code.optimal<-as.numeric(result.lrm.matching.code.optimal)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.lrm.optimal/sum(cross.classification.lrm.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    }
    # EXPORT OF PLOT FOR LRM MODEL
    pdf(file = "result_LRM_FourfoldPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.lrm.optimal/sum(cross.classification.lrm.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.lrm.optimal, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    dev.off()
    
    ### Optimal susceptibility classes identification 
    if(enable_probability_optimal_classification==TRUE)
    {
      lrm.unexplained.errors<-round(1-(lrm.probability.classification.optimal$tpr[lrm.optimal.index]+lrm.probability.classification.optimal$tpr[lrm.optimal.index])/2,2)
      
      if(type_probability_optimal_classification=="proportional")
      {
        lrm.unexplained.errors.partition<-1-(lrm.unexplained.errors/((length(breaks.histogram.values)/2)))*(1:((length(breaks.histogram.values)/2)))
        
        for(count_part in 1:(length(lrm.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-lrm.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((lrm.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          lrm.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((lrm.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          lrm.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        lrm.breaks.histogram.values.optimal<-c(0,lrm.probability.classification.optimal$prob_thres[which(lrm.probability.classification.optimal$breaks_sel==TRUE)],1)
        #lrm.probability.classification.optimal[which(lrm.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if(type_probability_optimal_classification=="fixed")
      {
        
        step.lrm.unexplained.fixed<-0.1
        if(lrm.unexplained.errors<=0.1) step.lrm.unexplained.fixed<-0.05
        if(lrm.unexplained.errors<=0.05) step.lrm.unexplained.fixed<-0.025
        lrm.unexplained.errors.partition<-seq(step.lrm.unexplained.fixed,1-step.lrm.unexplained.fixed,step.lrm.unexplained.fixed)[seq(step.lrm.unexplained.fixed,1-step.lrm.unexplained.fixed,step.lrm.unexplained.fixed)>(1-lrm.unexplained.errors)]
        
        for(count_part in 1:(length(lrm.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-lrm.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((lrm.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          lrm.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((lrm.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          lrm.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        lrm.breaks.histogram.values.optimal<-c(0,lrm.probability.classification.optimal$prob_thres[which(lrm.probability.classification.optimal$breaks_sel==TRUE)],1)
        #lrm.probability.classification.optimal[which(lrm.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if (enable_screen_plotting==TRUE)
      {
        dev.new()
        hist(predict(result.lrm, type="response"), breaks=lrm.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal LRM susceptibility", col=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1))
        plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("LRM OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
        mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
        mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
        for (count in 1:(length(lrm.breaks.histogram.values.optimal)-1))
        {
          #count=1
          polygon(c(lrm.breaks.histogram.values.optimal[count:(count+1)],rev(lrm.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1)[count])  
        }
        polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
        lines(lrm.probability.classification.optimal$prob_thres,lrm.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
        lines(lrm.probability.classification.optimal$prob_thres,lrm.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
        index_points_plot<-c(1,which(lrm.probability.classification.optimal$breaks_sel==TRUE),dim(lrm.probability.classification.optimal)[1])
        points(lrm.probability.classification.optimal[index_points_plot[1:floor(length(lrm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
        text(lrm.probability.classification.optimal[index_points_plot[1:floor(length(lrm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(lrm.probability.classification.optimal[index_points_plot[1:floor(length(lrm.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(lrm.breaks.histogram.values.optimal)/2)-1)))
        points(lrm.probability.classification.optimal[index_points_plot[ceiling(length(lrm.breaks.histogram.values.optimal)/2):length(lrm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
        text(lrm.probability.classification.optimal[index_points_plot[ceiling(length(lrm.breaks.histogram.values.optimal)/2):length(lrm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(lrm.probability.classification.optimal[index_points_plot[ceiling(length(lrm.breaks.histogram.values.optimal)/2):length(lrm.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(lrm.breaks.histogram.values.optimal)-ceiling(length(lrm.breaks.histogram.values.optimal)/2)),3))
      }
      
      pdf(file = "result_LRM_Histogram_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      hist(predict(result.lrm, type="response"), breaks=lrm.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal LRM susceptibility", col=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1))
      dev.off()
      
      pdf(file = "result_LRM_ModelEvaluationPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("LRM OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
      mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
      mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
      for (count in 1:(length(lrm.breaks.histogram.values.optimal)-1))
      {
        #count=1
        polygon(c(lrm.breaks.histogram.values.optimal[count:(count+1)],rev(lrm.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1)[count])  
      }
      polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
      lines(lrm.probability.classification.optimal$prob_thres,lrm.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
      lines(lrm.probability.classification.optimal$prob_thres,lrm.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
      index_points_plot<-c(1,which(lrm.probability.classification.optimal$breaks_sel==TRUE),dim(lrm.probability.classification.optimal)[1])
      points(lrm.probability.classification.optimal[index_points_plot[1:floor(length(lrm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(lrm.probability.classification.optimal[index_points_plot[1:floor(length(lrm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(lrm.probability.classification.optimal[index_points_plot[1:floor(length(lrm.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(lrm.breaks.histogram.values.optimal)/2)-1)))
      points(lrm.probability.classification.optimal[index_points_plot[1+ceiling(length(lrm.breaks.histogram.values.optimal)/2):length(lrm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
      text(lrm.probability.classification.optimal[index_points_plot[1+ceiling(length(lrm.breaks.histogram.values.optimal)/2):length(lrm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(lrm.probability.classification.optimal[index_points_plot[1+ceiling(length(lrm.breaks.histogram.values.optimal)/2):length(lrm.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(lrm.breaks.histogram.values.optimal)-1-ceiling(length(lrm.breaks.histogram.values.optimal)/2)),3))
      points(lrm.probability.classification.optimal[lrm.optimal.index,c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(lrm.probability.classification.optimal[lrm.optimal.index,c("prob_thres","tpr")],labels=with(round(lrm.probability.classification.optimal[lrm.optimal.index,c("prob_thres","tpr")],3),paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(4))
      dev.off()
    }
  }
  
  ###########################################
  
  
  contingency.table.matrix.lrm<-matrix(nrow=dim(roc.plot.lrm.series$plot.data)[1],ncol=8)
  colnames(contingency.table.matrix.lrm)<-c("Threshold","TP","TN","FP","FN","TPR","FPR","COHEN_KAPPA")
  contingency.table.matrix.lrm[,1]<-roc.plot.lrm.series$plot.data[,1,1]
  contingency.table.matrix.lrm[,6]<-roc.plot.lrm.series$plot.data[,2,1]
  contingency.table.matrix.lrm[,7]<-roc.plot.lrm.series$plot.data[,3,1]
  values.observed<-training.table[,2]
  #values.predicted<-result.lrm$result$fitted.values
  values.predicted<-predict(result.lrm,type="response")
  for (count.threshold.series in 1:dim(roc.plot.lrm.series$plot.data)[1])
  {
    value.threshold<-contingency.table.matrix.lrm[count.threshold.series,1]
    values.probability.reclassified<-NULL
    values.probability.reclassified<-as.numeric(values.predicted>value.threshold) 
    #sum(values.probability.reclassified-round(values.predicted)) # Check sum: It has to be 0 if threshold is equal to 1
    series.pasted<-paste(values.observed,values.probability.reclassified,sep="")
    series.pasted<-gsub("00","1",series.pasted)
    series.pasted<-gsub("01","2",series.pasted)
    series.pasted<-gsub("10","3",series.pasted)
    series.pasted<-gsub("11","4",series.pasted)
    series.pasted<-as.numeric(series.pasted)
    TP<-as.numeric(sum(series.pasted>=4)) # True Positive
    FN<-as.numeric(sum(series.pasted>=3 & series.pasted<4)) # False Negative
    FP<-as.numeric(sum(series.pasted>=2 & series.pasted<3)) # False Positive
    TN<-as.numeric(sum(series.pasted>=1 & series.pasted<2)) # True Negative              
    #TPR<-TP/(TP+FN) # Hit Rate or True Positive Rate or Sensitivity - Assigned before the for cicle using rocplot data
    #FPR<-FP/(FP+TN) # False Alarm Rate or False Positive Rate or 1-Specificity
    # Cohen's Kappa = (agreement-chance)/(1-chance)  where agreement=(TP+TN)/(TP+TN+FP+FN) and chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    agreement=(TP+TN)/(TP+TN+FP+FN)
    chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    cohen.kappa.value<-(agreement-chance)/(1-chance)
    #Other
    #library(vcd)
    #cohen.kappa.value<-Kappa(cross.classification.table)
    contingency.table.matrix.lrm[count.threshold.series,2]<-TP
    contingency.table.matrix.lrm[count.threshold.series,3]<-TN
    contingency.table.matrix.lrm[count.threshold.series,4]<-FP
    contingency.table.matrix.lrm[count.threshold.series,5]<-FN
    contingency.table.matrix.lrm[count.threshold.series,8]<-cohen.kappa.value
  }
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    plot(roc.plot.lrm.series$plot.data[,1,1],roc.plot.lrm.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="LRM MODEL EVALUATION PLOT")
    lines(roc.plot.lrm.series$plot.data[,1,1],1-roc.plot.lrm.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
    lines(roc.plot.lrm.series$plot.data[,1,1], contingency.table.matrix.lrm[,8],col="blue",lty=1,lwd=1)
    mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
    mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
    mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  }
  pdf(file = "result_LRM_ModelEvaluationPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(roc.plot.lrm.series$plot.data[,1,1],roc.plot.lrm.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="LRM MODEL EVALUATION PLOT")
  lines(roc.plot.lrm.series$plot.data[,1,1],1-roc.plot.lrm.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.lrm.series$plot.data[,1,1], contingency.table.matrix.lrm[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  dev.off()
  
  ## VALIDATION OF LRM MODEL (Matching LRM posterior probability results and validation grouping variable)
  
  # Result Predicted
  #test.explanatory.variables.validation.lrm<-setx(result.lrm, data=validation.table[,2:dim(validation.table)[2]],fn=NULL)
  #predict.result.lrm.validation.posterior<-sim(result.lrm, x=test.explanatory.variables.validation.lrm)$result$fitted.values
  #predict.result.lrm.validation.class<-as.numeric(round(sim(result.lrm, x=test.explanatory.variables.validation.lrm)$result$fitted.values))
  predict.result.lrm.validation.posterior<-predict(result.lrm,newdata=validation.variables,type="response")
  predict.result.lrm.validation.class<-as.numeric(round(predict(result.lrm,newdata=validation.variables,type="response")))
  
  cross.classification.validation.lrm<-table(validation.table[,2],predict.result.lrm.validation.class,dnn=c("Observed","Predicted"))
  rownames(cross.classification.validation.lrm)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.validation.lrm)<-list("No Landslide","Landslide") # Predicted
  str(cross.classification.validation.lrm)
  #cross.classification.validation.lrm<-table(validation.grouping.variable,round(result.lrm$result$fitted.values),dnn=c("Observed","Predicted"))
  
  
  #Elaboration of Coefficient of association for contingency table
  #load package (vcd)
  library(vcd)
  
  #help(package=vcd)
  contingency.table.validation.lrm<-table2d_summary(cross.classification.validation.lrm)
  test.table.validation.lrm<-assocstats(cross.classification.validation.lrm)
  
  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.validation.lrm/sum(cross.classification.validation.lrm)*100,2), std="margin", main="VALIDATION LRM MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.lrm, std="margin", main="VALIDATION LRM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #load package (verification)
  library(verification)
  
  # 2nd method using verify function
  verification.validation.lrm<-verify(validation.table[,2],predict.result.lrm.validation.posterior, frcst.type="prob", obs.type="binary")
  #summary(verification.validation.lrm)
  
  # showing confidence intervals.  MAY BE SLOW
  area.under.roc.curve.validation.lrm<-roc.area(validation.table[,2],predict.result.lrm.validation.posterior)
  
  if (cross.classification.validation.lrm[1,2]==0 | cross.classification.validation.lrm[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.validation.lrm, main = "ROC PLOT: VALIDATION LRM MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[3] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.validation.lrm$A,2),";  Sample size = ",area.under.roc.curve.validation.lrm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[3], sep=""), side=3, col="red", cex=0.8)
  }
  
  # EXPORT OF PLOT FOR VALIDATION OF LRM MODEL
  
  pdf(file = "result_LRM_FourfoldPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.validation.lrm/sum(cross.classification.validation.lrm)*100,2), std="margin", main="VALIDATION LRM MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  #fourfold(cross.classification.validation.lrm, std="margin", main="VALIDATION LRM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  dev.off()
  
  #pdf(file = "result_LRM_ROCPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.validation.lrm, main = "ROC PLOT: VALIDATION LRM MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.validation.lrm<-roc.area(verification.table[,2],result.lrm.validation$fitted.values)
  #dev.off()
  
  pdf(file = "result_LRM_ROCPlot_bootstrap_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.validation.lrm, main = "ROC PLOT: VALIDATION LRM MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[3] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.validation.lrm$A,2),";  Sample size = ",area.under.roc.curve.validation.lrm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[3], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  # Assignation of a matching code between observed and predicted values calculated using the validation dataset
  validation.lrm.matching.code<-paste(validation.grouping.variable,round(predict.result.lrm.validation.posterior),sep="")
  validation.lrm.matching.code<-gsub("00","1",validation.lrm.matching.code)
  validation.lrm.matching.code<-gsub("01","2",validation.lrm.matching.code)
  validation.lrm.matching.code<-gsub("10","3",validation.lrm.matching.code)
  validation.lrm.matching.code<-gsub("11","4",validation.lrm.matching.code)
  validation.lrm.matching.code<-as.numeric(validation.lrm.matching.code)
  
  ##########################################
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    cross.classification.validation.lrm.optimal<-table(validation.grouping.variable,as.numeric(predict.result.lrm.validation.posterior>lrm.probability.optimal.binary.threshold),dnn=c("Observed","Predicted"))
    rownames(cross.classification.validation.lrm.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.validation.lrm.optimal)<-list("No Landslide","Landslide") # Predicted
    
    #Different plots for contingency table
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.validation.lrm.optimal/sum(cross.classification.validation.lrm.optimal)*100,2), std="margin", main="VALIDATION LRM MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
      #fourfold(cross.classification.validation.lrm, std="margin", main="VALIDATION LRM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    }
    if (cross.classification.validation.lrm.optimal[1,2]==0 | cross.classification.validation.lrm.optimal[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
    pdf(file = "result_LRM_FourfoldPlot_Validation_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.validation.lrm.optimal/sum(cross.classification.validation.lrm.optimal)*100,2), std="margin", main="VALIDATION LRM MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.lrm, std="margin", main="VALIDATION LRM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    dev.off()
    
    
    # Assignation of a optimal matching code between observed and predicted values calculated using the validation dataset
    validation.lrm.matching.code.optimal<-paste(validation.grouping.variable,as.numeric(predict.result.lrm.validation.posterior>lrm.probability.optimal.binary.threshold),sep="")
    validation.lrm.matching.code.optimal<-gsub("00","1",validation.lrm.matching.code.optimal)
    validation.lrm.matching.code.optimal<-gsub("01","2",validation.lrm.matching.code.optimal)
    validation.lrm.matching.code.optimal<-gsub("10","3",validation.lrm.matching.code.optimal)
    validation.lrm.matching.code.optimal<-gsub("11","4",validation.lrm.matching.code.optimal)
    validation.lrm.matching.code.optimal<-as.numeric(validation.lrm.matching.code.optimal)
  }
  #########################################
  
  
  # EXPORT OF LRM MODEL RESULTS
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.lrm$table[,1,],contingency.table.lrm$table[,2,],contingency.table.lrm$table[,3,])),file="result_LRM_slu.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind("rocarea",round(area.under.roc.curve.lrm$A,4)),file="result_LRM_slu.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("RESULTS OF LOGISTIC REGRESSION MODEL",file="result_LRM.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("LRM MODEL OUTPUTS",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  #write.table("Logistic Regression coefficients",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  ##Scaling coefficients
  ##write.table(cbind(names(result.lrm$result$coefficients),result.lrm$result$coefficients),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  #write.table(cbind(names(result.lrm$coefficients),result.lrm$coefficients),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("LRM results summary table",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(t(c("variable",colnames(summary(result.lrm)$coefficients))),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(data.frame(cbind(names(result.lrm$coefficients),summary(result.lrm)$coefficients),stringsAsFactors=FALSE),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE MODEL RESULT",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.lrm$table[,1,],contingency.table.lrm$table[,2,],contingency.table.lrm$table[,3,])),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE VALIDATION",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.validation.lrm$table[,1,],contingency.table.validation.lrm$table[,2,],contingency.table.validation.lrm$table[,3,])),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("MATCHING CODE DEFINITION",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(c("","OBSERVED NO LANDSLIDES: 0","OBSERVED LANDSLIDES: 1"), c("PREDICTED NO LANDSLIDES: 0","00 -> Code 1","10 -> Code 3"), c("PREDICTED LANDSLIDES: 1","01 -> Code 2","11 -> Code 4")),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  
  
  ########
  if(enable_probability_optimal_binary_classification==FALSE) 
  {
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","LRM MODEL POSTERIOR PROBABILITY","LRM MODEL CLASSIFICATION","LRM MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict(result.lrm, type="response"),round(predict(result.lrm, type="response")),result.lrm.matching.code)),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","LRM VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.lrm.validation.posterior,round(predict.result.lrm.validation.posterior),validation.lrm.matching.code)),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } else
  {
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    if(enable_probability_optimal_classification==TRUE) 
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY PARTITION -> Method: ",type_probability_optimal_classification,sep=""),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(lrm.probability.classification.optimal[index_points_plot,c("prob_thres","tnr","tpr")]),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
      #ivan
      write.table(paste("OPTIMAL probability binary threshold: ",lrm.probability.optimal.binary.threshold,sep=""),file="optimalthreshold_LRM2.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      ##fineivan
    } else
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY BINARY PARTITION",sep=""),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(lrm.probability.classification.optimal[lrm.optimal.index,c("prob_thres","tnr","tpr")]),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    }
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","LRM MODEL POSTERIOR PROBABILITY","LRM MODEL CLASSIFICATION","LRM MODEL RESULT MATCHING CODE","LRM OPTIMAL MODEL CLASSIFICATION","LRM OPTIMAL MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict(result.lrm, type="response"),round(predict(result.lrm, type="response")),result.lrm.matching.code,as.numeric(predict(result.lrm, type="response")>lrm.probability.optimal.binary.threshold),result.lrm.matching.code.optimal)),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","LRM VALIDATION MATCHING CODE","OPTIMAL VALIDATION CLASSIFICATION","OPTIMAL LRM VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.lrm.validation.posterior,round(predict.result.lrm.validation.posterior),validation.lrm.matching.code,as.numeric(predict.result.lrm.validation.posterior>lrm.probability.optimal.binary.threshold),validation.lrm.matching.code.optimal)),file="result_LRM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  }
  ########
  
  # PLOT AND EXPORT OF LRM MODEL MAPS
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POLYGONS")
  {
    ################################################
    ### Prima di cancellare testare con polygoni
    ##############################################
    #shape_training_lrm<-shape_training
    ##result_training_lrm_shape<-cbind(identification.value,result.lrm$result$y,result.lrm$result$fitted.values,round(result.lrm$result$fitted.values),result.lrm.matching.code)
    #result_training_lrm_shape<-cbind(identification.value,training.table[,2],predict(result.lrm, type="response"),round(predict(result.lrm, type="response")),result.lrm.matching.code)
    #colnames(result_training_lrm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH")
    #shape_training_lrm@data <- merge(x=shape_training_lrm@data,y=result_training_lrm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    ##writeOGR(shape_training_lrm,dsn="result_LRM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_training_lrm,dsn="result_LRM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    #shape_validation_lrm<-shape_validation
    #result_validation_lrm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.lrm.validation.posterior,round(predict.result.lrm.validation.posterior),validation.lrm.matching.code)
    #colnames(result_validation_lrm_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    #shape_validation_lrm@data <- merge(x=shape_validation_lrm@data,y=result_validation_lrm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    ##writeOGR(shape_validation_lrm,dsn="result_LRM_validation.shp",layer="validation",driver="ESRI Shapefile")  # Version of rgdal older than 2.13.1
    #writeOGR(shape_validation_lrm,dsn="result_LRM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ### DA qui Nuovo
    ##################################################################   
    #### aggiungere a poligoni colonna incertezza ed esportazione pdf incertezza
    shape_training_lrm<-shape_training
    #result_training_lrm_shape<-cbind(identification.value,result.lrm$result$y,result.lrm$result$fitted.values,round(result.lrm$result$fitted.values),result.lrm.matching.code)
    result_training_lrm_shape<-cbind(identification.value,training.table[,2],predict(result.lrm, type="response"),round(predict(result.lrm, type="response")),result.lrm.matching.code,ID.bootstrap.model.lrm.count,bootstrap.model.lrm.probability.mean,bootstrap.model.lrm.probability.sd,bootstrap.model.lrm.probability.min,bootstrap.model.lrm.probability.max,bootstrap.model.lrm.probability.sderror,t(bootstrap.model.lrm.probability.quantiles),bootstrap.model.lrm.prediction.mean,bootstrap.model.lrm.prediction.sd,bootstrap.model.lrm.prediction.min,bootstrap.model.lrm.prediction.max,bootstrap.model.lrm.prediction.sderror,t(bootstrap.model.lrm.prediction.quantiles))
    colnames(result_training_lrm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LRM_SAMP","LRM_PMean","LRM_PSd","LRM_PMin","LRM_PMax","LRM_PSder","LRM_PQ0","LRM_PQ_005","LRM_PQ_025","LRM_PQ05","LRM_PQ_075","LRM_PQ095","LRM_PQ1","LRM_PrMean","LRM_PrSd","LRM_PrMin","LRM_PrMax","LRM_PrSder","LRM_PrQ0","LRM_PrQ005","LRM_PrQ025","LRM_PrQ05","LRM_PrQ075","LRM_PrQ095","LRM_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_lrm_shape<-cbind(result_training_lrm_shape,as.numeric(predict(result.lrm, type="response")>lrm.probability.optimal.binary.threshold),result.lrm.matching.code.optimal)
      colnames(result_training_lrm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LRM_SAMP","LRM_PMean","LRM_PSd","LRM_PMin","LRM_PMax","LRM_PSder","LRM_PQ0","LRM_PQ_005","LRM_PQ_025","LRM_PQ05","LRM_PQ_075","LRM_PQ095","LRM_PQ1","LRM_PrMean","LRM_PrSd","LRM_PrMin","LRM_PrMax","LRM_PrSder","LRM_PrQ0","LRM_PrQ005","LRM_PrQ025","LRM_PrQ05","LRM_PrQ075","LRM_PrQ095","LRM_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    
    
    shape_training_lrm@data <- merge(x=shape_training_lrm@data,y=result_training_lrm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_lrm,dsn="result_LRM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_lrm,dsn="result_LRM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    # WARNING: The validation does't have the unvertainty estimation: probably this can be daone using the parabolic error function 
    shape_validation_lrm<-shape_validation
    result_validation_lrm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.lrm.validation.posterior,round(predict.result.lrm.validation.posterior),validation.lrm.matching.code)
    colnames(result_validation_lrm_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_lrm_shape<-cbind(result_validation_lrm_shape,as.numeric(predict.result.lrm.validation.posterior>lrm.probability.optimal.binary.threshold),validation.lrm.matching.code.optimal)
      colnames(result_validation_lrm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    
    
    shape_validation_lrm@data <- merge(x=shape_validation_lrm@data,y=result_validation_lrm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_lrm@data <- cbind(shape_validation_lrm@data,PROB_SDMOD=(coefficients(fit.parabola.probability.lrm)*(shape_validation_lrm@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.lrm)*shape_validation_lrm@data$VAL_PROB))
    #writeOGR(shape_validation_lrm,dsn="result_LRM_validation.shp",layer="validation",driver="ESRI Shapefile")  # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_lrm,dsn="result_LRM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ##################################################################
    
    
    # Plot and export of maps
    # LRM Susceptibility
    #dev.new()
    pdf(file = "result_LRM_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_lrm, zcol=c("MOD_PROB"), names.attr=c("LRM MODEL PROBABILITY"), main="LRM MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # LRM Model Matching Code
    #dev.new()
    pdf(file = "result_LRM_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_lrm, zcol=c("MOD_MATCH"), names.attr=c("LRM MODEL MATCHING CODE"), main="LRM MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
    
    # LRM Model Uncertainty
    #dev.new()
    pdf(file = "result_LRM_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_lrm, zcol=c("LRM_PrSd"), names.attr=c("LRM MODEL UNCERTAINTY"), main="LRM MODEL UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()
    
    # LRM Validation Susceptibility
    #dev.new()
    pdf(file = "result_LRM_Validation_Susceptibility_Map.pdf")
    print(spplot(obj=shape_validation_lrm, zcol=c("VAL_PROB"), names.attr=c("LRM VALIDATION PROBABILITY"), main="LRM VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # LRM Validaiton Matching Code
    #dev.new()
    pdf(file = "result_LRM_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_lrm, zcol=c("VAL_MATCH"), names.attr=c("LRM VALIDATION MATCHING CODE"), main="LRM VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
        
    # LRM Validation Uncertainty
    #dev.new()
    pdf(file = "result_LRM_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_lrm, zcol=c("PROB_SDMOD"), names.attr=c("LRM VALIDATION UNCERTAINTY"), main="LRM VALIDATION UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()

    ########### TBT
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      # LRM Model Matching Code
      #dev.new()
      pdf(file = "result_LRM_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_training_lrm, zcol=c("OPT_MATCH"), names.attr=c("LRM MODEL MATCHING CODE"), main="LRM MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      # LRM Model Matching Code
      #dev.new()
      pdf(file = "result_LRM_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_validation_lrm, zcol=c("OPT_MATCH"), names.attr=c("LRM VALIDATION MATCHING CODE"), main="LRM VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      if(enable_probability_optimal_classification==TRUE)
      {
        pdf(file = "result_LRM_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
        print(spplot(obj=shape_training_lrm, zcol=c("MOD_PROB"), names.attr=c("LRM MODEL PROBABILITY"), main="LRM MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=(lrm.breaks.histogram.values.optimal)+c(rep(0,(length(lrm.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1)))
        dev.off()
        
        pdf(file = "result_LRM_Validation_Susceptibility_Map_Optimal.pdf")
        print(spplot(obj=shape_validation_lrm, zcol=c("VAL_PROB"), names.attr=c("LRM VALIDATION PROBABILITY"), main="LRM VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=(lrm.breaks.histogram.values.optimal)+c(rep(0,(length(lrm.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1)))
        dev.off()
      }
    }
    ###########
    
    
    ### Success & prediction rate curve
    susc_values_lrm<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_lrm@data) %in% c(configuration.spatial.data.table[c(5,6)],"MOD_PROB"))
    shape_training_lrm@data<-shape_training_lrm@data[order(shape_training_lrm@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_training_lrm@data[,1]<-cumsum(shape_training_lrm@data[,1])/sum(shape_training_lrm@data[,1])*100
    shape_training_lrm@data[,2]<-cumsum(shape_training_lrm@data[,2])/sum(shape_training_lrm@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_lrm@data[,3],shape_training_lrm@data[,1])
    area_susc_values_lrm_mod<-c(0,ind_pro_fun_mod(susc_values_lrm[-c(1,length(susc_values_lrm))]),100)
    
    #dev.new()
    pdf(file = "result_LRM_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="LRM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_lrm)-1))
    {
      #count=1
      polygon(c(area_susc_values_lrm_mod[count:(count+1)],rev(area_susc_values_lrm_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_lrm@data[,1],shape_training_lrm@data[,2],col="black")
    dev.off()
    
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_lrm_mod_optimal<-c(0,ind_pro_fun_mod(rev(lrm.breaks.histogram.values.optimal)[-c(1,length(lrm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_LRM_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="LRM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(lrm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_lrm_mod_optimal[count:(count+1)],rev(area_susc_values_lrm_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_lrm@data[,1],shape_training_lrm@data[,2],col="black")
      dev.off()
    }
    ############
    
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_lrm@data) %in% c(configuration.spatial.data.table[c(5,6)],"VAL_PROB"))
    shape_validation_lrm@data<-shape_validation_lrm@data[order(shape_validation_lrm@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_validation_lrm@data[,1]<-cumsum(shape_validation_lrm@data[,1])/sum(shape_validation_lrm@data[,1])*100
    shape_validation_lrm@data[,2]<-cumsum(shape_validation_lrm@data[,2])/sum(shape_validation_lrm@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_lrm@data[,3],shape_validation_lrm@data[,1])
    area_susc_values_lrm_val<-c(0,ind_pro_fun_val(susc_values_lrm[-c(1,length(susc_values_lrm))]),100)
    
    #dev.new()
    pdf(file = "result_LRM_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="LRM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_lrm)-1))
    {
      #count=1
      polygon(c(area_susc_values_lrm_val[count:(count+1)],rev(area_susc_values_lrm_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_lrm@data[,1],shape_validation_lrm@data[,2],col="black")
    dev.off()
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_lrm_val_optimal<-c(0,ind_pro_fun_val(rev(lrm.breaks.histogram.values.optimal)[-c(1,length(lrm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_LRM_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="LRM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(lrm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_lrm_val_optimal[count:(count+1)],rev(area_susc_values_lrm_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_lrm@data[,1],shape_validation_lrm@data[,2],col="black")
      dev.off()
    }
    ############
    
  }
  
  
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POINTS")
  {
    shape_training_lrm<-shape_training
    #result_training_lrm_shape<-cbind(identification.value,training.table[,2],predict.result.lrm$posterior[,2],as.numeric(levels(predict.result.lrm$class))[predict.result.lrm$class],result.lrm.matching.code,ID.bootstrap.model.lrm.count,bootstrap.model.lrm.probability.mean,bootstrap.model.lrm.probability.sd,bootstrap.model.lrm.probability.min,bootstrap.model.lrm.probability.max,bootstrap.model.lrm.probability.sderror,t(bootstrap.model.lrm.probability.quantiles),bootstrap.model.lrm.prediction.mean,bootstrap.model.lrm.prediction.sd,bootstrap.model.lrm.prediction.min,bootstrap.model.lrm.prediction.max,bootstrap.model.lrm.prediction.sderror,t(bootstrap.model.lrm.prediction.quantiles))
    result_training_lrm_shape<-cbind(identification.value,training.table[,2],predict(result.lrm, type="response"),as.numeric(round(predict(result.lrm, type="response"))),result.lrm.matching.code,ID.bootstrap.model.lrm.count,bootstrap.model.lrm.probability.mean,bootstrap.model.lrm.probability.sd,bootstrap.model.lrm.probability.min,bootstrap.model.lrm.probability.max,bootstrap.model.lrm.probability.sderror,t(bootstrap.model.lrm.probability.quantiles),bootstrap.model.lrm.prediction.mean,bootstrap.model.lrm.prediction.sd,bootstrap.model.lrm.prediction.min,bootstrap.model.lrm.prediction.max,bootstrap.model.lrm.prediction.sderror,t(bootstrap.model.lrm.prediction.quantiles))
    colnames(result_training_lrm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LRM_SAMP","LRM_PMean","LRM_PSd","LRM_PMin","LRM_PMax","LRM_PSder","LRM_PQ0","LRM_PQ_005","LRM_PQ_025","LRM_PQ05","LRM_PQ_075","LRM_PQ095","LRM_PQ1","LRM_PrMean","LRM_PrSd","LRM_PrMin","LRM_PrMax","LRM_PrSder","LRM_PrQ0","LRM_PrQ005","LRM_PrQ025","LRM_PrQ05","LRM_PrQ075","LRM_PrQ095","LRM_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_lrm_shape<-cbind(result_training_lrm_shape,as.numeric(predict(result.lrm, type="response")>lrm.probability.optimal.binary.threshold),result.lrm.matching.code.optimal)
      colnames(result_training_lrm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","LRM_SAMP","LRM_PMean","LRM_PSd","LRM_PMin","LRM_PMax","LRM_PSder","LRM_PQ0","LRM_PQ_005","LRM_PQ_025","LRM_PQ05","LRM_PQ_075","LRM_PQ095","LRM_PQ1","LRM_PrMean","LRM_PrSd","LRM_PrMin","LRM_PrMax","LRM_PrSder","LRM_PrQ0","LRM_PrQ005","LRM_PrQ025","LRM_PrQ05","LRM_PrQ075","LRM_PrQ095","LRM_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    shape_training_lrm@data <- merge(x=shape_training_lrm@data,y=result_training_lrm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_lrm,dsn="result_LRM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    
    if(enable_detailed_data_export==TRUE) writeOGR(shape_training_lrm,dsn="result_LRM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    
    # WARNING: The validation does't have the unvertainty estimation: probabilty this can be daone using the parabolic error function 
    shape_validation_lrm<-shape_validation
    result_validation_lrm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.lrm.validation.posterior,round(predict.result.lrm.validation.posterior),validation.lrm.matching.code)
    colnames(result_validation_lrm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_lrm_shape<-cbind(result_validation_lrm_shape,as.numeric(predict.result.lrm.validation.posterior>lrm.probability.optimal.binary.threshold),validation.lrm.matching.code.optimal)
      colnames(result_validation_lrm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    shape_validation_lrm@data <- merge(x=shape_validation_lrm@data,y=result_validation_lrm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_lrm@data <- cbind(shape_validation_lrm@data,PROB_SDMOD=(coefficients(fit.parabola.probability.lrm)*(shape_validation_lrm@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.lrm)*shape_validation_lrm@data$VAL_PROB))
    #writeOGR(shape_validation_lrm,dsn="result_LRM_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    if(enable_detailed_data_export==TRUE) if(enable_detailed_data_export==TRUE) writeOGR(shape_validation_lrm,dsn="result_LRM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    require(raster)
    
    # Plot and export of maps
    # LRM Susceptibility
    #dev.new()
    pdf(file = "result_LRM_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_lrm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_LRM_Model_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_LRM_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1),breaks=round(lrm.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    if(enable_detailed_data_export==TRUE)
    {
    # LRM Model Matching Code
    #dev.new()
    pdf(file = "result_LRM_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_lrm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
    print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    
    writeRaster(layer_gridded_raster, filename="result_LRM_Model_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    }
    
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_LRM_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_lrm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
      print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_LRM_Model_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_LRM_Model_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_lrm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_LRM_Model_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
    }
    ############
    if(enable_detailed_data_export==TRUE)
    {
    # LRM Model uncertainity
    #dev.new()
    pdf(file = "result_LRM_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_lrm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"LRM_PrSd"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_LRM_Model_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    }
    
    # LRM Validation Susceptibility
    #dev.new()
    pdf(file = "result_LRM_Validation_Susceptibility_Map.pdf")
    layer_gridded<-shape_validation_lrm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_LRM_Validation_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_LRM_Validation_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1),breaks=round(lrm.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    
    # LRM Validation Matching Code
    #dev.new()
    pdf(file = "result_LRM_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_lrm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    #zoom(layer_gridded_raster)  	
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_LRM_Validation_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)

    ########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_LRM_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_lrm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      #zoom(layer_gridded_raster)  	
      dev.off()
      writeRaster(layer_gridded_raster, filename="result_LRM_Validation_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_LRM_Validation_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_lrm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_LRM_Validation_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
      
    }
    #######
    
    if(enable_detailed_data_export==TRUE) 
    {
    # LRM Model uncertainity
    #dev.new()
    pdf(file = "result_LRM_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_lrm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"PROB_SDMOD"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_LRM_Validation_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    }
    
    ### Success & prediction rate curve
    susc_values_lrm<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_lrm@data) %in% c("MOD_GROUP","MOD_PROB"))
    shape_training_lrm@data<-cbind(PIXEL_AREA=rep(configuration.spatial.data.table[c(8)]^2,dim(shape_training_lrm@data)[1]),shape_training_lrm@data[order(shape_training_lrm@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_training_lrm@data$MOD_GROUP<-shape_training_lrm@data$MOD_GROUP*configuration.spatial.data.table[c(8)]^2
    shape_training_lrm@data[,1]<-cumsum(shape_training_lrm@data[,1])/sum(shape_training_lrm@data[,1])*100
    shape_training_lrm@data[,2]<-cumsum(shape_training_lrm@data[,2])/sum(shape_training_lrm@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_lrm@data[,3],shape_training_lrm@data[,1])
    area_susc_values_lrm_mod<-c(0,ind_pro_fun_mod(susc_values_lrm[-c(1,length(susc_values_lrm))]),100)
    
    #dev.new()
    pdf(file = "result_LRM_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="LRM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_lrm)-1))
    {
      #count=1
      polygon(c(area_susc_values_lrm_mod[count:(count+1)],rev(area_susc_values_lrm_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_lrm@data[,1],shape_training_lrm@data[,2],col="black")
    dev.off()
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_lrm_mod_optimal<-c(0,ind_pro_fun_mod(rev(lrm.breaks.histogram.values.optimal)[-c(1,length(lrm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_LRM_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="LRM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(lrm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_lrm_mod_optimal[count:(count+1)],rev(area_susc_values_lrm_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_lrm@data[,1],shape_training_lrm@data[,2],col="black")
      dev.off()
    }
    ############
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_lrm@data) %in% c("VAL_GROUP","VAL_PROB"))
    shape_validation_lrm@data<-cbind(PIXEL_AREA=rep(as.numeric(configuration.spatial.data.table[c(8)])^2,dim(shape_validation_lrm@data)[1]),shape_validation_lrm@data[order(shape_validation_lrm@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_validation_lrm@data$VAL_GROUP<-shape_validation_lrm@data$VAL_GROUP*as.numeric(configuration.spatial.data.table[c(8)])^2
    shape_validation_lrm@data[,1]<-cumsum(shape_validation_lrm@data[,1])/sum(shape_validation_lrm@data[,1])*100
    shape_validation_lrm@data[,2]<-cumsum(shape_validation_lrm@data[,2])/sum(shape_validation_lrm@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_lrm@data[,3],shape_validation_lrm@data[,1])
    area_susc_values_lrm_val<-c(0,ind_pro_fun_val(susc_values_lrm[-c(1,length(susc_values_lrm))]),100)
    
    #dev.new()
    pdf(file = "result_LRM_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="LRM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_lrm)-1))
    {
      #count=1
      polygon(c(area_susc_values_lrm_val[count:(count+1)],rev(area_susc_values_lrm_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_lrm@data[,1],shape_validation_lrm@data[,2],col="black")
    dev.off()
    
    ###########
    
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_lrm_val_optimal<-c(0,ind_pro_fun_val(rev(lrm.breaks.histogram.values.optimal)[-c(1,length(lrm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_LRM_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="LRM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(lrm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_lrm_val_optimal[count:(count+1)],rev(area_susc_values_lrm_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(lrm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_lrm@data[,1],shape_validation_lrm@data[,2],col="black")
      dev.off()
    }
    ############
  }
}
 

#------------------ NEURAL NETWORK MODEL ANALISYS -------------------#

if(model.run.matrix[4] == "YES")
  {
  library(nnet)
  
  
  if (analysis.parameter.matrix[4] == "OPT")
  {
    # This "for" cicle search the best number of Weight Decay and Hidden Layer Nodes basing on the minimization of Sum of Squared Error (SSE) or on the mamixization of ROC area results  
    # INITIATE A NULL TABLE
    nnm.table <- NULL
    total.iteration.number<-3*length(round((length(explanatory.variables)/2)):(length(explanatory.variables)-2))*10
    iteration.number<-0
    
    
    # SEARCH FOR OPTIMAL WEIGHT DECAY WITH RANGE OF WEIGHT DECAYS SUGGESTED BY B. RIPLEY
    for (weight.decay in c(0.0001, 0.001, 0.01))
    {
      # SEARCH FOR OPTIMAL NUMBER OF HIDDEN UNITS
      for (n.nodes in round((length(explanatory.variables)/2)):(length(explanatory.variables)-2))
      {
        # UNITIATE A NULL VECTOR
        sse <- NULL
        # FOR EACH SETTING, RUN NEURAL NET MULTIPLE TIMES
        for (i.counts in 1:10)
        {
          # INITIATE THE RANDOM STATE FOR EACH NET
          set.seed(i.counts)
          # TRAIN NEURAL NETS
          result.nnm <- nnet(explanatory.variables, training.table[,2], size = n.nodes, rang = 0.5, maxit = 300, MaxNWts=10000, decay = weight.decay, skip = FALSE, trace = TRUE,reltol =0.00001) # Default reltol 10^-8: greater reltol -> fast convergence 
          # CALCULATE SSE (Sum of Squared Error) and ROC.area
          test.sse <- sum(((as.numeric(grouping.variable)-1) - round(predict(result.nnm)))^2)
          library(verification)
          test.ROC.area <- (roc.area((as.numeric(grouping.variable)-1),result.nnm$fitted.values))$A
          
          iteration.number<-iteration.number+1
          print(paste("Iteration",iteration.number,"of",total.iteration.number,"-",round((iteration.number/total.iteration.number*100),1),"%","completed"))
          
          # APPEND EACH SSE and ROC.area TO A VECTOR
          if (i.counts == 1) sse <- test.sse else sse <- rbind(sse, test.sse)
          if (i.counts == 1) ROC.area <- test.ROC.area else ROC.area <- rbind(ROC.area, test.ROC.area)
        }
        # APPEND AVERAGED SSE and AVERAGED ROC.area WITH RELATED PARAMETERS TO A TABLE
        nnm.table <- rbind(nnm.table, c(WEIGHT_DECAY = weight.decay, HYDDEN_LAYER_NODES = n.nodes, SUM_SQUARED_ERROR = mean(sse), ROC_AREA = mean(ROC.area)))
      }
    }
    # PRINT OUT THE RESULT
    print(nnm.table)
    
    # Extracting value of Weight Decay and Number of nodes that minimize the SSE
    pos.sse.min<-which.min(nnm.table[,3])
    weight.decay.sse.min<-nnm.table[pos.sse.min,1]
    n.nodes.sse.min<-nnm.table[pos.sse.min,2]
    sse.min<-min(nnm.table[,3])
    
    
    # Extracting value of Weight Decay and Number of nodes that maximize the ROC Area
    pos.roc.area.max<-which.max(nnm.table[,4])
    weight.decay.roc.area.max<-nnm.table[pos.roc.area.max,1]
    n.nodes.roc.area.max<-nnm.table[pos.roc.area.max,2]
    roc.area.max<-max(nnm.table[,4])
    
    n.nodes.selected <- n.nodes.sse.min
    weight.decay.selected <- weight.decay.sse.min
    #n.nodes.selected <- n.nodes.roc.area.max
    #weight.decay.selected <- weight.decay.roc.area.max
  }
  
  if (analysis.parameter.matrix[4] == "NOR")
  {
    # Manual choice 
    n.nodes.selected <-round(dim(explanatory.variables)[2]/2)
    weight.decay.selected <- 0 # default value nnet 0
  }
  
  
  
  if (class(try(nnet(explanatory.variables, training.table[,2], size = n.nodes.selected, rang = 0.5, maxit = 100, MaxNWts=10000, decay = weight.decay.selected, skip = FALSE, trace = TRUE, reltol =0.000001)))=="try-error") # Default reltol 10^-8: greater reltol -> fast convergence
  { 
    #nnet(explanatory.variables, training.table[,2], size = n.nodes.selected, rang = 0.5, maxit = 10000, MaxNWts=100, decay = weight.decay.selected, skip = FALSE, trace = TRUE)
    write.table("Analysis based on Neural Network Model was not completed",file="Error_NNM_Analysis.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="Error_NNM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("Error LOG",file="Error_NNM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind("Message",rev(1:length(as.vector(.Traceback)))," ->",as.vector(.Traceback)),file="Error_NNM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  }
  
  result.nnm<-NULL
  set.seed(seed.value)
  result.nnm<-nnet(explanatory.variables, training.table[,2], size = n.nodes.selected, rang = 0.5, maxit = 10000, MaxNWts=10000, decay = weight.decay.selected, skip = FALSE, trace = TRUE, reltol =0.000001) # Default reltol 10^-8: greater reltol -> fast convergence 
  #names(result.nnm)
  
  # Result Predicted
  predict.result.nnm<-predict(result.nnm)
  str(predict.result.nnm)
  
  # As predicted values also result.nnm$fitted.values can be considered because it corresponds to predict.result.nnm as the sum of the difference of the two vectors two is 0
  # sum(result.nnm$fitted.values-predict.result.nnm)
  
  
  cross.classification.nnm<-table((as.numeric(grouping.variable)-1),round(predict.result.nnm),dnn=c("Observed","Predicted"))
  rownames(cross.classification.nnm)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.nnm)<-list("No Landslide","Landslide") # Predicted    
  str(cross.classification.nnm)
  
  # Assignation of a matching code between observed and predicted values
  result.nnm.matching.code<-paste(grouping.variable,round(predict.result.nnm),sep="")
  result.nnm.matching.code<-gsub("00","1",result.nnm.matching.code)
  result.nnm.matching.code<-gsub("01","2",result.nnm.matching.code)
  result.nnm.matching.code<-gsub("10","3",result.nnm.matching.code)
  result.nnm.matching.code<-gsub("11","4",result.nnm.matching.code)
  result.nnm.matching.code<-as.numeric(result.nnm.matching.code)
  
  #Elaboration of Coefficient of association for contingency table 
  #load package (vcd)  
  library(vcd)
  
  #help(package=vcd)         
  contingency.table.nnm<-table2d_summary(cross.classification.nnm)
  test.table.nnm<-assocstats(cross.classification.nnm)
  
  
  #Different plots for contingency table 
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.nnm/sum(cross.classification.nnm)*100,2),std="margin", main="NEURAL NETWORK MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.nnm,std="margin", main="NEURAL NETWORK MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #A ROC curve plots the false alarm rate against the hit rate
  #for a probablistic forecast for a range of thresholds. 
  
  #load package (verification)  
  library(verification)
  
  #verify function
  #Based on the type of inputs, this function calculates a range of verification statistics and skill scores.
  #Additionally, it creates a verify class object that can be further analyzed.
  
  ##### ROC PLOT OBS - POSTERIOR PROBABILITY ASSOCIATED TO 1
  
  # Method using verify function
  verification.results.nnm<-verify((as.numeric(grouping.variable)-1),predict.result.nnm, frcst.type="prob", obs.type="binary")
  #summary(verification.results.nnm)
  
  #str(verification.results.nnm)
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(verification.results.nnm, main = "ROC PLOT: NEURAL NETWORK MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #}
  area.under.roc.curve.nnm<-roc.area((as.numeric(grouping.variable)-1),predict.result.nnm)
  
  ## showing confidence intervals.  MAY BE SLOW
  
  
  if (cross.classification.nnm[1,2]==0 | cross.classification.nnm[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.results.nnm, main = "ROC PLOT: NEURAL NETWORK MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[4] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.nnm$A,2),";  Sample size = ",area.under.roc.curve.nnm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[4], sep=""), side=3, col="red", cex=0.8)
    ## Histogram of posterior probability
    dev.new()                            
    hist(predict.result.nnm, breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Neural Network Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  }
  pdf(file = "result_NNM_Histogram.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  hist(predict.result.nnm, breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Neural Network Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  dev.off() 
  
  # EXPORT OF PLOT FOR NNM MODEL
  
  pdf(file = "result_NNM_FourfoldPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.nnm/sum(cross.classification.nnm)*100,2),std="margin", main="NEURAL NETWORK MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  #fourfold(cross.classification.nnm,std="margin", main="NEURAL NETWORK MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  dev.off()
  
  #pdf(file = "result_NNM_ROCPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.results.nnm, main = "ROC PLOT: NEURAL NETWORK MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #dev.off()
  
  pdf(file = "result_NNM_ROCPlot_bootstrap.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.results.nnm, main = "ROC PLOT: NEURAL NETWORK MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[4] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.nnm$A,2),";  Sample size = ",area.under.roc.curve.nnm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[4], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  ## BOOTSTRAP PROCEDURE FOR THE ESTIMATION OF MODEL PREDICTION VARIABILITY
  if(bootstrap.model.variability[4] == "YES")
  {
    bootstrap.sample.model.nnm<-bootstrap.sample.model[4]
    
    matrix.bootstrap.model.nnm<-matrix(data=NA, nrow=dim(training.table)[1], ncol=(bootstrap.sample.model.nnm*3)+1)
    colnames(matrix.bootstrap.model.nnm)<-rep("na",(bootstrap.sample.model.nnm*3)+1)
    matrix.bootstrap.model.nnm[,1]<-identification.value
    colnames(matrix.bootstrap.model.nnm)[1]<-"ID"
    name.sel.run<-paste(rep("ID_Selection_Run",bootstrap.sample.model.nnm),1:bootstrap.sample.model.nnm,sep="_")
    colnames(matrix.bootstrap.model.nnm)[seq(2,(bootstrap.sample.model.nnm*3)-1,3)]<-name.sel.run
    name.prob.run<-paste(rep("Probability_Run",bootstrap.sample.model.nnm),1:bootstrap.sample.model.nnm,sep="_")
    colnames(matrix.bootstrap.model.nnm)[seq(3,(bootstrap.sample.model.nnm*3),3)]<-name.prob.run
    name.pred.run<-paste(rep("Prediction_Run",bootstrap.sample.model.nnm),1:bootstrap.sample.model.nnm,sep="_")
    colnames(matrix.bootstrap.model.nnm)[seq(4,(bootstrap.sample.model.nnm*3)+1,3)]<-name.pred.run
    
    selection.index<-NULL
    library(nnet)
    #Bootstrap procedure
    for (count.boot in 1:bootstrap.sample.model.nnm)
    {
      selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
      matrix.bootstrap.model.nnm[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
      explanatory.variables.bootstrap.model.nnm<-training.table[selection.index,3:dim(training.table)[2]]
      grouping.variable.bootstrap.model.nnm<-training.table[selection.index,2]
      #result.bootstrap.model.nnm<-nnet(explanatory.variables.bootstrap.model.nnm, grouping.variable.bootstrap.model.nnm, size = n.nodes.selected, rang = 0.5, maxit = 10000, MaxNWts=10000, decay = weight.decay.selected, skip = FALSE, trace = TRUE)
      while(inherits(try(result.bootstrap.model.nnm<-nnet(explanatory.variables.bootstrap.model.nnm, grouping.variable.bootstrap.model.nnm, size = n.nodes.selected, rang = 0.5, maxit = 10000, MaxNWts=10000, decay = weight.decay.selected, skip = FALSE, trace = TRUE, reltol =0.000001),silent=TRUE),what="try-error"))  # # Default reltol 10^-8: greater reltol -> fast convergence 
      {
        print(paste("Count boot: ",count.boot," - Boostrap while resampling",sep=""))
        selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
        matrix.bootstrap.model.nnm[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
        explanatory.variables.bootstrap.model.nnm<-training.table[selection.index,3:dim(training.table)[2]]
        if(bootstrap_constant_correction==TRUE)
        {
          print("Performing bootstrap 0 value correction")
          indexbootstrapcosntant<-which(explanatory.variables.bootstrap.model.nnm==0,arr.ind=TRUE)
          explanatory.variables.bootstrap.model.nnm[indexbootstrapcosntant]<-runif(length(indexbootstrapcosntant)/2, min = 0.00001, max = 0.01)
        }
        grouping.variable.bootstrap.model.nnm<-as.factor(training.table[selection.index,2])
      }
      matrix.bootstrap.model.nnm[,(count.boot*3)+1]<-predict(result.bootstrap.model.nnm,newdata=explanatory.variables)
      matrix.bootstrap.model.nnm[as.numeric(names(table(selection.index))),(count.boot*3)]<-matrix.bootstrap.model.nnm[as.numeric(names(table(selection.index))),(count.boot*3)+1]
    }
    # Export of bootstrap sample
    write.table(matrix.bootstrap.model.nnm,file="result_NNM_BootstrapSamples.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    
    ID.bootstrap.model.nnm.count<-numeric(length=dim(training.table)[1])
    #Probability (selected values)
    bootstrap.model.nnm.probability.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.probability.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.probability.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.probability.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.probability.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.probability.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #Prediction (all values)
    bootstrap.model.nnm.prediction.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.prediction.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.prediction.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.prediction.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.prediction.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.nnm.prediction.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #    for (count.row.variability in 1:dim(training.table)[1])
    #        {
    #        # Statistics on boostrapped probability
    #        ID.bootstrap.model.nnm.count[count.row.variability]<-length(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(2,(bootstrap.sample.model.nnm*3)-1,3)]))
    #        bootstrap.model.nnm.probability.mean[count.row.variability]<-mean(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(3,(bootstrap.sample.model.nnm*3),3)]))
    #        bootstrap.model.nnm.probability.sd[count.row.variability]<-sd(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(3,(bootstrap.sample.model.nnm*3),3)]))
    #        bootstrap.model.nnm.probability.min[count.row.variability]<-min(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(3,(bootstrap.sample.model.nnm*3),3)]))
    #        bootstrap.model.nnm.probability.max[count.row.variability]<-max(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(3,(bootstrap.sample.model.nnm*3),3)]))
    #        bootstrap.model.nnm.probability.sderror[count.row.variability]<-bootstrap.model.nnm.probability.sd[count.row.variability]/ID.bootstrap.model.nnm.count[count.row.variability]
    #        bootstrap.model.nnm.probability.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(3,(bootstrap.sample.model.nnm*3),3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        # Statistics on boostrapped prediction
    #        bootstrap.model.nnm.prediction.mean[count.row.variability]<-mean(matrix.bootstrap.model.nnm[count.row.variability,seq(4,(bootstrap.sample.model.nnm*3)+1,3)])
    #        bootstrap.model.nnm.prediction.sd[count.row.variability]<-sd(matrix.bootstrap.model.nnm[count.row.variability,seq(4,(bootstrap.sample.model.nnm*3)+1,3)])
    #        bootstrap.model.nnm.prediction.min[count.row.variability]<-min(matrix.bootstrap.model.nnm[count.row.variability,seq(4,(bootstrap.sample.model.nnm*3)+1,3)])
    #        bootstrap.model.nnm.prediction.max[count.row.variability]<-max(matrix.bootstrap.model.nnm[count.row.variability,seq(4,(bootstrap.sample.model.nnm*3)+1,3)])
    #        bootstrap.model.nnm.prediction.sderror[count.row.variability]<-bootstrap.model.nnm.prediction.sd[count.row.variability]/bootstrap.sample.model.nnm
    #        bootstrap.model.nnm.prediction.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.nnm[count.row.variability,seq(4,(bootstrap.sample.model.nnm*3)+1,3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        }
    
    fun_length<-function(x) {length_var<-length(x[which(is.finite(x))]); return(length_var)}
    fun_quantile<-function(x) {quantile_var<-t(quantile(x[which(is.finite(x))],probs=c(0,0.05,0.25,0.5,0.75,0.95,1))); return(quantile_var)}
    ID.bootstrap.model.nnm.count<-apply(matrix.bootstrap.model.nnm[,grep("ID_Selection",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=fun_length)
    bootstrap.model.nnm.probability.mean<-apply(matrix.bootstrap.model.nnm[,grep("Probability",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.nnm.probability.sd<-apply(matrix.bootstrap.model.nnm[,grep("Probability",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.nnm.probability.min<-apply(matrix.bootstrap.model.nnm[,grep("Probability",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.nnm.probability.max<-apply(matrix.bootstrap.model.nnm[,grep("Probability",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.nnm.probability.sderror<-bootstrap.model.nnm.probability.sd/bootstrap.sample.model.nnm
    bootstrap.model.nnm.probability.quantiles<-apply(matrix.bootstrap.model.nnm[,grep("Probability",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=fun_quantile)
    bootstrap.model.nnm.prediction.mean<-apply(matrix.bootstrap.model.nnm[,grep("Prediction",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.nnm.prediction.sd<-apply(matrix.bootstrap.model.nnm[,grep("Prediction",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.nnm.prediction.min<-apply(matrix.bootstrap.model.nnm[,grep("Prediction",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.nnm.prediction.max<-apply(matrix.bootstrap.model.nnm[,grep("Prediction",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.nnm.prediction.sderror<-bootstrap.model.nnm.prediction.sd/bootstrap.sample.model.nnm
    bootstrap.model.nnm.prediction.quantiles<-apply(matrix.bootstrap.model.nnm[,grep("Prediction",colnames(matrix.bootstrap.model.nnm))],MARGIN=1,FUN=fun_quantile)
    
    
    # Export of bootstrap sample statistics
    write.table(cbind("ID","NNM_NumberSelectedSamples","NNM_Probability_Mean","NNM_Probability_Sd","NNM_Probability_Min","NNM_Probability_Max","NNM_Probability_Sderror","NNM_Probability_Quantiles_0","NNM_Probability_Quantiles_0.05","NNM_Probability_Quantiles_0.25","NNM_Probability_Quantiles_0.5","NNM_Probability_Quantiles_0.75","NNM_Probability_Quantiles_0.95","NNM_Probability_Quantiles_1","NNM_Prediction_Mean","NNM_Prediction_Sd","NNM_Prediction_Min","NNM_Prediction_Max","NNM_Prediction_Sderror","NNM_Prediction_Quantiles_0","NNM_Prediction_Quantiles_0.05","NNM_Prediction_Quantiles_0.25","NNM_Prediction_Quantiles_0.5","NNM_Prediction_Quantiles_0.75","NNM_Prediction_Quantiles_0.95","NNM_Prediction_Quantiles_1"),file="result_NNM_BootstrapStatistics.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind(identification.value,ID.bootstrap.model.nnm.count,bootstrap.model.nnm.probability.mean,bootstrap.model.nnm.probability.sd,bootstrap.model.nnm.probability.min,bootstrap.model.nnm.probability.max,bootstrap.model.nnm.probability.sderror,t(bootstrap.model.nnm.probability.quantiles),bootstrap.model.nnm.prediction.mean,bootstrap.model.nnm.prediction.sd,bootstrap.model.nnm.prediction.min,bootstrap.model.nnm.prediction.max,bootstrap.model.nnm.prediction.sderror,t(bootstrap.model.nnm.prediction.quantiles)),file="result_NNM_BootstrapStatistics.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    
    if (enable_screen_plotting==TRUE)
    {
      #dev.new()
      #double.sd.histogram.variability<-hist(bootstrap.model.nnm.probability.sd*2,breaks=seq(0,1,0.05),labels=TRUE)
      #plot(double.sd.histogram.variability$counts, seq(0,0.95,0.05), type="S",ylim=c(0,1), labels=TRUE)
      dev.new()
      plot(bootstrap.model.nnm.probability.mean,bootstrap.model.nnm.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="NNM BOOTSTRAP: Mean Probability vs Mean Prediction")
      abline(a=0,b=1,col="red",lty=1,lwd=1)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.nnm,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    }
    
    pdf(file = "result_NNM_BootstrapMeansComparison.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(bootstrap.model.nnm.probability.mean,bootstrap.model.nnm.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="NNM BOOTSTRAP: Mean Probability vs Mean Prediction")
    abline(a=0,b=1,col="red",lty=1,lwd=1)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.nnm,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    dev.off()
    
    
    # BOOTSTRAPPED PROBABILITY - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.probability.nnm<-cbind(bootstrap.model.nnm.probability.mean,2*bootstrap.model.nnm.probability.sd)
    parabola.probability.nnm<-na.omit(parabola.probability.nnm[order(parabola.probability.nnm[,1]),])
    colnames(parabola.probability.nnm)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.probability.nnm <- nls(parabola.probability.nnm[,"ordinate"] ~ coeff.a*(parabola.probability.nnm[,"abscissa"]^2) + (-1)*coeff.a*parabola.probability.nnm[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.probability.nnm<-predict(fit.parabola.probability.nnm)
    #coef(fit.parabola.probability.nnm)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.probability.nnm[,"abscissa"],parabola.probability.nnm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="NNM Model Probability Variability (Bootstrap)")
      lines(parabola.probability.nnm[,"abscissa"],value.parabola.probability.nnm,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.nnm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.nnm),3),coeff.b= -round(coef(fit.parabola.probability.nnm),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_NNM_BootstrapProbabilityVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.probability.nnm[,"abscissa"],parabola.probability.nnm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="NNM Model Probability Variability (Bootstrap)")
    lines(parabola.probability.nnm[,"abscissa"],value.parabola.probability.nnm,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.nnm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.nnm),3),coeff.b= -round(coef(fit.parabola.probability.nnm),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
    
    # BOOTSTRAPPED PREDICTION - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.prediction.nnm<-cbind(bootstrap.model.nnm.prediction.mean,2*bootstrap.model.nnm.prediction.sd)
    parabola.prediction.nnm<-parabola.prediction.nnm[order(parabola.prediction.nnm[,1]),]
    colnames(parabola.prediction.nnm)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.prediction.nnm <- nls(parabola.prediction.nnm[,"ordinate"] ~ coeff.a*(parabola.prediction.nnm[,"abscissa"]^2) + (-1)*coeff.a*parabola.prediction.nnm[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.prediction.nnm<-predict(fit.parabola.prediction.nnm)
    #coef(fit.parabola.prediction.nnm)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.prediction.nnm[,"abscissa"],parabola.prediction.nnm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="NNM Model Prediction Variability (Bootstrap)")
      lines(parabola.prediction.nnm[,"abscissa"],value.parabola.prediction.nnm,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.nnm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.nnm),3),coeff.b= -round(coef(fit.parabola.prediction.nnm),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    
    pdf(file = "result_NNM_BootstrapPredictionVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.prediction.nnm[,"abscissa"],parabola.prediction.nnm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="NNM Model Prediction Variability (Bootstrap)")
    lines(parabola.prediction.nnm[,"abscissa"],value.parabola.prediction.nnm,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.nnm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.nnm),3),coeff.b= -round(coef(fit.parabola.prediction.nnm),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
  }
  
  
  ## Sensitivity, Specificity, Cohens kappa plot
  dev.new()
  roc.plot.nnm.series<-roc.plot(verification.results.nnm,binormal=bin,plot=FALSE,show.thres=plot_thres,thresholds=threshold_series)
  dev.off()
  #str(roc.plot.nnm.series)
  #roc.plot.nnm.series$plot.data
  #str(roc.plot.nnm.series$plot.data)
  
  ###########################################
  # min(abs(TPR - (1-FPR)))
  if(enable_probability_optimal_binary_classification==FALSE)
  {
    dev.new()
    roc.plot.nnm.series<-roc.plot(verification.results.nnm,binormal=bin,plot=FALSE,show.thres=FALSE,thresholds=threshold_series)
    dev.off()
  } else
  {
    dev.new()
    roc.plot.nnm.series<-roc.plot(verification.results.nnm,binormal=bin,plot=FALSE,show.thres=FALSE)
    dev.off()  
  }
  
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    nnm.probability.classification.optimal<-data.frame(prob_thres=roc.plot.nnm.series$plot.data[,1,1],tpr=roc.plot.nnm.series$plot.data[,2,1],fpr=roc.plot.nnm.series$plot.data[,3,1],tnr=(1-roc.plot.nnm.series$plot.data[,3,1]),diff_abs_tpr_tnr=abs(roc.plot.nnm.series$plot.data[,2,1]-(1-roc.plot.nnm.series$plot.data[,3,1])),optimal_sel=NA,breaks_sel=NA)
    index.nnm.filter<-which(nnm.probability.classification.optimal$prob_thres>0 & nnm.probability.classification.optimal$prob_thres<1) # removing strnge thresh values
    nnm.probability.classification.optimal<-rbind(c(0,1,1,0,1,NA,NA),nnm.probability.classification.optimal[index.nnm.filter,],c(1,0,0,1,1,NA,NA))
    nnm.optimal.index<-which(nnm.probability.classification.optimal$diff_abs_tpr_tnr==min(nnm.probability.classification.optimal$diff_abs_tpr_tnr))
    nnm.probability.classification.optimal$optimal_sel[nnm.optimal.index]<-TRUE
    nnm.probability.optimal.binary.threshold<-nnm.probability.classification.optimal$prob_thres[nnm.optimal.index]
    
    ### Generating the optimal fourfold plot
    cross.classification.nnm.optimal<-table(grouping.variable,predict.result.nnm>nnm.probability.optimal.binary.threshold,dnn=c("Observed","Predicted"))
    rownames(cross.classification.nnm.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.nnm.optimal)<-list("No Landslide","Landslide") # Predicted    
    str(cross.classification.nnm.optimal)
    # Assignation of a matching code between observed and predicted values
    result.nnm.matching.code.optimal<-paste(grouping.variable,as.numeric(predict.result.nnm>nnm.probability.optimal.binary.threshold),sep="")
    result.nnm.matching.code.optimal<-gsub("00","1",result.nnm.matching.code.optimal)
    result.nnm.matching.code.optimal<-gsub("01","2",result.nnm.matching.code.optimal)
    result.nnm.matching.code.optimal<-gsub("10","3",result.nnm.matching.code.optimal)
    result.nnm.matching.code.optimal<-gsub("11","4",result.nnm.matching.code.optimal)
    result.nnm.matching.code.optimal<-as.numeric(result.nnm.matching.code.optimal)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.nnm.optimal/sum(cross.classification.nnm.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    }
    # EXPORT OF PLOT FOR NNM MODEL
    pdf(file = "result_NNM_FourfoldPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.nnm.optimal/sum(cross.classification.nnm.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.nnm.optimal, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    dev.off()
    
    ### Optimal susceptibility classes identification 
    if(enable_probability_optimal_classification==TRUE)
    {
      nnm.unexplained.errors<-round(1-(nnm.probability.classification.optimal$tpr[nnm.optimal.index]+nnm.probability.classification.optimal$tpr[nnm.optimal.index])/2,2)
      
      if(type_probability_optimal_classification=="proportional")
      {
        nnm.unexplained.errors.partition<-1-(nnm.unexplained.errors/((length(breaks.histogram.values)/2)))*(1:((length(breaks.histogram.values)/2)))
        
        for(count_part in 1:(length(nnm.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-nnm.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((nnm.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          nnm.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((nnm.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          nnm.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        nnm.breaks.histogram.values.optimal<-c(0,nnm.probability.classification.optimal$prob_thres[which(nnm.probability.classification.optimal$breaks_sel==TRUE)],1)
        #nnm.probability.classification.optimal[which(nnm.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if(type_probability_optimal_classification=="fixed")
      {
        step.nnm.unexplained.fixed<-0.1
        if(nnm.unexplained.errors<=0.1) step.nnm.unexplained.fixed<-0.05
        if(nnm.unexplained.errors<=0.05) step.nnm.unexplained.fixed<-0.025
        nnm.unexplained.errors.partition<-seq(step.nnm.unexplained.fixed,1-step.nnm.unexplained.fixed,step.nnm.unexplained.fixed)[seq(step.nnm.unexplained.fixed,1-step.nnm.unexplained.fixed,step.nnm.unexplained.fixed)>(1-nnm.unexplained.errors)]
        
        for(count_part in 1:(length(nnm.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-nnm.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((nnm.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          nnm.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((nnm.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          nnm.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        nnm.breaks.histogram.values.optimal<-c(0,nnm.probability.classification.optimal$prob_thres[which(nnm.probability.classification.optimal$breaks_sel==TRUE)],1)
        #nnm.probability.classification.optimal[which(nnm.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if (enable_screen_plotting==TRUE)
      {
        dev.new()
        hist(predict.result.nnm, breaks=nnm.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal NNM susceptibility", col=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1))
        plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("NNM OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
        mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
        mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
        for (count in 1:(length(nnm.breaks.histogram.values.optimal)-1))
        {
          #count=1
          polygon(c(nnm.breaks.histogram.values.optimal[count:(count+1)],rev(nnm.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1)[count])  
        }
        polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
        lines(nnm.probability.classification.optimal$prob_thres,nnm.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
        lines(nnm.probability.classification.optimal$prob_thres,nnm.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
        index_points_plot<-c(1,which(nnm.probability.classification.optimal$breaks_sel==TRUE),dim(nnm.probability.classification.optimal)[1])
        points(nnm.probability.classification.optimal[index_points_plot[1:floor(length(nnm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
        text(nnm.probability.classification.optimal[index_points_plot[1:floor(length(nnm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(nnm.probability.classification.optimal[index_points_plot[1:floor(length(nnm.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(nnm.breaks.histogram.values.optimal)/2)-1)))
        points(nnm.probability.classification.optimal[index_points_plot[ceiling(length(nnm.breaks.histogram.values.optimal)/2):length(nnm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
        text(nnm.probability.classification.optimal[index_points_plot[ceiling(length(nnm.breaks.histogram.values.optimal)/2):length(nnm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(nnm.probability.classification.optimal[index_points_plot[ceiling(length(nnm.breaks.histogram.values.optimal)/2):length(nnm.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(nnm.breaks.histogram.values.optimal)-ceiling(length(nnm.breaks.histogram.values.optimal)/2)),3))
      }
      
      pdf(file = "result_NNM_Histogram_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      hist(predict.result.nnm, breaks=nnm.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal NNM susceptibility", col=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1))
      dev.off()
      
      pdf(file = "result_NNM_ModelEvaluationPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("NNM OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
      mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
      mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
      for (count in 1:(length(nnm.breaks.histogram.values.optimal)-1))
      {
        #count=1
        polygon(c(nnm.breaks.histogram.values.optimal[count:(count+1)],rev(nnm.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1)[count])  
      }
      polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
      lines(nnm.probability.classification.optimal$prob_thres,nnm.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
      lines(nnm.probability.classification.optimal$prob_thres,nnm.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
      index_points_plot<-c(1,which(nnm.probability.classification.optimal$breaks_sel==TRUE),dim(nnm.probability.classification.optimal)[1])
      points(nnm.probability.classification.optimal[index_points_plot[1:floor(length(nnm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(nnm.probability.classification.optimal[index_points_plot[1:floor(length(nnm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(nnm.probability.classification.optimal[index_points_plot[1:floor(length(nnm.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(nnm.breaks.histogram.values.optimal)/2)-1)))
      points(nnm.probability.classification.optimal[index_points_plot[1+ceiling(length(nnm.breaks.histogram.values.optimal)/2):length(nnm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
      text(nnm.probability.classification.optimal[index_points_plot[1+ceiling(length(nnm.breaks.histogram.values.optimal)/2):length(nnm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(nnm.probability.classification.optimal[index_points_plot[1+ceiling(length(nnm.breaks.histogram.values.optimal)/2):length(nnm.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(nnm.breaks.histogram.values.optimal)-1-ceiling(length(nnm.breaks.histogram.values.optimal)/2)),3))
      points(nnm.probability.classification.optimal[nnm.optimal.index,c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(nnm.probability.classification.optimal[nnm.optimal.index,c("prob_thres","tpr")],labels=with(round(nnm.probability.classification.optimal[nnm.optimal.index,c("prob_thres","tpr")],3),paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(4))
      dev.off()
    }
  }
  
  ###########################################
  
  
  contingency.table.matrix.nnm<-matrix(nrow=dim(roc.plot.nnm.series$plot.data)[1],ncol=8)
  colnames(contingency.table.matrix.nnm)<-c("Threshold","TP","TN","FP","FN","TPR","FPR","COHEN_KAPPA")
  contingency.table.matrix.nnm[,1]<-roc.plot.nnm.series$plot.data[,1,1]
  contingency.table.matrix.nnm[,6]<-roc.plot.nnm.series$plot.data[,2,1]
  contingency.table.matrix.nnm[,7]<-roc.plot.nnm.series$plot.data[,3,1]
  values.observed<-training.table[,2]
  values.predicted<-predict.result.nnm
  for (count.threshold.series in 1:dim(roc.plot.nnm.series$plot.data)[1])
  {
    value.threshold<-contingency.table.matrix.nnm[count.threshold.series,1]
    values.probability.reclassified<-NULL
    values.probability.reclassified<-as.numeric(values.predicted>value.threshold) 
    #sum(values.probability.reclassified-round(values.predicted)) # Check sum: It has to be 0 if threshold is equal to 1
    series.pasted<-paste(values.observed,values.probability.reclassified,sep="")
    series.pasted<-gsub("00","1",series.pasted)
    series.pasted<-gsub("01","2",series.pasted)
    series.pasted<-gsub("10","3",series.pasted)
    series.pasted<-gsub("11","4",series.pasted)
    series.pasted<-as.numeric(series.pasted)
    TP<-as.numeric(sum(series.pasted>=4)) # True Positive
    FN<-as.numeric(sum(series.pasted>=3 & series.pasted<4)) # False Negative
    FP<-as.numeric(sum(series.pasted>=2 & series.pasted<3)) # False Positive
    TN<-as.numeric(sum(series.pasted>=1 & series.pasted<2)) # True Negative              
    #TPR<-TP/(TP+FN) # Hit Rate or True Positive Rate or Sensitivity - Assigned before the for cicle using rocplot data
    #FPR<-FP/(FP+TN) # False Alarm Rate or False Positive Rate or 1-Specificity
    # Cohen's Kappa = (agreement-chance)/(1-chance)  where agreement=(TP+TN)/(TP+TN+FP+FN) and chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    agreement=(TP+TN)/(TP+TN+FP+FN)
    chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    cohen.kappa.value<-(agreement-chance)/(1-chance)
    #Other
    #library(vcd)
    #cohen.kappa.value<-Kappa(cross.classification.table)
    contingency.table.matrix.nnm[count.threshold.series,2]<-TP
    contingency.table.matrix.nnm[count.threshold.series,3]<-TN
    contingency.table.matrix.nnm[count.threshold.series,4]<-FP
    contingency.table.matrix.nnm[count.threshold.series,5]<-FN
    contingency.table.matrix.nnm[count.threshold.series,8]<-cohen.kappa.value
  }
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    plot(roc.plot.nnm.series$plot.data[,1,1],roc.plot.nnm.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="NNM MODEL EVALUATION PLOT")
    lines(roc.plot.nnm.series$plot.data[,1,1],1-roc.plot.nnm.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
    lines(roc.plot.nnm.series$plot.data[,1,1], contingency.table.matrix.nnm[,8],col="blue",lty=1,lwd=1)
    mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
    mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
    mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  }
  
  pdf(file = "result_NNM_ModelEvaluationPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(roc.plot.nnm.series$plot.data[,1,1],roc.plot.nnm.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="NNM MODEL EVALUATION PLOT")
  lines(roc.plot.nnm.series$plot.data[,1,1],1-roc.plot.nnm.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.nnm.series$plot.data[,1,1], contingency.table.matrix.nnm[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  dev.off()
  
  ## VALIDATION OF NNM MODEL (Matching NNM posterior probability results and validation grouping variable)
  
  # Predict NNM
  # Result Predicted
  
  predict.result.nnm.validation<-predict(result.nnm,newdata=validation.explanatory.variables)
  #str(predict.result.nnm.validation)
  #summary(predict.result.nnm.validation)
  cross.classification.validation.nnm<-table(validation.grouping.variable,round(predict.result.nnm.validation),dnn=c("Observed","Predicted"))
  rownames(cross.classification.validation.nnm)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.validation.nnm)<-list("No Landslide","Landslide") # Predicted
  #str(cross.classification.validation.nnm)
  #cross.classification.validation.nnm<-table(validation.grouping.variable,round(predict.result.nnm),dnn=c("Observed","Predicted"))
  
  
  #Elaboration of Coefficient of association for contingency table
  #load package (vcd)
  library(vcd)
  
  #help(package=vcd)
  contingency.table.validation.nnm<-table2d_summary(cross.classification.validation.nnm)
  test.table.validation.nnm<-assocstats(cross.classification.validation.nnm)
  
  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.validation.nnm/sum(cross.classification.validation.nnm)*100,2), std="margin", main="VALIDATION NNM MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.nnm, std="margin", main="VALIDATION NNM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #load package (verification)
  library(verification)
  
  # 2nd method using verify function
  verification.validation.nnm<-verify(validation.table[,2],predict.result.nnm.validation, frcst.type="prob", obs.type="binary")
  #summary(verification.validation.lrm)
  
  # showing confidence intervals.  MAY BE SLOW
  area.under.roc.curve.validation.nnm<-roc.area(validation.table[,2],predict.result.nnm.validation)
  
  if (cross.classification.validation.nnm[1,2]==0 | cross.classification.validation.nnm[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.validation.nnm, main = "ROC PLOT: VALIDATION NNM MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[4] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.validation.nnm$A,2),";  Sample size = ",area.under.roc.curve.validation.nnm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[4], sep=""), side=3, col="red", cex=0.8)
  }
  
  # EXPORT OF PLOT FOR VALIDATION OF NNM MODEL
  
  pdf(file = "result_NNM_FourfoldPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.validation.nnm/sum(cross.classification.validation.nnm)*100,2), std="margin", main="VALIDATION NNM MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  #fourfold(cross.classification.validation.nnm, std="margin", main="VALIDATION NNM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  dev.off()
  
  #pdf(file = "result_NNM_ROCPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.validation.nnm, main = "ROC PLOT: VALIDATION NNM MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.validation.nnm<-roc.area(verification.table[,2],predict.result.nnm.validation)
  #dev.off()
  
  pdf(file = "result_NNM_ROCPlot_bootstrap_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.validation.nnm, main = "ROC PLOT: VALIDATION NNM MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[4] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.validation.nnm$A,2),";  Sample size = ",area.under.roc.curve.validation.nnm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[4], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  # Assignation of a matching code between observed and predicted values calculated using the validation dataset
  validation.nnm.matching.code<-paste(validation.grouping.variable,round(predict.result.nnm.validation),sep="")
  validation.nnm.matching.code<-gsub("00","1",validation.nnm.matching.code)
  validation.nnm.matching.code<-gsub("01","2",validation.nnm.matching.code)
  validation.nnm.matching.code<-gsub("10","3",validation.nnm.matching.code)
  validation.nnm.matching.code<-gsub("11","4",validation.nnm.matching.code)
  validation.nnm.matching.code<-as.numeric(validation.nnm.matching.code)
  
  ##########################################
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    cross.classification.validation.nnm.optimal<-table(validation.grouping.variable,as.numeric(predict.result.nnm.validation>nnm.probability.optimal.binary.threshold),dnn=c("Observed","Predicted"))
    rownames(cross.classification.validation.nnm.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.validation.nnm.optimal)<-list("No Landslide","Landslide") # Predicted
    
    #Different plots for contingency table
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.validation.nnm.optimal/sum(cross.classification.validation.nnm.optimal)*100,2), std="margin", main="VALIDATION NNM MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
      #fourfold(cross.classification.validation.nnm, std="margin", main="VALIDATION NNM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    }
    if (cross.classification.validation.nnm.optimal[1,2]==0 | cross.classification.validation.nnm.optimal[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
    pdf(file = "result_NNM_FourfoldPlot_Validation_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.validation.nnm.optimal/sum(cross.classification.validation.nnm.optimal)*100,2), std="margin", main="VALIDATION NNM MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.nnm, std="margin", main="VALIDATION NNM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    dev.off()
    
    
    # Assignation of a optimal matching code between observed and predicted values calculated using the validation dataset
    validation.nnm.matching.code.optimal<-paste(validation.grouping.variable,as.numeric(predict.result.nnm.validation>nnm.probability.optimal.binary.threshold),sep="")
    validation.nnm.matching.code.optimal<-gsub("00","1",validation.nnm.matching.code.optimal)
    validation.nnm.matching.code.optimal<-gsub("01","2",validation.nnm.matching.code.optimal)
    validation.nnm.matching.code.optimal<-gsub("10","3",validation.nnm.matching.code.optimal)
    validation.nnm.matching.code.optimal<-gsub("11","4",validation.nnm.matching.code.optimal)
    validation.nnm.matching.code.optimal<-as.numeric(validation.nnm.matching.code.optimal)
  }
  #########################################
  
  
  # EXPORT OF NNM MODEL RESULTS
  write.table("RESULTS OF NEURAL NETWORK MODEL",file="result_NNM.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("NNM MODEL OUTPUTS",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  if (analysis.parameter.matrix[4] == "OPT")
  {
    write.table("Selection of Neural Network Structure",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("WEIGHT DECAY","HIDDEN LAYER NODES","SUM SQUARED ERROR","ROC AREA"),nnm.table),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  }
  write.table("Neural Network Structure Selected",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("INPUT NODES","HIDDEN LAYER NODES","OUTPUT NODES"),cbind(result.nnm$n[1],result.nnm$n[2],result.nnm$n[3])),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Value of Weight Decay Term Selected",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(weight.decay.selected,file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Value of Fitting Criterion Plus Weight Decay Term",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(result.nnm$value,file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Best Set of Weights Found",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(result.nnm$wts,file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE MODEL RESULT",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.nnm$table[,1,],contingency.table.nnm$table[,2,],contingency.table.nnm$table[,3,])),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE VALIDATION",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.validation.nnm$table[,1,],contingency.table.validation.nnm$table[,2,],contingency.table.validation.nnm$table[,3,])),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("MATCHING CODE DEFINITION",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(c("","OBSERVED NO LANDSLIDES: 0","OBSERVED LANDSLIDES: 1"), c("PREDICTED NO LANDSLIDES: 0","00 -> Code 1","10 -> Code 3"), c("PREDICTED LANDSLIDES: 1","01 -> Code 2","11 -> Code 4")),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  
  ########
  if(enable_probability_optimal_binary_classification==FALSE) 
  {
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","NNM MODEL POSTERIOR PROBABILITY","NNM MODEL CLASSIFICATION","NNM MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict.result.nnm,round(predict.result.nnm),result.nnm.matching.code)),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","NNM VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.nnm.validation,round(predict.result.nnm.validation),validation.nnm.matching.code)),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } else
  {
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    if(enable_probability_optimal_classification==TRUE) 
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY PARTITION -> Method: ",type_probability_optimal_classification,sep=""),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(nnm.probability.classification.optimal[index_points_plot,c("prob_thres","tnr","tpr")]),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    } else
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY BINARY PARTITION",sep=""),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(nnm.probability.classification.optimal[nnm.optimal.index,c("prob_thres","tnr","tpr")]),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    }
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","NNM MODEL POSTERIOR PROBABILITY","NNM MODEL CLASSIFICATION","NNM MODEL RESULT MATCHING CODE","NNM OPTIMAL MODEL CLASSIFICATION","NNM OPTIMAL MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict.result.nnm,round(predict.result.nnm),result.nnm.matching.code,as.numeric(predict.result.nnm>nnm.probability.optimal.binary.threshold),result.nnm.matching.code.optimal)),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","NNM VALIDATION MATCHING CODE","OPTIMAL VALIDATION CLASSIFICATION","OPTIMAL NNM VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.nnm.validation,round(predict.result.nnm.validation),validation.nnm.matching.code,as.numeric(predict.result.nnm.validation>nnm.probability.optimal.binary.threshold),validation.nnm.matching.code.optimal)),file="result_NNM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  }
  ########
  
  
  
  # PLOT AND EXPORT OF NNM MODEL MAPS
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POLYGONS")
  {
    ################################################
    ### Prima di cancellare testare con polygoni
    ##############################################
    #shape_training_nnm<-shape_training
    #result_training_nnm_shape<-cbind(identification.value,as.numeric(grouping.variable)-1,predict.result.nnm,round(predict.result.nnm),result.nnm.matching.code)
    #colnames(result_training_nnm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH")
    #shape_training_nnm@data <- merge(x=shape_training_nnm@data,y=result_training_nnm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)  
    ##writeOGR(shape_training_nnm,dsn="result_NNM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_training_nnm,dsn="result_NNM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    #shape_validation_nnm<-shape_validation
    #result_validation_nnm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.nnm.validation,round(predict.result.nnm.validation),validation.nnm.matching.code)
    #colnames(result_validation_nnm_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    #shape_validation_nnm@data <- merge(x=shape_validation_nnm@data,y=result_validation_nnm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    ##writeOGR(shape_validation_nnm,dsn="result_NNM_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_validation_nnm,dsn="result_NNM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ### DA qui Nuovo
    ##################################################################   
    #### aggiungere a poligoni colonna incertezza ed esportazione pdf incertezza
    shape_training_nnm<-shape_training
    result_training_nnm_shape<-cbind(identification.value,as.numeric(grouping.variable)-1,predict.result.nnm,round(predict.result.nnm),result.nnm.matching.code,ID.bootstrap.model.nnm.count,bootstrap.model.nnm.probability.mean,bootstrap.model.nnm.probability.sd,bootstrap.model.nnm.probability.min,bootstrap.model.nnm.probability.max,bootstrap.model.nnm.probability.sderror,t(bootstrap.model.nnm.probability.quantiles),bootstrap.model.nnm.prediction.mean,bootstrap.model.nnm.prediction.sd,bootstrap.model.nnm.prediction.min,bootstrap.model.nnm.prediction.max,bootstrap.model.nnm.prediction.sderror,t(bootstrap.model.nnm.prediction.quantiles))
    colnames(result_training_nnm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","NNM_SAMP","NNM_PMean","NNM_PSd","NNM_PMin","NNM_PMax","NNM_PSder","NNM_PQ0","NNM_PQ_005","NNM_PQ_025","NNM_PQ05","NNM_PQ_075","NNM_PQ095","NNM_PQ1","NNM_PrMean","NNM_PrSd","NNM_PrMin","NNM_PrMax","NNM_PrSder","NNM_PrQ0","NNM_PrQ005","NNM_PrQ025","NNM_PrQ05","NNM_PrQ075","NNM_PrQ095","NNM_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_nnm_shape<-cbind(result_training_nnm_shape,as.numeric(predict.result.nnm>nnm.probability.optimal.binary.threshold),result.nnm.matching.code.optimal)
      colnames(result_training_nnm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","NNM_SAMP","NNM_PMean","NNM_PSd","NNM_PMin","NNM_PMax","NNM_PSder","NNM_PQ0","NNM_PQ_005","NNM_PQ_025","NNM_PQ05","NNM_PQ_075","NNM_PQ095","NNM_PQ1","NNM_PrMean","NNM_PrSd","NNM_PrMin","NNM_PrMax","NNM_PrSder","NNM_PrQ0","NNM_PrQ005","NNM_PrQ025","NNM_PrQ05","NNM_PrQ075","NNM_PrQ095","NNM_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    
    shape_training_nnm@data <- merge(x=shape_training_nnm@data,y=result_training_nnm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)  
    #writeOGR(shape_training_nnm,dsn="result_NNM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_nnm,dsn="result_NNM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    # WARNING: The validation does't have the unvertainty estimation: probably this can be daone using the parabolic error function 
    shape_validation_nnm<-shape_validation
    result_validation_nnm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.nnm.validation,round(predict.result.nnm.validation),validation.nnm.matching.code)
    colnames(result_validation_nnm_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_nnm_shape<-cbind(result_validation_nnm_shape,as.numeric(predict.result.nnm.validation>nnm.probability.optimal.binary.threshold),validation.nnm.matching.code.optimal)
      colnames(result_validation_nnm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    
    shape_validation_nnm@data <- merge(x=shape_validation_nnm@data,y=result_validation_nnm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_nnm@data <- cbind(shape_validation_nnm@data,PROB_SDMOD=(coefficients(fit.parabola.probability.nnm)*(shape_validation_nnm@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.nnm)*shape_validation_nnm@data$VAL_PROB))
    #writeOGR(shape_validation_nnm,dsn="result_NNM_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_nnm,dsn="result_NNM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ##################################################################
    
    
    # Plot and export of maps
    # NNM Susceptibility
    #dev.new()
    pdf(file = "result_NNM_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_nnm, zcol=c("MOD_PROB"), names.attr=c("NNM MODEL PROBABILITY"), main="NNM MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # NNM Model Matching Code
    #dev.new()
    pdf(file = "result_NNM_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_nnm, zcol=c("MOD_MATCH"), names.attr=c("NNM MODEL MATCHING CODE"), main="NNM MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
    
    # NNM Model Uncertainty
    #dev.new()
    pdf(file = "result_NNM_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_nnm, zcol=c("NNM_PrSd"), names.attr=c("NNM MODEL UNCERTAINTY"), main="NNM MODEL UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()
    
    # NNM Validation Susceptibility
    #dev.new()
    pdf(file = "result_NNM_Validation_Susceptibility_Map.pdf")
    print(spplot(obj=shape_validation_nnm, zcol=c("VAL_PROB"), names.attr=c("NNM VALIDATION PROBABILITY"), main="NNM VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # NNM Validation Matching Code
    #dev.new()
    pdf(file = "result_NNM_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_nnm, zcol=c("VAL_MATCH"), names.attr=c("NNM VALIDATION MATCHING CODE"), main="NNM VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
        
    # NNM Validation Uncertainty
    #dev.new()
    pdf(file = "result_NNM_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_nnm, zcol=c("PROB_SDMOD"), names.attr=c("NNM VALIDATION UNCERTAINTY"), main="NNM VALIDATION UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()

    ########### TBT
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      # NNM Model Matching Code
      #dev.new()
      pdf(file = "result_NNM_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_training_nnm, zcol=c("OPT_MATCH"), names.attr=c("NNM MODEL MATCHING CODE"), main="NNM MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      # NNM Model Matching Code
      #dev.new()
      pdf(file = "result_NNM_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_validation_nnm, zcol=c("OPT_MATCH"), names.attr=c("NNM VALIDATION MATCHING CODE"), main="NNM VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      if(enable_probability_optimal_classification==TRUE)
      {
        pdf(file = "result_NNM_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
        print(spplot(obj=shape_training_nnm, zcol=c("MOD_PROB"), names.attr=c("NNM MODEL PROBABILITY"), main="NNM MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=(nnm.breaks.histogram.values.optimal)+c(rep(0,(length(nnm.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1)))
        dev.off()
        
        pdf(file = "result_NNM_Validation_Susceptibility_Map_Optimal.pdf")
        print(spplot(obj=shape_validation_nnm, zcol=c("VAL_PROB"), names.attr=c("NNM VALIDATION PROBABILITY"), main="NNM VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=(nnm.breaks.histogram.values.optimal)+c(rep(0,(length(nnm.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1)))
        dev.off()
      }
    }
    ###########
    
    
    ### Success & prediction rate curve
    susc_values_nnm<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_nnm@data) %in% c(configuration.spatial.data.table[c(5,6)],"MOD_PROB"))
    shape_training_nnm@data<-shape_training_nnm@data[order(shape_training_nnm@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_training_nnm@data[,1]<-cumsum(shape_training_nnm@data[,1])/sum(shape_training_nnm@data[,1])*100
    shape_training_nnm@data[,2]<-cumsum(shape_training_nnm@data[,2])/sum(shape_training_nnm@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_nnm@data[,3],shape_training_nnm@data[,1])
    area_susc_values_nnm_mod<-c(0,ind_pro_fun_mod(susc_values_nnm[-c(1,length(susc_values_nnm))]),100)
    
    #dev.new()
    pdf(file = "result_NNM_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="NNM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_nnm)-1))
    {
      #count=1
      polygon(c(area_susc_values_nnm_mod[count:(count+1)],rev(area_susc_values_nnm_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_nnm@data[,1],shape_training_nnm@data[,2],col="black")
    dev.off()
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_nnm_mod_optimal<-c(0,ind_pro_fun_mod(rev(nnm.breaks.histogram.values.optimal)[-c(1,length(nnm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_NNM_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="NNM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(nnm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_nnm_mod_optimal[count:(count+1)],rev(area_susc_values_nnm_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_nnm@data[,1],shape_training_nnm@data[,2],col="black")
      dev.off()
    }
    ############
    
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_nnm@data) %in% c(configuration.spatial.data.table[c(5,6)],"VAL_PROB"))
    shape_validation_nnm@data<-shape_validation_nnm@data[order(shape_validation_nnm@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_validation_nnm@data[,1]<-cumsum(shape_validation_nnm@data[,1])/sum(shape_validation_nnm@data[,1])*100
    shape_validation_nnm@data[,2]<-cumsum(shape_validation_nnm@data[,2])/sum(shape_validation_nnm@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_nnm@data[,3],shape_validation_nnm@data[,1])
    area_susc_values_nnm_val<-c(0,ind_pro_fun_val(susc_values_nnm[-c(1,length(susc_values_nnm))]),100)
    
    #dev.new()
    pdf(file = "result_NNM_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="NNM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_nnm)-1))
    {
      #count=1
      polygon(c(area_susc_values_nnm_val[count:(count+1)],rev(area_susc_values_nnm_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_nnm@data[,1],shape_validation_nnm@data[,2],col="black")
    dev.off()
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_nnm_val_optimal<-c(0,ind_pro_fun_val(rev(nnm.breaks.histogram.values.optimal)[-c(1,length(nnm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_NNM_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="NNM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(nnm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_nnm_val_optimal[count:(count+1)],rev(area_susc_values_nnm_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_nnm@data[,1],shape_validation_nnm@data[,2],col="black")
      dev.off()
    }
    ############
    
  }
  
  
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POINTS")
  {
    shape_training_nnm<-shape_training
    result_training_nnm_shape<-cbind(identification.value,training.table[,2],predict.result.nnm$posterior[,2],as.numeric(levels(predict.result.nnm$class))[predict.result.nnm$class],result.nnm.matching.code,ID.bootstrap.model.nnm.count,bootstrap.model.nnm.probability.mean,bootstrap.model.nnm.probability.sd,bootstrap.model.nnm.probability.min,bootstrap.model.nnm.probability.max,bootstrap.model.nnm.probability.sderror,t(bootstrap.model.nnm.probability.quantiles),bootstrap.model.nnm.prediction.mean,bootstrap.model.nnm.prediction.sd,bootstrap.model.nnm.prediction.min,bootstrap.model.nnm.prediction.max,bootstrap.model.nnm.prediction.sderror,t(bootstrap.model.nnm.prediction.quantiles))
    colnames(result_training_nnm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","NNM_SAMP","NNM_PMean","NNM_PSd","NNM_PMin","NNM_PMax","NNM_PSder","NNM_PQ0","NNM_PQ_005","NNM_PQ_025","NNM_PQ05","NNM_PQ_075","NNM_PQ095","NNM_PQ1","NNM_PrMean","NNM_PrSd","NNM_PrMin","NNM_PrMax","NNM_PrSder","NNM_PrQ0","NNM_PrQ005","NNM_PrQ025","NNM_PrQ05","NNM_PrQ075","NNM_PrQ095","NNM_PrQ1")
    colnames(result_training_nnm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","NNM_SAMP","NNM_PMean","NNM_PSd","NNM_PMin","NNM_PMax","NNM_PSder","NNM_PQ0","NNM_PQ_005","NNM_PQ_025","NNM_PQ05","NNM_PQ_075","NNM_PQ095","NNM_PQ1","NNM_PrMean","NNM_PrSd","NNM_PrMin","NNM_PrMax","NNM_PrSder","NNM_PrQ0","NNM_PrQ005","NNM_PrQ025","NNM_PrQ05","NNM_PrQ075","NNM_PrQ095","NNM_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_nnm_shape<-cbind(result_training_nnm_shape,as.numeric(predict.result.nnm>nnm.probability.optimal.binary.threshold),result.nnm.matching.code.optimal)
      colnames(result_training_nnm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","NNM_SAMP","NNM_PMean","NNM_PSd","NNM_PMin","NNM_PMax","NNM_PSder","NNM_PQ0","NNM_PQ_005","NNM_PQ_025","NNM_PQ05","NNM_PQ_075","NNM_PQ095","NNM_PQ1","NNM_PrMean","NNM_PrSd","NNM_PrMin","NNM_PrMax","NNM_PrSder","NNM_PrQ0","NNM_PrQ005","NNM_PrQ025","NNM_PrQ05","NNM_PrQ075","NNM_PrQ095","NNM_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    shape_training_nnm@data <- merge(x=shape_training_nnm@data,y=result_training_nnm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_nnm,dsn="result_NNM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_nnm,dsn="result_NNM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    
    # WARNING: The validation does't have the unvertainty estimation: probabilty this can be daone using the parabolic error function 
    shape_validation_nnm<-shape_validation
    result_validation_nnm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.nnm.validation$posterior[,2],as.numeric(levels(predict.result.nnm.validation$class))[predict.result.nnm.validation$class],validation.nnm.matching.code)
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_nnm_shape<-cbind(result_validation_nnm_shape,as.numeric(predict.result.nnm.validation>nnm.probability.optimal.binary.threshold),validation.nnm.matching.code.optimal)
      colnames(result_validation_nnm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    colnames(result_validation_nnm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
    shape_validation_nnm@data <- merge(x=shape_validation_nnm@data,y=result_validation_nnm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_nnm@data <- cbind(shape_validation_nnm@data,PROB_SDMOD=(coefficients(fit.parabola.probability.nnm)*(shape_validation_nnm@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.nnm)*shape_validation_nnm@data$VAL_PROB))
    #writeOGR(shape_validation_nnm,dsn="result_NNM_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_nnm,dsn="result_NNM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    require(raster)
    
    # Plot and export of maps
    # NNM Susceptibility
    #dev.new()
    pdf(file = "result_NNM_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_nnm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_NNM_Model_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_NNM_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1),breaks=round(nnm.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    
    # NNM Model Matching Code
    #dev.new()
    pdf(file = "result_NNM_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_nnm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
    print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    
    writeRaster(layer_gridded_raster, filename="result_NNM_Model_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_NNM_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_nnm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
      print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_NNM_Model_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_NNM_Model_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_nnm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_NNM_Model_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
    }
    ############
    
    # NNM Model uncertainity
    #dev.new()
    pdf(file = "result_NNM_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_nnm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"NNM_PrSd"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_NNM_Model_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    
    
    # NNM Validation Susceptibility
    #dev.new()
    pdf(file = "result_NNM_Validation_Susceptibility_Map.pdf")
    layer_gridded<-shape_validation_nnm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_NNM_Validation_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_NNM_Validation_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1),breaks=round(nnm.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    # NNM Validation Matching Code
    #dev.new()
    pdf(file = "result_NNM_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_nnm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    #zoom(layer_gridded_raster)  	
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_NNM_Validation_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    
    ########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_NNM_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_nnm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      #zoom(layer_gridded_raster)  	
      dev.off()
      writeRaster(layer_gridded_raster, filename="result_NNM_Validation_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_NNM_Validation_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_nnm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_NNM_Validation_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
      
    }
    #######
    
    # NNM Model uncertainity
    #dev.new()
    pdf(file = "result_NNM_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_nnm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"PROB_SDMOD"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_NNM_Validation_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    
    ### Success & prediction rate curve
    susc_values_nnm<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_nnm@data) %in% c("MOD_GROUP","MOD_PROB"))
    shape_training_nnm@data<-cbind(PIXEL_AREA=rep(configuration.spatial.data.table[c(8)]^2,dim(shape_training_nnm@data)[1]),shape_training_nnm@data[order(shape_training_nnm@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_training_nnm@data$MOD_GROUP<-shape_training_nnm@data$MOD_GROUP*configuration.spatial.data.table[c(8)]^2
    shape_training_nnm@data[,1]<-cumsum(shape_training_nnm@data[,1])/sum(shape_training_nnm@data[,1])*100
    shape_training_nnm@data[,2]<-cumsum(shape_training_nnm@data[,2])/sum(shape_training_nnm@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_nnm@data[,3],shape_training_nnm@data[,1])
    area_susc_values_nnm_mod<-c(0,ind_pro_fun_mod(susc_values_nnm[-c(1,length(susc_values_nnm))]),100)
    
    #dev.new()
    pdf(file = "result_NNM_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="NNM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_nnm)-1))
    {
      #count=1
      polygon(c(area_susc_values_nnm_mod[count:(count+1)],rev(area_susc_values_nnm_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_nnm@data[,1],shape_training_nnm@data[,2],col="black")
    dev.off()
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_nnm_mod_optimal<-c(0,ind_pro_fun_mod(rev(nnm.breaks.histogram.values.optimal)[-c(1,length(nnm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_NNM_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="NNM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(nnm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_nnm_mod_optimal[count:(count+1)],rev(area_susc_values_nnm_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_nnm@data[,1],shape_training_nnm@data[,2],col="black")
      dev.off()
    }
    ############
    
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_nnm@data) %in% c("VAL_GROUP","VAL_PROB"))
    shape_validation_nnm@data<-cbind(PIXEL_AREA=rep(as.numeric(configuration.spatial.data.table[c(8)])^2,dim(shape_validation_nnm@data)[1]),shape_validation_nnm@data[order(shape_validation_nnm@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_validation_nnm@data$VAL_GROUP<-shape_validation_nnm@data$VAL_GROUP*as.numeric(configuration.spatial.data.table[c(8)])^2
    shape_validation_nnm@data[,1]<-cumsum(shape_validation_nnm@data[,1])/sum(shape_validation_nnm@data[,1])*100
    shape_validation_nnm@data[,2]<-cumsum(shape_validation_nnm@data[,2])/sum(shape_validation_nnm@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_nnm@data[,3],shape_validation_nnm@data[,1])
    area_susc_values_nnm_val<-c(0,ind_pro_fun_val(susc_values_nnm[-c(1,length(susc_values_nnm))]),100)
    
    #dev.new()
    pdf(file = "result_NNM_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="NNM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_nnm)-1))
    {
      #count=1
      polygon(c(area_susc_values_nnm_val[count:(count+1)],rev(area_susc_values_nnm_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_nnm@data[,1],shape_validation_nnm@data[,2],col="black")
    dev.off()
    
    ###########
    
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_nnm_val_optimal<-c(0,ind_pro_fun_val(rev(nnm.breaks.histogram.values.optimal)[-c(1,length(nnm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_NNM_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="NNM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(nnm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_nnm_val_optimal[count:(count+1)],rev(area_susc_values_nnm_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(nnm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_nnm@data[,1],shape_validation_nnm@data[,2],col="black")
      dev.off()
    }
    ############
    
  }
  
  
}

  

#-------------------- FORECAST COMBINATION MODEL --------------------#
#####  FORECAST COMBINATION USING A LOGISTIC REGRESSION MODEL
#####(A constrained Ordinary Least Squared estimation can also be adopted)

if(model.run.matrix[5] == "YES")
  {
  
  #library(Zelig) 
  #forecasting.combined.variables<-as.data.frame(cbind(data.variables[,1],predict.result.lda$posterior[,2],predict.result.qda$posterior[,2],result.cfm$result$fitted.values,predict.result.nnm))
  #colnames(forecasting.combined.variables)<-c("FRAX","resultlda","resultqda","resultcfm","resultnnm")   # Names of column mustn't have points
  
  forecasting.combined.variables<-as.data.frame(cbind(FRAX=data.variables[,1]))
  
  if(model.run.matrix[1]=="YES") { forecasting.combined.variables<-cbind(forecasting.combined.variables,resultlda=predict.result.lda$posterior[,2])}
  if(model.run.matrix[2]=="YES") { forecasting.combined.variables<-cbind(forecasting.combined.variables,resultqda=predict.result.qda$posterior[,2])}
  if(model.run.matrix[3]=="YES") { forecasting.combined.variables<-cbind(forecasting.combined.variables,resultcfm=predict(result.lrm,type="response"))}
  if(model.run.matrix[4]=="YES") { forecasting.combined.variables<-cbind(forecasting.combined.variables,resultnnm=predict.result.nnm)}
  
  forecasting.combined.variables.validation<-as.data.frame(cbind(FRAX=validation.table[,2]))
  
  if(model.run.matrix[1]=="YES") { forecasting.combined.variables.validation<-cbind(forecasting.combined.variables.validation,resultlda=predict.result.lda.validation$posterior[,2])}
  if(model.run.matrix[2]=="YES") { forecasting.combined.variables.validation<-cbind(forecasting.combined.variables.validation,resultqda=predict.result.qda.validation$posterior[,2])}
  if(model.run.matrix[3]=="YES") { forecasting.combined.variables.validation<-cbind(forecasting.combined.variables.validation,resultcfm=predict.result.lrm.validation.posterior)}
  if(model.run.matrix[4]=="YES") { forecasting.combined.variables.validation<-cbind(forecasting.combined.variables.validation,resultnnm=predict.result.nnm.validation)}
  
  #if (class(try(zelig(as.formula(paste(names(forecasting.combined.variables)[1],"~",paste(names(forecasting.combined.variables)[2:dim(forecasting.combined.variables)[2]],collapse= "+"))), data=forecasting.combined.variables, model="logit",cite=FALSE)))=="try-error")  
  if (class(try(glm(as.formula(paste(names(forecasting.combined.variables)[1],"~",paste(names(forecasting.combined.variables)[2:dim(forecasting.combined.variables)[2]],collapse= "+"))), data=forecasting.combined.variables,family=binomial())))=="try-error")
    
  { 
    #zelig(as.formula(paste(names(forecasting.combined.variables)[1],"~",paste(names(forecasting.combined.variables)[2:dim(forecasting.combined.variables)[2]],collapse= "+"))), data=forecasting.combined.variables, model="logit")
    write.table("The combination of forecast using Logistic Regression Model was not completed",file="Error_CFM_Analysis.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="Error_CFM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("Error LOG",file="Error_CFM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind("Message",rev(1:length(as.vector(.Traceback)))," ->",as.vector(.Traceback)),file="Error_CFM_Analysis.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } 
  
  result.cfm<-NULL
  #result.cfm<-zelig(as.formula(paste(names(forecasting.combined.variables)[1],"~",paste(names(forecasting.combined.variables)[2:dim(forecasting.combined.variables)[2]],collapse= "+"))), data=forecasting.combined.variables, model="logit",cite=FALSE)
  result.cfm<-glm(as.formula(paste(names(forecasting.combined.variables)[1],"~",paste(names(forecasting.combined.variables)[2:dim(forecasting.combined.variables)[2]],collapse= "+"))), data=forecasting.combined.variables,family=binomial())
  #str(result.cfm)    
  #names(result.cfm)
  
  #for predicted value (posterior probablity calculated with model) result.cfm$result$fitted.values was considered
  
  #cross.classification.cfm<-table(as.numeric(result.cfm$result$y),round(result.cfm$result$fitted.values),dnn=c("Observed","Predicted"))
  cross.classification.cfm<-table(grouping.variable,round(predict(result.cfm, type="response")),dnn=c("Observed","Predicted"))
  rownames(cross.classification.cfm)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.cfm)<-list("No Landslide","Landslide") # Predicted    
  #str(cross.classification.cfm)
  
  # Assignation of a matching code between observed and predicted values
  #result.cfm.matching.code<-paste(grouping.variable,round(result.cfm$result$fitted.values),sep="")
  result.cfm.matching.code<-paste(grouping.variable,round(predict(result.cfm, type="response")),sep="")
  result.cfm.matching.code<-gsub("00","1",result.cfm.matching.code)
  result.cfm.matching.code<-gsub("01","2",result.cfm.matching.code)
  result.cfm.matching.code<-gsub("10","3",result.cfm.matching.code)
  result.cfm.matching.code<-gsub("11","4",result.cfm.matching.code)
  result.cfm.matching.code<-as.numeric(result.cfm.matching.code)
  
  #Elaboration of Coefficient of association for contingency table
  #load package (vcd)  
  library(vcd)
  
  #help(package=vcd)         
  contingency.table.cfm<-table2d_summary(cross.classification.cfm)
  test.table.cfm<-assocstats(cross.classification.cfm)
  
  
  #Different plots for contingency table 
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.cfm/sum(cross.classification.cfm)*100,2),std="margin", main="COMBINATION LOGISTIC REGRESSION MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.cfm,std="margin", main="COMBINATION LOGISTIC REGRESSION MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #A ROC curve plots the false alarm rate against the hit rate
  #for a probablistic forecast for a range of thresholds. 
  
  #load package (verification)  
  library(verification)
  
  #verify function
  #Based on the type of inputs, this function calculates a range of verification statistics and skill scores.
  #Additionally, it creates a verify class object that can be further analyzed.
  
  
  
  ##### ROC PLOT OBS - POSTERIOR PROBABILITY ASSOCIATED TO 1                                                                                 
  
  # Method using verify function
  #verification.results.cfm<-verify(result.cfm$result$y,result.cfm$result$fitted.values, frcst.type="prob", obs.type="binary")
  verification.results.cfm<-verify(training.table[,2],predict(result.cfm, type="response"), frcst.type="prob", obs.type="binary")
  
  #str(verification.results.cfm)
  #if (enable_screen_plotting==TRUE)
  #{
  #dev.new()
  #roc.plot(verification.results.cfm, main = "ROC PLOT: COMBINATION LOGISTIC REGRESSION", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.cfm<-roc.area(result.cfm$result$y,result.cfm$result$fitted.values)
  #}
  area.under.roc.curve.cfm<-roc.area(training.table[,2],predict(result.cfm, type="response"))
  
  
  ## showing confidence intervals.  MAY BE SLOW
  
  if (cross.classification.cfm[1,2]==0 | cross.classification.cfm[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()                                                              
    roc.plot(verification.results.cfm, main = "ROC PLOT: COMBINATION LOGISTIC REGRESSION", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[5] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.cfm$A,2),";  Sample size = ",area.under.roc.curve.cfm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[5], sep=""), side=3, col="red", cex=0.8)
    ## Histogram of posterior probability
    dev.new()                            
    #hist(result.cfm$result$fitted.values, breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Combination Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
    hist(predict(result.cfm, type="response"), breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Combination Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  }
  pdf(file = "result_CFM_Histogram.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #hist(result.cfm$result$fitted.values, breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Combination Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  hist(predict(result.cfm, type="response"), breaks=breaks.histogram.values, freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of Combination Logistic Regression Model susceptibility", col=color_ramp_palette_fun(length(breaks.histogram.values)-1))
  dev.off() 
  
  
  # EXPORT OF PLOT FOR CFM MODEL
  
  pdf(file = "result_CFM_FourfoldPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.cfm/sum(cross.classification.cfm)*100,2),std="margin", main="COMBINATION LOGISTIC REGRESSION MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  #fourfold(cross.classification.cfm,std="margin", main="COMBINATION LOGISTIC REGRESSION MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
  dev.off()
  
  #pdf(file = "result_CFM_ROCPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.results.cfm, main = "ROC PLOT: COMBINATION LOGISTIC REGRESSION", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #dev.off()
  
  pdf(file = "result_CFM_ROCPlot_bootstrap.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.results.cfm, main = "ROC PLOT: COMBINATION LOGISTIC REGRESSION", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[5] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.cfm$A,2),";  Sample size = ",area.under.roc.curve.cfm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[5], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  ## BOOTSTRAP PROCEDURE FOR THE ESTIMATION OF MODEL PREDICTION VARIABILITY
  if(bootstrap.model.variability[5] == "YES")
  {
    bootstrap.sample.model.cfm<-bootstrap.sample.model[5]
    
    matrix.bootstrap.model.cfm<-matrix(data=NA, nrow=dim(training.table)[1], ncol=(bootstrap.sample.model.cfm*3)+1)
    colnames(matrix.bootstrap.model.cfm)<-rep("na",(bootstrap.sample.model.cfm*3)+1)
    matrix.bootstrap.model.cfm[,1]<-identification.value
    colnames(matrix.bootstrap.model.cfm)[1]<-"ID"
    name.sel.run<-paste(rep("ID_Selection_Run",bootstrap.sample.model.cfm),1:bootstrap.sample.model.cfm,sep="_")
    colnames(matrix.bootstrap.model.cfm)[seq(2,(bootstrap.sample.model.cfm*3)-1,3)]<-name.sel.run
    name.prob.run<-paste(rep("Probability_Run",bootstrap.sample.model.cfm),1:bootstrap.sample.model.cfm,sep="_")
    colnames(matrix.bootstrap.model.cfm)[seq(3,(bootstrap.sample.model.cfm*3),3)]<-name.prob.run
    name.pred.run<-paste(rep("Prediction_Run",bootstrap.sample.model.cfm),1:bootstrap.sample.model.cfm,sep="_")
    colnames(matrix.bootstrap.model.cfm)[seq(4,(bootstrap.sample.model.cfm*3)+1,3)]<-name.pred.run
    
    selection.index<-NULL
    #library(Zelig)
    #Bootstrap procedure
    for (count.boot in 1:bootstrap.sample.model.cfm)
    {
      selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
      matrix.bootstrap.model.cfm[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
      data.variables.bootstrap.model.cfm<-forecasting.combined.variables[selection.index,]
      explanatory.variables.bootstrap.model.cfm<-forecasting.combined.variables[selection.index,2:dim(forecasting.combined.variables)[2]]
      grouping.variable.bootstrap.model.cfm<-as.factor(forecasting.combined.variables[selection.index,1])
      #result.bootstrap.model.cfm<-zelig(as.formula(paste(names(data.variables.bootstrap.model.cfm)[1],"~",paste(names(data.variables.bootstrap.model.cfm[,2:dim(data.variables.bootstrap.model.cfm)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.cfm, model="logit",cite=FALSE)
      #result.bootstrap.model.cfm<-glm(as.formula(paste(names(data.variables.bootstrap.model.cfm)[1],"~",paste(names(data.variables.bootstrap.model.cfm[,2:dim(data.variables.bootstrap.model.cfm)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.cfm,family=binomial())
      while(inherits(try(result.bootstrap.model.cfm<-glm(as.formula(paste(names(data.variables.bootstrap.model.cfm)[1],"~",paste(names(data.variables.bootstrap.model.cfm[,2:dim(data.variables.bootstrap.model.cfm)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.cfm,family=binomial()),silent=TRUE),what="try-error"))
      {
        print(paste("Count boot: ",count.boot," - Boostrap while resampling",sep=""))
        selection.index<-sample(1:dim(training.table)[1], replace=TRUE, prob=NULL)
        matrix.bootstrap.model.cfm[as.numeric(names(table(selection.index))),(count.boot*3)-1]<-table(selection.index)
        explanatory.variables.bootstrap.model.cfm<-training.table[selection.index,3:dim(training.table)[2]]
        if(bootstrap_constant_correction==TRUE)
        {
          print("Performing bootstrap 0 value correction")
          indexbootstrapcosntant<-which(explanatory.variables.bootstrap.model.cfm==0,arr.ind=TRUE)
          explanatory.variables.bootstrap.model.cfm[indexbootstrapcosntant]<-runif(length(indexbootstrapcosntant)/2, min = 0.00001, max = 0.01)
        }
        grouping.variable.bootstrap.model.cfm<-as.factor(training.table[selection.index,2])
      }
      excluded.variables.bootstrap.model.cfm<-which(match(result.bootstrap.model.cfm$coefficients,NA)==1)
      if (length(excluded.variables.bootstrap.model.cfm) != 0)
      {
        data.variables.bootstrap.model.cfm.selected<-data.variables.bootstrap.model.cfm[,-excluded.variables.bootstrap.model.cfm]
        #setx.data.prediction<-forecasting.combined.variables[,-excluded.variables.bootstrap.model.cfm]
      } else
      {
        data.variables.bootstrap.model.cfm.selected<-data.variables.bootstrap.model.cfm
        #setx.data.prediction<-forecasting.combined.variables
      }
      #result.bootstrap.model.cfm.selected<-zelig(as.formula(paste(names(data.variables.bootstrap.model.cfm.selected)[1],"~",paste(names(data.variables.bootstrap.model.cfm.selected[,2:dim(data.variables.bootstrap.model.cfm.selected)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.cfm.selected, model="logit",cite=FALSE)
      result.bootstrap.model.cfm.selected<-glm(as.formula(paste(names(data.variables.bootstrap.model.cfm.selected)[1],"~",paste(names(data.variables.bootstrap.model.cfm.selected[,2:dim(data.variables.bootstrap.model.cfm.selected)[2]]),collapse= "+"))), data=data.variables.bootstrap.model.cfm.selected,family=binomial())
      #x.result.bootstrap.model.cfm.selected.prediction<-setx(result.bootstrap.model.cfm.selected,data=setx.data.prediction,fn=NULL)
      #matrix.bootstrap.model.cfm[,(count.boot*3)+1]<-sim(result.bootstrap.model.cfm.selected,x=x.result.bootstrap.model.cfm.selected.prediction)$result$fitted.values
      #matrix.bootstrap.model.cfm[as.numeric(names(table(selection.index))),(count.boot*3)]<- matrix.bootstrap.model.cfm[as.numeric(names(table(selection.index))),(count.boot*3)+1]
      matrix.bootstrap.model.cfm[,(count.boot*3)+1]<-predict(result.bootstrap.model.cfm.selected,newdata=forecasting.combined.variables[,-1],type="response")
      matrix.bootstrap.model.cfm[as.numeric(names(table(selection.index))),(count.boot*3)]<-matrix.bootstrap.model.cfm[as.numeric(names(table(selection.index))),(count.boot*3)+1]
    }
    
    # Export of bootstrap sample
    write.table(matrix.bootstrap.model.cfm,file="result_CFM_BootstrapSamples.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    
    ID.bootstrap.model.cfm.count<-numeric(length=dim(training.table)[1])
    #Probability (selected values)
    bootstrap.model.cfm.probability.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.probability.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.probability.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.probability.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.probability.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.probability.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #Prediction (all values)
    bootstrap.model.cfm.prediction.mean<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.prediction.sd<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.prediction.min<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.prediction.max<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.prediction.sderror<-numeric(length=dim(training.table)[1])
    bootstrap.model.cfm.prediction.quantiles<-matrix(nrow=dim(training.table)[1],ncol=7)
    
    #    for (count.row.variability in 1:dim(training.table)[1])
    #        {
    #        # Statistics on boostrapped probability
    #        ID.bootstrap.model.cfm.count[count.row.variability]<-length(na.omit(matrix.bootstrap.model.cfm[count.row.variability,seq(2,(bootstrap.sample.model.cfm*3)-1,3)]))
    #        bootstrap.model.cfm.probability.mean[count.row.variability]<-mean(na.omit(matrix.bootstrap.model.cfm[count.row.variability,seq(3,(bootstrap.sample.model.cfm*3),3)]))
    #        bootstrap.model.cfm.probability.sd[count.row.variability]<-sd(na.omit(matrix.bootstrap.model.cfm[count.row.variability,seq(3,(bootstrap.sample.model.cfm*3),3)]))
    #        bootstrap.model.cfm.probability.min[count.row.variability]<-min(na.omit(matrix.bootstrap.model.cfm[count.row.variability,seq(3,(bootstrap.sample.model.cfm*3),3)]))
    #        bootstrap.model.cfm.probability.max[count.row.variability]<-max(na.omit(matrix.bootstrap.model.cfm[count.row.variability,seq(3,(bootstrap.sample.model.cfm*3),3)]))
    #        bootstrap.model.cfm.probability.sderror[count.row.variability]<-bootstrap.model.cfm.probability.sd[count.row.variability]/ID.bootstrap.model.cfm.count[count.row.variability]
    #        bootstrap.model.cfm.probability.quantiles[count.row.variability,]<-quantile(na.omit(matrix.bootstrap.model.cfm[count.row.variability,seq(3,(bootstrap.sample.model.cfm*3),3)]),probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        # Statistics on boostrapped prediction
    #        bootstrap.model.cfm.prediction.mean[count.row.variability]<-mean(matrix.bootstrap.model.cfm[count.row.variability,seq(4,(bootstrap.sample.model.cfm*3)+1,3)])
    #        bootstrap.model.cfm.prediction.sd[count.row.variability]<-sd(matrix.bootstrap.model.cfm[count.row.variability,seq(4,(bootstrap.sample.model.cfm*3)+1,3)])
    #        bootstrap.model.cfm.prediction.min[count.row.variability]<-min(matrix.bootstrap.model.cfm[count.row.variability,seq(4,(bootstrap.sample.model.cfm*3)+1,3)])
    #        bootstrap.model.cfm.prediction.max[count.row.variability]<-max(matrix.bootstrap.model.cfm[count.row.variability,seq(4,(bootstrap.sample.model.cfm*3)+1,3)])
    #        bootstrap.model.cfm.prediction.sderror[count.row.variability]<-bootstrap.model.cfm.prediction.sd[count.row.variability]/bootstrap.sample.model.cfm
    #        bootstrap.model.cfm.prediction.quantiles[count.row.variability,]<-quantile(matrix.bootstrap.model.cfm[count.row.variability,seq(4,(bootstrap.sample.model.cfm*3)+1,3)],probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
    #        }
    
    fun_length<-function(x) {length_var<-length(x[which(is.finite(x))]); return(length_var)}
    fun_quantile<-function(x) {quantile_var<-t(quantile(x[which(is.finite(x))],probs=c(0,0.05,0.25,0.5,0.75,0.95,1))); return(quantile_var)}
    ID.bootstrap.model.cfm.count<-apply(matrix.bootstrap.model.cfm[,grep("ID_Selection",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=fun_length)
    bootstrap.model.cfm.probability.mean<-apply(matrix.bootstrap.model.cfm[,grep("Probability",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.cfm.probability.sd<-apply(matrix.bootstrap.model.cfm[,grep("Probability",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.cfm.probability.min<-apply(matrix.bootstrap.model.cfm[,grep("Probability",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.cfm.probability.max<-apply(matrix.bootstrap.model.cfm[,grep("Probability",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.cfm.probability.sderror<-bootstrap.model.cfm.probability.sd/bootstrap.sample.model.cfm
    bootstrap.model.cfm.probability.quantiles<-apply(matrix.bootstrap.model.cfm[,grep("Probability",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=fun_quantile)
    bootstrap.model.cfm.prediction.mean<-apply(matrix.bootstrap.model.cfm[,grep("Prediction",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=mean,na.rm = TRUE)
    bootstrap.model.cfm.prediction.sd<-apply(matrix.bootstrap.model.cfm[,grep("Prediction",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=sd,na.rm = TRUE)
    bootstrap.model.cfm.prediction.min<-apply(matrix.bootstrap.model.cfm[,grep("Prediction",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=min,na.rm = TRUE)
    bootstrap.model.cfm.prediction.max<-apply(matrix.bootstrap.model.cfm[,grep("Prediction",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=max,na.rm = TRUE)
    bootstrap.model.cfm.prediction.sderror<-bootstrap.model.cfm.prediction.sd/bootstrap.sample.model.cfm
    bootstrap.model.cfm.prediction.quantiles<-apply(matrix.bootstrap.model.cfm[,grep("Prediction",colnames(matrix.bootstrap.model.cfm))],MARGIN=1,FUN=fun_quantile)
    
    # Export of bootstrap sample statistics
    write.table(cbind("ID","CFM_NumberSelectedSamples","CFM_Probability_Mean","CFM_Probability_Sd","CFM_Probability_Min","CFM_Probability_Max","CFM_Probability_Sderror","CFM_Probability_Quantiles_0","CFM_Probability_Quantiles_0.05","CFM_Probability_Quantiles_0.25","CFM_Probability_Quantiles_0.5","CFM_Probability_Quantiles_0.75","CFM_Probability_Quantiles_0.95","CFM_Probability_Quantiles_1","CFM_Prediction_Mean","CFM_Prediction_Sd","CFM_Prediction_Min","CFM_Prediction_Max","CFM_Prediction_Sderror","CFM_Prediction_Quantiles_0","CFM_Prediction_Quantiles_0.05","CFM_Prediction_Quantiles_0.25","CFM_Prediction_Quantiles_0.5","CFM_Prediction_Quantiles_0.75","CFM_Prediction_Quantiles_0.95","CFM_Prediction_Quantiles_1"),file="result_CFM_BootstrapStatistics.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(cbind(identification.value,ID.bootstrap.model.cfm.count,bootstrap.model.cfm.probability.mean,bootstrap.model.cfm.probability.sd,bootstrap.model.cfm.probability.min,bootstrap.model.cfm.probability.max,bootstrap.model.cfm.probability.sderror,t(bootstrap.model.cfm.probability.quantiles),bootstrap.model.cfm.prediction.mean,bootstrap.model.cfm.prediction.sd,bootstrap.model.cfm.prediction.min,bootstrap.model.cfm.prediction.max,bootstrap.model.cfm.prediction.sderror,t(bootstrap.model.cfm.prediction.quantiles)),file="result_CFM_BootstrapStatistics.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    
    if (enable_screen_plotting==TRUE)
    {
      #dev.new()
      #double.sd.histogram.variability<-hist(bootstrap.model.cfm.probability.sd*2,breaks=seq(0,1,0.05),labels=TRUE)
      #plot(double.sd.histogram.variability$counts, seq(0,0.95,0.05), type="S",ylim=c(0,1), labels=TRUE)
      dev.new()
      plot(bootstrap.model.cfm.probability.mean,bootstrap.model.cfm.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="CFM BOOTSTRAP: Mean Probability vs Mean Prediction")
      abline(a=0,b=1,col="red",lty=1,lwd=1)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.cfm,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    }
    pdf(file = "result_CFM_BootstrapMeansComparison.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(bootstrap.model.cfm.probability.mean,bootstrap.model.cfm.prediction.mean,xlab="Probability mean",ylab="Prediction mean", type="p",main="CFM BOOTSTRAP: Mean Probability vs Mean Prediction")
    abline(a=0,b=1,col="red",lty=1,lwd=1)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.cfm,sep=""),side=3, padj=-0.5, adj=0.5, col="red",cex=0.8)
    dev.off()
    
    
    # BOOTSTRAPPED PROBABILITY - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.probability.cfm<-cbind(bootstrap.model.cfm.probability.mean,2*bootstrap.model.cfm.probability.sd)
    parabola.probability.cfm<-na.omit(parabola.probability.cfm[order(parabola.probability.cfm[,1]),])
    colnames(parabola.probability.cfm)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.probability.cfm <- nls(parabola.probability.cfm[,"ordinate"] ~ coeff.a*(parabola.probability.cfm[,"abscissa"]^2) + (-1)*coeff.a*parabola.probability.cfm[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.probability.cfm<-predict(fit.parabola.probability.cfm)
    #coef(fit.parabola.probability.cfm)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.probability.cfm[,"abscissa"],parabola.probability.cfm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="CFM Model Probability Variability (Bootstrap)")
      lines(parabola.probability.cfm[,"abscissa"],value.parabola.probability.cfm,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.cfm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.cfm),3),coeff.b= -round(coef(fit.parabola.probability.cfm),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    pdf(file = "result_CFM_BootstrapProbabilityVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.probability.cfm[,"abscissa"],parabola.probability.cfm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped probability mean",ylab="2 Standard Deviations", type="p",main="CFM Model Probability Variability (Bootstrap)")
    lines(parabola.probability.cfm[,"abscissa"],value.parabola.probability.cfm,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.cfm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.probability.cfm),3),coeff.b= -round(coef(fit.parabola.probability.cfm),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
    
    # BOOTSTRAPPED PREDICTION - Fit parabola 3 parameter y = ax^2 + bx + c
    parabola.prediction.cfm<-cbind(bootstrap.model.cfm.prediction.mean,2*bootstrap.model.cfm.prediction.sd)
    parabola.prediction.cfm<-parabola.prediction.cfm[order(parabola.prediction.cfm[,1]),]
    colnames(parabola.prediction.cfm)<-c("abscissa","ordinate")
    
    #If y has to be 0 in x=0 and x=1, this means that c=0 and a+b=0, so in our case since a<0, a has to be equal to -b
    fit.parabola.prediction.cfm <- nls(parabola.prediction.cfm[,"ordinate"] ~ coeff.a*(parabola.prediction.cfm[,"abscissa"]^2) + (-1)*coeff.a*parabola.prediction.cfm[,"abscissa"], start = c("coeff.a"=-1), control=list(maxiter=1000))
    value.parabola.prediction.cfm<-predict(fit.parabola.prediction.cfm)
    #coef(fit.parabola.prediction.cfm)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      plot(parabola.prediction.cfm[,"abscissa"],parabola.prediction.cfm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="CFM Model Prediction Variability (Bootstrap)")
      lines(parabola.prediction.cfm[,"abscissa"],value.parabola.prediction.cfm,col="red",lwd=1.5)
      mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.cfm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
      espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
      list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.cfm),3),coeff.b= -round(coef(fit.parabola.prediction.cfm),3))
      as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
      mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    }
    pdf(file = "result_CFM_BootstrapPredictionVariability.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(parabola.prediction.cfm[,"abscissa"],parabola.prediction.cfm[,"ordinate"],xlim=c(0,1),ylim=c(0,1),xlab="Bootstrapped prediction mean",ylab="2 Standard Deviations", type="p",main="CFM Model Prediction Variability (Bootstrap)")
    lines(parabola.prediction.cfm[,"abscissa"],value.parabola.prediction.cfm,col="red",lwd=1.5)
    mtext(paste("Number of bootstrap samples: ",bootstrap.sample.model.cfm,sep=""),side=3, padj=-0.5, adj=0.5, col="blue",cex=1)
    espr <- expression(Y == coeff.a %*% X ^2 + coeff.b %*% X)
    list.espr.subs <- list(coeff.a = round(coef(fit.parabola.prediction.cfm),3),coeff.b= -round(coef(fit.parabola.prediction.cfm),3))
    as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]])
    mtext(as.expression(do.call(substitute, list(as.call(espr), list.espr.subs))[[1]]),side=1, padj=-1.5, adj=0.5,col="red",cex=1)
    dev.off()
  }
  
  ## Sensitivity, Specificity, Cohens kappa plot
  dev.new()
  roc.plot.cfm.series<-roc.plot(verification.results.cfm,binormal=bin,plot=FALSE,show.thres=plot_thres,thresholds=threshold_series)
  dev.off()
  #str(roc.plot.cfm.series)
  #roc.plot.cfm.series$plot.data
  #str(roc.plot.cfm.series$plot.data)
  
  ###########################################
  # min(abs(TPR - (1-FPR)))
  if(enable_probability_optimal_binary_classification==FALSE)
  {
    dev.new()
    roc.plot.cfm.series<-roc.plot(verification.results.cfm,binormal=bin,plot=FALSE,show.thres=FALSE,thresholds=threshold_series)
    dev.off()
  } else
  {
    dev.new()
    roc.plot.cfm.series<-roc.plot(verification.results.cfm,binormal=bin,plot=FALSE,show.thres=FALSE)
    dev.off()  
  }
  
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    cfm.probability.classification.optimal<-data.frame(prob_thres=roc.plot.cfm.series$plot.data[,1,1],tpr=roc.plot.cfm.series$plot.data[,2,1],fpr=roc.plot.cfm.series$plot.data[,3,1],tnr=(1-roc.plot.cfm.series$plot.data[,3,1]),diff_abs_tpr_tnr=abs(roc.plot.cfm.series$plot.data[,2,1]-(1-roc.plot.cfm.series$plot.data[,3,1])),optimal_sel=NA,breaks_sel=NA)
    index.cfm.filter<-which(cfm.probability.classification.optimal$prob_thres>0 & cfm.probability.classification.optimal$prob_thres<1) # removing strnge thresh values
    cfm.probability.classification.optimal<-rbind(c(0,1,1,0,1,NA,NA),cfm.probability.classification.optimal[index.cfm.filter,],c(1,0,0,1,1,NA,NA))
    cfm.optimal.index<-which(cfm.probability.classification.optimal$diff_abs_tpr_tnr==min(cfm.probability.classification.optimal$diff_abs_tpr_tnr))
    cfm.probability.classification.optimal$optimal_sel[cfm.optimal.index]<-TRUE
    cfm.probability.optimal.binary.threshold<-cfm.probability.classification.optimal$prob_thres[cfm.optimal.index]
    
    ### Generating the optimal fourfold plot
    cross.classification.cfm.optimal<-table(grouping.variable,predict(result.cfm, type="response")>cfm.probability.optimal.binary.threshold,dnn=c("Observed","Predicted"))
    rownames(cross.classification.cfm.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.cfm.optimal)<-list("No Landslide","Landslide") # Predicted    
    str(cross.classification.cfm.optimal)
    # Assignation of a matching code between observed and predicted values
    result.cfm.matching.code.optimal<-paste(grouping.variable,as.numeric(predict(result.cfm, type="response")>cfm.probability.optimal.binary.threshold),sep="")
    result.cfm.matching.code.optimal<-gsub("00","1",result.cfm.matching.code.optimal)
    result.cfm.matching.code.optimal<-gsub("01","2",result.cfm.matching.code.optimal)
    result.cfm.matching.code.optimal<-gsub("10","3",result.cfm.matching.code.optimal)
    result.cfm.matching.code.optimal<-gsub("11","4",result.cfm.matching.code.optimal)
    result.cfm.matching.code.optimal<-as.numeric(result.cfm.matching.code.optimal)
    
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.cfm.optimal/sum(cross.classification.cfm.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    }
    # EXPORT OF PLOT FOR CFM MODEL
    pdf(file = "result_CFM_FourfoldPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.cfm.optimal/sum(cross.classification.cfm.optimal)*100,2), std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    #fourfold(cross.classification.cfm.optimal, std="margin", main="LINEAR DISCRIMINANT ANALYSIS MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(150,220,255,max=255), rgb(0,0,128,max=255)))
    dev.off()
    
    ### Optimal susceptibility classes identification 
    if(enable_probability_optimal_classification==TRUE)
    {
      cfm.unexplained.errors<-round(1-(cfm.probability.classification.optimal$tpr[cfm.optimal.index]+cfm.probability.classification.optimal$tpr[cfm.optimal.index])/2,2)
      
      if(type_probability_optimal_classification=="proportional")
      {
        cfm.unexplained.errors.partition<-1-(cfm.unexplained.errors/((length(breaks.histogram.values)/2)))*(1:((length(breaks.histogram.values)/2)))
        
        for(count_part in 1:(length(cfm.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-cfm.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((cfm.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          cfm.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((cfm.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          cfm.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        cfm.breaks.histogram.values.optimal<-c(0,cfm.probability.classification.optimal$prob_thres[which(cfm.probability.classification.optimal$breaks_sel==TRUE)],1)
        #cfm.probability.classification.optimal[which(cfm.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if(type_probability_optimal_classification=="fixed")
      {
        
        step.cfm.unexplained.fixed<-0.1
        if(cfm.unexplained.errors<=0.1) step.cfm.unexplained.fixed<-0.05
        if(cfm.unexplained.errors<=0.05) step.cfm.unexplained.fixed<-0.025
        cfm.unexplained.errors.partition<-seq(step.cfm.unexplained.fixed,1-step.cfm.unexplained.fixed,step.cfm.unexplained.fixed)[seq(step.cfm.unexplained.fixed,1-step.cfm.unexplained.fixed,step.cfm.unexplained.fixed)>(1-cfm.unexplained.errors)]
        
        for(count_part in 1:(length(cfm.unexplained.errors.partition)-1))
        {
          #count_part<-1
          unexplained.errors.partition.sel<-cfm.unexplained.errors.partition[count_part]
          index_tpr_sel<-max(which((cfm.probability.classification.optimal$tpr>=unexplained.errors.partition.sel)))
          cfm.probability.classification.optimal$breaks_sel[index_tpr_sel]<-TRUE
          index_tnr_sel<-min(which((cfm.probability.classification.optimal$tnr>=unexplained.errors.partition.sel)))
          cfm.probability.classification.optimal$breaks_sel[index_tnr_sel]<-TRUE
        }
        cfm.breaks.histogram.values.optimal<-c(0,cfm.probability.classification.optimal$prob_thres[which(cfm.probability.classification.optimal$breaks_sel==TRUE)],1)
        #cfm.probability.classification.optimal[which(cfm.probability.classification.optimal$breaks_sel==TRUE),]
      }
      
      if (enable_screen_plotting==TRUE)
      {
        dev.new()
        hist(predict(result.cfm, type="response"), breaks=cfm.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal CFM susceptibility", col=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1))
        plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("CFM OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
        mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
        mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
        for (count in 1:(length(cfm.breaks.histogram.values.optimal)-1))
        {
          #count=1
          polygon(c(cfm.breaks.histogram.values.optimal[count:(count+1)],rev(cfm.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1)[count])  
        }
        polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
        lines(cfm.probability.classification.optimal$prob_thres,cfm.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
        lines(cfm.probability.classification.optimal$prob_thres,cfm.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
        index_points_plot<-c(1,which(cfm.probability.classification.optimal$breaks_sel==TRUE),dim(cfm.probability.classification.optimal)[1])
        points(cfm.probability.classification.optimal[index_points_plot[1:floor(length(cfm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
        text(cfm.probability.classification.optimal[index_points_plot[1:floor(length(cfm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(cfm.probability.classification.optimal[index_points_plot[1:floor(length(cfm.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(cfm.breaks.histogram.values.optimal)/2)-1)))
        points(cfm.probability.classification.optimal[index_points_plot[ceiling(length(cfm.breaks.histogram.values.optimal)/2):length(cfm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
        text(cfm.probability.classification.optimal[index_points_plot[ceiling(length(cfm.breaks.histogram.values.optimal)/2):length(cfm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(cfm.probability.classification.optimal[index_points_plot[ceiling(length(cfm.breaks.histogram.values.optimal)/2):length(cfm.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(cfm.breaks.histogram.values.optimal)-ceiling(length(cfm.breaks.histogram.values.optimal)/2)),3))
      }
      
      pdf(file = "result_CFM_Histogram_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      hist(predict(result.cfm, type="response"), breaks=cfm.breaks.histogram.values.optimal,freq=TRUE, xlab="Susceptibility Class", ylab="Frequency", main="Histogram of optimal CFM susceptibility", col=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1))
      dev.off()
      
      pdf(file = "result_CFM_ModelEvaluationPlot_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main=paste("CFM OPTIMAL MODEL EVALUATION PLOT: ",type_probability_optimal_classification,sep=""))
      mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="dark red",cex=0.8)
      mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="navy blue",cex=0.8)
      for (count in 1:(length(cfm.breaks.histogram.values.optimal)-1))
      {
        #count=1
        polygon(c(cfm.breaks.histogram.values.optimal[count:(count+1)],rev(cfm.breaks.histogram.values.optimal[count:(count+1)])),c(0,0,1,1),border="darkgray",lty="dotted",lwd=0.5,col=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1)[count])  
      }
      polygon(c(0,1,1,0),c(0,0,1,1),border="black",lty="solid",lwd=1,col=NULL)
      lines(cfm.probability.classification.optimal$prob_thres,cfm.probability.classification.optimal$tpr,lty=1,lwd=2,col="dark red")
      lines(cfm.probability.classification.optimal$prob_thres,cfm.probability.classification.optimal$tnr,lty=1,lwd=2,col="navy blue")
      index_points_plot<-c(1,which(cfm.probability.classification.optimal$breaks_sel==TRUE),dim(cfm.probability.classification.optimal)[1])
      points(cfm.probability.classification.optimal[index_points_plot[1:floor(length(cfm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(cfm.probability.classification.optimal[index_points_plot[1:floor(length(cfm.breaks.histogram.values.optimal)/2)],c("prob_thres","tpr")],labels=with(round(cfm.probability.classification.optimal[index_points_plot[1:floor(length(cfm.breaks.histogram.values.optimal)/2)],],3), paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(3,rep(2,floor(length(cfm.breaks.histogram.values.optimal)/2)-1)))
      points(cfm.probability.classification.optimal[index_points_plot[1+ceiling(length(cfm.breaks.histogram.values.optimal)/2):length(cfm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],pch=19,cex=1,col="black")
      text(cfm.probability.classification.optimal[index_points_plot[1+ceiling(length(cfm.breaks.histogram.values.optimal)/2):length(cfm.breaks.histogram.values.optimal)],c("prob_thres","tnr")],labels=with(round(cfm.probability.classification.optimal[index_points_plot[1+ceiling(length(cfm.breaks.histogram.values.optimal)/2):length(cfm.breaks.histogram.values.optimal)],],3),paste("(",prob_thres,";",tnr,")",sep="")),cex=0.7,pos=c(rep(4,length(cfm.breaks.histogram.values.optimal)-1-ceiling(length(cfm.breaks.histogram.values.optimal)/2)),3))
      points(cfm.probability.classification.optimal[cfm.optimal.index,c("prob_thres","tpr")],pch=19,cex=1,col="black")
      text(cfm.probability.classification.optimal[cfm.optimal.index,c("prob_thres","tpr")],labels=with(round(cfm.probability.classification.optimal[cfm.optimal.index,c("prob_thres","tpr")],3),paste("(",prob_thres,";",tpr,")",sep="")),cex=0.7,pos=c(4))
      dev.off()
    }
  }
  
  ###########################################
  
  
  
  contingency.table.matrix.cfm<-matrix(nrow=dim(roc.plot.cfm.series$plot.data)[1],ncol=8)
  colnames(contingency.table.matrix.cfm)<-c("Threshold","TP","TN","FP","FN","TPR","FPR","COHEN_KAPPA")
  contingency.table.matrix.cfm[,1]<-roc.plot.cfm.series$plot.data[,1,1]
  contingency.table.matrix.cfm[,6]<-roc.plot.cfm.series$plot.data[,2,1]
  contingency.table.matrix.cfm[,7]<-roc.plot.cfm.series$plot.data[,3,1]
  values.observed<-training.table[,2]
  #values.predicted<-result.cfm$result$fitted.values
  values.predicted<-predict(result.cfm,type="response")
  
  for (count.threshold.series in 1:dim(roc.plot.cfm.series$plot.data)[1])
  {
    value.threshold<-contingency.table.matrix.cfm[count.threshold.series,1]
    values.probability.reclassified<-NULL
    values.probability.reclassified<-as.numeric(values.predicted>value.threshold) 
    #sum(values.probability.reclassified-round(values.predicted)) # Check sum: It has to be 0 if threshold is equal to 1
    series.pasted<-paste(values.observed,values.probability.reclassified,sep="")
    series.pasted<-gsub("00","1",series.pasted)
    series.pasted<-gsub("01","2",series.pasted)
    series.pasted<-gsub("10","3",series.pasted)
    series.pasted<-gsub("11","4",series.pasted)
    series.pasted<-as.numeric(series.pasted)
    TP<-as.numeric(sum(series.pasted>=4)) # True Positive
    FN<-as.numeric(sum(series.pasted>=3 & series.pasted<4)) # False Negative
    FP<-as.numeric(sum(series.pasted>=2 & series.pasted<3)) # False Positive
    TN<-as.numeric(sum(series.pasted>=1 & series.pasted<2)) # True Negative              
    #TPR<-TP/(TP+FN) # Hit Rate or True Positive Rate or Sensitivity - Assigned before the for cicle using rocplot data
    #FPR<-FP/(FP+TN) # False Alarm Rate or False Positive Rate or 1-Specificity
    # Cohen's Kappa = (agreement-chance)/(1-chance)  where agreement=(TP+TN)/(TP+TN+FP+FN) and chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    agreement=(TP+TN)/(TP+TN+FP+FN)
    chance=((((TN+FN)*(TN+FP))/(TP+TN+FP+FN))+(((TP+FP)*(TP+FN))/(TP+TN+FP+FN)))/(TP+TN+FP+FN)
    cohen.kappa.value<-(agreement-chance)/(1-chance)
    #Other
    #library(vcd)
    #cohen.kappa.value<-Kappa(cross.classification.table)
    contingency.table.matrix.cfm[count.threshold.series,2]<-TP
    contingency.table.matrix.cfm[count.threshold.series,3]<-TN
    contingency.table.matrix.cfm[count.threshold.series,4]<-FP
    contingency.table.matrix.cfm[count.threshold.series,5]<-FN
    contingency.table.matrix.cfm[count.threshold.series,8]<-cohen.kappa.value
  }
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    plot(roc.plot.cfm.series$plot.data[,1,1],roc.plot.cfm.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="CFM MODEL EVALUATION PLOT")
    lines(roc.plot.cfm.series$plot.data[,1,1],1-roc.plot.cfm.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
    lines(roc.plot.cfm.series$plot.data[,1,1], contingency.table.matrix.cfm[,8],col="blue",lty=1,lwd=1)
    mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
    mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
    mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  }
  pdf(file = "result_CFM_ModelEvaluationPlot.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(roc.plot.cfm.series$plot.data[,1,1],roc.plot.cfm.series$plot.data[,2,1],type="l",lty=1,lwd=1,col="red",xlim=c(0,1),ylim=c(0,1),xlab="Probability threshold",ylab="Evaluation parameter", main="CFM MODEL EVALUATION PLOT")
  lines(roc.plot.cfm.series$plot.data[,1,1],1-roc.plot.cfm.series$plot.data[,3,1],col="dark green",lty=1,lwd=1)
  lines(roc.plot.cfm.series$plot.data[,1,1], contingency.table.matrix.cfm[,8],col="blue",lty=1,lwd=1)
  mtext("SENSITIVITY",side=3, padj=-0.5, adj=0.01, col="red",cex=0.8)
  mtext("COHEN'S KAPPA",side=3, padj=-0.5, adj=0.5, col="blue",cex=0.8)
  mtext("SPECIFICITY",side=3, padj=-0.5, adj=0.99, col="dark green",cex=0.8)
  dev.off()
  
  ## VALIDATION OF CFM MODEL (Matching CFM posterior probability results and validation grouping variable)
  
  # Result Predicted
  #test.explanatory.variables.validation.cfm<-setx(result.cfm, data=forecasting.combined.variables.validation,fn=NULL)
  #predict.result.cfm.validation.posterior<-sim(result.cfm, x=test.explanatory.variables.validation.cfm)$result$fitted.values
  #predict.result.cfm.validation.class<-as.numeric(round(sim(result.cfm, x=test.explanatory.variables.validation.cfm)$result$fitted.values))
  predict.result.cfm.validation.posterior<-predict(result.cfm,newdata=forecasting.combined.variables.validation,type="response")
  predict.result.cfm.validation.class<-as.numeric(round(predict(result.cfm,newdata=forecasting.combined.variables.validation,type="response")))
  
  cross.classification.validation.cfm<-table(validation.table[,2],predict.result.cfm.validation.class,dnn=c("Observed","Predicted"))
  rownames(cross.classification.validation.cfm)<-list("No Landslide","Landslide") # Observed
  colnames(cross.classification.validation.cfm)<-list("No Landslide","Landslide") # Predicted
  #str(cross.classification.validation.cfm)
  #cross.classification.validation.cfm<-table(validation.grouping.variable,round(predict.result.cfm.validation.posterior),dnn=c("Observed","Predicted"))
  
  #Elaboration of Coefficient of association for contingency table
  #load package (vcd)
  library(vcd)
  
  #help(package=vcd)
  contingency.table.validation.cfm<-table2d_summary(cross.classification.validation.cfm)
  test.table.validation.cfm<-assocstats(cross.classification.validation.cfm)
  
  #Different plots for contingency table
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    fourfold(round(cross.classification.validation.cfm/sum(cross.classification.validation.cfm)*100,2), std="margin", main="VALIDATION CFM MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.cfm, std="margin", main="VALIDATION CFM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  }
  
  #Receiver Operating Characteristic (ROC) plots for one or more models.
  #load package (verification)
  library(verification)
  
  # 2nd method using verify function
  verification.validation.cfm<-verify(validation.table[,2],predict.result.cfm.validation.posterior, frcst.type="prob", obs.type="binary")
  #summary(verification.validation.cfm)
  
  # showing confidence intervals.  MAY BE SLOW
  area.under.roc.curve.validation.cfm<-roc.area(validation.table[,2],predict.result.cfm.validation.posterior)
  
  if (cross.classification.validation.cfm[1,2]==0 | cross.classification.validation.cfm[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
  
  if (enable_screen_plotting==TRUE)
  {
    dev.new()
    roc.plot(verification.validation.cfm, main = "ROC PLOT: VALIDATION CFM MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[5] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
    mtext(paste("ROC area = ",round(area.under.roc.curve.validation.cfm$A,2),";  Sample size = ",area.under.roc.curve.validation.cfm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[5], sep=""), side=3, col="red", cex=0.8)
  }
  
  # EXPORT OF PLOT FOR VALIDATION OF CFM MODEL
  
  pdf(file = "result_CFM_FourfoldPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  fourfold(round(cross.classification.validation.cfm/sum(cross.classification.validation.cfm)*100,2), std="margin", main="VALIDATION CFM MODEL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  #fourfold(cross.classification.validation.cfm, std="margin", main="VALIDATION CFM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
  dev.off()
  
  #pdf(file = "result_CFM_ROCPlot_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  #roc.plot(verification.validation.cfm, main = "ROC PLOT: VALIDATION CFM MODEL", binormal = TRUE, plot = "both", extra=TRUE, legend=TRUE)
  #area.under.roc.curve.validation.cfm<-roc.area(verification.table[,2],result.cfm.validation$fitted.values)
  #dev.off()
  
  pdf(file = "result_CFM_ROCPlot_bootstrap_Validation.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  roc.plot(verification.validation.cfm, main = "ROC PLOT: VALIDATION CFM MODEL", binormal=bin, plot=bin_plot, CI=enable_rocplot_confidenceinterval, n.boot=bootstrap.sample.values[5] , alpha = 0.05, extra=TRUE, legend=TRUE,show.thres=plot_thres,thresholds=threshold_series)
  mtext(paste("ROC area = ",round(area.under.roc.curve.validation.cfm$A,2),";  Sample size = ",area.under.roc.curve.validation.cfm$n.total,";  Bootstrap samples = ",bootstrap.sample.values[5], sep=""), side=3, col="red", cex=0.8)
  dev.off()
  
  # Assignation of a matching code between observed and predicted values calculated using the validation dataset
  validation.cfm.matching.code<-paste(validation.grouping.variable,round(predict.result.cfm.validation.posterior),sep="")
  validation.cfm.matching.code<-gsub("00","1",validation.cfm.matching.code)
  validation.cfm.matching.code<-gsub("01","2",validation.cfm.matching.code)
  validation.cfm.matching.code<-gsub("10","3",validation.cfm.matching.code)
  validation.cfm.matching.code<-gsub("11","4",validation.cfm.matching.code)
  validation.cfm.matching.code<-as.numeric(validation.cfm.matching.code)
  
  ##########################################
  if(enable_probability_optimal_binary_classification==TRUE)
  {
    cross.classification.validation.cfm.optimal<-table(validation.grouping.variable,as.numeric(predict.result.cfm.validation.posterior>cfm.probability.optimal.binary.threshold),dnn=c("Observed","Predicted"))
    rownames(cross.classification.validation.cfm.optimal)<-list("No Landslide","Landslide") # Observed
    colnames(cross.classification.validation.cfm.optimal)<-list("No Landslide","Landslide") # Predicted
    
    #Different plots for contingency table
    if (enable_screen_plotting==TRUE)
    {
      dev.new()
      fourfold(round(cross.classification.validation.cfm.optimal/sum(cross.classification.validation.cfm.optimal)*100,2), std="margin", main="VALIDATION CFM MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
      #fourfold(cross.classification.validation.cfm, std="margin", main="VALIDATION CFM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    }
    if (cross.classification.validation.cfm.optimal[1,2]==0 | cross.classification.validation.cfm.optimal[2,1]==0) {bin=FALSE; bin_plot="emp"; plot_thres=FALSE} else {bin=TRUE; bin_plot="both"; plot_thres=TRUE}
    pdf(file = "result_CFM_FourfoldPlot_Validation_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    fourfold(round(cross.classification.validation.cfm.optimal/sum(cross.classification.validation.cfm.optimal)*100,2), std="margin", main="VALIDATION CFM MODEL OPTIMAL", conf_level=0.95, extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255), rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    #fourfold(cross.classification.validation.cfm, std="margin", main="VALIDATION CFM MODEL", extended=TRUE, space = 0.2, margin=1, color = c(rgb(255,0,0,max=255), rgb(255,128,0,max=255), rgb(56,168,0,max=255), rgb(170,255,0,max=255),  rgb(170,135,210,max=255), rgb(115,70,155,max=255)))
    dev.off()
    
    
    # Assignation of a optimal matching code between observed and predicted values calculated using the validation dataset
    validation.cfm.matching.code.optimal<-paste(validation.grouping.variable,as.numeric(predict.result.cfm.validation.posterior>cfm.probability.optimal.binary.threshold),sep="")
    validation.cfm.matching.code.optimal<-gsub("00","1",validation.cfm.matching.code.optimal)
    validation.cfm.matching.code.optimal<-gsub("01","2",validation.cfm.matching.code.optimal)
    validation.cfm.matching.code.optimal<-gsub("10","3",validation.cfm.matching.code.optimal)
    validation.cfm.matching.code.optimal<-gsub("11","4",validation.cfm.matching.code.optimal)
    validation.cfm.matching.code.optimal<-as.numeric(validation.cfm.matching.code.optimal)
  }
  #########################################
  
  
  # EXPORT OF CFM MODEL RESULTS
  
  write.table("RESULTS OF COMBINATION FORECAST LOGISTIC REGRESSION MODEL",file="result_CFM.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CFM MODEL OUTPUTS",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Logistic Regression coefficients",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  #Scaling coefficients
  #write.table(cbind(names(result.cfm$results$coefficients),result.cfm$results$coefficients),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(names(result.cfm$coefficients),result.cfm$coefficients),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE MODEL RESULT",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.cfm$table[,1,],contingency.table.cfm$table[,2,],contingency.table.cfm$table[,3,])),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("CONTINGENCY TABLE VALIDATION",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(c("","No Landslide Predicted","Landslide Predicted","Total"),cbind(c("No Landslide Observed","Landslide Observed","Total"),contingency.table.validation.cfm$table[,1,],contingency.table.validation.cfm$table[,2,],contingency.table.validation.cfm$table[,3,])),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("MATCHING CODE DEFINITION",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(cbind(c("","OBSERVED NO LANDSLIDES: 0","OBSERVED LANDSLIDES: 1"), c("PREDICTED NO LANDSLIDES: 0","00 -> Code 1","10 -> Code 3"), c("PREDICTED LANDSLIDES: 1","01 -> Code 2","11 -> Code 4")),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  
  ########
  if(enable_probability_optimal_binary_classification==FALSE) 
  {
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","CFM MODEL POSTERIOR PROBABILITY","CFM MODEL CLASSIFICATION","CFM MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict(result.cfm, type="response"),round(predict(result.cfm, type="response")),result.cfm.matching.code)),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","CFM VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.cfm.validation.posterior,round(predict.result.cfm.validation.posterior),validation.cfm.matching.code)),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  } else
  {
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    if(enable_probability_optimal_classification==TRUE) 
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY PARTITION -> Method: ",type_probability_optimal_classification,sep=""),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(cfm.probability.classification.optimal[index_points_plot,c("prob_thres","tnr","tpr")]),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    } else
    {
      write.table(paste("OPTIMAL SUSCEPTIBILITY BINARY PARTITION",sep=""),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
      write.table(data.frame(cfm.probability.classification.optimal[cfm.optimal.index,c("prob_thres","tnr","tpr")]),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=TRUE)
    }
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","GROUPING VARIABLE","CFM MODEL POSTERIOR PROBABILITY","CFM MODEL CLASSIFICATION","CFM MODEL RESULT MATCHING CODE","CFM OPTIMAL MODEL CLASSIFICATION","CFM OPTIMAL MODEL RESULT MATCHING CODE"),cbind(identification.value,training.table[,2],predict(result.cfm, type="response"),round(predict(result.cfm, type="response")),result.cfm.matching.code,as.numeric(predict(result.cfm, type="response")>cfm.probability.optimal.binary.threshold),result.cfm.matching.code.optimal)),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table("FINAL RESULTS VALIDATION",file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
    write.table(rbind(c("ID","VALIDATION GROUPING VARIABLE","VALIDATION POSTERIOR PROBABILITY","VALIDATION CLASSIFICATION","CFM VALIDATION MATCHING CODE","OPTIMAL VALIDATION CLASSIFICATION","OPTIMAL CFM VALIDATION MATCHING CODE"),cbind(validation.table[,1],validation.table[,2],predict.result.cfm.validation.posterior,round(predict.result.cfm.validation.posterior),validation.cfm.matching.code,as.numeric(predict.result.cfm.validation.posterior>cfm.probability.optimal.binary.threshold),validation.cfm.matching.code.optimal)),file="result_CFM.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  }
  ########
  
  # PLOT AND EXPORT OF CFM MODEL MAPS
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POLYGONS")
  {
    ################################################
    ### Prima di cancellare testare con polygoni
    ##############################################
    #shape_training_cfm<-shape_training
    ##result_training_cfm_shape<-cbind(identification.value,result.cfm$result$y,result.cfm$result$fitted.values,round(result.cfm$result$fitted.values),result.cfm.matching.code)
    #result_training_cfm_shape<-cbind(identification.value,training.table[,2],predict(result.cfm, type="response"),round(predict(result.cfm, type="response")),result.cfm.matching.code)
    #colnames(result_training_cfm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH")
    #shape_training_cfm@data <- merge(x=shape_training_cfm@data,y=result_training_cfm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)  
    ## writeOGR(shape_training_cfm,dsn="result_CFM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_training_cfm,dsn="result_CFM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    #shape_validation_cfm<-shape_validation
    #result_validation_cfm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.cfm.validation.posterior,round(predict.result.cfm.validation.posterior),validation.cfm.matching.code)
    #colnames(result_validation_cfm_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    #shape_validation_cfm@data <- merge(x=shape_validation_cfm@data,y=result_validation_cfm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    ##writeOGR(shape_validation_cfm,dsn="result_CFM_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    #writeOGR(shape_validation_cfm,dsn="result_CFM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ### DA qui Nuovo
    ##################################################################   
    #### aggiungere a poligoni colonna incertezza ed esportazione pdf incertezza
    shape_training_cfm<-shape_training
    #result_training_cfm_shape<-cbind(identification.value,result.cfm$result$y,result.cfm$result$fitted.values,round(result.cfm$result$fitted.values),result.cfm.matching.code)
    result_training_cfm_shape<-cbind(identification.value,training.table[,2],predict(result.cfm, type="response"),round(predict(result.cfm, type="response")),result.cfm.matching.code,ID.bootstrap.model.cfm.count,bootstrap.model.cfm.probability.mean,bootstrap.model.cfm.probability.sd,bootstrap.model.cfm.probability.min,bootstrap.model.cfm.probability.max,bootstrap.model.cfm.probability.sderror,t(bootstrap.model.cfm.probability.quantiles),bootstrap.model.cfm.prediction.mean,bootstrap.model.cfm.prediction.sd,bootstrap.model.cfm.prediction.min,bootstrap.model.cfm.prediction.max,bootstrap.model.cfm.prediction.sderror,t(bootstrap.model.cfm.prediction.quantiles))
    colnames(result_training_cfm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","CFM_SAMP","CFM_PMean","CFM_PSd","CFM_PMin","CFM_PMax","CFM_PSder","CFM_PQ0","CFM_PQ_005","CFM_PQ_025","CFM_PQ05","CFM_PQ_075","CFM_PQ095","CFM_PQ1","CFM_PrMean","CFM_PrSd","CFM_PrMin","CFM_PrMax","CFM_PrSder","CFM_PrQ0","CFM_PrQ005","CFM_PrQ025","CFM_PrQ05","CFM_PrQ075","CFM_PrQ095","CFM_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_cfm_shape<-cbind(result_training_cfm_shape,as.numeric(predict(result.cfm, type="response")>cfm.probability.optimal.binary.threshold),result.cfm.matching.code.optimal)
      colnames(result_training_cfm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","CFM_SAMP","CFM_PMean","CFM_PSd","CFM_PMin","CFM_PMax","CFM_PSder","CFM_PQ0","CFM_PQ_005","CFM_PQ_025","CFM_PQ05","CFM_PQ_075","CFM_PQ095","CFM_PQ1","CFM_PrMean","CFM_PrSd","CFM_PrMin","CFM_PrMax","CFM_PrSder","CFM_PrQ0","CFM_PrQ005","CFM_PrQ025","CFM_PrQ05","CFM_PrQ075","CFM_PrQ095","CFM_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    shape_training_cfm@data <- merge(x=shape_training_cfm@data,y=result_training_cfm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_cfm,dsn="result_CFM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_cfm,dsn="result_CFM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    # WARNING: The validation does't have the unvertainty estimation: probably this can be daone using the parabolic error function 
    shape_validation_cfm<-shape_validation
    result_validation_cfm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.cfm.validation.posterior,round(predict.result.cfm.validation.posterior),validation.cfm.matching.code)
    colnames(result_validation_cfm_shape)<-c("ID","GROUP_VAR","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_cfm_shape<-cbind(result_validation_cfm_shape,as.numeric(predict.result.cfm.validation.posterior>cfm.probability.optimal.binary.threshold),validation.cfm.matching.code.optimal)
      colnames(result_validation_cfm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    shape_validation_cfm@data <- merge(x=shape_validation_cfm@data,y=result_validation_cfm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_cfm@data <- cbind(shape_validation_cfm@data,PROB_SDMOD=(coefficients(fit.parabola.probability.cfm)*(shape_validation_cfm@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.cfm)*shape_validation_cfm@data$VAL_PROB))
    #writeOGR(shape_validation_cfm,dsn="result_CFM_validation.shp",layer="validation",driver="ESRI Shapefile")  # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_cfm,dsn="result_CFM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    ##################################################################
    ##################################################################
    
    
    # Plot and export of maps
    # CFM Susceptibility
    #dev.new()
    pdf(file = "result_CFM_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_cfm, zcol=c("MOD_PROB"), names.attr=c("CFM MODEL PROBABILITY"), main="CFM MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # CFM Model Matching Code
    #dev.new()
    pdf(file = "result_CFM_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_cfm, zcol=c("MOD_MATCH"), names.attr=c("CFM MODEL MATCHING CODE"), main="CFM MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
    
    # CFM Model Uncertainty
    #dev.new()
    pdf(file = "result_CFM_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_training_cfm, zcol=c("CFM_PrSd"), names.attr=c("CFM MODEL UNCERTAINTY"), main="CFM MODEL UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()
    
    # CFM Validation Susceptibility
    #dev.new()
    pdf(file = "result_CFM_Validation_Susceptibility_Map.pdf")
    print(spplot(obj=shape_validation_cfm, zcol=c("VAL_PROB"), names.attr=c("CFM VALIDATION PROBABILITY"), main="CFM VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.susceptibility, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.susceptibility))
    dev.off()
    
    # CFM Validation Matching Code
    #dev.new()
    pdf(file = "result_CFM_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_cfm, zcol=c("VAL_MATCH"), names.attr=c("CFM VALIDATION MATCHING CODE"), main="CFM VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
    dev.off()
        
    # CFM Validation Uncertainty
    #dev.new()
    pdf(file = "result_CFM_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    print(spplot(obj=shape_validation_cfm, zcol=c("PROB_SDMOD"), names.attr=c("CFM VALIDATION UNCERTAINTY"), main="CFM VALIDATION UNCERTAINTY",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.uncertainty, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.uncertainty))
    dev.off()

    
    ########### TBT
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      # CFM Model Matching Code
      #dev.new()
      pdf(file = "result_CFM_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_training_cfm, zcol=c("OPT_MATCH"), names.attr=c("CFM MODEL MATCHING CODE"), main="CFM MODEL MATCHING CODE",  sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      # CFM Model Matching Code
      #dev.new()
      pdf(file = "result_CFM_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(spplot(obj=shape_validation_cfm, zcol=c("OPT_MATCH"), names.attr=c("CFM VALIDATION MATCHING CODE"), main="CFM VALIDATION MATCHING CODE",  sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=breaks.map.matching.code, regions=TRUE, colorkey=list(space="bottom"), col.regions=color.vector.matching))
      dev.off()
      
      if(enable_probability_optimal_classification==TRUE)
      {
        pdf(file = "result_CFM_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
        print(spplot(obj=shape_training_cfm, zcol=c("MOD_PROB"), names.attr=c("CFM MODEL PROBABILITY"), main="CFM MODEL PROBABILITY", sp.layout=list(arrow.training,scale.training,text.scale1.training,text.scale2.training), scales = list(draw = TRUE), at=(cfm.breaks.histogram.values.optimal)+c(rep(0,(length(cfm.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1)))
        dev.off()
        
        pdf(file = "result_CFM_Validation_Susceptibility_Map_Optimal.pdf")
        print(spplot(obj=shape_validation_cfm, zcol=c("VAL_PROB"), names.attr=c("CFM VALIDATION PROBABILITY"), main="CFM VALIDATION PROBABILITY", sp.layout=list(arrow.validation,scale.validation,text.scale1.validation,text.scale2.validation), scales = list(draw = TRUE), at=(cfm.breaks.histogram.values.optimal)+c(rep(0,(length(cfm.breaks.histogram.values.optimal)-1)),0.0001), regions=TRUE, colorkey=list(space="bottom"), col.regions=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1)))
        dev.off()
      }
    }
    ###########
    
    ### Success & prediction rate curve
    susc_values_cfm<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_cfm@data) %in% c(configuration.spatial.data.table[c(5,6)],"MOD_PROB"))
    shape_training_cfm@data<-shape_training_cfm@data[order(shape_training_cfm@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_training_cfm@data[,1]<-cumsum(shape_training_cfm@data[,1])/sum(shape_training_cfm@data[,1])*100
    shape_training_cfm@data[,2]<-cumsum(shape_training_cfm@data[,2])/sum(shape_training_cfm@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_cfm@data[,3],shape_training_cfm@data[,1])
    area_susc_values_cfm_mod<-c(0,ind_pro_fun_mod(susc_values_cfm[-c(1,length(susc_values_cfm))]),100)
    
    #dev.new()
    pdf(file = "result_CFM_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="CFM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_cfm)-1))
    {
      #count=1
      polygon(c(area_susc_values_cfm_mod[count:(count+1)],rev(area_susc_values_cfm_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_cfm@data[,1],shape_training_cfm@data[,2],col="black")
    dev.off()
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_cfm_mod_optimal<-c(0,ind_pro_fun_mod(rev(cfm.breaks.histogram.values.optimal)[-c(1,length(cfm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_CFM_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="CFM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(cfm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_cfm_mod_optimal[count:(count+1)],rev(area_susc_values_cfm_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_cfm@data[,1],shape_training_cfm@data[,2],col="black")
      dev.off()
    }
    ############
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_cfm@data) %in% c(configuration.spatial.data.table[c(5,6)],"VAL_PROB"))
    shape_validation_cfm@data<-shape_validation_cfm@data[order(shape_validation_cfm@data[,ind_col_succ_pred_rate[3]],decreasing=TRUE),ind_col_succ_pred_rate]
    shape_validation_cfm@data[,1]<-cumsum(shape_validation_cfm@data[,1])/sum(shape_validation_cfm@data[,1])*100
    shape_validation_cfm@data[,2]<-cumsum(shape_validation_cfm@data[,2])/sum(shape_validation_cfm@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_cfm@data[,3],shape_validation_cfm@data[,1])
    area_susc_values_cfm_val<-c(0,ind_pro_fun_val(susc_values_cfm[-c(1,length(susc_values_cfm))]),100)
    
    #dev.new()
    pdf(file = "result_CFM_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="CFM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_cfm)-1))
    {
      #count=1
      polygon(c(area_susc_values_cfm_val[count:(count+1)],rev(area_susc_values_cfm_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_cfm@data[,1],shape_validation_cfm@data[,2],col="black")
    dev.off()
    
    ########### TBT
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_cfm_val_optimal<-c(0,ind_pro_fun_val(rev(cfm.breaks.histogram.values.optimal)[-c(1,length(cfm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_CFM_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="CFM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(cfm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_cfm_val_optimal[count:(count+1)],rev(area_susc_values_cfm_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_cfm@data[,1],shape_validation_cfm@data[,2],col="black")
      dev.off()
    }
    ############
    
  }
  
  
  if(configuration.spatial.data.table$PRESENCE == "YES" & configuration.spatial.data.table$GEOMETRY=="POINTS")
  {
    shape_training_cfm<-shape_training
    #result_training_cfm_shape<-cbind(identification.value,training.table[,2],predict.result.cfm$posterior[,2],as.numeric(levels(predict.result.cfm$class))[predict.result.cfm$class],result.cfm.matching.code,ID.bootstrap.model.cfm.count,bootstrap.model.cfm.probability.mean,bootstrap.model.cfm.probability.sd,bootstrap.model.cfm.probability.min,bootstrap.model.cfm.probability.max,bootstrap.model.cfm.probability.sderror,t(bootstrap.model.cfm.probability.quantiles),bootstrap.model.cfm.prediction.mean,bootstrap.model.cfm.prediction.sd,bootstrap.model.cfm.prediction.min,bootstrap.model.cfm.prediction.max,bootstrap.model.cfm.prediction.sderror,t(bootstrap.model.cfm.prediction.quantiles))
    result_training_cfm_shape<-cbind(identification.value,training.table[,2],predict(result.cfm, type="response"),as.numeric(round(predict(result.cfm, type="response"))),result.cfm.matching.code,ID.bootstrap.model.cfm.count,bootstrap.model.cfm.probability.mean,bootstrap.model.cfm.probability.sd,bootstrap.model.cfm.probability.min,bootstrap.model.cfm.probability.max,bootstrap.model.cfm.probability.sderror,t(bootstrap.model.cfm.probability.quantiles),bootstrap.model.cfm.prediction.mean,bootstrap.model.cfm.prediction.sd,bootstrap.model.cfm.prediction.min,bootstrap.model.cfm.prediction.max,bootstrap.model.cfm.prediction.sderror,t(bootstrap.model.cfm.prediction.quantiles))
    colnames(result_training_cfm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","CFM_SAMP","CFM_PMean","CFM_PSd","CFM_PMin","CFM_PMax","CFM_PSder","CFM_PQ0","CFM_PQ_005","CFM_PQ_025","CFM_PQ05","CFM_PQ_075","CFM_PQ095","CFM_PQ1","CFM_PrMean","CFM_PrSd","CFM_PrMin","CFM_PrMax","CFM_PrSder","CFM_PrQ0","CFM_PrQ005","CFM_PrQ025","CFM_PrQ05","CFM_PrQ075","CFM_PrQ095","CFM_PrQ1")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_training_cfm_shape<-cbind(result_training_cfm_shape,as.numeric(predict(result.cfm, type="response")>cfm.probability.optimal.binary.threshold),result.cfm.matching.code.optimal)
      colnames(result_training_cfm_shape)<-c("ID","MOD_GROUP","MOD_PROB","MOD_CLASS","MOD_MATCH","CFM_SAMP","CFM_PMean","CFM_PSd","CFM_PMin","CFM_PMax","CFM_PSder","CFM_PQ0","CFM_PQ_005","CFM_PQ_025","CFM_PQ05","CFM_PQ_075","CFM_PQ095","CFM_PQ1","CFM_PrMean","CFM_PrSd","CFM_PrMin","CFM_PrMax","CFM_PrSder","CFM_PrQ0","CFM_PrQ005","CFM_PrQ025","CFM_PrQ05","CFM_PrQ075","CFM_PrQ095","CFM_PrQ1","OPT_CLASS","OPT_MATCH")
    }
    #############
    shape_training_cfm@data <- merge(x=shape_training_cfm@data,y=result_training_cfm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    #writeOGR(shape_training_cfm,dsn="result_CFM_training.shp",layer="training",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_training_cfm,dsn="result_CFM_training",layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    
    # WARNING: The validation does't have the unvertainty estimation: probabilty this can be daone using the parabolic error function 
    shape_validation_cfm<-shape_validation
    #result_validation_cfm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.cfm.validation$posterior[,2],as.numeric(levels(predict.result.cfm.validation$class))[predict.result.cfm.validation$class],validation.cfm.matching.code)
    result_validation_cfm_shape<-cbind(validation.table[,1],validation.table[,2],predict.result.cfm.validation.posterior,as.numeric(round(predict.result.cfm.validation.posterior)),validation.cfm.matching.code)
    colnames(result_validation_cfm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH")
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      result_validation_cfm_shape<-cbind(result_validation_cfm_shape,as.numeric(predict.result.cfm.validation.posterior>cfm.probability.optimal.binary.threshold),validation.cfm.matching.code.optimal)
      colnames(result_validation_cfm_shape)<-c("ID","VAL_GROUP","VAL_PROB","VAL_CLASS","VAL_MATCH","OPT_CLASS","OPT_MATCH")
    }
    ############
    shape_validation_cfm@data <- merge(x=shape_validation_cfm@data,y=result_validation_cfm_shape,by.x=shape_merge_field, by.y="ID", all.x=T, sort=F)	
    shape_validation_cfm@data <- cbind(shape_validation_cfm@data,PROB_SDMOD=(coefficients(fit.parabola.probability.cfm)*(shape_validation_cfm@data$VAL_PROB^2)) + ((-1)*coefficients(fit.parabola.probability.cfm)*shape_validation_cfm@data$VAL_PROB))
    #writeOGR(shape_validation_cfm,dsn="result_CFM_validation.shp",layer="validation",driver="ESRI Shapefile") # Version of rgdal older than 2.13.1
    writeOGR(shape_validation_cfm,dsn="result_CFM_validation",layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
    
    require(raster)
    
    # Plot and export of maps
    # CFM Susceptibility
    #dev.new()
    pdf(file = "result_CFM_Model_Susceptibility_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_cfm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_CFM_Model_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_CFM_Model_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1),breaks=round(cfm.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    
    # CFM Model Matching Code
    #dev.new()
    pdf(file = "result_CFM_Model_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_cfm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"MOD_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
    print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_CFM_Model_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_CFM_Model_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_cfm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      index_col_macthing<-as.numeric(names(table(layer_gridded_raster@data@values)))
      print(plot(layer_gridded_raster,col=color.vector.matching[index_col_macthing],legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_CFM_Model_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_CFM_Model_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_training_cfm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_CFM_Model_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
    }
    ############
    
    # CFM Model uncertainity
    #dev.new()
    pdf(file = "result_CFM_Model_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_training_cfm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"CFM_PrSd"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_CFM_Model_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    
    
    # CFM Validation Susceptibility
    #dev.new()
    pdf(file = "result_CFM_Validation_Susceptibility_Map.pdf")
    layer_gridded<-shape_validation_cfm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_PROB"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.susceptibility,breaks=round(breaks.map.susceptibility,2)))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_CFM_Validation_Susceptibility_Map.tif", format="GTiff", overwrite=TRUE)
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      pdf(file = "result_CFM_Validation_Susceptibility_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1),breaks=round(cfm.breaks.histogram.values.optimal,3)))
      dev.off()
    }
    ############
    
    # CFM Validation Matching Code
    #dev.new()
    pdf(file = "result_CFM_Validation_MatchingCode_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_cfm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"VAL_MATCH"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
    legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
    #zoom(layer_gridded_raster)  	
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_CFM_Validation_MatchingCode_Map.tif", format="GTiff", overwrite=TRUE)
    
    ########
    if(enable_probability_optimal_binary_classification==TRUE) 
    {
      #dev.new()
      pdf(file = "result_CFM_Validation_MatchingCode_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_cfm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_MATCH"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color.vector.matching,round(breaks.map.matching.code),legend=FALSE))
      legend("topright", legend = c("1: TN","2: FP","3: FN","4: TP"), cex=0.8,fill = color.vector.matching,xjust=1.1,bg="transparent",box.col="transparent")
      #zoom(layer_gridded_raster)    
      dev.off()
      writeRaster(layer_gridded_raster, filename="result_CFM_Validation_MatchingCode_Map_Optimal.tif", format="GTiff", overwrite=TRUE)

      #dev.new()
      pdf(file = "result_CFM_Validation_SusceptibilityBinary_Map_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      layer_gridded<-shape_validation_cfm
      layer_gridded@data<-as.data.frame(layer_gridded@data[,"OPT_CLASS"])
      gridded(layer_gridded)<-TRUE
      layer_gridded_raster<-raster(layer_gridded)
      #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
      res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
      print(plot(layer_gridded_raster,col=color_ramp_palette_fun(2),legend=FALSE))
      legend("topright", legend = c("0: Not susceptible","1: Susceptible"), cex=0.8,fill = color_ramp_palette_fun(2),xjust=1.1,bg="transparent",box.col="transparent")
      dev.off()
      
      writeRaster(layer_gridded_raster, filename="result_CFM_Validation_SusceptibilityBinary_Map_Optimal.tif", format="GTiff", overwrite=TRUE)
      
      
    }
    #######
    
    # CFM Model uncertainity
    #dev.new()
    pdf(file = "result_CFM_Validation_Uncertainty_Map.pdf",onefile = TRUE, pagecentre=TRUE)
    layer_gridded<-shape_validation_cfm
    layer_gridded@data<-as.data.frame(layer_gridded@data[,"PROB_SDMOD"])
    gridded(layer_gridded)<-TRUE
    layer_gridded_raster<-raster(layer_gridded)
    #res(layer_gridded_raster)<-gridparameters(layer_gridded)[1,2]
    res(layer_gridded_raster)<-as.numeric(configuration.spatial.data.table[c(8)])
    print(plot(layer_gridded_raster,col=color.vector.uncertainty,breaks.map.uncertainty))
    #zoom(layer_gridded_raster)		
    dev.off()
    
    writeRaster(layer_gridded_raster, filename="result_CFM_Validation_Uncertainty_Map.tif", format="GTiff", overwrite=TRUE)
    
    ### Success & prediction rate curve
    susc_values_cfm<-c(1,0.8,0.55,0.45,0.2,0)
    
    # ordering data for susceptibility
    ind_col_succ_pred_rate<-which(colnames(shape_training_cfm@data) %in% c("MOD_GROUP","MOD_PROB"))
    shape_training_cfm@data<-cbind(PIXEL_AREA=rep(configuration.spatial.data.table[c(8)]^2,dim(shape_training_cfm@data)[1]),shape_training_cfm@data[order(shape_training_cfm@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_training_cfm@data$MOD_GROUP<-shape_training_cfm@data$MOD_GROUP*configuration.spatial.data.table[c(8)]^2
    shape_training_cfm@data[,1]<-cumsum(shape_training_cfm@data[,1])/sum(shape_training_cfm@data[,1])*100
    shape_training_cfm@data[,2]<-cumsum(shape_training_cfm@data[,2])/sum(shape_training_cfm@data[,2])*100
    
    ind_pro_fun_mod<-approxfun(shape_training_cfm@data[,3],shape_training_cfm@data[,1])
    area_susc_values_cfm_mod<-c(0,ind_pro_fun_mod(susc_values_cfm[-c(1,length(susc_values_cfm))]),100)
    
    #dev.new()
    pdf(file = "result_CFM_SuccessRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="CFM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_cfm)-1))
    {
      #count=1
      polygon(c(area_susc_values_cfm_mod[count:(count+1)],rev(area_susc_values_cfm_mod[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_training_cfm@data[,1],shape_training_cfm@data[,2],col="black")
    dev.off()
    
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_cfm_mod_optimal<-c(0,ind_pro_fun_mod(rev(cfm.breaks.histogram.values.optimal)[-c(1,length(cfm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_CFM_SuccessRateCurve_Optimal.pdf",onefile = TRUE, pagecentre=TRUE)
      plot(0,0,col="transparent",main="CFM SUCCESS RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(cfm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_cfm_mod_optimal[count:(count+1)],rev(area_susc_values_cfm_mod_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_training_cfm@data[,1],shape_training_cfm@data[,2],col="black")
      dev.off()
    }
    ############
    
    ind_col_succ_pred_rate<-which(colnames(shape_validation_cfm@data) %in% c("VAL_GROUP","VAL_PROB"))
    shape_validation_cfm@data<-cbind(PIXEL_AREA=rep(as.numeric(configuration.spatial.data.table[c(8)])^2,dim(shape_validation_cfm@data)[1]),shape_validation_cfm@data[order(shape_validation_cfm@data[,ind_col_succ_pred_rate[2]],decreasing=TRUE),ind_col_succ_pred_rate])
    shape_validation_cfm@data$VAL_GROUP<-shape_validation_cfm@data$VAL_GROUP*as.numeric(configuration.spatial.data.table[c(8)])^2
    shape_validation_cfm@data[,1]<-cumsum(shape_validation_cfm@data[,1])/sum(shape_validation_cfm@data[,1])*100
    shape_validation_cfm@data[,2]<-cumsum(shape_validation_cfm@data[,2])/sum(shape_validation_cfm@data[,2])*100
    
    ind_pro_fun_val<-approxfun(shape_validation_cfm@data[,3],shape_validation_cfm@data[,1])
    area_susc_values_cfm_val<-c(0,ind_pro_fun_val(susc_values_cfm[-c(1,length(susc_values_cfm))]),100)
    
    #dev.new()
    pdf(file = "result_CFM_PredictionRateCurve.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
    plot(0,0,col="transparent",main="CFM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
    for (count in 1:(length(susc_values_cfm)-1))
    {
      #count=1
      polygon(c(area_susc_values_cfm_val[count:(count+1)],rev(area_susc_values_cfm_val[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color.vector.susceptibility)[count])  
    }
    polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
    lines(shape_validation_cfm@data[,1],shape_validation_cfm@data[,2],col="black")
    dev.off()
    ###########
    if(enable_probability_optimal_classification==TRUE) 
    {
      area_susc_values_cfm_val_optimal<-c(0,ind_pro_fun_val(rev(cfm.breaks.histogram.values.optimal)[-c(1,length(cfm.breaks.histogram.values.optimal))]),100)
      pdf(file = "result_CFM_PredictionRateCurve_Optimal.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
      plot(0,0,col="transparent",main="CFM PREDICTION RATE CURVE",xlab="Cumulative percentage study area ",ylab="Cumulative percentage landslide area",xlim=c(0,100),ylim=c(0,100))
      for (count in 1:(length(rev(cfm.breaks.histogram.values.optimal))-1))
      {
        #count=1
        polygon(c(area_susc_values_cfm_val_optimal[count:(count+1)],rev(area_susc_values_cfm_val_optimal[count:(count+1)])),c(0,0,100,100),border="darkgray",lty="dotted",lwd=0.5,col=rev(color_ramp_palette_fun(length(cfm.breaks.histogram.values.optimal)-1))[count])  
      }
      polygon(c(0,100,100,0),c(0,0,100,100),border="black",lty="solid",lwd=1,col=NULL)
      lines(shape_validation_cfm@data[,1],shape_validation_cfm@data[,2],col="black")
      dev.off()
    }
    ############
    
    
  }
}
  
                                                                      #
#------------------- MODEL PROBABILITY COMPARISON --------------------#


if(model.run.matrix[1] == "YES" & model.run.matrix[2] == "YES")
  {
  # LDA - QDA
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.lda$posterior[,2],predict.result.qda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LDA Model Probability",ylab="QDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_LDA_QDA.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.lda$posterior[,2],predict.result.qda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LDA Model Probability",ylab="QDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }
  
if(model.run.matrix[1] == "YES" & model.run.matrix[3] == "YES")
  {  
  # LDA - LRM
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.lda$posterior[,2],predict(result.lrm, type="response"),type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LDA Model Probability",ylab="LRM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_LDA_LRM.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.lda$posterior[,2],predict(result.lrm, type="response"),type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LDA Model Probability",ylab="LRM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[1] == "YES" & model.run.matrix[4] == "YES")
  {
  # LDA - NNM
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.lda$posterior[,2],predict.result.nnm,type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LDA Model Probability",ylab="NNM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_LDA_NNM.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.lda$posterior[,2],predict.result.nnm,type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LDA Model Probability",ylab="NNM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[2] == "YES" & model.run.matrix[1] == "YES")
  {
  # QDA - LDA
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.qda$posterior[,2],predict.result.lda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="QDA Model Probability",ylab="LDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_QDA_LDA.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.qda$posterior[,2],predict.result.lda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="QDA Model Probability",ylab="LDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[2] == "YES" & model.run.matrix[3] == "YES")
  {
  # QDA - LRM
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.qda$posterior[,2],predict(result.lrm, type="response"),type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="QDA Model Probability",ylab="LRM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_QDA_LRM.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.qda$posterior[,2],predict(result.lrm, type="response"),type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="QDA Model Probability",ylab="LRM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[2] == "YES" & model.run.matrix[4] == "YES")
  {  
  # QDA - NNM
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.qda$posterior[,2],predict.result.nnm,type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="QDA Model Probability",ylab="NNM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_QDA_NNM.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.qda$posterior[,2],predict.result.nnm,type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="QDA Model Probability",ylab="NNM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[3] == "YES" & model.run.matrix[1] == "YES")
  {
  # LRM - LDA
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict(result.lrm, type="response"),predict.result.lda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LRM Model Probability",ylab="LDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_LRM_LDA.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict(result.lrm, type="response"),predict.result.lda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LRM Model Probability",ylab="LDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[3] == "YES" & model.run.matrix[2] == "YES")
  {
  # LRM - QDA
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict(result.lrm, type="response"),predict.result.qda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LRM Model Probability",ylab="QDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_LRM_QDA.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict(result.lrm, type="response"),predict.result.qda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LRM Model Probability",ylab="QDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[3] == "YES" & model.run.matrix[4] == "YES")
  {
  # LRM - NNM
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict(result.lrm, type="response"),predict.result.nnm,type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LRM Model Probability",ylab="NNM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_LRM_NNM.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict(result.lrm, type="response"),predict.result.nnm,type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="LRM Model Probability",ylab="NNM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[4] == "YES" & model.run.matrix[1] == "YES")
  {
  # NNM - LDA
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.nnm,predict.result.lda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="NNM Model Probability",ylab="LDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_NNM_LDA.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.nnm,predict.result.lda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="NNM Model Probability",ylab="LDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[4] == "YES" & model.run.matrix[2] == "YES")
  {
  # NNM - QDA
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.nnm,predict.result.qda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="NNM Model Probability",ylab="QDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_NNM_QDA.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.nnm,predict.result.qda$posterior[,2],type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="NNM Model Probability",ylab="QDA Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

if(model.run.matrix[4] == "YES" & model.run.matrix[3] == "YES")
  {
  # NNM - LRM
  if (enable_screen_plotting==TRUE)
  {
  dev.new()
  plot(predict.result.nnm,predict(result.lrm, type="response"),type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="NNM Model Probability",ylab="LRM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  }
  pdf(file = "result_ModelComparison_NNM_LRM.pdf", width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
  plot(predict.result.nnm,predict(result.lrm, type="response"),type="p",pch=1,cex=0.85,col="dark blue",xlim=c(0,1),ylim=c(0,1),xlab="NNM Model Probability",ylab="LRM Model Probability", main="MODEL COMPARISON")
  abline(a=0,b=1,col="red",lty=1,lwd=1)
  dev.off()
  }

save.image(file=paste("Susceptibility_",unlist(strsplit(rdata_file,"/"))[length(unlist(strsplit(rdata_file,"/")))],sep=""))

time_end_calculation<-Sys.time()
total_calculation_time<-difftime(time_end_calculation, time_start_calculation,units="hours")
print(paste("Total calculation time: ",total_calculation_time," hours",sep=""))



