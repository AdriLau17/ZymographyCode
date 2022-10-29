# Main code

library(EBImage)
library("R.matlab")
library(mnormt)    
library(mixtools)

############################### Image processing ###############################

# Load the image (Gel.jpg file) and the masks information (Gel.mat file) 
# In the main folder of this repository there are two sample files: Gel.jpg and Gel.mat.
# IMPORTANT: To upload the information with the following instructions, the files 
# must be saved in your documents folder or you must specify the working 
# directory where the files are located, for example:
# readImage("C:/User/Data/Images/Gel.jpg")
# Note that you must use forward slash "/" 

im <- readImage("Gel.jpg")
plot(im)

Informacion <- readMat("Gel.mat")
NumMask <- length(Informacion$puntajes)

colorMode(im) = Grayscale

############################## Background analysis #############################
BTIL <- BTIL.f(im, Informacion, NumMask)
BTIL

######################## Automatization of Um parameter ########################
Um <- P.Umbral(im, Informacion, NumMask)
Um

######################## Automatization of B parameter #########################
B <- P.B(im, Informacion, NamMask)
B

################################# Methodology ##################################
Results <- Methodology(im, Informacion, NumMask, B, Um)
Results
