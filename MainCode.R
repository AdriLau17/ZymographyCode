# Main code

library(EBImage)
library("R.matlab")
library(mnormt)    
library(mixtools)

############################### Image processing ###############################

# Load the image (.jpg file) and the masks information (.mat file)
# IMPORTANT: To upload the information with the following instructions, the files 
# must be saved in your documents folder or you must specify the working 
# directory where the files are located, for example:
# readImage("C:/User/Data/Images/Gel 2 (recortado) 24-06-19.jpg")
# Note that you must use forward slash "/" 

im <- readImage("Gel 2 (recortado) 24-06-19.jpg")
plot(im)

Informacion <- readMat("Gel 2 (recortado) 24-06-19.mat")
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