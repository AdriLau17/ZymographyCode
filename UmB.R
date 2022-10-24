# Functions to obtain the threshold value "Um" and the "B" value
#Inputs:
#   im - Image
#   Informacion - Information about the image masks
#   NumMask - Number of masks 


################################################################################
# Function to obtain the threshold value "Um"
#Output: Um value

P.Umbral <- function(im, Informacion, NumMask){
  im <- im
  Informacion <- Informacion 
  NumMask <- NumMask
  
  mask <- t(Informacion$mascaras[,, 1])
  for(w in 2:NumMask){
    mask <- mask + t(Informacion$mascaras[,, w])
  }
  
  # Select the positions within the masks
  pos.pixels <- which(mask==1, arr.ind=TRUE)
  
  # Grayscale tone array
  M <- imageData(im)[,,1]
  #M <- imageData(im)
  
  dimp <- dim(pos.pixels)[1]
  Mpos <- NULL
  for(u in 1:dimp){
    Mpos[u] <- M[pos.pixels[u,1], pos.pixels[u,2]]
  }
  
  # Split by tones (threshold)
  l_espacio <- (max(Mpos)-min(Mpos))/5
  Umbral <- min(Mpos) + l_espacio*4
  Umbral <- round(Umbral,2)
  return(Umbral)
}

##########################################################################
# Funtion to obtain the B value
#Output: B Value

P.B <-function(im, Informacion, NamMask){
  im <- im 
  Informacion <- Informacion
  NumMask <- NumMask
  
  # Dimensions of the masks
  resul <- dim.mask(im, Informacion, NumMask)
  lap <- resul$lap
  Esquinas <- resul$ESQUINAS
  
  # Avarages 
  p.lap <- c(mean(lap[,1]), mean(lap[,2]), mean(lap[,3]))
  
  # Difference of dimensions with respect to the average 
  dif.lap <- sweep(lap,2,p.lap)
  
  # Bigger masks than the average
  indnoneg <- which((dif.lap[,3] > 0)==TRUE)
  noneg2 <- lap[indnoneg, 3]
  
  # Avarage of the biggest masks
  mnoneg2 <- mean(noneg2)
  
  dnoneg2 <- NULL
  for(w in 1:length(noneg2)){
    dnoneg2[w] <- abs(noneg2[w]-mnoneg2)
  }
  posB2 <- which.min(dnoneg2)
  
  pos2 <- indnoneg[posB2]
  
  limB2 <- lap[pos2,3]
  B <- sum((lap[,3]>=limB2)==TRUE)
  return(B)
}
