# Functions to obtain the Umbral value "Um" and the B value
# im - Image
# Informacion - Information about the image masks
# NumMask - Number of masks 

################################################################################
# Um value
P.Umbral <- function(im, Informacion, NumMask){
  im <- im
  Informacion <- Informacion 
  NumMask <- NumMask
  
  mask <- t(Informacion$mascaras[,, 1])
  for(w in 2:NumMask){
    mask <- mask + t(Informacion$mascaras[,, w])
  }
  
  #Seleccionar las posiciones dentro de las máscaras
  pos.pixels <- which(mask==1, arr.ind=TRUE)
  
  #Matriz de tonos de grises
  M <- imageData(im)[,,1]
  #M <- imageData(im)
  
  dimp <- dim(pos.pixels)[1]
  Mpos <- NULL
  for(u in 1:dimp){
    Mpos[u] <- M[pos.pixels[u,1], pos.pixels[u,2]]
  }
  
  #División por tonos con posible umbral
  l_espacio <- (max(Mpos)-min(Mpos))/5
  Umbral <- min(Mpos) + l_espacio*4
  Umbral <- round(Umbral,2)
  return(Umbral)
}

##########################################################################
# B value
P.B <-function(im, Informacion, NamMask){
  im <- im 
  Informacion <- Informacion
  NumMask <- NumMask
  
  #Dimensiones de la máscara
  resul <- dim.mask(im, Informacion, NumMask)
  lap <- resul$lap
  Esquinas <- resul$ESQUINAS
  
  #Promedio de dimensiones 
  p.lap <- c(mean(lap[,1]), mean(lap[,2]), mean(lap[,3]))
  
  #Diferencia de dimensiones con respecto al promedio 
  dif.lap <- sweep(lap,2,p.lap)
  
  #considerar las máscaras más grandes que el valor promedio
  indnoneg <- which((dif.lap[,3] > 0)==TRUE)
  noneg2 <- lap[indnoneg, 3]
  
  #Promedio de éstas
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
