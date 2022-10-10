# Background analysis
# Background tone interval length

BTIL.f <- function(im, Informacion, NumMask){
  im <- im
  Ejey <- dim(im)[2]
  Informacion <- Informacion
  NumMask <- NumMask
  
  Dim <- dim(im)
  
  #Cambiar a tonalidades de gris
  colorMode(im) = Grayscale
  
  #Cortar la parte de la zimografía que es de interés
  im_cut = im[1:Dim[1], 1:Ejey,1]
  
  #Matriz de tonos de grises
  M <- imageData(im)[,,1]
  
  #Encontrar los pixeles en las máscaras
  mask1 <- t(Informacion$mascaras[1:Ejey, 1:Dim[1], 1])
  for(w in 2:NumMask){
    mask <- t(Informacion$mascaras[1:Ejey, 1:Dim[1], w])
    mask1 <- mask1 + mask
  }
  #Quitar los repetidos
  if(max(mask1)==2){
    val2 <- which(mask1 == 2, arr.ind=TRUE)
    mask1[val2]<-0
  }
  #Pixeles del fondo
  p_fondo <- which(mask1 == 0, arr.ind=TRUE)
  
  #Pixeles de las máscaras
  p_mask <- which(mask1 == 1, arr.ind=TRUE)
  
  #Tonos de gris de los pixeles de fondo
  dimpf <- dim(p_fondo)[1]
  Mpos <- NULL
  for(u in 1:dimpf){
    Mpos[u] <- M[p_fondo[u,1], p_fondo[u,2]]
  }
  
  #Juntar la posición con su correspondiente tono en M	
  FONDO <- cbind(p_fondo, Mpos)
  
  #Para las esquinas de los pixeles de las máscaras
  mM <- c(min(p_mask[,1]), max(p_mask[,1]), min(p_mask[,2]), max(p_mask[,2]))
  
  #Discriminación de fondo
  fondo2 <- NULL
  fondo3 <- NULL
  indices <- which(FONDO[,1]>=mM[1] & FONDO[,1] <= mM[2])
  fondo2 <- FONDO[indices,]
  indices <- which(fondo2[,2]>=mM[3] & fondo2[,2]<=mM[4])
  fondo3 <- fondo2[indices,]

  intervalo <- max(fondo3[,3])-min(fondo3[,3])
  return(intervalo)
}