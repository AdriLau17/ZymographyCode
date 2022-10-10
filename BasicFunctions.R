# BASIC FUNCTIONS

# General function (Process) 
# In:
# im - image
#	Informacion - matlab information
#	NumMask - Number of masks
#	B - Number of bigger masks
#	Um - Umbral

# Out:    TO - Table 
#			| Proportions | Width | Lenght | Pixels | Class
#		Graphics (considered pixels and classes of the masks)

Methodology <- function(im, Informacion, NumMask, B, Um){
  im <- im
  Informacion <- Informacion 
  NumMask <- NumMask
  B <- B
  Um <- Um
  
  #Para las etiquetas de identificaci�n
  Letters <- c(LETTERS,letters)[1:NumMask]
  
  #Para calcular los puntos centrales de las m�scaras
  medias<-matrix(ncol=2, nrow=NumMask)
  for(w in 1:NumMask){
    mask <- t(Informacion$mascaras[,, w])
    xy <- which(mask==1, arr.ind = TRUE)
    mx <- (max(xy[,1])+min(xy[,1]))/2
    my <- (max(xy[,2])+min(xy[,2]))/2
    medias[w,]<-c(mx,my)
  }
  Perfil.mascaras(im, Informacion, NumMask)	
  text(x=medias[,1],y=medias[,2], labels=Letters, col="blue", cex=0.85)
  print("perfiles")
  #x11()
  #Matriz lp de dimensiones de cada m�scara  | largo | ancho | num. pixeles |
  LAP <- dim.mask(im, Informacion, NumMask)
  lp <- LAP$lap
  
  if(B!=0){
    #Identificar las B m�scaras m�s grandes
    OrdenPMask <- sort(lp[,3], decreasing = TRUE)
    Bmask <- OrdenPMask[1:B]
    BMask <- NULL
    for(u in 1:length(Bmask)){
      BMask[u] <- which(lp[,3]==Bmask[u])
    }
    BMask2 <- BMask
    BMask <- sort(BMask)
    letras <- Letters[BMask2]
    Letters <- Letters[-BMask2]
    
    #Encontrar los pixeles en las m�scaras por debajo del umbral marcado
    #sin considerar las m�scaras m�s grandes
    pixe <- Pos.MaskUmbral(im, Um, NumMask, Informacion, BMask)
    print("Identification of the interest pixels set P complete.")
    #Ajuste de las Gaussianas sin considerar las B m�scaras m�s grandes
    Resultados <- Ajuste.Umbral(pixe, NumMask, Informacion, "Si", BMask)
    Resultados
    print("Gaussian fitting with the GI algorithm complete.")
    #Ordenar las m�scaras con respecto a los valores de pr
    IM <- Identificar.Mascara(Informacion, NumMask, Resultados, BMask)
    print("Ordering")
    #Obtener tabla con proporciones, ancho, largo, pixeles y etiqueta ordenados
    lp <- lp[-BMask,]
    TO <- Tabla.Ordenada(lp, IM, NumMask, BMask)
    letras2 <- Letters[IM$mascaras.or]		
    
    #Agregar a la tabla las m�scaras ya clasificadas
    if(B==1){
      to <- c(1,LAP$lap[BMask2,], as.numeric(unlist(Informacion$puntajes[BMask2])))
      to <- matrix(to, nrow=1)
      to <- as.data.frame(to)
    }else{
      to <- cbind(rep(1,B), LAP$lap[BMask2,], as.numeric(unlist(Informacion$puntajes[BMask2])))
      to <- as.data.frame(to)
    }
    names(to)<- c("Proporciones", "Ancho", "Largo", "Pixeles", "Clase")
    TO <- rbind(to, TO)
    
    #Normalizar las proporciones de la tabla TO
    TO[,1]<-TO$Proporciones / sum(TO$Proporciones)
    
    #Agregar etiquetas para identificaci�n
    Etiquetas <- c(letras, letras2)
    TO <- cbind(Etiquetas,TO)
    
    #Tabla resultado
    TO2 <- cbind(TO[1],TO[2]*100, TO[6])
    names(TO2)<- c("Spot", "Proportions", "Class")
    
  }else{
    BMask <- NULL
    
    #Encontrar los pixeles en las m�scaras por debajo del umbral marcado
    pixe <- Pos.MaskUmbral(im, Um, NumMask, Informacion, BMask)
    
    #Ajuste de las Gaussianas sin considerar las B m�scaras m�s grandes
    Resultados <- Ajuste.Umbral(pixe, NumMask, Informacion, "Si", BMask)
    Resultados
    
    #Ordenar las m�scaras con respecto a los valores de pr
    IM <- Identificar.Mascara(Informacion, NumMask, Resultados, BMask)
    
    #Obtener tabla con proporciones, ancho, largo, pixeles y etiqueta ordenados
    TO <- Tabla.Ordenada(lp, IM, NumMask, BMask)
    
    #Agregar etiquetas para identificaci�n
    Etiquetas <- Letters[IM$mascaras.or]
    TO <- cbind(Etiquetas,TO)
    
    #Tabla resultado
    TO2 <- cbind(TO[1],TO[2]*100, TO[6])
    names(TO2)<- c("Spot", "Proportions", "Class")
  }
  return(TO2)
}

#############################################################################
# Masks

Perfil.mascaras <- function(im, Informacion, NumMask){
  im <- im
  Informacion <- Informacion 
  NumMask <- NumMask
  
  dIm <- dim(im)
  
  #Matriz de largo, ancho y n�mero de pixeles
  lap <- matrix(nrow=NumMask, ncol=3)
  ESQUINAS <- matrix(nrow=NumMask, ncol=8)
  for(r in 1:NumMask){
    mask <- t(Informacion$mascaras[1:dIm[2], 1:dIm[1], r])
    xy <- which(mask==1, arr.ind = TRUE)
    
    ejex <- c(min(xy[,1]), max(xy[,1]))
    ejey <- c(min(xy[,2]), max(xy[,2]))
    
    esquinas <- matrix(nrow=4,ncol=2)
    esquinas[,1] <- c(ejex[1], ejex[1], ejex[2], ejex[2]) 
    esquinas[,2] <- c(ejey[1], ejey[2], ejey[1], ejey[2])
    
    ESQUINAS[r,] <- c(esquinas[1,], esquinas[2,], esquinas[3,], esquinas[4,])
  }
  
  #plot(im)
  for(w in 1:NumMask){
    segments(ESQUINAS[w,1], ESQUINAS[w,2],ESQUINAS[w,3],ESQUINAS[w,4])
    segments(ESQUINAS[w,5], ESQUINAS[w,6],ESQUINAS[w,7],ESQUINAS[w,8])
    segments(ESQUINAS[w,1], ESQUINAS[w,2],ESQUINAS[w,5],ESQUINAS[w,6])
    segments(ESQUINAS[w,3], ESQUINAS[w,4],ESQUINAS[w,7],ESQUINAS[w,8])
  }
}

##############################################################################
#Find the dimensions and the number of pixels of each mask
#In: 
#	im - image
#	Informacion - matlab information
#	NumMask - Number of masks
#Out: lap matrix with the dimensions of each mask
#    			|  Lenght  | Width | Number of Pixels |
#	  ESQUINAS matrix with the points that are corners of each mask
#		|x1|y1|x2|y2|x3|y3|x4|y4|

dim.mask <- function(im, Informacion, NumMask){
  im <- im
  Informacion <- Informacion 
  NumMask <- NumMask
  
  dIm <- dim(im)
  #plot(im)
  
  #Matriz de largo, ancho y n�mero de pixeles
  lap <- matrix(nrow=NumMask, ncol=3)
  ESQUINAS <- matrix(nrow=NumMask, ncol=8)
  for(r in 1:NumMask){
    mask <- t(Informacion$mascaras[1:dIm[2], 1:dIm[1], r])
    xy <- which(mask==1, arr.ind = TRUE)
    
    ejex <- c(min(xy[,1]), max(xy[,1]))
    ejey <- c(min(xy[,2]), max(xy[,2]))
    
    esquinas <- matrix(nrow=4,ncol=2)
    esquinas[,1] <- c(ejex[1], ejex[1], ejex[2], ejex[2]) 
    esquinas[,2] <- c(ejey[1], ejey[2], ejey[1], ejey[2])
    
    lap[r,] <- c(dist(ejex),  dist(ejey), dim(xy)[1])
    
    #points(esquinas[,1], esquinas[,2], col="blue", pch=20)
    
    ESQUINAS[r,] <- c(esquinas[1,], esquinas[2,], esquinas[3,], esquinas[4,])
  }
  lista <- list("lap" = lap, "ESQUINAS" = ESQUINAS)
  return(lista)
}

###############################################################################
# Pixels positions that are obtained after applying the umbral (Graphic)

Pos.MaskUmbral <- function(im, Um, NumMask, Informacion, BMask){
  im <- im 
  Um <- Um
  NumMask <- NumMask
  Informacion <- Informacion
  BMask <- BMask 
  
  #Matriz de tonos de grises
  M <- imageData(im)[,,1]
  
  #No seleccionar las posiciones de los pixeles de las m�scaras m�s grandes
  if(length(BMask)!=0){
    qmask <- matrix(0, ncol=dim(im)[2], nrow=dim(im)[1])
    for(w in 1:length(BMask)){
      qmask <- qmask + t(Informacion$mascaras[,,BMask[w]])
    }
    q.pixels <- which(qmask==1, arr.ind=TRUE)
    
    #Segunda opci�n
    x <- 1:NumMask
    x <- x[-BMask]
    
    mask <- matrix(0, ncol=dim(im)[2], nrow=dim(im)[1])
    for(w in x){
      mask <- mask + t(Informacion$mascaras[,, w])
    }
    
  }else{
    mask <- t(Informacion$mascaras[,, 1])
    for(w in 2:NumMask){
      mask <- mask + t(Informacion$mascaras[,, w])
    }
    
  }
  pos.pixels <- which(mask==1,arr.ind=TRUE)	
  
  dimp <- dim(pos.pixels)[1]
  Mpos <- NULL
  for(u in 1:dimp){
    Mpos[u] <- M[pos.pixels[u,1], pos.pixels[u,2]]
  }
  #Juntar la posici�n con su correspondiente tono en M	
  datos <- cbind(pos.pixels, Mpos)
  
  #Seleccionar los pixeles de los tonos menores al umbral
  pix <- which(Mpos<Um)
  datos2 <- datos[pix,]
  pix <- NULL
  pix <- datos2[,1:2]
  
  #Para calcular los puntos centrales de las m�scaras
  medias<-matrix(ncol=2, nrow=NumMask)
  for(w in 1:NumMask){
    mask <- t(Informacion$mascaras[,, w])
    xy <- which(mask==1, arr.ind = TRUE)
    mx <- (max(xy[,1])+min(xy[,1]))/2
    my <- (max(xy[,2])+min(xy[,2]))/2
    medias[w,]<-c(mx,my)
  }
  
  #Gr�ficas
  #par(mfrow=c(1,2)) 
  #plot(im)
  #labelsm <- unlist(Informacion$puntajes)
  #text(x=medias[,1],y=medias[,2], labels=labelsm, col="red", cex=1.5)
  
  #dIm <- dim(M)
  #if(dIm[1]>dIm[2]){
  #  x <- seq(1:dIm[1])
  #  y <- c(seq(1:dIm[2]),seq(1:(dIm[1]-dIm[2])))
  #  plot(x,-y, col="white",main=paste("Menor que ", Um),asp=1)
  #}else{if(dIm[1]==dIm[2]){
  #  x <- seq(1:dIm[1])
  #  plot(x,-x,col="white",main=paste("Menor que ", Um),asp=1)
  #}else{
  #  x <- c(seq(1:dIm[1]),seq(1:(dIm[2]-dIm[1])))
  #  y <- seq(1:dIm[2])
  #  plot(x,-y, col="white",main=paste("Menor que ", Um),asp=1)
  #}
  #}
  points(pix[,1], -pix[,2], col="black", pch=20)
  if(length(BMask)!=0){
    points(q.pixels[,1], -q.pixels[,2], col="gray", pch=20)
  }
  labelsm <- unlist(Informacion$puntajes)
  text(x=medias[,1],y=-medias[,2], labels=labelsm, col="red", cex=1.5)	
  return(pix)
}

###############################################################################
#Adjust NumMask gaussians in the data set pixe
#In: 	
#   pixe - posicion matrix
#		NumMask - Number of masks
#		Informacion - matlab information
#		Usomu - Si o No (Use of means)
#Out:	
#   mus - means matrix
#		sigmas - standard deviations matrix
#		pr - mixture proportions

Ajuste.Umbral <- function(pixe, NumMask, Informacion, Usomu, BMask){
  pixe <- pixe
  NumMask <- NumMask
  Informacion <- Informacion
  Usomu <- Usomu
  BMask <- BMask
  
  medias<-matrix(ncol=2, nrow=NumMask)
  nc <- length(BMask)	
  
  if(Usomu=="Si"){
    for(w in 1:NumMask){
      mask <- t(Informacion$mascaras[,, w])
      xy <- which(mask==1, arr.ind = TRUE)
      mx <- (max(xy[,1])+min(xy[,1]))/2
      my <- (max(xy[,2])+min(xy[,2]))/2
      medias[w,]<-c(mx,my)
    }
    if(length(BMask)!=0){
      medias <- medias[-BMask,]
    }
  }else{
    maximx <- dim(Informacion$mascara)[1]
    maximy <- 375	
    mx <- sample(1:maximx, (NumMask-nc), replace=FALSE)
    my <- sample(1:maximy, (NumMask-nc), replace=FALSE)
    medias <- cbind(mx,my)
  }
  
  s <- NumMask-nc
  mui <- medias
  x <- 4:8
  sig.in<-sample(x,(s*4), replace=TRUE)
  sigmai<-matrix(sig.in,s,2)
  pri <- rep(1/s,s)
  
  #print("Valores mu, sigma, pr")
  Salida <- MIGB(s, pixe, mui, sigmai, pri)
  R <- Salida$R
  pr <- Salida$pr
  print("Obtenci�n de resultados")
  mus <- R[1:(length(R)/2)]
  mus <- matrix(mus, s,2)
  sigmas <- R[(length(R)/2 + 1): length(R)]
  sigmas <- matrix(sigmas, s, 2)
  
  lista <- list("mus"=mus, "sigmas"=sigmas, "pr"=pr)
  return(lista)
}

################################################################################
#Identify masks with the ordered proportions
#In:
#	Informacion - matlab information
#	NumMask - Number of masks
#	Resultados - Vector of results

Identificar.Mascara <- function(Informacion, NumMask, Resultados, BMask){
  Informacion <- Informacion
  NumMask <- NumMask
  Resultados <- Resultados
  BMask <- BMask
  
  nc <- length(BMask)
  mus <- Resultados$mus
  sigmas <- Resultados$sigmas
  pr <- Resultados$pr
  
  pr.ordenado <- sort(pr, index.return = TRUE, decreasing=TRUE)
  pos <- pr.ordenado$ix
  
  pr.or <- pr.ordenado$x
  mus.or <- matrix(ncol=2, nrow=(NumMask-nc))
  for(w in 1:(NumMask-nc)){
    mus.or[w,] <- mus[pos[w],]
  }
  
  sigmas.or <- matrix(ncol=2, nrow=(NumMask-nc))
  for(w in 1:(NumMask-nc)){
    sigmas.or[w,] <- sigmas[pos[w],]
  }
  
  #Identificar m�scara con ayuda de las medias
  #Puntos centrales de las m�scaras
  medias<-matrix(ncol=2, nrow=NumMask)
  for(w in 1:NumMask){
    mask <- t(Informacion$mascaras[,, w])
    xy <- which(mask==1, arr.ind = TRUE)
    mx <- (max(xy[,1])+min(xy[,1]))/2
    my <- (max(xy[,2])+min(xy[,2]))/2
    medias[w,]<-c(mx,my)
  }
  if(length(BMask)!=0){
    medias <- medias[-BMask,]
  }
  MascaraRelacionada <- matrix(nrow=1, ncol=(NumMask-nc))
  medias2 <- medias
  for(u in 1:(NumMask-nc)){
    mu.con <- mus.or[u,]
    dmume <- NULL
    for(w in 1:(NumMask-nc)){
      dmume[w] <- dist(rbind(mu.con,medias2[w,]))
    }
    numas <- which(dmume<20)
    if(length(numas)==1){
      MascaraRelacionada[1,u] <- numas
      medias2[numas,] <- c(10000, 10000)
    }
    if(length(numas)==0){
      pos.min <- which.min(dmume)
      MascaraRelacionada[1,u] <- pos.min
      medias2[pos.min,]<-c(10000,10000)
    }
    if(length(numas)>0){
      p.re <- which.min(dmume[numas])
      MascaraRelacionada[1,u] <- numas[p.re]
      medias2[numas[p.re],] <- c(10000, 10000)
    }
  }
  
  etiquetas <- unlist(Informacion$puntajes)
  if(length(BMask)!=0){
    etiquetas <- etiquetas[-BMask]
  }
  etiquetas.or <- NULL
  for(w in 1:(NumMask-nc)){
    etiquetas.or[w] <- etiquetas[MascaraRelacionada[1,w]]
  }
  
  MascaraRelacionada <- as.vector(MascaraRelacionada)
  
  lista <- list("pr.or"=pr.or, "etiquetas.or"=etiquetas.or, "mascaras.or"=MascaraRelacionada, "mus.or"=mus.or, "sigmas.or"=sigmas.or)
  return(lista)
}

################################################################################

# Table with ordered proportions, wide, length, pixels and classes
Tabla.Ordenada <- function(lp, IM, NumMask, BMask){
  lp <- lp
  IM <- IM
  NumMask <- NumMask
  BMask <- BMask
  nc <- length(BMask)
  
  pr.or <- as.matrix(IM$pr.or, ncol=1)
  etiquetas.or <- as.matrix(as.numeric(IM$etiquetas.or),ncol=1)
  orden <- IM$mascaras.or
  lp2 <- matrix(ncol=3, nrow=(NumMask-nc))
  for(w in 1:(NumMask-nc)){
    lp2[w,] <- lp[orden[w],]
  }
  TABLA <- cbind(pr.or, lp2, etiquetas.or)
  TABLA <- as.data.frame(TABLA)
  names(TABLA)<- c("Proporciones", "Ancho", "Largo", "Pixeles", "Clase")
  return(TABLA)
}