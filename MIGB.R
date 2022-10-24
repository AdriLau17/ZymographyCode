# MIGB 
# GI algorithm for a 2-dimensional data set
#Inputs:
#	  s - Number of adjusted gaussians
#	  datos - dataframe in dimension 2
#	  mui - initial values matrix for the means |x|y| (s rows)
#	  sigmai - initial values matrix for the standard deviation (idem mu)
#	  pri - initial values vector for the mixture proportions
#Outputs:
#   R - vector with the calculated parameters  (s vectors for the means, s vectors for the covariance matrices)
#   pr - vector with s values for the mixture proportions

MIGB <- function(s, datos, mui, sigmai, pri){
  s <- s
  datos <- datos
  mu <- mui
  sigma <- sigmai
  pr <- pri
  
  if(dim(mu)[1]==s && dim(sigma)[1]==s && length(pr)==s){
    muc <- mu + 1
    sigmac <- sigma + 1
    prc <- pr
    cont <- 1
    
    val <- length(pr)
    valc <- length(pr)
    val <- c(as.vector(mu), as.vector(sigma))
    valc <- c(as.vector(muc), as.vector(sigmac))
    
    while(as.numeric(dist(rbind(val, valc)))>1e-06 && sum(sigma!=0)==(dim(sigma)[1]*dim(sigma)[2]) && cont!=100){
      muc <- mu
      sigmac <- sigma
      prc <- pr
      valc <- c(as.vector(muc), as.vector(sigmac))
      
      #Process in each column
      G<-NULL
      for(w in 1:2){
        dat <- NULL
        dat <- datos[,w]
        RJ <- matrix(ncol=s, nrow=length(dat))
        for(f in 1:s){	
          RJ[,f] <- exp((-(dat-mu[f,w])^2)/(2*(sigma[f,w]^2)))
        }
        
        # New means and standard deviations
        for(f in 1:s){
          mu[f,w] <- sum(dat*RJ[,f]) / sum(RJ[,f])
          sigma[f,w] <- sqrt(sum(((dat-mu[f,w])^2)*RJ[,f]) / sum(RJ[,f]))
        }
        
        # Mixture proportions
        g <- matrix(ncol=s, nrow=length(dat))
        for(f in 1:s){
          g[,f] <- dnorm(dat,mu[f,w],sigma[f,w])
        }
        G<-cbind(G,g)	
      }
      
      # P(x|Gi) probabilities
      PG <- NULL
      RG <- NULL
      for(r in 1:s){
        PG<-G[,r]*G[,r+s]
        RG<-cbind(RG,PG)
      }
      
      # Numerator
      num <- matrix(ncol=s, nrow=length(dat))
      for(i in 1:s){
        num[,i] <-RG[,i]*pr[i]	
      }
      
      # Denominator
      den <- NULL
      for(i in 1:length(dat)){
        den[i] <- sum(num[i,])
      }
      
      # P(Gi|x)
      prx <- matrix(ncol=s, nrow=length(dat))
      for(i in 1:length(dat)){
        prx[i,]<-num[i,]/den[i]
      }
      prx[is.nan(prx)]<-0
      
      #Class of each element
      cl <- NULL
      for(i in 1:length(dat)){
        cl[i]<-which.max(as.vector(prx[i,]))
      }
      
      # New pr value
      for(i in 1:s){
        pr[i]<-sum(as.numeric(cl==i))/length(dat)
      }		
      cont <- cont + 1 
      valc <- c(as.vector(mu), as.vector(sigma))
    }
    lista<-list("R"=valc, "pr"=pr)
    return(lista)
    print(mu)
    print(pr)
  }else{print("Different length in initial parameters")
  }
}
