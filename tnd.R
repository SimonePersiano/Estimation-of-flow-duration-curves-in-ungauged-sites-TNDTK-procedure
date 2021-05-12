#
# Do remember to upload "pracma" package
# library(pracma)
#
# IMPORTANT: fdc has to be already standardized (divided by MAF value)
#
fdc.tnd<-function(fdc,norm,maxd,logq){
  zeros<-which(fdc<0.00003)
  if (length(zeros)!=0){
    fdc<-fdc[-zeros]
    d<-fdc[-zeros]
  }
  l<-length(fdc)
  d<-1:l/(l+1) ##### weibull
  #
  d1<-d[which(fdc<=1)[1]]  
  #
  if (is.null(logq)){ ##default
    logq=FALSE
  }
  #
  if (logq==FALSE){
    if (norm==TRUE){
      f<-approxfun(qnorm(d),fdc,rule=2)  
      tnd<-(qnorm(maxd)-qnorm(d1))-quad(f,xa=qnorm(d1),xb=qnorm(maxd))
    }
    
    else { 
      f<-approxfun(d,fdc,rule=2)
      tnd<-(maxd-d1)-quad(f,xa=d1,xb=maxd)
    }
  }
  
  if (logq==TRUE){
    
    if (norm==TRUE){
      f<-approxfun(qnorm(d),log(fdc),rule=2)  
      tnd <- -quad(f,xa=qnorm(d1),xb=qnorm(maxd))
    }
    
    else { 
      f<-approxfun(d,log(fdc),rule=2)
      tnd<-(maxd-d1)-quad(f,xa=d1,xb=maxd)
    }
    
  }
  return(tnd)
}
