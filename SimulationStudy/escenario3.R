library(hdrcde)
library(nor1mix)
library(multimode)
library(np)
library(lpme)
library(nor1mix)
library(pracma)
library(multimode)
# función de estimación de la reg modal
source("moderegwe.R")
# estimación de la ventana
source("moderegbw_we.R")
# rutinas de C modificadas necesarias para la estimación y la ventana
source("LCFITmodereg1.R")
source("CVmode_LCfit1.R")
source("binormal.dep.x.R")
conditional_density=function(x0,y0,X,Y,w,bwx,bwy){
  n=length(X)
  a=dnorm(x0-X, sd=bwx)
  b=dnorm(y0-Y, sd=bwy)
  pyx=(1/bwy)*(sum(a*b*w)/sum(a*w))
}
extraerimpares=function(x){
  n1=length(x)
  x1=vector()
  for(i in 1:n1){
    if(i%%2==1)
      x1=c(x1,x[i])
  }
  return(x1)
}
pobs<-function(x,tipo.mis){
  if(tipo.mis==1){pobs<-0.6+0.3*cos(pi*x)}
  if(tipo.mis==2){pobs<-0.6+0.3*cos(2*pi*x)}
  if(tipo.mis==3){pobs<-0.7+0.3*cos(2*pi*x^2)} 
  if(tipo.mis==4){pobs<-0.75} 
  pobs
}
simulacion<-function(nr,nsamp,tipo.mis,a1,a2,b1,b2,p1){
  ARCHIVOSALIDA<-paste("BI","_M",tipo.mis,"b1b2",b1,b2,p1,"p.txt",sep="")  
  # p1 está al final porque no me cogía el nombre tan largo del archivo"
  n=nsamp
  ve<-c(n,tipo.mis, a1,a2,b1,b2,p1)
  write(t(ve),file=ARCHIVOSALIDA,ncolumns=length(ve),append=T)
  # selección de la ventana
  #vent="boostrap"
  vent="CV-mode"
  #funciones
  m1=function(s){
    2*sin(2*pi*s)
  }
  y=rep(0,n)
  eps=rep(0,n)
  ## COMIENZO E LA SIMULACION
  for (l in 1:nr){
    set.seed(999+l)
    print(l)
    x<-sort(runif(n))
     eps <- error.dep.x(x,mu=c(a1,a2), sig=c(b1,b2),w=c(p1,1-p1))
    y<-m1(x)+eps
    # valor real de la moda conocida para este modelo 
    xgrid = seq(0, 1, l = 200)
    real1=cbind(xgrid,m1(xgrid)+a1)
    real2=cbind(xgrid,m1(xgrid)+a2)
    real=rbind(real1,real2)
    #####missing data
    u<-runif(n,min=0,max=1)
    deltay=rep(1,n)
    py=rep(1,n)
    for (i in 1: n){
      py[i]<-pobs(x[i],tipo.mis)
      if (u[i]>py[i]){deltay[i]=0}
    }
    #  plot(x,y)
    #  plot(x,py)
    mean(deltay)
    ycomp=y[deltay==1]
    xcomp=x[deltay==1]
    ymis=y[deltay==0]
    xmis=x[deltay==0]
    pycomp=py[deltay==1]
    # de la comparativa de missing
    # plot(y~x, main="TODA LA MUESTRA")
    # plot(ycomp~xcomp,main="MUESTRA OBSERVADA")
    # points(xmis,ymis,col="red")
    ## Estimación con datos COMPLETOS
    # Estimación del parámetro por validación cruzada o bootstrap  Zhou y Huang(2019)
    hhxy = moderegbw(y, x, method=vent, p.order=0,h1=NULL, h2=NULL)$bw
    # Estimación de la regresión modal. Suponemos nstar=4  
    fit_c=modereg(y, W=x, xgrid=xgrid, bw=hhxy, nstart=4,sig=NULL,
                   p.order=0,PLOT=FALSE)
    ## Estimación SIMPLIFICADA
    # Estimación del parámetro por validación cruzada o bootstrap  Zhou y Huang(2019)
    hhxy_s = moderegbw.w(Y=ycomp, X=xcomp, method=vent, p.order=0,h1=NULL,
                         h2=NULL,weight=rep(1,length(ycomp)))$bw
    # Estimación de la regresión modal. Suponemos nstar=4  
    fit_s=modereg.w(Y=ycomp, W=xcomp, xgrid=xgrid, bw=hhxy_s, nstart=4,
                    sig=NULL, p.order=0, weight=rep(1,length(ycomp)),PLOT=FALSE)
    ## Estimación Ponderada
    # Estimación del parámetro por validación cruzada o bootstrap  Zhou y Huang(2019)
    hhxy = moderegbw.w(Y=ycomp, X=xcomp, method=vent, p.order=0,h1=NULL, 
                       h2=NULL,weight=1/pycomp)$bw
    # Estimación de la regresión modal. Suponemos nstar=4  
    fit_w=modereg.w(Y=ycomp, W=xcomp, xgrid=xgrid, bw=hhxy, nstart=4,sig=NULL,
                    p.order=0, weight=1/pycomp,PLOT=FALSE)
    ## Estimación Imputada
    fit_prev=modereg.w(ycomp, W=xcomp, xgrid=xmis, bw=hhxy_s, nstart=4,
                       sig=NULL, p.order=0, weight=rep(1,length(ycomp)),PLOT=FALSE)
    ### IMPUTACIÓN
    ###########################
    # Estimación de la regresión modal. Suponemos nstar=4 
    # Una forma de imputación podría ser imputar por un solo valor, el más probable.
    # Calculamos la densidad condicional y en el candidato a moda donde alcance el máximo.
    xhat=xcomp
    yhat=ycomp
### repite a función ata que non da erro!!!
    n.try <- 1
    bw <- try(cde.bandwidths(x=xcomp,y=ycomp,deg=0,method=1),silent=TRUE)
    while ('try-error' %in% class(bw) & n.try < 50) {
          cat("\t número de intentos ",n.try,"\n")
         	n.try <- n.try +1
     	    bw <- try(cde.bandwidths(x=xcomp,y=ycomp,deg=0,method=1),silent=TRUE)
          }
    for(i in 1:length(xmis)){
      #primero calculo las modas
      ind1=4*i-3
      ind2=4*i
      modas=fit_prev$result$fitted[c(ind1:ind2),2]
      x0=fit_prev$result$fitted[c(ind1:ind2),1]
      pyx=rep(0,length(x0))
      for (i1 in 1:length(modas)){
        pyx[i1]=conditional_density(x0[i1],y0=modas[i1],X=xcomp,Y=ycomp,
                                    w=rep(1,length(xcomp)),bwx=bw$a,bwy=bw$b) 
        }
      moda=modas[which(pyx==max(pyx))]
      xhat=c(xhat,xmis[i])
      yhat=c(yhat,moda)
    }
    #ordeno las x
    dat1=cbind(xhat,yhat)
    dat=dat1[order(xhat),]
    xhat=dat[,1]
    yhat=dat[,2]
    n1=length(xhat)
    hhxy = moderegbw.w(yhat, xhat, method=vent, p.order=0,h1=NULL, 
                       h2=NULL,weight=rep(1,n1),)$bw
    fit_imput1=modereg.w(yhat, W=xhat, xgrid=xgrid, bw=hhxy, nstart=4,
                         sig=NULL, p.order=0, weight=rep(1,n1),PLOT=FALSE)
    # IMPUTACIÓN MÚLTIPLE
    # Supongamos nstar=4
    fit_imput2=list()
    for(nimput in 1:3){
      xhat=xcomp
      yhat=ycomp
      for(i in 1:length(xmis)){
        #primero calculo las modas
        ind1=4*i-4
        ind2=4*i
        modas=fit_prev$result$fitted[c(ind1:ind2),2]
        x0=fit_prev$result$fitted[c(ind1:ind2),1]
        pyx=rep(0,length(x0))
        for (i1 in 1:length(modas)){
          pyx[i1]=conditional_density(x0[i1],y0=modas[i1],X=xcomp,Y=ycomp,
                                      w=rep(1,length(xcomp)),bwx=bw$a,bwy=bw$b)
          }
        pro=pyx/sum(pyx)
        moda=sample(modas,size=1,prob=pro)
        xhat=c(xhat,xmis[i])
        yhat=c(yhat,moda)
      }
      #ordeno las x
      dat1=cbind(xhat,yhat)
      dat=dat1[order(xhat),]
      xhat=dat[,1]
      yhat=dat[,2]
      n1=length(xhat)
      hhxy_i = moderegbw.w(yhat, xhat, method=vent, p.order=0,h1=NULL, 
                           h2=NULL,weight=rep(1,n1),)$bw
      fit_imput2[[nimput]]=modereg.w(yhat, W=xhat, xgrid=xgrid, bw=hhxy_i, 
                      nstart=4,sig=NULL, p.order=0, weight=rep(1,n1),PLOT=FALSE)
    } 
    a_imput2=list()
    for(i in 1:length(xgrid)){
      #primero calculo las modas
      ind1=4*i-3
      ind2=4*i
      modas=c(fit_imput2[[1]]$result$fitted[c(ind1:ind2),2],
              fit_imput2[[2]]$result$fitted[c(ind1:ind2),2],
              fit_imput2[[3]]$result$fitted[c(ind1:ind2),2])
      x0=fit_imput2[[1]]$result$fitted[ind1,1]
      # aquí uso a función para estimar o número de modas e logo a que estima as modas 
      #do paquete de Ameijeiras e rosa multimode
      nmodas=nmodes(modas,bw=hhxy[2])
      res=locmodes(modas,mod0=nmodas,display=FALSE)
      mod=extraerimpares(res$locations)
      r=expand.grid(x0,mod)
      a_imput2=rbind(a_imput2,r)
    }
    ## cálculo del error  MEDIDA DE HAUSDORFF 
    ## medida hausdorf del ERROR
    a_c=fit_c$result$fitted
    a_s=fit_s$result$fitted
    a_p=fit_w$result$fitted
    a_imput1=fit_imput1$result$fitted
    # a_imput2=calculado antes
    dis=rep(0, length(xgrid));
    for (j in 1:length(xgrid)){
      ymhat=subset( a_c, a_c[,1]==xgrid[j])
      ym=subset( real, real[,1]==xgrid[j])
      dis[j] = hausdorff_dist(ymhat,ym)^2
    } 
    ISE_c=sum( dis )/length(dis)*abs(xgrid[2]-xgrid[1]);
    dis=rep(0, length(xgrid));
    for (j in 1:length(xgrid)){
      ymhat=subset( a_s, a_s[,1]==xgrid[j])
      ym=subset( real, real[,1]==xgrid[j])
      dis[j] = hausdorff_dist(ymhat,ym)^2
    } 
    ISE_s=sum( dis )/length(dis)*abs(xgrid[2]-xgrid[1]);
    dis=rep(0, length(xgrid));
    for (j in 1:length(xgrid)){
      ymhat=subset( a_p, a_p[,1]==xgrid[j])
      ym=subset( real, real[,1]==xgrid[j])
      dis[j] = hausdorff_dist(ymhat,ym)^2
    } 
    ISE_p=sum( dis )/length(dis)*abs(xgrid[2]-xgrid[1]);
    dis=rep(0, length(xgrid));
    for (j in 1:length(xgrid)){
      ymhat=subset(a_imput1, a_imput1[,1]==xgrid[j])
      ym=subset( real, real[,1]==xgrid[j])
      dis[j] = hausdorff_dist(ymhat,ym)^2
    } 
    ISE_imput1=sum( dis )/length(dis)*abs(xgrid[2]-xgrid[1]);
    dis=rep(0, length(xgrid));
    for (j in 1:length(xgrid)){
      ymhat=as.matrix(subset(a_imput2, a_imput2[,1]==xgrid[j]))
      ym=subset( real, real[,1]==xgrid[j])
      dis[j] = hausdorff_dist(ymhat,ym)^2
    } 
    ISE_imput2=sum( dis )/length(dis)*abs(xgrid[2]-xgrid[1]);
    vek<-cbind(l,n, mean(deltay),ISE_c,ISE_s,ISE_p,ISE_imput1,ISE_imput2)
    lvec=length(vek)
    write(t(vek),file=ARCHIVOSALIDA,ncolumns=lvec,append=T)    
  }  
}