################################################################################
## REGRESIÓN MODAL CON DATOS MISSING NA COVARIABLE (unidimensional)
##           aplicación a un caso real
################################################################################


# Paquetes necesarios
pkg<-c("hdrcde", "nor1mix","multimode", "np", "lpme", "pracma", 
       "RcppArmadillo", "inline", "rgl", "RColorBrewer")
for(p in pkg){
	if(!require(p, character.only=TRUE)){
 		install.packages(p,dependencies=TRUE)
 		library(p, character.only=TRUE)
	}
}
  # función de estimación de la reg modal
  source("moderegwe.R")
  # estimación de la ventana
  source("moderegbw_we.R")
  # rutinas de C modificadas necesarias para la estimación y la ventana
  source("LCFITmodereg1.R")
  source("CVmode_LCfit1.R")


######################################
## Funcións auxiliares
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


######################################
## Estimación simplificada - datos completos
#####################################
#x=dat[,c("age")];y=dat$cd496;method="CV-mode"; p.order=0; nstart=4; h1=NULL;h2=NULL;weight=rep(1,length(dat$cd496))
SimpModeReg <- function(x, y, method="CV-mode", H, p.order=0, nstart=4,
                         weight=rep(1,length(y[complete.cases(y)])), xnew)
{
 if (length(y[complete.cases(y)])==length(y)) stop("No missing data in the response")
 x <- x[complete.cases(y)]
 y <- y[complete.cases(y)]

  if(missing(xnew)) {
   if (ncol(as.data.frame(x))==1) xnew <- seq(min(x),max(x), length=25)
   if (ncol(as.data.frame(x))==2) stop("Pendente")
   }
  # Estimación del parámetro por validación cruzada o bootstrap  Zhou y Huang(2019)
    # modificado da función lpme::moderegbw
  if(missing(H))  ## H vector c(h1,h2)
  H <- moderegbw.w(X=x, Y=y, method=method, p.order=p.order, 
                       weight=weight)$bw
  # Estimación de la regresión modal. Suponemos nstart=4  
#  Y=y; W<- x; xgrid=xnew; bw=H; p.order=0; mesh=NULL
  fit_s=modereg.w(W=x, Y=y, xgrid=xnew, bw=H, nstart=nstart, 
                  sig=NULL, p.order=0, weight=weight,PLOT=FALSE)
  fit_s$H <- H
  return(fit_s)
}


######################################
## Estimación Ponderada
#####################################

WeigModeReg <- function(x, y, method="CV-mode", H, p.order=0, nstart=4, 
                         weight, xnew)
{
  if (missing(xnew)) {
   if (ncol(as.data.frame(x))==1) xnew <- seq(min(x),max(x), length=25)
   if (ncol(as.data.frame(x))==2) stop("Pendente")
   }
  if (missing(weight)||is.null(weight)) {          # Axuste GLS
      delta <- integer(length=length(y))
      delta[!is.na(y)] <- 1
      pi <- glm(delta ~ x, data=data.frame(x=x,delta=delta), family=binomial)
#        summary(pi)
        weight <- 1/pi$fitted.values
#      plot(x, pi$fitted.values, type='p', ylim=c(0,1)); points(x,delta, pch=19,col="blue")
      }
  # Estimación del parámetro por validación cruzada o bootstrap  Zhou y Huang(2019)
    # modificado da función lpme::moderegbw
      x <- x[complete.cases(y)]
      y <- y[complete.cases(y)]
  if(missing(H))  ## H vector c(h1,h2)
    H <-  moderegbw.w(Y=y, X=x, method=method, p.order=p.order,
                     weight=weight[complete.cases(y)])$bw
    # Estimación de la regresión modal. Suponemos nstar=4  
    fit_w=modereg.w(Y=y, W=x, xgrid=xnew, bw=H, nstart=nstart,sig=NULL, p.order=0,
                    weight=weight[complete.cases(y)],PLOT=FALSE)
   fit_w$H <- H
    return(fit_w) 
 }
 
 
######################################
## Estimación Imputación Simple
##################################### 
SImpModeReg <- function(x, y, method=method, H, nstart=nstart, xnew)
{
 dat1 <- data.frame(x=x[complete.cases(y)], y=y[complete.cases(y)])
 dat0 <- data.frame(x=x[!complete.cases(y)],y=y[!complete.cases(y)])

# Axuste previo --> simplificado
  if(missing(H))  ## H vector c(h1,h2)
    H <- moderegbw.w(X=dat1$x, Y=dat1$y, method="CV-mode", p.order=0, 
                       weight=rep(1,nrow(dat1)))$bw
  
  est_0 <- modereg.w(Y=dat1$y, W=dat1$x, xgrid=dat0$x, bw=H, nstart=nstart, 
                      sig=NULL, p.order=0, weight=rep(1,nrow(dat1)),PLOT=FALSE) 
  # Calculamos la densidad condicional y en el candidato a moda donde alcance el máximo.
  xhat=dat1$x
  yhat=dat1$y
  
  ### repite a función ata que non da erro!!!
  n.try <- 1
  bw <- try(hdrcde::cde.bandwidths(x=dat1$x,y=dat1$y,deg=0,method=1),silent=TRUE)
    while ('try-error' %in% class(bw)) {
          cat("\t número de intentos ",n.try,"\n")
         	n.try <- n.try +1
     	    bw <- try(cde.bandwidths(x=dat1$x,y=dat1$y,deg=0,method=1),silent=TRUE)
          }

    for(i in 1:length(dat0$x)){
      #primero calculo las modas
      ind1=nstart*i-nstart # 4*i-3
      ind2=nstart*i
      modas=est_0$result$fitted[c(ind1:ind2),2]
      x0=est_0$result$fitted[c(ind1:ind2),1]
      pyx=rep(0,length(x0))
      for (i1 in 1:length(modas)){
  #      pyx[i1]=conditional_density(x0[i1],y0=modas[i1],X=xcomp,Y=ycomp,w=rep(1,length(xcomp)),bwx=hhxy[1],bwy=hhxy[2]) 
        pyx[i1]=conditional_density(x0[i1],y0=modas[i1],X=dat1$x,Y=dat1$y,w=rep(1,length(dat1$x)),bwx=bw$a,bwy=bw$b) 
        }
      moda=modas[which(pyx==max(pyx))]
      xhat=c(xhat,dat1$x[i])
      yhat=c(yhat,moda)
    }
    #ordeno las x
    dat <- cbind(xhat,yhat)
      dat <- dat[order(xhat),]
      n1=nrow(dat)

    hhxy = moderegbw.w(Y=dat[,2], X=dat[,1], method=method, p.order=0,h1=NULL, h2=NULL, 
                       weight=rep(1,n1),)$bw
    simp.fit <- modereg.w(Y=dat[,2], W=dat[,1], xgrid=xnew, bw=hhxy, nstart=nstart,
                             sig=NULL, p.order=0, weight=rep(1,n1),PLOT=FALSE)
                             
    simp.fit$H <- H
    return(simp.fit)    
 }


######################################
## Estimación Imputación Múltiple
##################################### 
MImpModeReg <- function(x, y, method=method, H, nstart=4, xnew, fit0, modes.n="cts")
{
 dat1 <- data.frame(x=x[complete.cases(y)], y=y[complete.cases(y)])
 dat0 <- data.frame(x=x[!complete.cases(y)],y=y[!complete.cases(y)])

  if(missing(H)){  ## H vector c(h1,h2)
#    H <- moderegbw.w(X=dat1$x, Y=dat1$y, method="CV-mode", p.order=0, 
#                       weight=rep(1,nrow(dat1)))$bw
  ### repite a función ata que non da erro!!!
  n.try <- 1
  H <- try(hdrcde::cde.bandwidths(x=dat1$x,y=dat1$y,deg=0,method=1),silent=TRUE)
    while ('try-error' %in% class(H)) {
          cat("\t número de intentos ",n.try,"\n")
         	n.try <- n.try +1
     	    H <- try(cde.bandwidths(x=dat1$x,y=dat1$y,deg=0,method=1),silent=TRUE)
          }
  H <- c(H$a,H$b)
}
# Axuste previo --> simplificado
  if(missing(fit0))  ## H vector c(h1,h2)
    fit0 <- modereg(Y=dat1$y, W=dat1$x, xgrid=dat0$x, bw=H, nstart=nstart, 
                      sig=NULL, p.order=0,PLOT=FALSE)

#    fit0 <- modereg.w(Y=dat1$y, W=dat1$x, xgrid=dat0$x, bw=H, nstart=nstart, 
#                      sig=NULL, p.order=0, weight=rep(1,nrow(dat1)),PLOT=FALSE)

    fit_imput2=list()
    for(nimput in 1:3){
      xhat=dat1$x
      yhat=dat1$y
      for(i in 1:length(dat0$x)){
        #primero calculo las modas
        ind1=nstart*i-(nstart-1)
        ind2=nstart*i
        modas=fit0$result$fitted[c(ind1:ind2),2]
        x0=fit0$result$fitted[c(ind1:ind2),1]
        pyx=rep(0,length(x0))
        for (i1 in 1:length(modas)){
          pyx[i1]=conditional_density(x0[i1], y0=modas[i1], X=dat1$x, Y=dat1$y,
                                      w=rep(1,length(dat1$x)), bwx=H[1], bwy=H[2])
          }
        pro=pyx/sum(pyx)
        moda <- ifelse(length(modas)==1,
               modas,
               sample(modas,size=1,prob=pro)
               )
        xhat=c(xhat,dat0$x[i])
        yhat=c(yhat,moda)
      }
      #ordeno las x
      dat1 <- data.frame(x=xhat,y=yhat)
      dat1 <- dat1[order(xhat),]
      n1=nrow(dat1)
      hhxy_i = moderegbw(Y=dat1$y, X=dat1$x, method=method, p.order=0)$bw
      fit_imput2[[nimput]]=modereg(Y=dat1$y, W=dat1$x, xgrid=xnew, bw=hhxy_i,
                      nstart=nstart,sig=NULL, p.order=0, PLOT=FALSE)

#      hhxy_i = moderegbw.w(Y=dat1$y, X=dat1$x, method=method, p.order=0,
#                           weight=rep(1,n1))$bw
#      fit_imput2[[nimput]]=modereg.w(Y=dat1$y, W=dat1$x, xgrid=xnew, bw=hhxy_i,
#                      nstart=nstart,sig=NULL, p.order=0, weight=rep(1,n1),PLOT=FALSE)
    }
    
    
    a_imput2=list()
    nmodas <- nstart
    for(i in 1:length(xnew)){
      #primero calculo las modas
      ind1=nstart*i-(nstart-1)
      ind2=nstart*i
      modas=c(fit_imput2[[1]]$result$fitted[c(ind1:ind2),2],fit_imput2[[2]]$result$fitted[c(ind1:ind2),2],fit_imput2[[3]]$result$fitted[c(ind1:ind2),2])
      x0=fit_imput2[[1]]$result$fitted[ind1,1]
      # aquí uso a función para estimar o número de modas e logo a que estima as modas 
      #do paquete de Ameijeiras e rosa multimode
      if (modes.n!="cts")
        nmodas <- multimode::nmodes(modas,bw=hhxy_i[2])# hhxy[2]) nas simulacións debería ser o obtido pola imputación simple!!!!
      res <- multimode::locmodes(modas,mod0=nmodas,display=FALSE)
      mod <- extraerimpares(res$locations)
      r <- expand.grid(x0,mod)
      a_imput2 <- rbind(a_imput2,r)
    }
    # Saída --> list()
    xx <- vector("list",length(xnew))    
    for (j in 1:length(xnew))
      xx[[j]] <- subset(a_imput2, a_imput2[,1]==xnew[j])[,2]

      # Saída: x: valores para obter as modas
      #        nmodes: número de modas para cada valor de xç
      #        result: modas
  return(list(x=xnew,nmodes=sapply(xx,length),result=xx, H=H))
}

################################################################################
## Todas as estimacións xuntas
################################################################################
modal.fit.missing <- function(x, y, nstart=2, H, p.order=0, xnew=NULL,
                  i.plot=TRUE, modes.n="cts") {
# p.order <- 0; nstart=2

require(sm)
if(is.null(xnew)) xnew <- seq(min(x), max(x), length=100)

if(sum(is.na(x))>0) stop("Só missing na resposta")
if(sum(is.na(y))==0) stop("Sin missing na resposta")

dat <- data.frame(x=x,y=y)
dat1 <- dat[complete.cases(dat),]

## Regresión Non paramétrica en media -- para comparar
s.fit.mean <- sm.regression(x=dat1$x,y=dat1$y, display="none") #h=1.5,

## Caso simplificado:
s.fit <- SimpModeReg(x=dat$x,y=dat$y, method="CV-mode", H, p.order=p.order,
                         nstart=nstart, weight=rep(1,nrow(dat1)), xnew=xnew)
## Caso Ponderado:
w.fit <- WeigModeReg(x=dat$x, y=dat$y, method="CV-mode", H, p.order=p.order,
                         nstart=nstart, xnew=xnew, weight=NULL)
                        
## Estimación Imputada Simple
simp.fit <- SImpModeReg(x=dat$x, y=dat$y, method="CV-mode", H, nstart=nstart, 
                      xnew=xnew)

## Estimación Imputada Múltiple
mimp.fit <- MImpModeReg(x=dat$x, y=dat$y, method="CV-mode", H, nstart=nstart, 
                      xnew=xnew, modes.n=modes.n)
# objeto Modal Regression with Missing Data 
obj.MRMD <- list(x=x, y=y, nstart=nstart, xnew=xnew, 
                 mean.reg = s.fit.mean, modal.r=s.fit,
                 imp.modal.r=w.fit, simp.modal.r = simp.fit,
                 mimp.modal.r = mimp.fit)
if (i.plot) plot.MRMD(obj.MRMD)
return(obj.MRMD)
}


################################################################################
## Rutinas para debuxar
################################################################################
# Ordea a columna Y e Densidade en función da densidade!
# o df debe ter o mesmo número de modas para cada valor de X
ord.modes <- function(df.modes, nmodes=2) {
  df.modes2 <- df.modes
  if (nmodes==1) stop("Only for modes number greater than 1")
  i <- seq(1,nrow(df.modes), by=nmodes)
  for (j in 1:(nrow(df.modes2)/nmodes-1)) #j <-1
    df.modes2[i[j]:(i[j+1]-1),2:3] <- df.modes2[i[j]-1+order(df.modes2[i[j]:(i[j+1]-1),3], decreasing = TRUE),2:3]
#  cat("\n ",j,"\t",order(df.modes2[i[j]:(i[j+1]-1),3], decreasing = TRUE))
  return(df.modes2)
  }
# Non vale para o obxeto de imputación múltiple!!!!!
graf.modal <- function(dat, obj.modal.est, nmodes=2, mimp=FALSE){
  ifelse(!mimp,
         df.modes <- data.frame(x=obj.modal.est$result$fitted[,1], y=obj.modal.est$result$fitted[,2]),
         df.modes <- data.frame(x=rep(obj.modal.est$x,obj.modal.est$nmodes),y=unlist(obj.modal.est$result))
        )
  n.try <- 1
  H <- try(hdrcde::cde.bandwidths(x=dat$x,y=dat$y,deg=0,method=1),silent=TRUE)
    while ('try-error' %in% class(H)) {
          cat("\t número de intentos ",n.try,"\n")
         	n.try <- n.try +1
     	    H <- try(cde.bandwidths(x=dat$x,y=dat$y,deg=0,method=1),silent=TRUE)
          }
  H <- c(H$a,H$b)
H <- H * 1.1     # Así evito 6 NAN ao principio
 for (i in 1:nrow(df.modes))
     df.modes$cond.density[i]=conditional_density(df.modes$x[i],y0=df.modes$y[i],X=dat$x,Y=dat$y,w=rep(1,length(dat$x)),bwx=H[1],bwy=H[2]) 
 
  if (nmodes>1) df.modes <- ord.modes(df.modes,nmodes=nmodes)
  attr(df.modes, "H") <- H
  return(df.modes)                       
}



## Representación gráfica
plot.MRMD <- function(obj.MRMD, legend=TRUE) {
  datos <- data.frame(x=obj.MRMD$x,y=obj.MRMD$y)
  dat <- datos[complete.cases(datos),]
  n.modas <- obj.MRMD$nstart
  s.fit <- obj.MRMD$modal.r
  s.fit.mean <- obj.MRMD$mean.reg
  w.fit <- obj.MRMD$imp.modal.r
  simp.fit <- obj.MRMD$simp.modal.r
  mimp.fit <- obj.MRMD$mimp.modal.r

# png("2021 Boia de Cies - datos diarios.png", width = 640, height = 480,)
# png("cd496-modes-start4.png", width = 640, height = 480,)
colores <- brewer.pal(n=5, name="Oranges")

par(mfrow=c(2,2),mar = c(2.5,2.25,2.5,0.5))
  # Simplificado
  xx <- graf.modal(dat,s.fit,nmodes=n.modas)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="Simplified estimator")
  if (legend) legend("topleft", legend=levels(cortes), bty = "n", fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
  # Ponderado
  xx <- graf.modal(dat,w.fit,nmodes=n.modas)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="Inversely Probability Weighting estimator")
  if (legend) legend("topleft", legend=levels(cortes), bty = "n", fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")   
   rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
 # Imputación simple  
  xx <- graf.modal(dat,simp.fit,nmodes=n.modas)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="Simple imputation estimator")
  if (legend) legend("topleft", legend=levels(cortes), bty = "n", fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
  # Imputación Múltiple
  xx <- graf.modal(dat,mimp.fit,nmodes=n.modas, mimp=TRUE)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="Multiple imputation estimator")
  if (legend) legend("topleft", legend=levels(cortes), bty = "n", fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
par(mfrow=c(1,1))
# dev.off()
}
## Representación gráfica
plot.MRMD2 <- function(obj.MRMD) {
  datos <- data.frame(x=obj.MRMD$x,y=obj.MRMD$y)
  dat <- datos[complete.cases(datos),]
  n.modas <- obj.MRMD$nstart
  s.fit <- obj.MRMD$modal.r
  s.fit.mean <- obj.MRMD$mean.reg
  w.fit <- obj.MRMD$imp.modal.r
  simp.fit <- obj.MRMD$simp.modal.r
  mimp.fit <- obj.MRMD$mimp.modal.r
  
# png("cd496-modes-start4.png", width = 640, height = 480,)
par(mfrow=c(2,2))
  # Simplificado
  plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="x",
       ylab="y", main="Simplified estimator")
    points(rep(s.fit$xgrid,s.fit$x.num), s.fit$mode, pch=as.character(1:n.modas),
          col="3",lwd="2", cex=0.75) #c("1","2","3","4")
    lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2)
    rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
  # Ponderado
  plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="x",
       ylab="y", main="Inversely Probability Weighting estimator")
    points(rep(w.fit$xgrid,w.fit$x.num), w.fit$mode, pch=as.character(1:n.modas),
          col="3",lwd="2", cex=0.75) #c("1","2","3","4")
    lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2)
    rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
 # Imputación simple  
  plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="x",
       ylab="y", main="Simple imputation estimator")
    points(rep(simp.fit$xgrid,simp.fit$x.num), simp.fit$mode, pch=as.character(1:n.modas),
          col="3",lwd="2", cex=0.75) #c("1","2","3","4")
    lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2)
    rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
  # Imputación Múltiple
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="x",
       ylab="x", main="Multiple imputation estimator")
    points(rep(mimp.fit$x,mimp.fit$nmodes), unlist(mimp.fit$result), pch=19,
          col="3",lwd="2", cex=0.75) #c("1","2","3","4")
    lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2)
    rug(jitter(obj.MRMD$x[!is.na(obj.MRMD$x)], amount = 0.01), side = 1, col = "light blue")
par(mfrow=c(1,1))
# dev.off()
}
