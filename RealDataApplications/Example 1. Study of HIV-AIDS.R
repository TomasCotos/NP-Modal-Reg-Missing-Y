#setwd("C:/users/tcotos/dropbox/Modal Regression/PROGRAMAS_MODAL/datos-reais")
source("funcions-datos-reais.r")


library(BART)
data(ACTG175)
  datos <- ACTG175[ACTG175$arms==0,c("age","cd496")]

ll <- vector("list",length=4)
for (i in 1:length(ll)) {
  nstart <- i
  ll[[nstart]] <- modal.fit.missing(x=datos$age, y= datos$cd496, nstart=nstart, p.order=0,
                            xnew=seq(min(datos$age), max(datos$age), length=100),
                            i.plot=FALSE)
# png(paste("cd496-modes-start",nstart,".png",sep=""), width = 640, height = 480,)
   plot.MRMD(ll[[nstart]], legend=FALSE)
#plot.MRMD2(cd496)
# dev.off()
}

# png("cd496-modes-MI.png", width = 640, height = 480)
colores <- brewer.pal(n=5, name="Oranges")
par(mfrow=c(2,2),mar = c(2.5,2.25,1,0.5))
  # Imputación Múltiple
for (i in 1:length(ll)) {
  datos <- data.frame(x=ll[[i]]$x,y=ll[[i]]$y)
  dat <- datos[complete.cases(datos),]
  s.fit.mean <- ll[[i]]$mean.reg
  n.modas <- ll[[i]]$nstart
  xx <- graf.modal(dat, ll[[i]]$mimp.modal.r,nmodes=n.modas, mimp=TRUE)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main=titulo[i])
#  if (legend) legend("topleft", legend=levels(cortes), fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(ll[[i]]$x[!is.na(ll[[i]]$x)], amount = 0.01), side = 1, col = "light blue")
}
par(mfrow=c(1,1))
# dev.off()

# png("cd496-modes-S.png", width = 640, height = 480)
par(mfrow=c(2,2),mar = c(2.5,2.25,1,0.5))
  # Imputación Múltiple
for (i in 1:length(ll)) {
  datos <- data.frame(x=ll[[i]]$x,y=ll[[i]]$y)
  dat <- datos[complete.cases(datos),]
  s.fit.mean <- ll[[i]]$mean.reg
  n.modas <- ll[[i]]$nstart
  xx <- graf.modal(dat, ll[[i]]$modal.r,nmodes=n.modas)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="")
#  if (legend) legend("topleft", legend=levels(cortes), fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(ll[[i]]$x[!is.na(ll[[i]]$x)], amount = 0.01), side = 1, col = "light blue")
}
par(mfrow=c(1,1))
# dev.off()

# png("cd496-modes-W.png", width = 640, height = 480)
par(mfrow=c(2,2),mar = c(2.5,2.25,1,0.5))
  # Imputación Múltiple
for (i in 1:length(ll)) {
  datos <- data.frame(x=ll[[i]]$x,y=ll[[i]]$y)
  dat <- datos[complete.cases(datos),]
  s.fit.mean <- ll[[i]]$mean.reg
  n.modas <- ll[[i]]$nstart
  xx <- graf.modal(dat, ll[[i]]$imp.modal.r,nmodes=n.modas)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="")
#  if (legend) legend("topleft", legend=levels(cortes), fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(ll[[i]]$x[!is.na(ll[[i]]$x)], amount = 0.01), side = 1, col = "light blue")
}
par(mfrow=c(1,1))
# dev.off()


# png("cd496-modes-SI.png", width = 640, height = 480)
par(mfrow=c(2,2),mar = c(2.5,2.25,1,0.5))
  # Imputación Múltiple
for (i in 1:length(ll)) {
  datos <- data.frame(x=ll[[i]]$x,y=ll[[i]]$y)
  dat <- datos[complete.cases(datos),]
  s.fit.mean <- ll[[i]]$mean.reg
  n.modas <- ll[[i]]$nstart
  xx <- graf.modal(dat, ll[[i]]$simp.modal.r,nmodes=n.modas)
   cortes <- cut(xx$cond.density, breaks=5, labels=NULL, dig.lab = 1)
   plot(x=dat$x,y=dat$y, col="blue", pch=19, cex=0.95, xlab="",
     ylab="", main="")
#  if (legend) legend("topleft", legend=levels(cortes), fill=colores)
   points(x=xx$x,y=xx$y, pch=19, col=colores[cortes],lwd="2", cex=0.75)
       lines(x=s.fit.mean$eval.points, y=s.fit.mean$estimate, lwd=2, col="black")
   rug(jitter(ll[[i]]$x[!is.na(ll[[i]]$x)], amount = 0.01), side = 1, col = "light blue")
}
par(mfrow=c(1,1))
# dev.off()



################################################################################
##
##  Figure 9: Histogram of (age, cd496) (left) 
##  and conditional histogram cd496 given age (right)
##
################################################################################
library(ggplot2)
library(gridExtra)
library(reshape2)

datos <- ACTG175[ACTG175$arms==0,]
summary(datos)
# cd496: 211 datos missing de 532
dat <- datos[complete.cases(datos),]

  p <- ggplot(dat,aes(x=age,y=cd496))
  p + geom_point() + 
    stat_density2d(aes(alpha=..density..), #h=c(5,100),
                   geom="tile",contour=FALSE)


data2<-dat[,c("age","cd496")]
  data2$group.x<-cut(data2$age, breaks = 11)
  data2$group.y<-cut(data2$cd496, breaks = 11)
data2 <- data2[,c("group.x","group.y")]
dd <- table(data2)
dd1 <- melt(dd)

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#1º leyend
dd1$group<-cut(dd1$value, breaks = 11)
levels(dd1$group) <- c("0-0.060","0.06,0.12" ,"0.12-0.18","0.18-0.24","0.24-0.30",
       "0.30-0.36","0.36-0.42","0.42-0.48","0.48-0.54","0.54-0.60","0.60-0.66" )
p1 <- ggplot(dd1,aes(x=group.x, y=group.y, fill=group))+geom_tile()+
  theme(legend.position = "right") +
  xlab("age")+ylab("cd496") +
  theme(
#        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
#       axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()     
        ) +
  scale_fill_manual(breaks = levels(dd1$group),
                    values = grey(10:0/11))
mylegend<-g_legend(p1)


## Bidimensional histogram
dd1$group<-cut(dd1$value, breaks = 11)
levels(dd1$group) <- c("0-0.060","0.06,0.12" ,"0.12-0.18","0.18-0.24","0.24-0.30",
       "0.30-0.36","0.36-0.42","0.42-0.48","0.48-0.54","0.54-0.60","0.60-0.66" )
p1 <- ggplot(dd1,aes(x=group.x, y=group.y, fill=group))+geom_tile()+
  theme(legend.position = "none") +
  xlab("age")+ylab("cd496") +
  theme(
#        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
#       axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()     
        ) +
  scale_fill_manual(breaks = levels(dd1$group),
                    values = grey(10:0/11))

dd_all <- dd1

## Conditional histogram
dd.cond <- dd/rowSums(dd)
dd1 <- melt(dd.cond)

dd1$group<-cut(dd1$value, breaks = 11)
levels(dd1$group) <- c("0-0.060","0.06,0.12" ,"0.12-0.18","0.18-0.24","0.24-0.30",
       "0.30-0.36","0.36-0.42","0.42-0.48","0.48-0.54","0.54-0.60","0.60-0.66" )
p2 <- ggplot(dd1,aes(x=group.x, y=group.y, fill=group))+geom_tile()+
  theme(legend.position = "none") +
  xlab("age")+ylab("cd496") +
  theme(
#        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
#       axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()     
        ) +
  scale_fill_manual(breaks = levels(dd1$group),
                    values = grey(10:0/11))
#png("conditional.density.png")
#  grid.arrange(p1, p2, nrow = 1, ncol=2)
#  grid.arrange(p1, p2, mylegend, nrow = 1, ncol=3)
#dev.off()
dd_all <- rbind(dd_all,dd1)


dd_all$Variable.name <- rep(c("bidimensional","conditional"),each=nrow(dd1))
names(dd_all)[4] <- "density"
# png("conditional.density2.png")
ggplot(dd_all,aes(x=group.x, y=group.y, fill=density))+geom_tile()+
facet_wrap(~Variable.name,scales="free_x") + 
  xlab("age")+ylab("cd496") +
  theme(
#        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
#       axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()     
        ) +
  scale_fill_manual(breaks = levels(dd1$group),
                    values = grey(10:0/11))
# dev.off()