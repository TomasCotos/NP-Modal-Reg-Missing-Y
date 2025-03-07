#setwd("C:/users/tcotos/dropbox/Modal Regression/PROGRAMAS_MODAL/datos-reais")
source("funcions-datos-reais.r")


arquivo <-  "2021 Boia de Cíes - datos diarios.txt"   # 23.01% datos missing
datos <- read.table(arquivo, sep="\t", dec=",")


# VARIABLES CONSULTADAS:
# VAR 0:Velocidade do Vento(m/s)	
# VAR 1:Temperatura Media do aire(ºC)	
# VAR 2:Temperatura Máxima do Aire(ºC)	
# VAR 3:Temperatura Mínima do Aire(ºC)	
# VAR 4:Humidade Relativa Media(%)	
# VAR 5:Refacho(m/s)	
# VAR 6:Humidade Relativa Máxima(%)	
# VAR 7:Humidade Relativa Mínima(%)	
# VAR 8:Dirección do Refacho(graos)	
# VAR 9:Temperatura de Orballo(ºC)	
# VAR 10:Dirección do vento predominante(º)	
# VAR 11:Osíxeno(mL/L)	
# VAR 12:Máximo Osíxeno Disuelto(mL/L)	
# VAR 13:Mínimo Osíxeno Disuelto(mL/L)	
# VAR 14:Temperatura (Superficie)(ºC)	
# VAR 15:Conductividade (Superficie)(mS/cm)	
# VAR 16:Salinidade (Superficie)()	
# VAR 17:Anomalía da Densidade (Superficie)(kg/m3)	
# VAR 18:Temperatura Máxima (Superficie)(ºC)	
# VAR 19:Conductividade Máxima (Superficie)(mS/cm)	
# VAR 20:Salinidade Máxima (superficie)()	
# VAR 21:Anomalía da Densidade Máxima (Superficie)(kg/m3)	
# VAR 22:Temperatura Mínima (Superficie)(ºC)	
# VAR 23:Conductividade Mínima (Superficie)(mS/cm)	
# VAR 24:Salinidade Mínima (superficie)()	
# VAR 25:Anomalía da Densidade Mínima (Superficie)(kg/m3)	
nomes <- paste("Var",0:25,sep="")
colnames(datos) <- c("data",nomes)                    
# Transformar -9999 en missing
datos[datos==-9999] <- NA


## Lag 3
ii <- seq(1,nrow(datos),by=3)
df <- data.frame(x=datos$Var1[ii], y=datos$Var18[ii])
  df <- df[complete.cases(df$x),]
  summary(df)
  dim(df)


mm <- modal.fit.missing(x=df$x,y=df$y,nstart=2, p.order=0,
                        H= c(1,1),
                        xnew=seq(min(df$x), max(df$x), length=100) )
                        i.plot=FALSE)
# png("2021 Boia de Cies - datos diarios.png", width = 640, height = 480)
plot.MRMD(mm)
# dev.off()