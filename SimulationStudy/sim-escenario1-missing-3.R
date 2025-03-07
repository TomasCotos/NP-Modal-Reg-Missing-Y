source("escenario1.R")

M<-500
n <- 200
a1=-1.5  # media 
a2=1.5   # media  
b1=0.5  # desv. tipica
b2=0.5  # desv. tipica

################
  tipo.mis<-3
  p1=0.5
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
  p1=0.75
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
  p1=0.85
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
