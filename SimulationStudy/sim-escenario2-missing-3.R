source("escenario2.R")

M<-500
n <- 200
# a1=-1.5  # media 
a2=1.5   # media  
b1=0.5  # desv. tipica
b2=0.5  # desv. tipica

################
  tipo.mis<-3

p1=0.75


a1=-1.5  # media 
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
a1=-1.0  # media   
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
a1=-.5  # media   
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
a1=0  # media   
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
a1=0.5  # media   
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
a1=1.0  # media   
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
a1=1.5  # media   
  simulacion(M,n,tipo.mis,a1,a2,b1,b2,p1)
