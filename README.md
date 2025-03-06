#  Supplementary material for Nonparametric modal regression with missing observations in response (NP-Modal-Reg-Missing-Y) 
Supplementary codes and data used in paper


Please cite this paper as:

```
@article{,
  title={Nonparametric modal regression with missing observations in response},
  author={Ana Pérez-González , Tomás R. Cotos-Yañez and Rosa M. Crujeiras},
  journal={},
  year={2025}
}
```

# Installation  
In order to use paper implementation and run all files (numerical and real examples), the following prerequisites (packages) are needed:


## 1. packages

```{r , eval = FALSE}
paquetes <- c("nor1mix", "pracma", "np", "hdrcde", "multimode", "lpme")
install.packages(paquetes)
```

In order to compile C code is necessary:
```{r , eval = FALSE}
paquetes.2 <- c("RcppArmadillo", "inline")
install.packages(paquetes)
```

# Simulation (Numerical Studies)

+ `simulation.R`: Code for all simulations.

  - Linear smooth (LS)
  - Linear non-smooth (LNS)
  - Nonlinear smooth (NLS)
  - Nonlinear non-smooth (NLNS)
  
+ `sim--missing--1.R`, `sim--missing--2.R`, `sim--missing--3.R`, `sim--missing--4.R`, scripts with the parameters for each simulation scenario (Scenarios 1--4).

# Real Data Applications

Consult a detailed documentation of the data examples and R code of used.

## 1. Study of HIV-AIDS
Our first example is the `ACTG 175` dataset available from the R-package `BART`. To develop our analysis,
conditional modes of `CD4` count at $96\pm 5$ weeks (response variable) were estimated given `age` variable. The number of missing observations in the response is 211 (about $39.7\%$ of the sample).

```{r , eval =T,message=F,warning=F}
library(BART)
data(ACTG175)
  datos <- ACTG175[ACTG175$arms==0,c("age","cd496")]
  summary(datos)
```

The corresponding plots are displayed in Table 4 of paper. 

```{r , eval =T,message=F,warning=F,echo=FALSE}
# incluir gráficos para número de modas == 4

```

Histogram of (`age`,`cd496`) (left) and conditional histogram `cd496` given `age` (rigth). Figure 8 of paper

```{r , eval =F,message=F,warning=F,echo=FALSE}
# Histograma bidi e conditional
```

The achievement of this work is the study of the behavior of missing data mechanisms in the response to modal regression. Several methodologies of missing data literature are adapted to this situation: complete case analysis, inverse probability weighting,  simple imputation or multiple imputation.

+ `/RealDataApplications/cd496.R`: Code for HIV-AIDS data example.

```{r , eval = F}
source("/RealDataApplications/cd496.R")
```

 
## 2. Maximum surface temperature

Now, we have analysed daily data from the ocean-meteorological buoy near the Cies Islands, located in the municipality of Vigo. (Spain) (lon/lat: -8°53.59’,42°10.69’). The floating device collects data from different measurements. One of them is the daily maximum surface temperature, which had a 23.01% percentage of missing data in 2021.

Our objective is to study the relation between the daily maximum surface temperature (`Y`) and the average air temperature (`X`). In order to eliminate the temporal dependence, data each 3 days are selected.

```{r , eval = T}
arquivo <-  "2021 Boia de Cíes - datos diarios.txt"
datos <- read.table(arquivo, sep="\t", dec=",")
nomes <- paste("Var",0:25,sep="")
colnames(datos) <- c("data",nomes)                    
# Transformar -9999 en missing
datos[datos==-9999] <- NA
ii <- seq(1,nrow(datos),by=3)
df <- data.frame(x=datos$Var1[ii], y=datos$Var18[ii])
  df <- df[complete.cases(df$x),]
  summary(df)
  dim(df)

# high computational time
# mm <- modal.fit.missing(x=df$x,y=df$y,nstart=2, p.order=0,
#                        H= c(1,1),
#                        xnew=seq(min(df$x), max(df$x), length=100))
```

The corresponding plot is displayed in Figure 9 of paper. 

## References

Sparapani, R., Spanbauer, C., McCulloch, R., 2021. Nonparametric machine learning and
efficient computation with Bayesian additive regression trees: The BART R package.
Journal of Statistical Software 97, 1–66.
