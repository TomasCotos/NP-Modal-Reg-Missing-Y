
# Supplementary material for Nonparametric modal regression with missing observations in response (NP-Modal-Reg-Missing-Y)

Supplementary codes and data used in paper Ana Pérez-González and
Crujeiras (2025)

Please cite this paper as:

    @article{PerezGonzalez2025,
      title={Nonparametric modal regression with missing observations in response},
      author={Ana Pérez-González, Tomás R. Cotos-Yañez and Rosa M. Crujeiras},
      journal={},
      year={2025}
    }

# R packages

In order to use paper implementation and run all files (numerical and
real examples), the following prerequisites (packages) are needed:

``` r
pkg<-c("hdrcde", "nor1mix","multimode", "np", "lpme", "pracma", "sm",
       "BART")
install.packages(pkg)
```

In order to compile C code is necessary:

``` r
pkg <- c("RcppArmadillo", "inline")
install.packages(pkg)
```

In order to obtain the graphics

``` r
pkg <- c("rgl", "RColorBrewer", "ggplot2", "gridExtra", "reshape2")
install.packages(pkg)
```

# Simulation (Numerical Studies)

- `escenario1.R`, `escenario2.R`, `escenario3.R`, `escenario4.R`: Code
  of scenario 1 – 4 for all simulations.

  - Simplified estimator (S)
  - Inversely Probability Weighting estimator (IPW)
  - Imputed simple estimator (IS)
  - Imputed multiple estimator (IM)

- `sim-escenario1-missing-1.R`, `sim-escenario1-missing-2.R`,
  `sim-escenario1-missing-3.R`, `sim-escenario1-missing-4.R`, scripts of
  scenario 1 with the parameters for each missing data mode given by:

  - M1: $p\left(x\right)=0.6+0.3\cos(\pi x)$
  - M2: $p\left(x\right)=0.6+0.3\cos(2 \pi x)$
  - M3 $p\left(x\right)=0.7+0.3\cos(2\pi x^2)$
  - M4 $p\left(x\right)=0.75$

- idem for all other scenarios.

# Real Data Applications

## 1. Study of HIV-AIDS

Our first example is the `ACTG 175` dataset available from the R-package
`BART` Sparapani, Spanbauer, and McCulloch (2021). To develop our
analysis, conditional modes of `CD4` count at $96\pm 5$ weeks (response
variable) were estimated given `age` variable. The number of missing
observations in the response is 211 (about $39.7\%$ of the sample).

``` r
library(BART)
data(ACTG175)
  datos <- ACTG175[ACTG175$arms==0,c("age","cd496")]
  summary(datos)
```

    ##       age            cd496      
    ##  Min.   :13.00   Min.   :  8.0  
    ##  1st Qu.:29.00   1st Qu.:163.0  
    ##  Median :34.00   Median :283.0  
    ##  Mean   :35.23   Mean   :287.6  
    ##  3rd Qu.:40.00   3rd Qu.:396.0  
    ##  Max.   :70.00   Max.   :857.0  
    ##                  NA's   :211

The achievement of this work is the study of the behavior of missing
data mechanisms in the response to modal regression. Several
methodologies of missing data literature are adapted to this situation:
complete case analysis, inverse probability weighting, simple imputation
or multiple imputation.

Smooth bidimensional estimator of (`age`,`cd496`). Similar to histogram
of Figure 9.

![](README-copy_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The corresponding plot are displayed in Figure 8 of paper.

![](RealDataApplications/cd496-modes-MI.png)

Figure 8: Conditional modes with multiple imputation of cd496 given age
with modes number 1,2,3,4 (from top left corner to bottom right corner)
respectively. Smooth regression estimates in black (with only complete
cases), and modal regression in orange. Orange light colors are shown
for low density values and orange dark colors for high density values.

![](RealDataApplications/cd496-modes-start1.png)

Conditional modes of $\texttt{cd496}$ given $\texttt{age}$ with mode
number 1. Smooth regression estimates in black (with only complete
cases), and modal regression in orange. Orange light colors are shown
for low density values and orange dark colors for high density values.

![](RealDataApplications/cd496-modes-start2.png)

Conditional modes of $\texttt{cd496}$ given $\texttt{age}$ with mode
number 2. Smooth regression estimates in black (with only complete
cases), and modal regression in orange. Orange light colors are shown
for low density values and orange dark colors for high density values.
![](RealDataApplications/cd496-modes-start3.png)

Conditional modes of $\texttt{cd496}$ given $\texttt{age}$ with mode
number 3. Smooth regression estimates in black (with only complete
cases), and modal regression in orange. Orange light colors are shown
for low density values and orange dark colors for high density values.

![](RealDataApplications/cd496-modes-start4.png)

Conditional modes of $\texttt{cd496}$ given $\texttt{age}$ with mode
number 4. Smooth regression estimates in black (with only complete
cases), and modal regression in orange. Orange light colors are shown
for low density values and orange dark colors for high density values.

- `/RealDataApplications/cd496.R`: Code for HIV-AIDS data example.

``` r
source("/RealDataApplications/Example 1. Study of HIV-AIDS.R")
```

## 2. Maximum surface temperature

Now, we have analysed daily data from the ocean-meteorological buoy near
the Cies Islands, located in the municipality of Vigo. (Spain) (lon/lat:
-8°53.59’,42°10.69’). The floating device collects data from different
measurements. One of them is the daily maximum surface temperature,
which had a 23.01% percentage of missing data in 2021.

Our objective is to study the relation between the daily maximum surface
temperature (`Y`) and the average air temperature (`X`). In order to
eliminate the temporal dependence, data each 3 days are selected.

``` r
arquivo <-  "RealDataApplications/2021 Boia de Cíes - datos diarios.txt"
datos <- read.table(arquivo, sep="\t", dec=",")
nomes <- paste("Var",0:25,sep="")
colnames(datos) <- c("data",nomes)                    
# Transformar -9999 en missing
datos[datos==-9999] <- NA
ii <- seq(1,nrow(datos),by=3)
df <- data.frame(x=datos$Var1[ii], y=datos$Var18[ii])
  df <- df[complete.cases(df$x),]
  summary(df)
```

    ##        x               y        
    ##  Min.   : 7.30   Min.   :12.50  
    ##  1st Qu.:13.45   1st Qu.:14.90  
    ##  Median :15.60   Median :16.20  
    ##  Mean   :15.38   Mean   :16.23  
    ##  3rd Qu.:17.30   3rd Qu.:17.30  
    ##  Max.   :23.20   Max.   :21.90  
    ##                  NA's   :26

``` r
  dim(df)
```

    ## [1] 119   2

``` r
# high computational time
# mm <- modal.fit.missing(x=df$x,y=df$y,nstart=2, p.order=0,
#                        H= c(1,1),
#                        xnew=seq(min(df$x), max(df$x), length=100))
```

The corresponding plot is displayed in Figure 9 of paper.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-PerezGonzalez2025" class="csl-entry">

Ana Pérez-González, Tomás R. Cotos-Yañez, and Rosa M. Crujeiras. 2025.
“Nonparametric Modal Regression with Missing Observations in Response.”

</div>

<div id="ref-JSSv097i01" class="csl-entry">

Sparapani, Rodney, Charles Spanbauer, and Robert McCulloch. 2021.
“Nonparametric Machine Learning and Efficient Computation with Bayesian
Additive Regression Trees: The BART r Package.” *Journal of Statistical
Software* 97 (1): 1–66. <https://doi.org/10.18637/jss.v097.i01>.

</div>

</div>
