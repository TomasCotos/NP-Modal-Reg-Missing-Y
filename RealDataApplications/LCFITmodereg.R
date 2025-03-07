library(RcppArmadillo)
library(inline)


code <- "
 
  // Transfer R variables into C++;
  NumericVector x(x_);
  NumericVector y(y_);
  IntegerVector yindx(yindx_);
  NumericVector X(X_);
  NumericVector Y(Y_);
   double h1 = as<double>(h1_);
  double h2 = as<double>(h2_);
   NumericVector w(w_);
  int max_iterations = as<int>(max_iterations_);
  double eps = as<double>(eps_);
  int nx = x.size();
  int ny = y.size();
  int n = X.size();
  
  // results to save 
  NumericVector ym(ny);
  
  // temp variables
  NumericMatrix Ku0ij(n,nx);
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;

  for(int i=0; i<n; ++i){
    for(int j=0; j<nx; ++j){
      Ku0ij(i,j) = exp( -0.5*std::pow(((X[i]-x[j])/h1), 2) );
    }
  }
  
  for(int j=0; j<nx; ++j){
    int ind1 = yindx[j];
    int ind2 = yindx[j+1]-1;
    for(int jj=ind1; jj<=ind2; ++jj){
      newy = y[jj];
      err_now=1e10;
      iter_now=0;
      while((iter_now < max_iterations)&&(err_now > eps)){
        YKG_tot=0;
        KG_tot=0;
        oldy = newy;
        for(int i=0; i<n; ++i){
          KGj = w(i) * Ku0ij(i,j)*exp( -0.5*std::pow(((newy-Y[i])/h2), 2) );
          KG_tot += KGj;
          YKG_tot += Y[i]*KGj;
        }
        if((KG_tot)<1e-10){
          newy=NA_REAL;
          break;
        }else{
          newy = YKG_tot/KG_tot;
          err_now = std::abs(newy-oldy);
          iter_now++;
        }
      }
      if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=NA_REAL;
      ym[jj]=newy;
    }
  }
  return List::create(Named(\"mode\")=ym);
"


LCfitModeReg1 <- cxxfunction(signature(x_="numeric", y_="numeric", yindx_="integer", X_="numeric", Y_="numeric",  h1_="numeric", h2_="numeric", w_="numeric",max_iterations_="integer", eps_="numeric"), code, plugin="RcppArmadillo")



# ahora hay que llamarlo como una funci?n desde R
##
##fit = LCfitModeReg1( x, meshy, yindx, W,Y, bw[1], bw[2],weight,
##                     maxiter, tol)$mode;


