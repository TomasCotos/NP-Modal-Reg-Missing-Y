library(RcppArmadillo)
library(inline)


code <- "

  // Transfer R variables into C++;
  NumericVector X(X_);
  NumericVector Y(Y_);
  NumericVector pX(pX_);
  NumericVector y(y_);
  IntegerVector yindx(yindx_);
   NumericVector h1(h1_);
  NumericVector h2(h2_);
   NumericVector w(w_);
  int max_iterations = as<int>(max_iterations_);
  double eps = as<double>(eps_);
  int n = X.size();
  int nh1 = h1.size();
  int nh2 = h2.size();
  
  // temp variable
  NumericMatrix CVh(nh1, nh2);
  NumericVector Ku0ijArray(n*n*nh1);
  arma::cube Ku0ij(Ku0ijArray.begin(), n, n, nh1, false);
  double KGj = 0;
  double KG_tot=0;
  double YKG_tot=0;
  double oldy=0;
  double newy=0;
  int iter_now;
  double err_now;
  
  for(int k=0; k<nh1; ++k){
    for(int j=0; j<n; ++j){
      for(int i=0; i<n; ++i){
        Ku0ij(i,j,k) = exp( -0.5*std::pow(((X[i]-X[j])/h1[k]), 2) );
      }
    }
  }
  
  for(int k1=0; k1<nh1; ++k1){
    for(int k2=0; k2<nh2; ++k2){
      R_CheckUserInterrupt();
      double cvj=0;
      for(int j=0; j<n; ++j){
        if(pX[j]>1e-10){
          int ind1 = yindx[j];
          int ind2 = yindx[j+1]-1;
          arma::vec ym(ind2-ind1+1);
          for(int jj=ind1; jj<=ind2; ++jj){
            newy = y[jj];
            err_now=1e10;
            iter_now=0;
            while((iter_now < max_iterations)&&(err_now > eps)){
              //R_CheckUserInterrupt();
              YKG_tot=0;
              KG_tot=0;
              oldy = newy;
              for(int i=0; i<n; ++i){
                if(i!=j){
                  KGj =w(i) * Ku0ij(i,j,k1)*exp( -0.5*std::pow(((newy-Y[i])/h2[k2]), 2) );
                  KG_tot += KGj;
                  YKG_tot += Y[i]*KGj;
                }
              }
              if((KG_tot)<1e-10){
                newy=99999;
                break;
              }else{
                newy = YKG_tot/KG_tot;
                err_now = std::abs(newy-oldy);
                iter_now++;
              }
            }
            if((iter_now==max_iterations)&&(err_now > (eps*10))) newy=99999;
            ym[jj-ind1]=newy;
          }
          int ii=0;
          int ym_num = ym.n_elem;
          while((abs(ym[ii]-99999)<1e-20)&(ii<ym_num)){
            ++ii;
          }
          if(ii==y.size()){
            ym.fill(99999);
          }else{
            for(int jj=0; jj<ym_num; ++jj){
              if(abs(ym[ii]-99999)<1e-20) ym[jj] = ym[ii];
            }
          }
          ym = arma::round(ym*100);
          arma::vec ymnew = arma::unique(ym)/100;
          cvj += std::pow(arma::min(arma::abs(ymnew-Y[j]))*(ymnew.n_elem), 2)*pX[j]*w[j];
        }
      }
      CVh(k1,k2) = cvj/(n+0.0);
    }
  }
  
  return List::create(Named(\"CV\")=CVh);
"






CVmode_LCfit1<- cxxfunction(signature(X_="numeric", Y_="numeric",pX_="numeric", y_="numeric",yindx_="integer",h1_="numeric", h2_="numeric", w_="numeric",max_iterations_="integer", eps_="numeric"), code, plugin="RcppArmadillo")





