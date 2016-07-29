#include <Rcpp.h>
using namespace Rcpp;


//' Rcpp function to compute constant used in \code{spectraldensity} evaluation
//'
//' @param rho parameter of limiting model corresponding to index of regular variation
//' @param alpha vector of Dirichlet allocations (must be a vector of integers)
//' 
//' @keywords internal
//' 
//' @return a constant
// [[Rcpp::export(.dens_const)]]
NumericVector dens_const(NumericVector alpha, NumericVector rho){
  if (rho[0]>=1 || rho[0]<0 || alpha[0]<rho[0] || alpha.size()!=1 || rho.size()!=1)
  Rcpp::stop("Constant: Invalid arguments");
  return exp(lgamma(alpha-rho)-lgamma(alpha)-lgamma(1-rho));
}


//' Rcpp function to compute constant factor used in \code{spectraldensity} evaluation
//'
//' @param rho parameter of limiting model corresponding to index of regular variation
//' @param index the summation index corresponding to current factor
//' 
//' @keywords internal
//' 
//' @return a constant
// [[Rcpp::export(.b)]]
NumericVector b_const(NumericVector index, NumericVector rho = NumericVector::create(1.0)){
  NumericVector result(1);
  if(index[0]==0){
    result[0]=0;
  }
  if(index[0]>=1){
    result[0]=rho[0];
  }
  if(index[0]>=2){
    for(int i=1; i<index[0]; i++){
      result[0]*=(i-rho[0]);
    }
  }
  return result;
}

//' Rcpp function for the spectral density of the Liouville EV model.
//' 
//' Function used to plot the spectral density function. Not well-defined at endpoints.
//'
//' @param w pseudo-angle between (0,1) at which to eveluate the spectraldensity of the CDA of C
//' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
//' @param rho parameter of limiting model corresponding to index of regular variation
//' 
//' 
//' @return Vector of values of spectral density at \code{w}
// [[Rcpp::export(.spectraldensity)]]
NumericVector spectraldensity(NumericVector w, NumericVector alphavec, NumericVector rho){
  //TO DO: Must fix the cases w=0, w=1. Currently handled by hbvevdliouv.
  if (alphavec.size()!=2) Rcpp::stop("Not implemented beyond bivariate case");
  //Redefine constants as NumericVector objects for convenience
  NumericVector a1(1);
  a1[0]= alphavec[0];
  NumericVector a2(1);
  a2[0]= alphavec[1];
  NumericVector k1=dens_const(a1, rho);
  NumericVector k2=dens_const(a2, rho);
  NumericVector result(w.size());
  //placeholders for partial results
  NumericVector result1(w.size());
  NumericVector result2(w.size());
  NumericVector result3(w.size());
  NumericVector ind(1);
  NumericVector ind_a(1);
  NumericVector ind_b(1);
  for(int count=0; count<w.size();count++){
    result1[count]=0;
    result2[count]=0;
    result3[count]=1;
    //First component
    if(alphavec[0]>1){
      for(int l=1;l<alphavec[0];l++){
        ind[0]=l;
        result1[count]=result1[count]+b_const(ind,rho)[0]/gamma(ind)[0]*exp((1/rho[0])*log(k1[0]*w[count])+
        (l/rho[0])*log(k2[0]*(1-w[count]))-
        (l+2)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))))*
        (l*exp((1/rho[0])*log(k1[0]*w[count]))-exp((1/rho[0])*log(k2[0]*(1-w[count]))));
      }
    }
    if(alphavec[1]>1){
      for(int j=1;j<alphavec[1];j++){
        ind[0]=j;
        result1[count]=result1[count]+b_const(ind,rho)[0]/gamma(ind)[0]*exp((j/rho[0])*log(k1[0]*w[count])+
        (1/rho[0])*log(k2[0]*(1-w[count]))-
        (j+2)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))))*
        (j*exp((1/rho[0])*log(k2[0]*(1-w[count])))-exp((1/rho[0])*log(k1[0]*w[count])));
      }
    }
    if(alphavec[0]>0 && alphavec[1]>0){
      for(int l=1;l<alphavec[0];l++){
        for(int j=1;j<alphavec[1];j++){
          ind[0]=j+l;
          ind_a[0]=l+1;
          ind_b[0]=j+1;
          result1[count]=result1[count]+b_const(ind,rho)[0]/(gamma(ind_a)[0]*gamma(ind_b)[0])*exp((j/rho[0])*
          log(k1[0]*w[count])+(l/rho[0])*log(k2[0]*(1-w[count]))-
          (l+j+2)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))))*
          (-(l+2*l*j+j)*exp((1/rho[0])*log(k1[0]*k2[0]*w[count]*(1-w[count])))+pow(j,2.0)*exp((2/rho[0])*
          log(k2[0]*(1-w[count])))+
          pow(l,2.0)*exp((2/rho[0])*log(k1[0]*(w[count]))));
        }
      }
    }
    result1[count]=result1[count]*0.5*exp(rho[0]*log(exp(-(1/rho[0])*log(k1[0]*w[count]))+
    exp(-(1/rho[0])*log(k2[0]*(1-w[count])))))/(pow(rho[0],2.0)*w[count]*(1-w[count]));
    //Second component
    if(alphavec[0]>1){
      for(int l=1;l<alphavec[0];l++){
        ind[0]=l;
        result2[count]=result2[count]+b_const(ind,rho)[0]/gamma(ind)[0]*exp((1/rho[0])*log(k1[0]*w[count])+
        (l/rho[0])*log(k2[0]*(1-w[count]))-
        (l+1)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))));
      }
    }
    if(alphavec[1]>1){
      for(int j=1;j<alphavec[1];j++){
        ind[0]=j;
        result2[count]=result2[count]-b_const(ind,rho)[0]/gamma(ind)[0]*exp((j/rho[0])*log(k1[0]*w[count])+
        (1/rho[0])*log(k2[0]*(1-w[count]))-
        (j+1)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))));
      }
    }
    if(alphavec[0]>0 && alphavec[1]>0){
      for(int l=1;l<alphavec[0];l++){
        for(int j=1;j<alphavec[1];j++){
          ind[0]=j+l;
          ind_a[0]=l+1;
          ind_b[0]=j+1;
          result2[count]=result2[count]+b_const(ind,rho)[0]/(gamma(ind_a)[0]*gamma(ind_b)[0])*exp((j/rho[0])*
          log(k1[0]*w[count])+(l/rho[0])*log(k2[0]*(1-w[count]))-
          (l+j+1)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))))*
          (l*exp((1/rho[0])*log(k1[0]*(w[count])))-j*exp((1/rho[0])*log(k2[0]*(1-w[count]))));
        }
      }
    }
    result2[count]=result2[count]*
    0.5*exp((rho[0]-1)*log(exp(-(1/rho[0])*log(k1[0]*w[count]))+
    exp(-(1/rho[0])*log(k2[0]*(1-w[count])))))/rho[0]*
    (-exp(-(1/rho[0])*log(k1[0])-(1/rho[0]+1)*log(w[count])-log(1-w[count]))+
    exp(-(1/rho[0])*log(k2[0])-(1/rho[0]+1)*log(1-w[count])-log(w[count])));
    //Third component
    if(alphavec[0]>0 && alphavec[1]>0){
      for(int l=0;l<alphavec[0];l++){
        for(int j=0;j<alphavec[1];j++){
          ind[0]=j+l;
          ind_a[0]=l+1;
          ind_b[0]=j+1;
          result3[count]=result3[count]-b_const(ind,rho)[0]/(gamma(ind_a)[0]*gamma(ind_b)[0])*exp((j/rho[0])*
          log(k1[0]*w[count])+(l/rho[0])*log(k2[0]*(1-w[count]))-
          (l+j)*log(exp((1/rho[0])*log(k1[0]*w[count]))+exp((1/rho[0])*log(k2[0]*(1-w[count])))));
        }
      }
    }
    result3[count]=result3[count]*0.5*(1-rho[0])/rho[0]*exp((rho[0]-2)*log(exp(-(1/rho[0])*log(k1[0]*
    w[count]))+exp(-(1/rho[0])*log(k2[0]*(1-w[count]))))-
    (1/rho[0])*log(k1[0]*k2[0])-(1/rho[0]+1)*log(w[count]*(1-w[count])));
    result[count]=-result1[count]+result2[count]+result3[count];
    if(result[count]<0){
      result[count]=0;	
    }
  }
  return result;
  //List::create(Named("r1")=result1,Named("r2")=result2, Named("r3")=result3,Named("result")=result);
}



