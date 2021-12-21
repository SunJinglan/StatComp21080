#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp in Chapter 9 homework
//' @description A Gibbs sampler using Rcpp in Chapter 9 homework
//' @param N the number of samples
//' @param burn the number of between-sample random numbers
//' @param a the parameter a in the density
//' @param b the parameter b in the density
//' @param n the parameter n in the density
//' @return a matrix with \code{n} rows and 2 columns
//' @examples
//' \dontrun{
//' rnC <- gibbsc(100,10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsc(int N,int burn,int a,int b,int n){
  NumericMatrix X(N,2);
  NumericMatrix XX(N-burn,2);
  float x,y;
  X(0,0)=1;
  X(0,1)=0.5;
  for(int i=1;i<N;i++){
    y=X(i-1,1);
    X(i,0)=rbinom(1,n,y)[0];
    x=X(i,0);
    X(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  for(int k=0;k<N-burn;k++){
    XX(k,0)=X(k+burn,0);
    XX(k,1)=X(k+burn,1);
  }
  return XX;
}


