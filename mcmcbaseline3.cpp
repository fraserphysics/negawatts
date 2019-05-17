#include <RcppCommon.h>
#include <iostream>
#include <cmath>
#include <future>
#include <vector>
using namespace std;
using namespace Rcpp;
class MCMCparam{
public:
  int n;
  int p;
  std::vector<double> tau;
  std::vector<double> pi;
  std::vector<int>  wj;
  MCMCparam(int nn, int pp)
    :n{nn}, p{pp}, tau(nn,0.0), pi(nn,0.0), wj(nn*pp)
  { }
  int& operator()(int i, int j) {
    return wj[i+j*n];
  }
};

namespace Rcpp {
  template <> SEXP wrap( const MCMCparam &posterior) ;
}
#include <Rcpp.h>

namespace Rcpp {
  template <> SEXP wrap( const MCMCparam &posterior) {
    Rcpp::List monres;
    NumericMatrix wjj(posterior.n, posterior.p);
    // Fill with value
    int xsize = wjj.nrow() * wjj.ncol();
    for (int i = 0; i < xsize; i++) {
      wjj[i] = posterior.wj[i];
    }
    return 
      Rcpp::wrap(
        Rcpp::DataFrame::create(
	  Rcpp::Named("wj")=wjj,
	  Rcpp::Named("tau")=posterior.tau,
	  Rcpp::Named("pi")=posterior.pi));
      };
    }

// [[Rcpp::export]]
List mcmcbaseline3(IntegerVector paramintt,
		   NumericMatrix ccc, NumericVector yyy,
		   NumericVector hyperparamm, 
		   NumericVector ttau0, 
		   NumericVector ppp0, 
		   IntegerVector wwj0) 
{
  // Random generator (from R);
  RNGScope scope;
  // declarations;
  // variables locales ;
  int ii, kk, jj, kks, sommewjp;
  double w0, w1,  sum0, sum1, RSS, probabilite;
  // parametres de dimension;
  int NN = paramintt["N"];
  int mini = paramintt["min"];
  int nsave = paramintt["nsave"];
  int nburnin = paramintt["nburnin"]+1;
  int niter = nsave+nburnin;
  int TT=ccc.ncol();
  int MM=ccc.nrow();
  // vecteurs de calcul;
  std::vector<int> sommewj(nsave,0);
  std::vector<double> sumA(TT,0);
  std::vector<double> sumB(TT,0);
  std::vector<double> sumC(TT,0);
  // hyperparametres
  double alpha = hyperparamm["alpha"];
  double beta = hyperparamm["beta"];
  double alphap = hyperparamm["alphap"];
  double betap = hyperparamm["betap"];
  // resultat;
  MCMCparam posterior(nsave,MM);
  // les sommes initiales;
  for (jj = 0; jj < MM; jj++){
    posterior(0,jj) = wwj0(jj);
    sommewj[0] += wwj0(jj);
    for (ii = 0; ii < TT; ii++){
      sumA[ii] +=  posterior(0,jj)*ccc(jj,ii);
    }
   }
  
  posterior.tau[0]=ttau0(0);
  posterior.pi[0]=ppp0(0);
  sommewj[0]=max(sommewj[0],mini);
  // boucle ;
  kks=1;
  sommewjp=sommewj[0];
  for (kk=1; kk < niter; kk++){
    for (jj=0; jj < MM; jj++){
      for (ii = 0; ii < TT; ii++){
        sumB[ii] = sumA[ii] - posterior(kks-1,jj)*ccc(jj,ii);
        sumC[ii] = sumB[ii] + ccc(jj,ii);
      }
      sum0=0.0;
      sum1=0.0;
      for (ii = 0; ii < TT; ii++){
        sum0 += pow(yyy(ii)-NN*sumB[ii]/(sommewjp - posterior(kks-1,jj)),2);
        sum1 += pow(yyy(ii)-NN*sumC[ii]/(sommewjp - posterior(kks-1,jj)+1),2);
      }
      w0=-posterior.tau[kks-1]*0.5*sum0;
      w1=-posterior.tau[kks-1]*0.5*sum1;
      probabilite = 1/(1+(1-posterior.pi[kks-1])/posterior.pi[kks-1]*exp(w0-w1));
      if (R::runif(0,1)<probabilite) {
        posterior(kks,jj)=1;
      } else {
        posterior(kks,jj)=0;
      }
      for (ii = 0; ii < TT; ii++){
        sumA[ii] += (posterior(kks,jj) -posterior(kks-1,jj)) * ccc(jj,ii);
      }
      sommewjp+=(posterior(kks,jj) -posterior(kks-1,jj));
      sommewjp=max(sommewjp,mini);
    }
    sommewj[kks]=sommewjp;
    RSS=0.0;
    for (ii = 0; ii < TT; ii++){
      RSS += pow(yyy(ii)-NN*sumA[ii]/sommewjp,2);
    }
    posterior.pi[kks] = R::rbeta(sommewjp+alphap,MM-sommewjp+betap);
    posterior.tau[kks] = R::rgamma(alpha+TT/2,1/(beta+0.5*RSS));
    /* si le compteur arrive au nombre de burnin
       ou si on depasse la capacite du tableau (nburnin >nsave) 
       on repart a 0 dans les tableaux de resultats */
    if (((kks==(nsave-1))&(kk!=(niter-1)))||(kk==nburnin)) {
      // remise a zero des parametres pour debut
      posterior.tau[0]=posterior.tau[kks];
      posterior.pi[0]=posterior.pi[kks];
      sommewj[0]=sommewj[kks];
      for (jj = 0; jj < MM; jj++){
        posterior(0,jj) = posterior(kks,jj);
      }
      // redemarrage
      kks=1;
    } else {
      kks++;   
    }
   }
  return Rcpp::wrap(posterior);
}

MCMCparam mcmcbaseline3mtl(IntegerVector paramintt,
		   NumericMatrix ccc, NumericVector yyy,
		   NumericVector hyperparamm, 
		   NumericVector ttau0, 
		   NumericVector ppp0, 
		   IntegerVector wwj0) 
{
  // declarations;
  // variables locales ;
  int ii, kk, jj, kks, sommewjp;
  double w0, w1,  sum0, sum1, RSS, probabilite;
  // parametres de dimension;
  int NN = paramintt["N"];
  int mini = paramintt["min"];
  int nsave = paramintt["nsave"];
  int nburnin = paramintt["nburnin"]+1;
  int niter = nsave+nburnin;
  int TT=ccc.ncol();
  int MM=ccc.nrow();
  // vecteurs de calcul;
  std::vector<int> sommewj(nsave,0);
  std::vector<double> sumA(TT,0);
  std::vector<double> sumB(TT,0);
  std::vector<double> sumC(TT,0);
  // hyperparametres
  double alpha = hyperparamm["alpha"];
  double beta = hyperparamm["beta"];
  double alphap = hyperparamm["alphap"];
  double betap = hyperparamm["betap"];
  // resultat;
  MCMCparam posterior(nsave,MM);
  // les sommes initiales;
  for (jj = 0; jj < MM; jj++){
    posterior(0,jj) = wwj0(jj);
    sommewj[0] += wwj0(jj);
    for (ii = 0; ii < TT; ii++){
      sumA[ii] +=  posterior(0,jj)*ccc(jj,ii);
    }
   }
  
  posterior.tau[0]=ttau0(0);
  posterior.pi[0]=ppp0(0);
  sommewj[0]=max(sommewj[0],mini);
  // boucle ;
  kks=1;
  sommewjp=sommewj[0];
  for (kk=1; kk < niter; kk++){
    for (jj=0; jj < MM; jj++){
      for (ii = 0; ii < TT; ii++){
        sumB[ii] = sumA[ii] - posterior(kks-1,jj)*ccc(jj,ii);
        sumC[ii] = sumB[ii] + ccc(jj,ii);
      }
      sum0=0.0;
      sum1=0.0;
      for (ii = 0; ii < TT; ii++){
        sum0 += pow(yyy(ii)-NN*sumB[ii]/(sommewjp - posterior(kks-1,jj)),2);
        sum1 += pow(yyy(ii)-NN*sumC[ii]/(sommewjp - posterior(kks-1,jj)+1),2);
      }
      w0=-posterior.tau[kks-1]*0.5*sum0;
      w1=-posterior.tau[kks-1]*0.5*sum1;
      probabilite = 1/(1+(1-posterior.pi[kks-1])/posterior.pi[kks-1]*exp(w0-w1));
      if (R::runif(0,1)<probabilite) {
        posterior(kks,jj)=1;
      } else {
        posterior(kks,jj)=0;
      }
      for (ii = 0; ii < TT; ii++){
        sumA[ii] += (posterior(kks,jj) -posterior(kks-1,jj)) * ccc(jj,ii);
      }
      sommewjp+=(posterior(kks,jj) -posterior(kks-1,jj));
      sommewjp=max(sommewjp,mini);
    }
    sommewj[kks]=sommewjp;
    RSS=0.0;
    for (ii = 0; ii < TT; ii++){
      RSS += pow(yyy(ii)-NN*sumA[ii]/sommewjp,2);
    }
    posterior.pi[kks] = R::rbeta(sommewjp+alphap,MM-sommewjp+betap);
    posterior.tau[kks] = R::rgamma(alpha+TT/2,1/(beta+0.5*RSS));
    /* si le compteur arrive au nombre de burnin
       ou si on depasse la capacite du tableau (nburnin >nsave) 
       on repart a 0 dans les tableaux de resultats */
    if (((kks==(nsave-1))&(kk!=(niter-1)))||(kk==nburnin)) {
      // remise a zero des parametres pour debut
      posterior.tau[0]=posterior.tau[kks];
      posterior.pi[0]=posterior.pi[kks];
      sommewj[0]=sommewj[kks];
      for (jj = 0; jj < MM; jj++){
        posterior(0,jj) = posterior(kks,jj);
      }
      // redemarrage
      kks=1;
    } else {
      kks++;   
    }
   }
  return posterior;
}

//[[Rcpp::export]]
  List mcmcbaseline3mt(IntegerVector paramintt,
		       NumericMatrix ccc, 
		       NumericVector yyy,
		       NumericVector hyperparamm, 
		       NumericVector ttau0, 
		       NumericVector ppp0, 
		       IntegerVector wwj0){

RNGScope scope;
int nchains=paramintt["nchains"];
int MM=ccc.nrow();
int nsave = paramintt["nsave"];

std::vector<MCMCparam> resu(nchains,MCMCparam(nsave,MM));
std::vector<future<MCMCparam> > tasklist; 
    for (int t=0; t<nchains; t++)
      tasklist.push_back( async(launch::async,
				mcmcbaseline3mtl,
				ref(paramintt), ref(ccc), 
ref(yyy), ref(hyperparamm), ref(ttau0), ref(ppp0), ref(wwj0)));

for (int t=0; t<nchains; t++) 
  resu[t]=tasklist[t].get();

return Rcpp::wrap(resu);
  }

//[[Rcpp::export]]
List ICYbaseline3(NumericMatrix ccc, IntegerMatrix wj, NumericVector tau, int NN) {
  RNGScope scope;
  int TT=ccc.ncol();
  int MM=ccc.nrow();
  int KK=wj.nrow();
  int sommewjk;
  NumericMatrix sumA(KK,TT);
  NumericMatrix Yp(KK,TT);
  List res;
  // pour chaque param simule de la posterior
  for (int k = 0; k < KK; k++){
    sommewjk=0;
    // calcul de la moyenne
    for (int ii = 0; ii < TT; ii++){
      for (int jj = 0; jj < MM; jj++){
	sumA(k,ii) +=  wj(k,jj)*ccc(jj,ii);
	if (ii==0) sommewjk += wj(k,jj) ;
      }
      sumA(k,ii) = NN*sumA(k,ii)/sommewjk;
      // random draw from posterior predictive
      Yp(k,ii)=R::rnorm(sumA(k,ii),1/std::sqrt(tau(k)));
    }
  }
  res["mu"]=sumA;
  res["Y"]=Yp;
  return res;
}
