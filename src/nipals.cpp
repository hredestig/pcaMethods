#include <math.h>
#include <Rcpp.h>
#include <R.h>

using namespace std;

double difference(vector<double>& vec1, vector<double>& vec2) {
  double diff = 0;
  double a;
  int len = vec1.size();
  for(int i = 0; i < len; i++) {
    a = vec1[i] - vec2[i];
    diff += a * a;
  }
  return(diff);
}

void norm(vector<double>& vec) {
  double siz = 0;
  int len = vec.size();
  for(int i = 0; i < len; i++) {
    siz += vec[i] * vec[i];
  }
  siz = sqrt(siz);
  for(int i = 0; i < len; i++) {
    vec[i] = vec[i] / siz;
  }
}

RcppExport SEXP Nipals(SEXP Mat, SEXP params) {
  try{
    bool cnt;
    int count = 0;
    double tsize;
    Rcpp::List rl = R_NilValue;
  
    Rcpp::List rparams(params);
    int maxSteps = Rcpp::as<int>(rparams["maxSteps"]);
    double eps = Rcpp::as<double>(rparams["threshold"]);
    int nPcs = Rcpp::as<int>(rparams["nPcs"]);
    double varLimit = Rcpp::as<double>(rparams["varLimit"]);
    Rcpp::NumericMatrix mat(Mat);
    Rcpp::NumericMatrix omat = Rcpp::clone<Rcpp::NumericMatrix>( Mat );
    int nr = mat.nrow();
    int nc = mat.ncol();
    Rcpp::NumericMatrix est_mat(nr, nc);
    Rcpp::NumericMatrix tt(nr, nPcs);
    Rcpp::NumericMatrix pp(nc, nPcs);
    vector<double> r2cum;
    vector<double> thold(nr);
    vector<double> th(nr);
    vector<double> phold(nc);
    vector<double> ph(nc);
    double tss = 0;
    double sse = 0;
    int np = 0;
    double anotherPc = true;
    for (int r = 0; r < nr; r++) {
      for (int c = 0; c < nc; c++) {
	if(!ISNAN(mat(r,c))) {
	  tss += mat(r,c) * mat(r,c);
	}
      }
    }

    while(anotherPc) {
      for(int r = 0; r < nr; r++) {
	th[r] = 0;
	if(!ISNAN(mat(r,0))) {
	  th[r] = mat(r,0);
	}
      }

      cnt = true;
      count = 0;
      while(cnt) {
	count++;
	for(int c = 0; c < nc; c++) {
	  ph[c] = 0;
	}
	tsize = 0;
	for(int r = 0; r < nr; r++) {
	  tsize += th[r] * th[r];
	}
	for(int r = 0; r < nr; r++) {
	  double ti = th[r] / tsize;
	  for(int c = 0; c < nc; c++) {
	    if(!ISNAN(mat(r,c))) {
	      ph[c] += mat(r,c) * ti;
	    }
	  }
	}
	
	norm(ph);
	thold = th;
	for(int r = 0; r < nr; r++) {
	  th[r] = 0;
	  for(int c = 0; c < nc; c++) {
	    if(!ISNAN(mat(r,c))) {
	      th[r] += mat(r,c) * ph[c];
	    } 
	  }
	}

	if(count > maxSteps)  {
	  throw 1;
	}
	if(difference(thold, th) <= eps) {
	  cnt = false;
	}
      }
      //deflate mat
      sse = 0;
      double mathat = 0;
      double err = 0;
      for(int r = 0; r < nr; r++) {
	for(int c = 0; c < nc; c++) {
	  if(!ISNAN(mat(r,c))) {
	    mathat = th[r] * ph[c];
	    est_mat(r, c) += mathat;
	    err = omat(r,c) - est_mat(r, c);
	    sse += err * err;
	    mat(r,c) -= mathat;
	  }
	}
      }
      r2cum.push_back(1 - (sse / tss));

      for(int r = 0; r < nr; r++) {
	tt(r,np) = th[r];
      }
      for(int c = 0; c < nc; c++) {
	pp(c,np) = ph[c];
      }
      if(fabs(varLimit - 1) > 1e-4) {
	if(r2cum[np] >= varLimit) {
	  anotherPc = false;
	}
      } 
      if (np + 1 >= nPcs){
	anotherPc = false;
      }
      np++;
    }
    
    if(np != nPcs) {
      Rcpp::NumericMatrix ttt(nr, np);
      Rcpp::NumericMatrix ppp(nc, np);
      for(int r = 0; r < nr; r++) {
	for(int p = 0; p < np; p++) {
	  ttt(r,p) = tt(r,p);
	}
      }
      for(int c = 0; c < nc; c++) {
	for(int p = 0; p < np; p++) {
	  ppp(c,p) = pp(c,p);
	}
      }
      rl["scores"] = ttt;
      rl["loadings"] = ppp;
    } else {
      rl["scores"] = tt;
      rl["loadings"] = pp;
    }
    rl["R2cum"] = r2cum;
    return rl;
  }catch(int e) {
    if(e == 1) {
      ::Rf_error("Too many iterations, quitting");
    }else {
      ::Rf_error("unknown error");
    }
  } catch(std::exception& ex) {
    forward_exception_to_r(ex); 
  } catch(...) {
    ::Rf_error("unknown error");
  }
  return R_NilValue;
}
