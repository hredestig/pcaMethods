#include <math.h>
#include <Rcpp.h>

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
  SEXP rl = R_NilValue;
  bool cnt;
  int count = 0;
  double tsize;

  char* exceptionMesg = NULL;
  try {
    RcppParams rparams(params);
    int maxSteps = rparams.getIntValue("maxSteps");
    double eps = rparams.getDoubleValue("threshold");
    int nPcs = rparams.getIntValue("nPcs");
    double varLimit = rparams.getDoubleValue("varLimit");
    RcppMatrix<double> mat(Mat);
    RcppMatrix<double> omat(Mat);
    int nr = mat.getDim1();
    int nc = mat.getDim2();
    RcppMatrix<double> est_mat(nr, nc);
    RcppMatrix<double> tt(nr, nPcs);
    RcppMatrix<double> pp(nc, nPcs);
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
	if(!isnan(mat(r,c))) {
	  tss += mat(r,c) * mat(r,c);
	}
      }
    }

    while(anotherPc) {
      for(int r = 0; r < nr; r++) {
	th[r] = 0;
	if(!isnan(mat(r,0))) {
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
	    if(!isnan(mat(r,c))) {
	      ph[c] += mat(r,c) * ti;
	    }
	  }
	}
	
	norm(ph);
	thold = th;
	for(int r = 0; r < nr; r++) {
	  th[r] = 0;
	  for(int c = 0; c < nc; c++) {
	    if(!isnan(mat(r,c))) {
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
	  if(!isnan(mat(r,c))) {
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
    
    RcppResultSet rs;
    if(np != nPcs) {
      RcppMatrix<double> ttt(nr, np);
      RcppMatrix<double> ppp(nc, np);
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
      rs.add("scores", ttt);
      rs.add("loadings", ppp);
    } else {
      rs.add("scores", tt);
      rs.add("loadings", pp);
    }
    rs.add("R2cum", r2cum);
    rl = rs.getReturnList();
  } catch(int e) {
    if(e == 1) {
      exceptionMesg = copyMessageToR("Too many iterations, quitting");
    }else {
      exceptionMesg = copyMessageToR("unknown error");
    }
  } catch(std::exception& ex) {
    exceptionMesg = copyMessageToR(ex.what());
  } catch(...) {
    exceptionMesg = copyMessageToR("unknown error");
  }

  if (exceptionMesg != NULL) Rf_error(exceptionMesg);
  return rl;
}

