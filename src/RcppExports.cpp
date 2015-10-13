#include <Rcpp.h>

using namespace Rcpp;

// Nipals
List Nipals(SEXP Mat, SEXP params);
RcppExport SEXP pcaMethods_Nipals(SEXP MatSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< SEXP >::type Mat(MatSEXP );
        Rcpp::traits::input_parameter< SEXP >::type params(paramsSEXP );
        List __result = Nipals(Mat, params);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
