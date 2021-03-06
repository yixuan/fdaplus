#include <Rcpp.h>

using namespace Rcpp;

/*
Assuming starting point is 0, scale = sqrt(2 / T), omega = 2 * pi / T
B1(x) = scale / sqrt(2)
B2(x) = scale * sin(omega * x)
B3(x) = scale * cos(omega * x)
B4(x) = scale * sin(omega * 2 * x)
B5(x) = scale * cos(omega * 2 * x)
...
B_{2k}(x) = scale * sin(omega * k * x)
B_{2k+1}(x) = scale * cos(omega * k * x)
*/
RcppExport SEXP fourier_feval(SEXP x, SEXP ind, SEXP period)
{
BEGIN_RCPP

    NumericVector xval(x);
    IntegerVector index(ind);
    double periodT = as<double>(period);
    double w = M_2PI / periodT;
    double scale = sqrt(2.0 / periodT);
    
    int p = index.length();
    int T = xval.length();
    NumericMatrix res(p, T);
    
    int k;
    for(int i = 0; i < p; i++)
    {
        k = index[i] / 2;
        if(k == 0)
        {
            res(i, _) = NumericVector(T, scale / sqrt(2.0));
            continue;
        }
        if(index[i] % 2 == 0)
        {
            res(i, _) = scale * sin((w * k) * xval);
        } else {
            res(i, _) = scale * cos((w * k) * xval);
        }
    }
    
    return res;
    
END_RCPP
}
