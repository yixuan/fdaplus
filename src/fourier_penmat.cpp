#include <Rcpp.h>
#include "fourier.h"

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


m-th derivative of f(x) = sin(w * x):
    sign * w^m * fun(w * x)
fun(x) = sin(x) if m is even, and fun(x) = cos(x) if m is odd.
sign has the cycle 1, -1, -1, 1, ...

m-th derivative of f(x) = cos(w * x):
    sign * w^m * fun(w * x)
fun(x) = cos(x) if m is even, and fun(x) = sin(x) if m is odd.
sign has the cycle -1, -1, 1, 1, ...
*/
RcppExport SEXP fourier_penmat(SEXP range, SEXP ind, SEXP period, SEXP penderiv)
{
BEGIN_RCPP

    NumericVector range_(range);
    double lb = range_[0];
    double ub = range_[1];
    
    IntegerVector index(ind);
    int p = index.length();
    
    double T = as<double>(period);
    double scale2 = 2.0 / T;
    
    int m = as<int>(penderiv);  // we can assume m > 0 here
    
    NumericVector coef(p);
    IntegerVector type(p);  // 0 - constant, 1 - sin, 2 - cos
    for(int i = 0; i < p; i++)
    {
        coef[i] = (index[i] / 2) * M_2PI / T;
        if(index[i] == 1)
            type[i] = 0;
        else if(index[i] % 2 == 0)
            type[i] = 1;
        else
            type[i] = 2;
    }
    
    NumericMatrix res(p, p);   
    double innerprod = 0.0, tmp;
    int msign = (m % 2 == 0 ? 1 : -1);  // -1 -- m is odd, +1 -- m is even
    for(int j = 0; j < p; j++)
    {
        for(int i = j; i < p; i++)
        {
            if(type[i] * type[j] == 0)  // one of them is constant
            {
                res(i, j) = res(j, i) = 0.0;
                continue;
            }
            tmp =  msign * scale2 * pow(coef[i] * coef[j], (double) m);
            switch(msign * (type[i] + 2 * type[j]))
            {
            case 3:  // sin - sin, m is even
            case -6: // cos - cos, m is odd
                innerprod =  tmp * integral_sinsin(coef[i], coef[j], lb, ub);
                break;
            case -3: // sin - sin, m is odd
            case 6:  // cos - cos, m is even
                innerprod = tmp * integral_coscos(coef[i], coef[j], lb, ub);
                break;
            case 5:  // sin - cos, m is even
            case -4: // cos - sin, m is odd
                innerprod = tmp * integral_sincos(coef[i], coef[j], lb, ub);
                break;
            case -5:  // sin - cos, m is odd
            case 4:   // cos - sin, m is even
                innerprod = tmp * integral_sincos(coef[j], coef[i], lb, ub);
                break;
            default:
                innerprod = 0.0;
            }
            res(i, j) = innerprod;
            res(j, i) = innerprod;
        }
    }
    
    return res;
    
END_RCPP
}
