#include <Rcpp.h>

using namespace Rcpp;

// Calculate
// \int_{lower}^{upper} sin(w * x) dx
//   / cos(w * lower) / w - cos(w * upper) / w      if w != 0
// = |
//   \ 0                                            if w == 0
inline double integral_sin(double w, double lower, double upper)
{
    return fabs(w) < 1e-16 ?
        0.0 :
        (cos(w * lower) - cos(w * upper)) / w;
}

// Calculate
// \int_{lower}^{upper} cos(w * x) dx
//   / sin(w * upper) / w - sin(w * lower) / w      if w != 0
// = |
//   \ upper - lower                                if w == 0
inline double integral_cos(double w, double lower, double upper)
{
    return fabs(w) < 1e-16 ?
        upper - lower :
        (sin(w * upper) - sin(w * lower)) / w;
}

// Calculate
// \int_{lower}^{upper} sin(w1 * x) * cos(w2 * x) dx
// = 0.5 * \int sin((w1 - w2) * x) + sin((w1 + w2) * x) dx
inline double integral_sincos(double w1, double w2,
                              double lower, double upper)
{
    return 0.5 * (integral_sin(w1 - w2, lower, upper) +
                  integral_sin(w1 + w2, lower, upper));
}

// Calculate
// \int_{lower}^{upper} sin(w1 * x) * sin(w2 * x) dx
// = 0.5 * \int cos((w1 - w2) * x) - cos((w1 + w2) * x) dx
inline double integral_sinsin(double w1, double w2,
                              double lower, double upper)
{
    return 0.5 * (integral_cos(w1 - w2, lower, upper) -
                  integral_cos(w1 + w2, lower, upper));
}

// Calculate
// \int_{lower}^{upper} cos(w1 * x) * cos(w2 * x) dx
// = 0.5 * \int cos((w1 - w2) * x) + cos((w1 + w2) * x) dx
inline double integral_coscos(double w1, double w2,
                              double lower, double upper)
{
    return 0.5 * (integral_cos(w1 - w2, lower, upper) +
                  integral_cos(w1 + w2, lower, upper));
}

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
RcppExport SEXP fourier_inprod(SEXP range, SEXP indx, SEXP indy,
                               SEXP periodx, SEXP periody)
{
BEGIN_RCPP

    NumericVector rangeX(range);
    double lb = rangeX[0];
    double ub = rangeX[1];
    
    IntegerVector indX(indx);
    IntegerVector indY(indy);
    
    int nX = indX.length();
    int nY = indY.length();
    
    double TX = as<double>(periodx);
    double TY = as<double>(periody);
    double scale = 2.0 / sqrt(TX * TY);
    
    NumericVector coefX(nX);
    NumericVector coefY(nY);
    IntegerVector typeX(nX);  // 0 - constant, 1 - sin, 2 - cos
    IntegerVector typeY(nY);  // 4 - constant, 8 - sin, 16 - cos
    
    for(int i = 0; i < nX; i++)
    {
        coefX[i] = (indX[i] / 2) * M_2PI / TX;
        if(indX[i] == 1)
            typeX[i] = 0;
        else if(indX[i] % 2 == 0)
            typeX[i] = 1;
        else
            typeX[i] = 2;
    }
    for(int i = 0; i < nY; i++)
    {
        coefY[i] = (indY[i] / 2) * M_2PI / TY;
        if(indY[i] == 1)
            typeY[i] = 4;
        else if(indY[i] % 2 == 0)
            typeY[i] = 8;
        else
            typeY[i] = 16;
    }
    
    NumericMatrix res(nX, nY);   
    double innerprod = 0.0;
    for(int j = 0; j < nY; j++)
    {
        for(int i = 0; i < nX; i++)
        {
            switch(typeX[i] + typeY[j])
            {
            case 4:  // constant - constant
                innerprod = scale * 0.5 * (ub - lb);
                break;
            case 8:  // constant - sin
                innerprod = scale / sqrt(2.0) * integral_sin(coefY[j], lb, ub);
                break;
            case 16: // constant - cos
                innerprod = scale / sqrt(2.0) * integral_cos(coefY[j], lb, ub);
                break;
            case 5:  // sin - constant
                innerprod = scale / sqrt(2.0) * integral_sin(coefX[i], lb, ub);
                break;
            case 9:  // sin - sin
                innerprod = scale * integral_sinsin(coefX[i], coefY[j], lb, ub);
                break;
            case 17: // sin - cos
                innerprod = scale * integral_sincos(coefX[i], coefY[j], lb, ub);
                break;
            case 6: // cos - constant
                innerprod = scale / sqrt(2.0) * integral_cos(coefX[i], lb, ub);
                break;
            case 10: // cos - sin
                innerprod = scale * integral_sincos(coefY[j], coefX[i], lb, ub);
                break;
            case 18: // cos - cos
                innerprod = scale * integral_coscos(coefX[i], coefY[j], lb, ub);
                break;
            default:
                innerprod = 0.0;
            }
            res(i, j) = innerprod;
        }
    }
    
    return res;
    
END_RCPP
}
