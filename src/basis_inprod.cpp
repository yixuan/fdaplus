#include <RcppEigen.h>
#include "integral.h"

class BasisInprodIntegrand : public VectorFunc
{
private:
    typedef Rcpp::RObject RObject;
    typedef Rcpp::Environment Environment;
    typedef Rcpp::Function Function;
    typedef Rcpp::NumericVector NumericVector;
    typedef Rcpp::NumericMatrix NumericMatrix;
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;

    RObject basisX;
    RObject basisY;
    Environment fdaplus;
    Function feval;

    void eval(const double &x, double *res)
    {
        MapMat r(res, fun_dim(), 1L);
        this->eval(&x, r);
    }

    void eval(const double *x, MapMat &res)
    {
        int npts = res.cols();

        NumericVector t(x, x + npts);
        NumericMatrix bmatX = feval(basisX, t);
        NumericMatrix bmatY = feval(basisY, t);

        MapMat matX(Rcpp::as<MapMat>(bmatX));
        MapMat matY(Rcpp::as<MapMat>(bmatY));
        for(int j = 0; j < npts; j++)
        {
            MapMat outer(&res(0, j), matX.rows(), matY.rows());
            outer = matX.col(j) * matY.col(j).transpose();
        }
    }
public:
    BasisInprodIntegrand(RObject basisX_, RObject basisY_,
                         int ncoefX_, int ncoefY_) :
        VectorFunc(ncoefX_ * ncoefY_),
        basisX(basisX_), basisY(basisY_),
        fdaplus("package:fdaplus"), feval(fdaplus["feval"])
    {}
};



using Rcpp::RObject;
using Rcpp::NumericVector;
using Rcpp::wrap;
using Rcpp::as;

RcppExport SEXP basis_inprod(SEXP x, SEXP y)
{
BEGIN_RCPP

    RObject basisX(x);
    RObject basisY(y);
    NumericVector intrange(basisX.slot("range"));
    int ncoefX = as<int>(basisX.slot("ncoef"));
    int ncoefY = as<int>(basisY.slot("ncoef"));

    BasisInprodIntegrand integr(x, y, ncoefX, ncoefY);
    VectorCubature cuba(&integr, intrange[0], intrange[1]);
    return wrap(cuba.values());

END_RCPP
}
