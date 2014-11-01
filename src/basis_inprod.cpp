#include <RcppEigen.h>
#include "integral.h"

using namespace Rcpp;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

class BasisInprod : public VectorIntegrandBatch
{
private:
    RObject basisX;
    RObject basisY;
    Environment fdaplus;
    Function feval;
public:
    BasisInprod(RObject basisX_, RObject basisY_,
                int ncoefX_, int ncoefY_) :
        VectorIntegrandBatch(ncoefX_ * ncoefY_),
        basisX(basisX_), basisY(basisY_),
        fdaplus("package:fdaplus"), feval(fdaplus["feval"])
    {}
    void eval(const double *x)
    {
        NumericVector t(x, x + npoints);
        NumericMatrix bmatX = feval(basisX, t);
        NumericMatrix bmatY = feval(basisY, t);
        
        MapMat matX(as<MapMat>(bmatX));
        MapMat matY(as<MapMat>(bmatY));
        for(int j = 0; j < npoints; j++)
        {
            MapMat outer(&result(0, j), matX.rows(), matY.rows());
            outer = matX.col(j) * matY.col(j).transpose();
        }
        nevals += npoints;
    }
};



RcppExport SEXP basis_inprod(SEXP x, SEXP y)
{
BEGIN_RCPP

    RObject basisX(x);
    RObject basisY(y);
    NumericVector intrange(basisX.slot("range"));
    int ncoefX = as<int>(basisX.slot("ncoef"));
    int ncoefY = as<int>(basisY.slot("ncoef"));
   
    BasisInprod integr(x, y, ncoefX, ncoefY);
    VectorCubatureBatch cuba(&integr, intrange[0], intrange[1]);
    return wrap(cuba.values());

END_RCPP
}
