#include <RcppEigen.h>
#include "integral.h"

using namespace Rcpp;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;

class BsplineIntegrand : public VectorIntegrandBatch
{
private:
    RObject basisX;
    RObject basisY;
    Environment fdaplus;
    Function feval;
public:
    BsplineIntegrand(int nfuns_, RObject &basisX_, RObject &basisY_) :
        VectorIntegrandBatch(nfuns_), basisX(basisX_), basisY(basisY_),
        fdaplus("package:fdaplus"), feval(fdaplus["feval"])
    {}
    void eval(const double *x)
    {
        NumericVector t(x, x + npoints);
        NumericMatrix bsmatX = feval(basisX, t);
        NumericMatrix bsmatY = feval(basisY, t);
        
        MapMat matX(as<MapMat>(bsmatX));
        MapMat matY(as<MapMat>(bsmatY));
        for(int j = 0; j < npoints; j++)
        {
            MapMat outer(&result(0, j), matX.cols(), matY.cols());
            outer = matX.row(j).transpose() * matY.row(j);
        }
        nevals += npoints;
    }
};



RcppExport SEXP bspline_inprod(SEXP x, SEXP y, SEXP rangeval)
{
BEGIN_RCPP

    RObject basisX(x);
    RObject basisY(y);
    NumericVector range(rangeval);
    
    BsplineIntegrand integr(as<int>(basisX.slot("nbasis")) *
                                as<int>(basisY.slot("nbasis")),
                            basisX, basisY);
    VectorCubatureBatch cuba(&integr, range[0], range[1]);
    
    return wrap(cuba.values());
    /*
    return List::create(Named("value") = wrap(cuba.values()),
                        Named("abs.error") = wrap(cuba.errors()),
                        Named("niter") = wrap(cuba.iterations()));
    */
END_RCPP
}
