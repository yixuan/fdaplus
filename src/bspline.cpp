#include <RcppEigen.h>
#include "integral.h"

using namespace Rcpp;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

extern "C" {

SEXP spline_basis(SEXP knots, SEXP order, SEXP xvals, SEXP derivs);

}

class BsplineIntegrand : public VectorIntegrandBatch
{
private:
    int nbasisX;
    int nbasisY;
    SEXP knotsX;
    SEXP knotsY;
    int orderX;
    int orderY;
    IntegerVector SorderX;
    IntegerVector SorderY;
public:
    BsplineIntegrand(int nbasisx_, int nbasisy_,
                     SEXP knotsx_, SEXP knotsy_,
                     int orderx_, int ordery_) :
        VectorIntegrandBatch(nbasisx_ * nbasisy_),
        nbasisX(nbasisx_), nbasisY(nbasisy_),
        knotsX(knotsx_), knotsY(knotsy_),
        orderX(orderx_), orderY(ordery_),
        SorderX(1L, orderX), SorderY(1L, orderY)
    {}
    void eval(const double *x)
    {
        NumericVector t(x, x + npoints);
        IntegerVector derivs(npoints);
        
        NumericMatrix bsmatX = spline_basis(knotsX, SorderX, t, derivs);
        NumericMatrix bsmatY = spline_basis(knotsY, SorderY, t, derivs);
        /*
        There is some complication here. bsmatX[i, j] is NOT
        the value of B_j(x_i), but a compressed version of that.
        This is because for each x0, there are at most "norder"
        basis functions B_j so that B_j(x0) != 0
        For example, if B_j(x_i) is
        
               [,1]    [,2]    [,3]    [,4]    [,5]  [,6]
         [1,] 1.000 0.00000 0.00000 0.00000 0.00000 0.000
         [2,] 0.343 0.54225 0.11025 0.00450 0.00000 0.000
         [3,] 0.064 0.55800 0.34200 0.03600 0.00000 0.000
         [4,] 0.001 0.33075 0.54675 0.12150 0.00000 0.000
         [5,] 0.000 0.12800 0.58800 0.28200 0.00200 0.000
         [6,] 0.000 0.03125 0.46875 0.46875 0.03125 0.000
         [7,] 0.000 0.00200 0.28200 0.58800 0.12800 0.000
         [8,] 0.000 0.00000 0.12150 0.54675 0.33075 0.001
         [9,] 0.000 0.00000 0.03600 0.34200 0.55800 0.064
        [10,] 0.000 0.00000 0.00450 0.11025 0.54225 0.343
        [11,] 0.000 0.00000 0.00000 0.00000 0.00000 1.000
        
        Then trasnpose(bsmatX) is
        
                 [,1]    [,2]    [,3]    [,4]
         [1,] 1.00000 0.00000 0.00000 0.00000
         [2,] 0.34300 0.54225 0.11025 0.00450
         [3,] 0.06400 0.55800 0.34200 0.03600
         [4,] 0.00100 0.33075 0.54675 0.12150
         [5,] 0.12800 0.58800 0.28200 0.00200
         [6,] 0.03125 0.46875 0.46875 0.03125
         [7,] 0.00200 0.28200 0.58800 0.12800
         [8,] 0.12150 0.54675 0.33075 0.00100
         [9,] 0.03600 0.34200 0.55800 0.06400
        [10,] 0.00450 0.11025 0.54225 0.34300
        [11,] 0.00000 0.00000 0.00000 1.00000
        attr(,"Offsets")
         [1] 0 0 0 0 1 1 1 2 2 2 2

        The offset indicates the start location of data in each row of
        B_j(x_i).
        
        Now for each k, outer[i, j] = Bx_i(x_k) * By_j(x_k)
        and outer has at most orderX * orderY nonzero elements.
        i goes from offsetX[k] to offsetX[k] + orderX
        j goes from offsetY[k] to offsetY[k] + orderY
        */
        NumericVector offsetX = bsmatX.attr("Offsets");
        NumericVector offsetY = bsmatY.attr("Offsets");
        result.setZero();
        for(int k = 0; k < npoints; k++)
        {
            MapMat outer(&result(0, k), nbasisX, nbasisY);
            int startX = offsetX[k], endX = offsetX[k] + orderX;
            int startY = offsetY[k], endY = offsetY[k] + orderY;
            double *dataX = &bsmatX(0, k);
            double *dataY = &bsmatY(0, k);
            for(int i = startX; i < endX; i++)
            {
                for(int j = startY; j < endY; j++)
                {
                    outer(i, j) = dataX[i - startX] * dataY[j - startY];
                }
            }
        }
        nevals += npoints;
    }
};



RcppExport SEXP bspline_inprod(SEXP x, SEXP y, SEXP allknotsx,
                               SEXP allknotsy)
{
BEGIN_RCPP

    RObject basisX(x);
    RObject basisY(y);
    NumericVector intrange(basisX.slot("range"));
    NumericVector knotsX(allknotsx);
    NumericVector knotsY(allknotsy);
    int nbasisX = as<int>(basisX.slot("nbasis"));
    int nbasisY = as<int>(basisY.slot("nbasis"));
    int orderX = as<int>(basisX.slot("degree")) + 1;
    int orderY = as<int>(basisY.slot("degree")) + 1;
   
    BsplineIntegrand integr(nbasisX, nbasisY,
                            knotsX, knotsY, orderX, orderY);
    VectorCubatureBatch cuba(&integr, intrange[0], intrange[1]);
    return wrap(cuba.values());
    /*
    return List::create(Named("value") = wrap(cuba.values()),
                        Named("abs.error") = wrap(cuba.errors()),
                        Named("niter") = wrap(cuba.iterations()));
    */
END_RCPP
}
