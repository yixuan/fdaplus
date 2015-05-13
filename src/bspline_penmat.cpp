#include <RcppEigen.h>
#include "integral.h"

extern "C" {

SEXP spline_basis(SEXP knots, SEXP order, SEXP xvals, SEXP derivs);

}

class BsplinePenmatIntegrand : public VectorFunc
{
private:
    typedef Rcpp::RObject RObject;
    typedef Rcpp::IntegerVector IntegerVector;
    typedef Rcpp::NumericVector NumericVector;
    typedef Rcpp::NumericMatrix NumericMatrix;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;

    int nbasis;
    RObject knots;
    int order;
    IntegerVector Sorder;
    int deriv;
    IntegerVector Sderiv;

    void eval(const double &x, double *res)
    {
        MapMat r(res, fun_dim(), 1);
        this->eval(&x, r);
    }

    void eval(const double *x, MapMat &res)
    {
        int npts = res.cols();

        NumericVector t(x, x + npts);
        IntegerVector derivs(npts, deriv);

        NumericMatrix bsmat = spline_basis(knots, Sorder, t, derivs);
        NumericVector offset = bsmat.attr("Offsets");
        res.setZero();
        for(int k = 0; k < npts; k++)
        {
            MapVec outer(&res(0, k), fun_dim());
            int start = offset[k], end = offset[k] + order;
            double *data = &bsmat(0, k);
            double *writer;
            for(int j = start; j < end; j++)
            {
            // For a triangular n by n matrix A, if we only store its
            // lower part including the diagonal elements, then the index
            // of the A[j, j] element is
            //             (2 * nbasis - j + 1) * j / 2
                writer = &outer[(2 * nbasis - j + 1) * j / 2];
                for(int i = j; i < end; i++, writer++)
                {
                    *writer = data[i - start] * data[j - start];
                }
            }
        }
    }

public:
    BsplinePenmatIntegrand(int nbasis_, SEXP knots_, int order_, int deriv_) :
        VectorFunc(nbasis_ * (nbasis_ + 1) / 2),
        nbasis(nbasis_),
        knots(knots_),
        order(order_), Sorder(1, order_),
        deriv(deriv_), Sderiv(1, deriv_)
    {}
};



using Rcpp::RObject;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::wrap;
using Rcpp::as;

RcppExport SEXP bspline_penmat(SEXP x, SEXP allknots, SEXP penderiv)
{
BEGIN_RCPP

    RObject basis(x);
    NumericVector intrange(basis.slot("range"));
    NumericVector knots(allknots);
    int nbasis = as<int>(basis.slot("nbasis"));
    int order = as<int>(basis.slot("degree")) + 1;
    int deriv = as<int>(penderiv);

    BsplinePenmatIntegrand integr(nbasis, knots, order, deriv);
    VectorCubature cuba(&integr, intrange[0], intrange[1]);
    std::vector<double> values = cuba.values();

    NumericMatrix res(nbasis, nbasis);
    for(int j = 0; j < nbasis; j++)
    {
        double *colstart = &values[(2 * nbasis - j + 1) * j / 2];
        res(j, j) = *colstart;
        colstart++;
        for(int i = j + 1; i < nbasis; i++, colstart++)
        {
            res(i, j) = *colstart;
            res(j, i) = *colstart;
        }
    }
    return res;

END_RCPP
}
