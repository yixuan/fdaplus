#include <RcppEigen.h>

using Rcpp::wrap;
using Rcpp::as;
using Rcpp::List;
using Rcpp::Named;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;

typedef Map<MatrixXd> MapMat;

RcppExport SEXP sqrtm(SEXP mat)
{
BEGIN_RCPP

    const MapMat M(as<MapMat>(mat));
    SelfAdjointEigenSolver<MatrixXd> es(M);

    return wrap(es.operatorSqrt());

END_RCPP
}

RcppExport SEXP sqrtInvm(SEXP mat)
{
BEGIN_RCPP

    const MapMat M(as<MapMat>(mat));
    SelfAdjointEigenSolver<MatrixXd> es(M);

    return wrap(es.operatorInverseSqrt());

END_RCPP
}

RcppExport SEXP sqrtBothm(SEXP mat)
{
BEGIN_RCPP

    const MapMat M(as<MapMat>(mat));
    SelfAdjointEigenSolver<MatrixXd> es(M);

    return Rcpp::List::create(
        Named("sqrt") = wrap(es.operatorSqrt()),
        Named("sqrtInv") = wrap(es.operatorInverseSqrt())
    );

END_RCPP
}

