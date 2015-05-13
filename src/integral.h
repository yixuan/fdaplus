#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <RcppEigen.h>
#include "cubature/cubature.h"

class VectorFunc
{
private:
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;

    const int fdim;
    int nevals;

    // f(x), x is a scalar, f returns a vector of length fdim
    // res[i] = f_i(x)
    virtual void eval(const double &x, double *res) = 0;
    // f(x), x contains npts points, f returns a matrix of size fdim*npts
    // result[i, j] = f_i(x_j)
    virtual void eval(const double *x, MapMat &res)
    {
        int npts = res.cols();
        for(int i = 0; i < npts; i++)
        {
            this->eval(x[i], &res(0, i));
        }
    }

public:
    VectorFunc(const int fdim_) :
        fdim(fdim_), nevals(0)
    {}

    virtual ~VectorFunc() {}

    void eval_and_count(const double *x, const int npts, double *fval)
    {
        MapMat res(fval, fdim, npts);
        this->eval(x, res);
        nevals += npts;
    }
    int fun_dim() { return fdim; }
    int fun_nevals() { return nevals; }
};

inline int vintegrand(unsigned ndim, size_t npts, const double *x,
                      void *fdata, unsigned fdim, double *fval)
{
    VectorFunc *integr = (VectorFunc *) fdata;
    integr->eval_and_count(x, npts, fval);

    return 0;
}

class VectorCubature
{
private:
    VectorFunc *integr;
    const double lower;
    const double upper;
    const size_t max_eval;
    const double abs_eps;
    const double rel_eps;
    std::vector<double> val;
    std::vector<double> err;
    int niter;
public:
    VectorCubature(VectorFunc *integrand_,
        double lower_,
        double upper_,
        size_t max_eval_ = 2000,
        double abs_eps_ = 1e-8,
        double rel_eps_ = 1e-8) :
    integr(integrand_), lower(lower_), upper(upper_),
    max_eval(max_eval_), abs_eps(abs_eps_), rel_eps(rel_eps_),
    val(integr->fun_dim()), err(integr->fun_dim())
    {
        hcubature_v(integr->fun_dim(), vintegrand, integr,
                    1L, &lower, &upper,
                    max_eval, abs_eps, rel_eps,
                    ERROR_INDIVIDUAL,
                    val.data(), err.data());
        niter = integr->fun_nevals();
    }
    std::vector<double> values() { return val; }
    std::vector<double> errors() { return err; }
    int iterations() { return niter; }
};



#endif // INTEGRAL_H
