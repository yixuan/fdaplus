#include <RcppEigen.h>
#include "cubature.h"

// f(x), x is a scalar, f returns a vector of length nfuns
// result[i] = f_i(x)
class VectorIntegrand
{
protected:
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    int nfuns;
    int nevals;
    MapVec result;
public:
    VectorIntegrand(int nfuns_) : nfuns(nfuns_), result(NULL, nfuns) {}
    void setOutput(double *target_)
    {
        new (&result) MapVec(target_, nfuns);
    }
    // given x, this function should assign f_i(x) to result[i]
    // and nevals++
    virtual void eval(const double *x) = 0;
    int funDim() { return nfuns; }
    int funEvals() { return nevals; }
};

// similar to Integrand, but f can evaluate on a vector of x,
// and return a matrix, with result[i, j] = f_i(x_j)
class VectorIntegrandBatch
{
protected:
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;
    int nfuns;
    int npoints;
    int nevals;
    MapMat result;
public:
    VectorIntegrandBatch(int nfuns_) :
        nfuns(nfuns_), npoints(1), result(NULL, nfuns, npoints) {}
    void setOutput(int npoints_, double *target_)
    {
        npoints = npoints_;
        new (&result) MapMat(target_, nfuns, npoints);
    }
    // given x, this function should assign f_i(x[j]) to result[i, j]
    // and nevals++
    virtual void eval(const double *x) = 0;
    int funDim() { return nfuns; }
    int funEvals() { return nevals; }
};


inline int vintegrand(unsigned ndim, const double *x, void *fdata,
                       unsigned fdim, double *fval)
{
    VectorIntegrand *integr = (VectorIntegrand *) fdata;
    integr->setOutput(fval);
    integr->eval(x);
    return 0;
}

class VectorCubature
{
private:
    VectorIntegrand *integr;
    double lower;
    double upper;
    size_t maxEval;
    double absEps;
    double relEps;
    std::vector<double> value;
    std::vector<double> error;
    int niter;
public:
    VectorCubature(VectorIntegrand *integrand_,
        double lower_,
        double upper_,
        size_t maxEval_ = 2000,
        double absEps_ = 1e-8,
        double relEps_ = 1e-8) :
    integr(integrand_), lower(lower_), upper(upper_),
    maxEval(maxEval_), absEps(absEps_), relEps(relEps_),
    value(integrand_->funDim()), error(integrand_->funDim())
    {
        hcubature(integr->funDim(), vintegrand, integr,
                  1L, &lower, &upper,
                  maxEval, absEps, relEps,
                  ERROR_INDIVIDUAL,
                  value.data(), error.data());
        niter = integr->funEvals();
    }
    std::vector<double> values() { return value; }
    std::vector<double> errors() { return error; }
    int iterations() { return niter; }
};

inline int vintegrand_batch(unsigned ndim, size_t npts, const double *x,
                     void *fdata, unsigned fdim, double *fval)
{
    VectorIntegrandBatch *integr = (VectorIntegrandBatch *) fdata;
    integr->setOutput(npts, fval);
    integr->eval(x);
    return 0;
}

class VectorCubatureBatch
{
private:
    VectorIntegrandBatch *integr;
    double lower;
    double upper;
    size_t maxEval;
    double absEps;
    double relEps;
    std::vector<double> value;
    std::vector<double> error;
    int niter;
public:
    VectorCubatureBatch(VectorIntegrandBatch *integrand_,
        double lower_,
        double upper_,
        size_t maxEval_ = 2000,
        double absEps_ = 1e-8,
        double relEps_ = 1e-8) :
    integr(integrand_), lower(lower_), upper(upper_),
    maxEval(maxEval_), absEps(absEps_), relEps(relEps_),
    value(integrand_->funDim()), error(integrand_->funDim())
    {
        hcubature_v(integr->funDim(), vintegrand_batch, integr,
                    1L, &lower, &upper,
                    maxEval, absEps, relEps,
                    ERROR_INDIVIDUAL,
                    value.data(), error.data());
        niter = integr->funEvals();
    }
    std::vector<double> values() { return value; }
    std::vector<double> errors() { return error; }
    int iterations() { return niter; }
};
