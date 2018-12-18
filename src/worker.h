#ifndef _WORKER_H_
#define _WORKER_H_


#include "utils.h"
#include "ilp-min.h"
#include "input_instance.h"

class Worker
{
public:
    Worker(const DoubleMatrix& FA,
           const DoubleMatrix& FB,
           const DoubleArray& bins,
           const int m,
           const int k,
           const int n,
           const int cmax,
           const int d,
           const double mu,
           const int base,
           const bool ampdel,
           const CNMap& cn,
           const int iterConvergence,
           const int maxIter,
           const int timeLimit,
           const int memoryLimit,
           const int nrThreads,
           const DoubleMatrix& M0,
           const int seedIndex,
           const VERBOSITY_t v);
    
    /// Solve for given seed M0
    double solve();
    
    const IntMatrix getCA() const
    {
        return _allCA.back();
    }
    
    const IntMatrix getCB() const
    {
        return _allCB.back();
    }
    
    const DoubleMatrix getU() const
    {
        return _allU.back();
    }
    
private:
    ///
    const DoubleMatrix& _FA;
    ///
    const DoubleMatrix& _FB;
    ///
    const DoubleArray& _bins;
    ///
    const int _m;
    ///
    const int _k;
    /// Number of clones
    const int _n;
    ///
    const int _cmax;
    ///
    const int _d;
    ///
    const double _mu;
    ///
    const int _base;
    ///
    const bool _ampdel;
    ///
    const CNMap& _cn;
    /// Number of iterations for checking convergence
    const unsigned int _iterConvergence;
    /// Maximum number of iterations for each seed
    const unsigned int _maxIter;
    /// Time limit (seconds)
    const int _timeLimit;
    /// Memory limit (MB)
    const int _memoryLimit;
    /// Number of threads
    const int _nrThreads;
    /// Initial M0
    const DoubleMatrix& _M0;
    /// All C matrices
    Int3Matrix _allCA;
    /// All C matrices
    Int3Matrix _allCB;
    /// All M matrices
    Double3Matrix _allU;
    /// All objective values for the C-steps
    DoubleArray _allObjC;
    /// All objective values for the M-steps
    DoubleArray _allObjU;
    /// Seed index
    const int _seedIndex;
    ///
    const VERBOSITY_t _v;
    
};


#endif //_WORKER_H_

