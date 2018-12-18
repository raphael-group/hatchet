#include "worker.h"


Worker::Worker(const DoubleMatrix& FA,
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
               const VERBOSITY_t v)
    : _FA(FA)
    , _FB(FB)
    , _bins(bins)
    , _m(m)
    , _k(k)
    , _n(n)
    , _cmax(cmax)
    , _d(d)
    , _mu(mu)
    , _base(base)
    , _ampdel(ampdel)
    , _cn(cn)
    , _iterConvergence(iterConvergence)
    , _maxIter(maxIter)
    , _timeLimit(timeLimit)
    , _memoryLimit(memoryLimit)
    , _nrThreads(nrThreads)
    , _M0(M0)
    , _allCA()
    , _allCB()
    , _allU()
    , _allObjC()
    , _allObjU()
    , _seedIndex(seedIndex)
    , _v(v)
{
}


double Worker::solve()
{
    unsigned int iter_convergence = 0;
    unsigned int iter = 0;
    bool first = true;
    
    std::pair<IntMatrix, IntMatrix> hotstart;
    hotstart = IlpSubset::firstHotStart(_FA, _FB, _m, _k, _n, _cmax, _d, _base, _ampdel, _cn);
    
    while((iter_convergence < _iterConvergence) && (iter < _maxIter))
    {
        g_mutex.lock();
        IlpSubset carch(_n, _m, _k, _cmax, _d, _mu, _base, _ampdel, _cn, _FA, _FB, _bins, _v);
        g_mutex.unlock();

        carch.fixU(_allU.empty() ? _M0 : _allU.back());
        carch.init();
        carch.hotStart(hotstart.first, hotstart.second);
        bool status = carch.solve(_timeLimit, _memoryLimit, _nrThreads);
        assert(status);
        
        _allObjC.push_back(carch.getObjs()[0]);
        _allCA.push_back(carch.getACNs()[0]);
        _allCB.push_back(carch.getBCNs()[0]);
        hotstart = std::make_pair(_allCA.back(), _allCB.back());
        //assert(first || _allObjC.back() - TOL <= carch.getObjs()[0]);
        
        if(_v >= VERBOSITY_t::VERBOSE)
        {
            std::lock_guard<std::mutex> lock(g_output_mutex);
            std::ofstream coordinate_out("coordinate_descence.tsv", std::ofstream::out | std::ofstream::app);
            char buf[1024];
            snprintf(buf, 2014, "%d\t%d\t%s\t%f\t%f\t%f", _seedIndex, iter, "C", _allObjC.back(), carch.gap(), carch.runtime());
            coordinate_out << buf << std::endl;
        }
        
        g_mutex.lock();
        IlpSubset uarch(_n, _m, _k, _cmax, _d, _mu, _base, _ampdel, _cn, _FA, _FB, _bins, _v);
        g_mutex.unlock();
        
        uarch.fixC(_allCA.back(), _allCB.back());
        uarch.init();
        status = uarch.solve(_timeLimit, _memoryLimit, _nrThreads);
        assert(status);
        
        _allObjU.push_back(uarch.getObjs()[0]);
        _allU.push_back(uarch.getUs()[0]);
        //assert(first || _allObjU.back() - TOL <= carch.getObjs()[0]);

        if(_allObjC.back() - TOL <= _allObjU.back() && _allObjU.back() <= _allObjC.back() + TOL)
        {
            ++iter_convergence;
        } else {
            iter_convergence = 0;
        }
        
        if(_v >= VERBOSITY_t::VERBOSE)
        {
            std::lock_guard<std::mutex> lock(g_output_mutex);
            std::ofstream coordinate_out("coordinate_descence.tsv", std::ofstream::out | std::ofstream::app);
            char buf[1024];
            snprintf(buf, 2014, "%d\t%d\t%s\t%f\t%f\t%f", _seedIndex, iter, "U", _allObjU.back(), 0.0, uarch.runtime());
            coordinate_out << buf << std::endl;
        }
        
        first = false;
        ++iter;
    }
    assert(iter_convergence >= _iterConvergence | iter == _maxIter);
    return _allObjU.back();
}

