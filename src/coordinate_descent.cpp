#include "coordinate_descent.h"

CoordinateDescent::CoordinateDescent(const DoubleMatrix& FA,
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
                                     const unsigned int iterConvergence,
                                     const unsigned int maxIter,
                                     const int nrSeeds,
                                     const int nrWorkers,
                                     const int nrILPthreads,
                                     const int timeLimit,
                                     const int memoryLimit,
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
    , _nrSeeds(nrSeeds)
    , _nrWorkers(nrWorkers)
    , _nrILPthreads(nrILPthreads)
    , _timeLimit(timeLimit)
    , _memoryLimit(memoryLimit)
    , _sizeBubbles(-1)
    , _seeds()
    , _seedidx(0)
    , _threads()
    , _mutexTasks()
    , _mutexResults()
    , _objs()
    , _CAs()
    , _CBs()
    , _Us()
    , _best(-1)
    , _v(v)
{
    if(mu <= 0.1) {
        _sizeBubbles = 10;
    } else if(mu <= 0.15) {
        _sizeBubbles = 6;
    } else if(mu <= 0.2) {
        _sizeBubbles = 5;
    } else {
        _sizeBubbles = 3;
    }
}


void CoordinateDescent::run()
{
    for (int i = 0; i < _nrSeeds; ++i)
    {
        _seeds.push_back(buildRandomU());
    }
    
    log("Detailed data about coordinate descence is going to be written in: coordinate_descence.tsv", VERBOSITY_t::DEBUGGING, _v);
    if(_v >= VERBOSITY_t::VERBOSE)
    {
        std::ofstream coordinate_out("coordinate_descence.tsv", std::ofstream::out);
        coordinate_out << "";
    }
    log("Coordinate Descence {", VERBOSITY_t::VERBOSE, _v);
    for(int i = 0; i < _nrWorkers; ++i)
    {
        _threads.emplace_back(&CoordinateDescent::runWorker, this);
    }

    for(auto& thread : _threads)
    {
        thread.join();
    }
    
    log("}", VERBOSITY_t::VERBOSE, _v);
    _best = std::distance(_objs.begin(), std::min_element(_objs.begin(), _objs.end()));
}


DoubleMatrix CoordinateDescent::buildRandomU()
{
    DoubleMatrix U;
    for(unsigned int p = 0; p < _k; ++p)
    {
        std::uniform_int_distribution<> uni(1, _n);
        int non_zeros = std::max(uni(g_rng), uni(g_rng));
        non_zeros = std::min(non_zeros, _sizeBubbles);
        U.push_back(buildPartitionVector(non_zeros));
    }
    DoubleMatrix T(transpose(U));
    log("SEED " + std::to_string(_seeds.size()) + ":\n" + toString(T), VERBOSITY_t::DEBUGGING, _v);
    return T;
}


DoubleArray CoordinateDescent::buildPartitionVector(const int num_parts)
{
    assert(num_parts <= _sizeBubbles);
    IntArray positions(_n);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), g_rng);
    
    DoubleArray result(_n, 0.0);
    if(num_parts > 1)
    {
        IntArray samples(_sizeBubbles - 1);
        std::iota(samples.begin(), samples.end(), 1);
        assert(samples.front() == 1 && samples.back() == _sizeBubbles - 1);
        std::shuffle(samples.begin(), samples.end(), g_rng);
        
        DoubleArray bubbles(num_parts - 1, 0.0);
        for(unsigned int i = 0; i < num_parts - 1; ++i)
        {
            bubbles[i] = (double)samples.back() / (double)_sizeBubbles;
            samples.pop_back();
        }
        std::sort(bubbles.begin(), bubbles.end());
        
        for(unsigned int i = 0; i < num_parts; ++i)
        {
            if(i == 0)
            {
                result[positions[i]] = bubbles[i];
            }
            else if(i == (num_parts-1))
            {
                result[positions[i]] = 1.0 - bubbles.back();
            }
            else
            {
                result[positions[i]] = bubbles[i] - bubbles[i-1];
            }
            if(_mu - TOL <= result[positions[i]] && result[positions[i]] < _mu)
            {
                result[positions[i]] = _mu;
            }
            assert(result[positions[i]] >= _mu);
        }
    }
    else {
        result[positions[0]] = 1.0;
    }
    return result;
}


void CoordinateDescent::runWorker()
{
    while(true)
    {
        int seed = -1;
        {
            std::lock_guard<std::mutex> lock(_mutexTasks);
            seed = _seedidx < _seeds.size() ? _seedidx++ : -1;
        }
        std::this_thread::sleep_for(std::chrono::seconds(2));
        if(seed > -1) {
            Worker worker(_FA, _FB, _bins, _m, _k, _n, _cmax, _d, _mu, _base, _ampdel, _cn, _iterConvergence, _maxIter,
                          _timeLimit, _memoryLimit, _nrILPthreads, _seeds[seed], seed, _v);
            double obj = worker.solve();
            
            {
                std::lock_guard<std::mutex> lock(_mutexResults);
                _objs.push_back(obj);
                _CAs.push_back(worker.getCA());
                _CBs.push_back(worker.getCB());
                _Us.push_back(worker.getU());
            }
            
            {
                std::lock_guard<std::mutex> lock(g_output_mutex);
                std::locale::global(std::locale("C"));
                log(std::to_string(obj) + "; ", VERBOSITY_t::VERBOSE, _v);
            }
        } else {
            break;
        }
    }
}

