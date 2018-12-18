#ifndef _INPUTINSTANCE_H_
#define _INPUTINSTANCE_H_

#include "utils.h"


class InputInstance
{
    
public:
    InputInstance();
    
    friend std::ostream& operator<<(std::ostream& out, const InputInstance& instance);
    friend std::istream& operator>>(std::istream& in, InputInstance& instance);
    
    int m() const
    {
        return _m;
    }
    
    int k() const
    {
        return _k;
    }
    
    const DoubleMatrix& R() const
    {
        return _R;
    }
    
    const DoubleMatrix& B() const
    {
        return _B;
    }
    
    const DoubleArray& w() const
    {
        return _w;
    }
    
    const std::set<std::string>& clusters() const
    {
        return _clusters;
    }
    
    const std::map<std::string, int>& clusterToIdx() const
    {
        return _clusterToIdx;
    }
    
    const std::map<int, std::string>& idxToCluster() const
    {
        return _idxToCluster;
    }
    
    const std::set<std::string>& samples() const
    {
        return _samples;
    }
    
    const std::map<std::string, int>& sampleToIdx() const
    {
        return _sampleToIdx;
    }
    
    const std::map<int, std::string>& idxToSample() const
    {
        return _idxToSample;
    }

    
private:
    /// Number of segmental clusters
    int _m;
    /// Number of samples
    int _k;
    /// Read-depth ratios F[s][p] for each segmental cluster s in sample p
    DoubleMatrix _R;
    /// B-Allele frequency B[s][p] for each segmental cluster s in sample p
    DoubleMatrix _B;
    /// Bin size w[s] for each segmental cluster s
    DoubleArray _w;
    /// IDs of input segmental clusters
    std::set<std::string> _clusters;
    /// Map cluster to corresponding index
    std::map<std::string, int> _clusterToIdx;
    /// Map index to corresponding cluster
    std::map<int, std::string> _idxToCluster;
    /// Names of input samples
    std::set<std::string> _samples;
    /// Map sample to corresponding index
    std::map<std::string, int> _sampleToIdx;
    /// Map index to corresponding sample
    std::map<int, std::string> _idxToSample;

};

std::ostream& operator<<(std::ostream& out, const InputInstance& instance);
std::istream& operator>>(std::istream& in, InputInstance& instance);

#endif // _INPUTINSTANCE_H_
