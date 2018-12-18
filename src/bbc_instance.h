#ifndef _BBCINSTANCE_H_
#define _BBCINSTANCE_H_

#include "utils.h"


class BBCInstance
{
    
public:
    BBCInstance();
    
    friend std::ostream& operator<<(std::ostream& out, const BBCInstance& instance);
    friend std::istream& operator>>(std::istream& in, BBCInstance& instance);
    
    int m() const
    {
        return _m;
    }
    
    int k() const
    {
        return _k;
    }
    
    const StringMap& chromosomes() const
    {
        return _chromosomes;
    }
    
    const StringSet& clusters() const
    {
        return _clusters;
    }
    
    const BinMap& bins() const
    {
        return _bins;
    }
    
    const StringSet& samples() const
    {
        return _samples;
    }
    
    const StringMap& mapSampleToIdx() const
    {
        return _mapSampleToIdx;
    }

    const DoubleBBRecord& binsToSampleRdBaf() const
    {
        return _binsToSampleRdBaf;
    }

    const IntBBRecord& binsToSampleSCAB() const
    {
        return _binsToSampleSCAB;
    }

    const ClusterBBRecord& binsToCluster() const
    {
        return _binsToCluster;
    }
    
private:
    
    ///
    int _m;
    ///
    int _k;
    ///
    StringMap _chromosomes;
    ///
    StringSet _clusters;
    ///
    BinMap _bins;
    ///
    StringSet _samples;
    ///
    StringMap _mapSampleToIdx;
    ///
    DoubleBBRecord _binsToSampleRdBaf;
    ///
    IntBBRecord _binsToSampleSCAB;
    ///
    ClusterBBRecord _binsToCluster;

};
std::ostream& operator<<(std::ostream& out, const BBCInstance& instance);
std::istream& operator>>(std::istream& in, BBCInstance& instance);

#endif // _INPUTINSTANCE_H_


