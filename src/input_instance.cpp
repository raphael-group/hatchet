#include "input_instance.h"


InputInstance::InputInstance()
    : _m(-1)
    , _k(-1)
    , _R()
    , _B()
    , _w()
    , _clusters()
    , _clusterToIdx()
    , _idxToCluster()
    , _samples()
    , _sampleToIdx()
    , _idxToSample()
{
}


std::ostream& operator<<(std::ostream& out, const InputInstance& instance)
{
    out << "#PARAMS" << std::endl;
    out << instance._m << "# number of segmental clusters" << std::endl;
    out << instance._k << "# number of samples" << std::endl;
    out << "#ID\tSAMPLE\t#BINS\tRD\tBAF" << std::endl;
    for(int s = 0; s < instance._m; ++s)
    {
        for(int p = 0; p < instance._k; ++p)
        {
            out << instance._idxToCluster.find(s)->second << "\t" << instance._idxToSample.find(p)->second << "\t" << instance._w[s] << "\t" << instance.R()[s][p] << "\t" << instance.B()[s][p] << std::endl;
        }
    }
    return out;
};



std::istream& operator>>(std::istream& in, InputInstance& instance)
{
    instance._R.clear();
    instance._B.clear();
    
    std::string line;
    std::vector<std::string> tokens;
    
    std::vector<std::pair<std::string, std::string>> keys;
    std::vector<std::pair<double, double>> values;
    
    while(!in.eof())
    {
        std::getline(in, line, '\n');
        if(!line.empty() && line[0] != '#')
        {
            std::replace(line.begin(), line.end(), '\t', ' ');
            tokens = split(rtrim(ltrim(line)), ' ');
            
            std::string cluster(tokens[0]);
            if(instance._clusters.find(cluster) == instance._clusters.end())
            {
                int idx_cluster = instance._clusters.size();
                instance._clusters.insert(cluster);
                instance._clusterToIdx[cluster] = idx_cluster;
                instance._idxToCluster[idx_cluster] = cluster;
                double w = std::atoi(tokens[2].c_str());
                instance._w.push_back(w);
                assert(idx_cluster + 1 == instance._w.size());
            }
            
            std::string sample(tokens[1]);
            if(instance._samples.find(sample) == instance._samples.end())
            {
                int idx_sample = instance._samples.size();
                instance._samples.insert(sample);
                instance._sampleToIdx[sample] = idx_sample;
                instance._idxToSample[idx_sample] = sample;
            }
            
            double rd = std::atof(tokens[3].c_str());
            double baf = std::atof(tokens[8].c_str());
            
            keys.push_back(std::make_pair(cluster, sample));
            values.push_back(std::make_pair(rd, baf));
        }
    }
    
    instance._m = instance._clusters.size();
    instance._k = instance._samples.size();
    
    instance._R = DoubleMatrix(instance._m, DoubleArray(instance._k, 0.0));
    instance._B = DoubleMatrix(instance._m, DoubleArray(instance._k, 0.0));
    
    for(int i = 0; i < keys.size(); ++i)
    {
        int s = instance._clusterToIdx[keys[i].first];
        int p = instance._sampleToIdx[keys[i].second];

        instance._R[s][p] = values[i].first;
        instance._B[s][p] = values[i].second;
    }
    
    return in;
};

