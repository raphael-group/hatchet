
#include "bbc_instance.h"


BBCInstance::BBCInstance()
    : _m()
    , _k()
    , _chromosomes()
    , _clusters()
    , _bins()
    , _samples()
    , _mapSampleToIdx()
    , _binsToSampleRdBaf()
    , _binsToSampleSCAB()
    , _binsToCluster()
{
}


std::ostream& operator<<(std::ostream& out, const BBCInstance& instance)
{
    out << "#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\tCLUSTER" << std::endl;
    for(auto const& chr : instance._chromosomes)
    {
        size_t idx_chr = chr.second;
        for(auto const& bin : instance._bins[idx_chr])
        {
            size_t idx_bin = bin.second;
            for(auto const& sample : instance._binsToSampleRdBaf[idx_chr][idx_bin])
            {
                auto const& scab = instance._binsToSampleSCAB[idx_chr][idx_bin].at(sample.first);
                out << chr.first << '\t' << bin.first.first << '\t' << bin.first.second << '\t' << sample.first << '\t' << sample.second.first << '\t' << std::get<0>(scab) << '\t' << std::get<1>(scab) << '\t' << std::get<2>(scab) << '\t' << std::get<3>(scab) << '\t' << instance._binsToCluster[idx_chr][idx_bin] << std::endl;
            }
        }
    }
    return out;
};


std::istream& operator>>(std::istream& in, BBCInstance& instance)
{
    instance._chromosomes.clear();
    instance._clusters.clear();
    instance._samples.clear();
    instance._bins.clear();
    instance._binsToSampleRdBaf.clear();
    instance._binsToSampleSCAB.clear();
    instance._binsToCluster.clear();

    std::string line;
    std::vector<std::string> tokens;

    while(!in.eof())
    {
        std::getline(in, line, '\n');
        if(!line.empty() && line[0] != '#')
        {
            std::replace(line.begin(), line.end(), '\t', ' ');
            tokens = split(rtrim(ltrim(line)), ' ');

            std::string chr(tokens[0]);
            size_t idx_chr = -1;
            if(instance._chromosomes.find(chr) == instance._chromosomes.end())
            {
                idx_chr = instance._chromosomes.size();
                instance._chromosomes[chr] = idx_chr;
                instance._bins.emplace_back();
                instance._binsToSampleRdBaf.emplace_back();
                instance._binsToSampleSCAB.emplace_back();
                instance._binsToCluster.emplace_back();
            } else {
                idx_chr = instance._chromosomes.at(chr);
            }

            std::string sample(tokens[3]);
            instance._samples.insert(sample);

            std::pair<long long, long long> bin = std::make_pair(std::stoll(tokens[1]), std::stoll(tokens[2]));
            size_t idx_bin = -1;
            if(instance._bins[idx_chr].find(bin) == instance._bins[idx_chr].end())
            {
                idx_bin = instance._bins[idx_chr].size();
                instance._bins[idx_chr][bin] = idx_bin;

                assert(instance._binsToSampleRdBaf[idx_chr].size() == idx_bin);
                instance._binsToSampleRdBaf[idx_chr].emplace_back();
                instance._binsToSampleRdBaf[idx_chr][idx_bin][sample] = std::make_pair(std::stod(tokens[4]), std::stod(tokens[9]));

                assert(instance._binsToSampleSCAB[idx_chr].size() == idx_bin);
                instance._binsToSampleSCAB[idx_chr].emplace_back();
                instance._binsToSampleSCAB[idx_chr][idx_bin][sample] = std::make_tuple(std::stoi(tokens[5]), std::stoi(tokens[6]), std::stoi(tokens[7]), std::stoi(tokens[8]));

                assert(instance._binsToCluster[idx_chr].size() == idx_bin);
                instance._binsToCluster[idx_chr].emplace_back(tokens[10]);
            } else {
                idx_bin = instance._bins[idx_chr].at(bin);

                if(instance._binsToSampleRdBaf[idx_chr][idx_bin].find(sample) != instance._binsToSampleRdBaf[idx_chr][idx_bin].end())
                {
                    throw std::string("Found a bin defined multiple times for the same sample!\n\t" + line);
                } else {
                    instance._binsToSampleRdBaf[idx_chr][idx_bin][sample] = std::make_pair(std::stod(tokens[4]), std::stod(tokens[9]));
                }

                assert(instance._binsToSampleSCAB[idx_chr][idx_bin].find(sample) == instance._binsToSampleSCAB[idx_chr][idx_bin].end());
                instance._binsToSampleSCAB[idx_chr][idx_bin][sample] = std::make_tuple(std::stoi(tokens[5]), std::stoi(tokens[6]), std::stoi(tokens[7]), std::stoi(tokens[8]));

                if(instance._binsToCluster[idx_chr].at(idx_bin) != tokens[10])
                {
                    throw std::string("Found a bin with different clusters in different samples!\n\t" + line);
                }
            }
        }
    }

    for(auto const& chr : instance._chromosomes)
    {
        size_t idx_chr = chr.second;
        for(auto const& bin : instance._bins[idx_chr])
        {
            size_t idx_bin = bin.second;
            StringSet samples;
            for(auto const& sample : instance._binsToSampleRdBaf[idx_chr][idx_bin])
            {
                samples.insert(sample.first);
            }
            if(samples != instance._samples)
            {
                throw std::string("Found a bin not defined for every sample:\tChromosome= " + chr.first + " interval= (" + std::to_string(bin.first.first) + ", " + std::to_string(bin.first.second) + ")");
            }
        }
    }

    return in;
};
