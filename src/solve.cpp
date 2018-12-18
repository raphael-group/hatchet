
#include "utils.h"
#include "argparse.h"
#include "input_instance.h"
#include "bbc_instance.h"
#include "ilp-min.h"
#include "coordinate_descent.h"


CNAssignCluster inferDiploid(const InputInstance& instance, const ArgParse& argparse);
DoubleMatrix scaleReadDepth(const DoubleMatrix &R, const CNAssignCluster &cn, const int m, const int k, const std::map<std::string, int> &clusterToIdx, const VERBOSITY_t verbose);
void splitFractions(const DoubleMatrix &B, const DoubleMatrix &F, const int m, const int k, DoubleMatrix &FA, DoubleMatrix &FB);
void printSolutions(const int k, const double obj, const IntMatrix& CA, const IntMatrix& CB, const DoubleMatrix& U);
void printSolutions(const int k, const int numSolutions, const DoubleArray& obj, const Int3Matrix& CA, const Int3Matrix& CB, const Double3Matrix& U);
void segmentation(const InputInstance& instance, const BBCInstance& bbc, const int n, const IntMatrix& ACN, const IntMatrix& BCN, const DoubleMatrix& U, const std::string& outprefix, const VERBOSITY_t& verbosity);
CNMap getCNMap(const ArgParse& argparse, const InputInstance& instance);
std::string printFractions(const DoubleMatrix &F, const InputInstance &instance);



int main(int argc, char** argv)
{
    try {
        ArgParse argparse(argc, argv);
        
        log("Parsing and checking input arguments", VERBOSITY_t::ESSENTIAL, argparse.v());
        log(argparse.toString(), VERBOSITY_t::VERBOSE, argparse.v());
        g_rng = argparse.rnd() > 0 ? std::mt19937_64(argparse.rnd()) : std::mt19937_64();
        
        log("Reading the input SEG file", VERBOSITY_t::ESSENTIAL, argparse.v());
        InputInstance instance;
        std::ifstream input(argparse.inputSEG());
        input >> instance;
        
        if(argparse.v() >= VERBOSITY_t::DEBUGGING) {std::cerr << instance;}
        
        
        if (argparse.cn().size() > 0)
        {
            log("Scale the read-depth ratios into fractional copy numbers using the provided copy numbers", VERBOSITY_t::ESSENTIAL, argparse.v());
        } else {
            argparse.setClonalCluster(inferDiploid(instance, argparse));
            log("Automatically selecting cluster " + argparse.scalingCN().begin()->first + " as diploid", VERBOSITY_t::ESSENTIAL, argparse.v());
            log("Scale the read-depth ratios into fractional copy numbers using the automaticall selected diploid cluster", VERBOSITY_t::ESSENTIAL, argparse.v());
        }
        DoubleMatrix F(scaleReadDepth(instance.R(), argparse.scalingCN(), instance.m(), instance.k(), instance.clusterToIdx(), argparse.v()));
        log("Fractional copy numbers:\n" + printFractions(F, instance), VERBOSITY_t::DEBUGGING, argparse.v());
        
        CNMap cn(getCNMap(argparse, instance));
        
        int cmax = argparse.cmax() > 0 ? argparse.cmax() : maxCeil(F);
        log("Inferred maximum copy number:" + std::to_string(cmax), VERBOSITY_t::DEBUGGING, argparse.v());
        
        log("Compute allele-specific fractional copy numbers using BAF", VERBOSITY_t::ESSENTIAL, argparse.v());
        DoubleMatrix FA;
        DoubleMatrix FB;
        splitFractions(instance.B(), F, instance.m(), instance.k(), FA, FB);
        log("A-Allele fractional copy numbers:\n" + printFractions(FA, instance), VERBOSITY_t::DEBUGGING, argparse.v());
        log("B-Allele fractional copy numbers:\n" + printFractions(FB, instance), VERBOSITY_t::DEBUGGING, argparse.v());
        
        double obj;
        IntMatrix ACN;
        IntMatrix BCN;
        DoubleMatrix U;
        
        if(argparse.M() == SOLVE_t::BOTH)
        {
            log("Starting coordinate descent algorithm on " + std::to_string(argparse.nrSeeds()) + " seeds", VERBOSITY_t::ESSENTIAL, argparse.v());
            CoordinateDescent cd(FA, FB, instance.w(), instance.m(), instance.k(), argparse.n(), cmax, argparse.d(), argparse.minprop(), argparse.base(),
                                 argparse.forceAMPDEL(), cn, 2, argparse.maxiter(), argparse.nrSeeds(), argparse.j(), 1, argparse.s(), argparse.m(), argparse.v());
            cd.run();
            log("Objective value from coordinate descent:\t" + std::to_string(cd.getObjValue()), VERBOSITY_t::ESSENTIAL, argparse.v());
            
            log("Starting exact ILP for improving solution found by coordinate descent", VERBOSITY_t::ESSENTIAL, argparse.v());
            IlpSubset ilp(argparse.n(), instance.m(), instance.k(), cmax, argparse.d(), argparse.minprop(), argparse.base(),
                          argparse.forceAMPDEL(), cn, FA, FB, instance.w(), argparse.v());
            ilp.init();
            ilp.hotStart(cd.getACN(), cd.getBCN());
            ilp.solve(argparse.s(), argparse.m(), argparse.j());
            
            log("Runtime:\t" + std::to_string(ilp.runtime()), VERBOSITY_t::ESSENTIAL, argparse.v());
            log("Gap:\t" + std::to_string(ilp.gap()), VERBOSITY_t::ESSENTIAL, argparse.v());
            log("Final objective:\t" + std::to_string(ilp.getObjs()[0]), VERBOSITY_t::ESSENTIAL, argparse.v());
            printSolutions(instance.k(), 1, ilp.getObjs(), ilp.getACNs(), ilp.getBCNs(), ilp.getUs());
            
            obj = ilp.getObjs()[0];
            ACN = ilp.getACNs()[0];
            BCN = ilp.getBCNs()[0];
            U = ilp.getUs()[0];

        } else if (argparse.M() == SOLVE_t::ILP)
        {
            log("Starting exact ILP", VERBOSITY_t::ESSENTIAL, argparse.v());
            IlpSubset ilp(argparse.n(), instance.m(), instance.k(), cmax, argparse.d(), argparse.minprop(), argparse.base(),
                          argparse.forceAMPDEL(), cn, FA, FB, instance.w(), argparse.v());
            ilp.init();
            std::pair<IntMatrix, IntMatrix> hotstart;
            hotstart = IlpSubset::firstHotStart(FA, FB, instance.m(), instance.k(), argparse.n(), cmax, argparse.d(), argparse.base(), argparse.forceAMPDEL(), cn);
            ilp.hotStart(hotstart.first, hotstart.second);
            ilp.solve(argparse.s(), argparse.m(), argparse.j());
            
            log("Runtime:\t" + std::to_string(ilp.runtime()), VERBOSITY_t::ESSENTIAL, argparse.v());
            log("Gap:\t" + std::to_string(ilp.gap()), VERBOSITY_t::ESSENTIAL, argparse.v());
            log("Final objective:\t" + std::to_string(ilp.getObjs()[0]), VERBOSITY_t::ESSENTIAL, argparse.v());
            printSolutions(instance.k(), 1, ilp.getObjs(), ilp.getACNs(), ilp.getBCNs(), ilp.getUs());
            
            obj = ilp.getObjs()[0];
            ACN = ilp.getACNs()[0];
            BCN = ilp.getBCNs()[0];
            U = ilp.getUs()[0];

        } else if (argparse.M() == SOLVE_t::CD)
        {
            log("Starting coordinate descent algorithm on " + std::to_string(argparse.nrSeeds()) + " seeds", VERBOSITY_t::ESSENTIAL, argparse.v());
            CoordinateDescent cd(FA, FB, instance.w(), instance.m(), instance.k(), argparse.n(), cmax, argparse.d(), argparse.minprop(), argparse.base(),
                                 argparse.forceAMPDEL(), cn, 2, argparse.maxiter(), argparse.nrSeeds(), argparse.j(), 1, argparse.s(), argparse.m(), argparse.v());
            cd.run();
            
            log("Final objective:\t" + std::to_string(cd.getObjValue()), VERBOSITY_t::ESSENTIAL, argparse.v());
            printSolutions(instance.k(), cd.getObjValue(), cd.getACN(), cd.getBCN(), cd.getU());
            
            obj = cd.getObjValue();
            ACN = cd.getACN();
            BCN = cd.getBCN();
            U = cd.getU();

        } else {
            throw "Unknown solving mode!";
        }
        
        log("Reading the input BBC file", VERBOSITY_t::ESSENTIAL, argparse.v());
        BBCInstance bbc;
        std::ifstream bbffile(argparse.inputBBC());
        bbffile >> bbc;

        log("Segmenting the genome with the inferred copy numbers", VERBOSITY_t::ESSENTIAL, argparse.v());
        segmentation(instance, bbc, argparse.n(), ACN, BCN, U, argparse.o(), argparse.v());

        log("KTHXBYE!", VERBOSITY_t::ESSENTIAL, argparse.v());
    } catch(const std::string &e) {
        std::cerr << "\033[91m" << e << "\033[0m" << std::endl;
        std::exit(EXIT_FAILURE);
    } catch(const char * e) {
        std::cerr << "\033[91m" << e << "\033[0m" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


CNAssignCluster inferDiploid(const InputInstance& instance, const ArgParse& argparse)
{
    double max = -1;
    std::string name;
    for(std::string const& cluster : instance.clusters())
    {
        int s = instance.clusterToIdx().at(cluster);
        bool flag = true;
        for(auto const& sampleR : instance.R()[s])
        {
            flag = flag && (sampleR >= (0.5 - argparse.diploidThreshold()));
        }
        
        if (flag && (instance.w()[s] > max))
        {
            max = instance.w()[s];
            name = cluster;
        }
    }
    
    if (max == -1) {
        throw "No diploid cluster found with a diploid-read-depth-ratio threshold of " + std::to_string(argparse.diploidThreshold());
    }
    
    CNAssignCluster result;
    result.emplace(name, CNState(1, 1));
    
    return result;
}


DoubleMatrix scaleReadDepth(const DoubleMatrix &R, const CNAssignCluster &cn, const int m, const int k, const std::map<std::string, int> &clusterToIdx, const VERBOSITY_t verbose)
{
    std::vector<double> cscale(k);
    
    if(cn.size() == 1)
    {
        log("Scaling using:\t" + cn.begin()->first + ":" + std::to_string(cn.begin()->second.first) + ":" + std::to_string(cn.begin()->second.second), VERBOSITY_t::DEBUGGING, verbose);
        
        for(int p = 0; p < k; ++p)
        {
            double denom = R[(*clusterToIdx.find(cn.begin()->first)).second][p];
            if (denom <= 0)
            {
                throw "Only clusters with a positive read-depth ratio can be selected as diploid segments!";
            }
            cscale[p] = 2.0 / denom;
        }
        
    } else if (cn.size() >= 2) {
        auto it = cn.begin();
        std::string logging = "Scaling using:\t";
        
        int s1 = (*clusterToIdx.find(it->first)).second;
        CNState state1 = it->second;
        int c1 = state1.first + state1.second;
        logging += it->first + ":" + std::to_string(state1.first) + ":" + std::to_string(state1.second) + "\t";
        
        std::advance(it, 1);
        
        int s2 = (*clusterToIdx.find(it->first)).second;
        CNState state2 = it->second;
        int c2 = state2.first + state2.second;
        logging += it->first + ":" + std::to_string(state2.first) + ":" + std::to_string(state2.second);
        
        log(logging, VERBOSITY_t::DEBUGGING, verbose);

        for(int p = 0; p < k; ++p)
        {
            double purity = (2.0 * R[s1][p] - 2.0 * R[s2][p]) / (2.0 * R[s1][p] - c2 * R[s1][p] - 2.0 * R[s2][p] + c1 * R[s2][p]);
            cscale[p] = (2.0 - 2.0 * purity + purity * c1) / R[s1][p];
            
            if(purity < 0 || purity > 1 || cscale[p] < 0)
            {
	      throw "The given copy-numbers do not allow to scale the read-depth ratios into corresponding fractional copy numbers, having RDRs of " + std::to_string(R[s1][p]) + " and " + std::to_string(R[s2][p]) + "!";
            }
        }
        
    } else {
        throw "At least one clonal copy number should be provided!";
    }
    
    log("Scaling factors:\n" + toString(cscale), VERBOSITY_t::DEBUGGING, verbose);
    
    DoubleMatrix F(m);
    for(int s = 0; s < m; ++s)
    {
        F[s] = DoubleArray(k);
        for(int p = 0; p < k; ++p)
        {
            F[s][p] = R[s][p] * cscale[p];
        }
    }
    
    return F;
}


void splitFractions(const DoubleMatrix &B, const DoubleMatrix &F, const int m, const int k, DoubleMatrix &FA, DoubleMatrix &FB)
{
    FA.clear();
    FB.clear();
    
    FA = DoubleMatrix(m);
    FB = DoubleMatrix(m);
    for(int s = 0; s < m; ++s)
    {
        FA[s] = DoubleArray(k);
        FB[s] = DoubleArray(k);
        for(int p = 0; p < k; ++p)
        {
            FB[s][p] = B[s][p] * F[s][p];
            FA[s][p] = F[s][p] - FB[s][p];
        }
    }
}


void printSolutions(const int k, const double obj, const IntMatrix& CA, const IntMatrix& CB, const DoubleMatrix& U)
{
    std::cout << "Obj : " << obj << std::endl;
    std::cout << std::endl;
    
    for(int i = 0; i < CA[0].size(); ++i)
    {
        std::cout << "Clone " << i << " : ";
        for(int s = 0; s < CA.size(); ++s)
        {
            std::cout << CA[s][i] << ":" << CB[s][i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    for(int p = 0; p < k; ++p)
    {
        std::cout << "Sample " << p << ": ";
        for(int i = 0; i < U.size(); ++i)
        {
            std::cout << U[i][p] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


void printSolutions(const int k, const int numSolutions, const DoubleArray& obj, const Int3Matrix& CA, const Int3Matrix& CB, const Double3Matrix& U)
{
    for(int sol = 0; sol < numSolutions; ++sol)
    {
        std::cout << "=== SOL: " << sol+1 << std::endl;
        printSolutions(k, obj[sol], CA[sol], CB[sol], U[sol]);
    }
}


void segmentation(const InputInstance& instance, const BBCInstance& bbc, const int n, const IntMatrix& ACN, const IntMatrix& BCN, const DoubleMatrix& U, const std::string& outprefix, const VERBOSITY_t& verbosity)
{
    std::ofstream outbbc((outprefix + ".bbc.ucn.tsv").c_str());
    outbbc << "#CHR\tSTART\tEND\tSAMPLE\tRD\t#SNPS\tCOV\tALPHA\tBETA\tBAF\tCLUSTER";
    std::ofstream outseg((outprefix + ".seg.ucn.tsv").c_str());
    outseg << "#CHR\tSTART\tEND\tSAMPLE";
    outbbc << "\tcn_normal" << "\tu_normal";
    outseg << "\tcn_normal" << "\tu_normal";
    for(int i = 1; i < n; ++i)
    {
        outbbc << "\tcn_clone" << i << "\tu_clone" << i;
        outseg << "\tcn_clone" << i << "\tu_clone" << i;
    }
    outbbc << std::endl;
    outseg << std::endl;

    for(auto const& chr : bbc.chromosomes())
    {
        size_t idx_chr = chr.second;
        
        std::vector<Bin> supp;
        for(auto const& mbin : bbc.bins()[idx_chr])
        {
            supp.emplace_back(mbin.first.first, mbin.first.second);
        }
        std::sort(supp.begin(), supp.end(), [](const Bin& l, const Bin& r) {return l.first < r.first;});

        long long start = supp[0].first;
        long long end = supp[0].second;
        IntArray previousA(n, -1);
        IntArray previousB(n, -1);
        bool first = true;
        bool sameseg = true;
        
        for(auto const& bin : supp)
        {
            size_t idx_bin = bbc.bins()[idx_chr].at(bin);
            for(auto const& sample : bbc.binsToSampleRdBaf()[idx_chr][idx_bin])
            {
                std::string cluster = bbc.binsToCluster()[idx_chr][idx_bin];
                auto const& scab = bbc.binsToSampleSCAB()[idx_chr][idx_bin].at(sample.first);
                outbbc << chr.first << '\t' << bin.first << '\t' << bin.second << '\t' << sample.first << '\t' << sample.second.first << '\t' << std::get<0>(scab) << '\t' << std::get<1>(scab) << '\t' << std::get<2>(scab) << '\t' << std::get<3>(scab) << '\t' << sample.second.second << '\t' << cluster;
                for(int i = 0; i < n; ++i)
                {
                    int s = instance.clusterToIdx().at(cluster);
                    outbbc << '\t' << ACN[s][i] << '|' << BCN[s][i];
                    outbbc << '\t' << U[i][instance.sampleToIdx().at(sample.first)];
                }
                outbbc << std::endl;
            }
            
            for(int i = 0; i < n; ++i)
            {
                std::string cluster = bbc.binsToCluster()[idx_chr][idx_bin];
                int s = instance.clusterToIdx().at(cluster);
                if(first)
                {
                    first = false;
                    start = bin.first;
                    end = bin.second;
                    previousA.assign(ACN[s].begin(), ACN[s].end());
                    previousB.assign(BCN[s].begin(), BCN[s].end());
                    break;
                } else if(ACN[s][i] != previousA[i] || BCN[s][i] != previousB[i] || bin.first != end) {
                    
                    sameseg = false;
                    for(auto const& sample : instance.samples())
                    {
                        outseg << chr.first << '\t' << start << '\t' << end << '\t' << sample;
                        for(int i = 0; i < n; ++i)
                        {
                            outseg << '\t' << previousA[i] << '|' << previousB[i];
                            outseg << '\t' << U[i][instance.sampleToIdx().at(sample)];
                        }
                        outseg << std::endl;
                    }
                    start = bin.first;
                    end = bin.second;
                    previousA.assign(ACN[s].begin(), ACN[s].end());
                    previousB.assign(BCN[s].begin(), BCN[s].end());
                    break;
                }
            }
            
            if(sameseg)
            {
                end = bin.second;
            } else {
                sameseg = true;
            }
        }
        
        for(auto const& sample : instance.samples())
        {
            outseg << chr.first << '\t' << start << '\t' << end << '\t' << sample;
            for(int i = 0; i < n; ++i)
            {
                outseg << '\t' << previousA[i] << '|' << previousB[i];
                outseg << '\t' << U[i][instance.sampleToIdx().at(sample)];
            }
            outseg << std::endl;
        }
    }
    
    outbbc.close();
    log("A BBC file with inferred copy numbers and proportion have written in: " + outprefix + ".bbc.ucn.tsv", VERBOSITY_t::ESSENTIAL, verbosity);
    outseg.close();
    log("The tumor-clone segmented genomes with corresponding copy numbers have been written in: " + outprefix + ".seg.ucn.tsv", VERBOSITY_t::ESSENTIAL, verbosity);
}


CNMap getCNMap(const ArgParse& argparse, const InputInstance& instance)
{
    CNMap cnmap;
    for(auto const& cn : argparse.cn())
    {
        std::string cluster = cn.first;
        if(instance.clusters().find(cluster) == instance.clusters().end())
        {
            throw "The cluster with ID " + cluster + " specified in the input arguments has not been found in the input seg file!";
        }
        
        int idx = instance.clusterToIdx().at(cluster);
        cnmap.emplace(idx, CNState(cn.second.first, cn.second.second));
    }
    return cnmap;
}

std::string printFractions(const DoubleMatrix &F, const InputInstance &instance)
{
    std::stringstream s;
    for(DoubleMatrix::const_iterator it = F.begin(); it != F.end(); ++it)
    {
        s << instance.idxToCluster().at(std::distance(F.begin(), it)) << " : ";
        for(DoubleArray::const_iterator itt = (*it).begin(); itt != (*it).end(); ++itt)
        {
            s << (*itt) << " ";
        }
        s << std::endl;
    }
    return s.str();
}

