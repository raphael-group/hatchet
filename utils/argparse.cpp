#include "argparse.h"


ArgParse::ArgParse(int argc, char** argv)
    : _argc(argc)
    , _argv(argv)
    , _input()
    , _inputSEG()
    , _inputBBC()
    , _cn()
    , _scalingcn()
    , _h(false)
    , _n()
    , _d(-1)
    , _e(-1)
    , _j(std::max((int)std::thread::hardware_concurrency(), 1))
    , _p(400)
    , _u(0.01)
    , _m(-1)
    , _s(-1)
    , _i(10)
    , _r(-1)
    , _M(SOLVE_t::BOTH)
    , _v(VERBOSITY_t::ESSENTIAL)
    , _o()
    , _t(0.1)
    , _base(1)
    , _f(false)
    , _args()
    , _required()
    , _optional()
{
    _args.insert("-h");
    _args.insert("-n");
    _args.insert("-d");
    _args.insert("-e");
    _args.insert("-j");
    _args.insert("-p");
    _args.insert("-u");
    _args.insert("-m");
    _args.insert("-s");
    _args.insert("-r");
    _args.insert("-i");
    _args.insert("-c");
    _args.insert("-M");
    _args.insert("-v");
    _args.insert("-o");
    _args.insert("-t");
    _args.insert("-f");

    _required.insert(std::make_pair("input", false));
    _required.insert(std::make_pair("-n", false));

    _optional["-h"] = false;
    _optional["-d"] = false;
    _optional["-e"] = false;
    _optional["-j"] = false;
    _optional["-p"] = false;
    _optional["-u"] = false;
    _optional["-m"] = false;
    _optional["-s"] = false;
    _optional["-r"] = false;
    _optional["-i"] = false;
    _optional["-c"] = false;
    _optional["-M"] = false;
    _optional["-v"] = false;
    _optional["-o"] = false;
    _optional["-t"] = false;
    _optional["-f"] = false;

    parse();
    
    if(!_optional.at("-o"))
    {
        _o = base_name(_input);
    }
    
    if(2 * _base > _e && _e != -1)
    {
        throw "Maximum copy number cannot be lower than " + std::to_string(2 * _base) + "!";
    }
}


std::string ArgParse::help()
{
    std::stringstream res;
    res << "\033[96m" << "\033[1m" << "USAGE:" << "\033[0m" << "\033[96m" << std::endl;
    res << "\t./solve input -n INT [-e INT] [-j INT] [-m INT] [-s INT] [-c COPYNUMBERS] [-v] [-h]" << std::endl;
    res << "\033[96m" << "\033[1m" << "WHERE:" << "\033[0m" << "\033[96m" << std::endl;
    res << "\tinput PREFIX-INPUT\n\t\tCommon prefix for the input files .seg and .bbc (required)" << std::endl;
    res << "\t-n INT-NUMCLONES\n\t\tNumber of tumor clones to infer (required)" << std::endl;
    res << "\t-c COPYNUBERS 'ID-CLUSTER:INT-COPYNUMBER [, ID-CLUSTER:INT-COPYNUMBER]'\n\t\tClonal copy numbers for one or two clusters in the input. When the copy number for only one cluster is provided, its copy number must be equal to 2 as it corresponds to the diploid cluster.\n\t\tThe format of copy numbers is 'ID-CLUSTER:INT-COPYNUMBER [, ID-CLUSTER:INT-COPYNUMBER]' (default: diploid cluster automatically infferred)" << std::endl;
    res << "\t-d INT-MAXCNSTATES\n\t\tMaximum number of copy number states (default: 3)" << std::endl;
    res << "\t-e INT-MAXCOPYNUMBER\n\t\tMaximum integer copy number (default: inferred from fractions)" << std::endl;
    res << "\t-j INT-JOBS\n\t\tNumber of parallel jobs (default: maximum number of available processors)" << std::endl;
    res << "\t-p SEEDS\n\t\tNumber of seeds for coordinate descent (default: 400)" << std::endl;
    res << "\t-u MIN-PROP\n\t\tMinimum clone proporion in each sample (default: 0.01)" << std::endl;
    res << "\t-m INT-MAXMEMORY\n\t\tMaximum resident memory used by Gurobi expressed in Mbytes (default: automatically inferred from the maximum of the machine)" << std::endl;
    res << "\t-s INT-TIMELIMIT\n\t\tTime limit in seconds for the ilp solving of Gurobi, after the time limit the best found solution is returned (default: no timelimit)" << std::endl;
    res << "\t-i INT-ITERATIONS\n\t\tMaximum number of iterations composed of C-step/U-step for each seed (default: 10)" << std::endl;
    res << "\t-r INT-RANDOMSEED\n\t\tRandom seed for generator (default: none)" << std::endl;
    res << "\t-M INT-SOLVINGMODE\n\t\tSolving mode among: Coordinate Descent + exact ILP (0), exact ILP only (1), and Coordinate-descent only (2) (default: 0)" << std::endl;
    res << "\t-v VERBOSITY\n\t\tLevel of verbosity among: none (0), essential (1), verbose (2), and debug (3) (default: 1)" << std::endl;
    res << "\t-o PREFIX-OUTPUT\n\t\tPrefix of the output files that will be generated (default: current directory with the basename of the input prefix)"<< std::endl;
    res << "\t-t DIPLOID\n\t\tMaximum BAF shift for diploid cluster used to automatically infer the diploid cluster (default: 0.1)" << std::endl;
    res << "\t-f FORCE-AMPDEL\n\t\tForce every mutated segment (not 2 in diploid tumors and not 4 in tetraploid tumors) to be either amplified or deleted in all the clones (default: not applied)" << "\033[0m" << std::endl;
    return res.str();
}


std::string ArgParse::toString()
{
    std::stringstream res;
    res << "\tInput prefix:  " << _input << std::endl;
    res << "\tInput SEG:  " << _inputSEG << std::endl;
    res << "\tInput BBC:  " << _inputBBC << std::endl;
    res << "\tNumber of clones:  " << _n << std::endl;
    res << "\tClonal copy numbers:  ";
    for(auto const& state : _cn)
    {
        res << "{ " << state.first << " [Cluster] : " << state.second.first << "|" << state.second.second << " [CN] } ";
    }
    res << std::endl;
    res << "\tHelp message:  " << _h << std::endl;
    res << "\tMaximum number of copy-number states:  " << _d << std::endl;
    res << "\tMaximum integer copy number:  " << _e << std::endl;
    res << "\tNumber of jobs:  " << _j << std::endl;
    res << "\tNumber of seeds:  " << _p << std::endl;
    res << "\tMinimum tumor-clone threshold:  " << _u << std::endl;
    res << "\tMaximum resident memory:  " << _m << std::endl;
    res << "\tTime limit:  " << _s << std::endl;
    res << "\tMaximum number of iteratios:  " << _i << std::endl;
    res << "\tRandom seed:  " << _r << std::endl;
    res << "\tSolving mode:  " << (_M == SOLVE_t::BOTH ? "Coordinate-descent + exact ILP" : _M == SOLVE_t::ILP ? "exact ILP only" : "Coordinate-descent only") << std::endl;
    res << "\tVerbose:  " << _v << std::endl;
    res << "\tOutput prefix:  " << _o << std::endl;
    res << "\tDiploid threshold:  " << _t << std::endl;
    res << "\tBase:  " << _base << std::endl;
    res << "\tForce amp-del:  " << _f;
    return res.str();
}


void ArgParse::setClonalCluster(const CNAssignCluster& inferred)
{
    _cn.clear();
    _scalingcn.clear();
    for(auto const& pair : inferred)
    {
        _cn.emplace(pair.first, pair.second);
	_scalingcn.emplace(pair.first, pair.second);
    }
}


void ArgParse::parse()
{
    for(int i = 0; i < _argc; ++i)                                                                                                                                                                                                                       
    {                                                                                                                                                                                                                                                    
        std::cerr << "Given parsing: " << _argv[i] << std::endl;                                                                                                                                                                                     
    }
    
    for(int idx = 1; idx < _argc; ++idx)
    {
        std::string arg = std::string(_argv[idx]);
        std::cerr << "Reading: " << arg << std::endl;
        if(_args.find(arg) != _args.end())
        {
            std::cerr << "Argument: " << arg << std::endl;
            if (arg == "-h") {
                parseOptional(arg, "");
            }else if (arg == "-f") {
                parseOptional(arg, "");
            } else if (++idx < _argc) {
                parseOptional(arg, std::string(_argv[idx]));
            } else {
                throw "Missing argument values!";
            }
        } else {
            std::cerr << "Input: " << arg << std::endl;
            if((*_required.find("input")).second)
            {
                throw "Multiple default arguments found! Please specify only one default argument as input file!";
            } else {
                _input = arg;
                (*_required.find("input")).second = true;
                
                std::string seg = arg + ".seg";
                std::ifstream inputSEG(seg.c_str());
                if(!inputSEG.good())
                {
                    throw "The specified input file \"" + seg + "\" does not exists or cannot be opened!";
                } else {
                    _inputSEG = seg;
                }
                inputSEG.close();
                
                std::string bbc = arg + ".bbc";
                std::ifstream inputBBC(bbc.c_str());
                if(!inputBBC.good())
                {
                    throw "The specified input file \"" + bbc + "\" does not exists or cannot be opened!";
                } else {
                    _inputBBC = bbc;
                }
                inputBBC.close();
            }
        }
    }
    
    if(_h)
    {
        std::cout << help();
        std::exit(EXIT_SUCCESS);
    } else {
        for(std::map<std::string, bool>::const_iterator it = _required.begin();
            it != _required.end();
            ++it)
        {
            if(!(*it).second)
            {
                throw "The required input argument \"" + (*it).first + "\" is missing!";
            }
        }
    }
}


void ArgParse::parseOptional(const std::string &arg, const std::string &value)
{
    std::cerr << "Parsing: " << arg << " with value: " << value << std::endl;
    if (arg == "-n")
    {
        std::cerr << "Already present: " << (*_required.find(arg)).second  << std::endl;
    } else {
        std::cerr << "Already present: " << (*_optional.find(arg)).second  << std::endl;
    }
    if(arg == "-h")
    {
        _h =  true;
        if((*_optional.find("-h")).second)
        {
            throw "Multiple -h argument found!";
        } else {
            (*_optional.find("-h")).second = true;
        }
    } else if (arg == "-n") {
        if((*_required.find("-n")).second)
        {
            throw "Multiple -n arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The number of clone argument -n must be a positive integer!";
            } else {
                (*_required.find("-n")).second = true;
                _n = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-d") {
        if((*_optional.find("-d")).second)
        {
            throw "Multiple -d arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The maximum-cn-states argument -d must be a positive integer!";
            } else {
                (*_optional.find("-d")).second = true;
                _d = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-e") {
        if((*_optional.find("-e")).second)
        {
            throw "Multiple -e arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The maximum-copy-number argument -e must be a positive integer!";
            } else {
                (*_optional.find("-e")).second = true;
                _e = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-j") {
        if((*_optional.find("-j")).second)
        {
            throw "Multiple -j arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The number of thread argument -j must be a positive integer!";
            } else {
                (*_optional.find("-j")).second = true;
                _j = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-p") {
        if((*_optional.find("-p")).second)
        {
            throw "Multiple -p arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The seed argument -p must be a positive integer!";
            } else {
                (*_optional.find("-p")).second = true;
                _p = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-u") {
        if((*_optional.find("-u")).second)
        {
            throw "Multiple -u arguments found!";
        } else {
            double cast = std::stod(value);
            if(cast < 0.0 || cast > 0.3)
            {
                throw "The minimum tumor-clone proportion argument -u must be in [0, 0.3]!";
            } else {
                (*_optional.find("-u")).second = true;
                _u = cast;
            }
        }
    } else if (arg == "-m") {
        if((*_optional.find("-m")).second)
        {
            throw "Multiple -m arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The number of thread argument -m must be a positive integer!";
            } else {
                (*_optional.find("-m")).second = true;
                _m = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-s") {
        if((*_optional.find("-s")).second)
        {
            throw "Multiple -s arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The number of thread argument -s must be a positive integer!";
            } else {
                (*_optional.find("-s")).second = true;
                _s = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-r") {
        if((*_optional.find("-r")).second)
        {
            throw "Multiple -r arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The random-seed argument -r must be a positive integer!";
            } else {
                (*_optional.find("-r")).second = true;
                _r = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-i") {
        if((*_optional.find("-i")).second)
        {
            throw "Multiple -i arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The number of max-iteration argument -i must be a positive integer!";
            } else {
                (*_optional.find("-i")).second = true;
                _i = std::atoi(value.c_str());
            }
        }
    } else if (arg == "-f") {
        if((*_optional.find("-f")).second)
        {
            throw "Multiple -f arguments found!";
        } else {
            (*_optional.find("-f")).second = true;
            _f = true;
        }
    } else if (arg == "-M") {
        if((*_optional.find("-M")).second)
        {
            throw "Multiple -M arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The solving-mode argument -M must be a positive integer!";
            } else {
                (*_optional.find("-M")).second = true;
                int M = std::atoi(value.c_str());
                switch (M) {
                    case 0: _M = SOLVE_t::BOTH;
                        break;
                    case 1: _M = SOLVE_t::ILP;
                        break;
                    case 2: _M = SOLVE_t::CD;
                        break;
                    default:
                        throw "The solving mode must be a value among {0, 1, 2}";
                        break;
                }
            }
        }
    } else if (arg == "-v") {
        if((*_optional.find("-v")).second)
        {
            throw "Multiple -v arguments found!";
        } else {
            if(!checkInt(value))
            {
                throw "The verbosity argument -v must be a positive integer!";
            } else {
                (*_optional.find("-v")).second = true;
                int v = std::atoi(value.c_str());
                switch (v) {
                    case 0: _v = VERBOSITY_t::NONE;
                        break;
                    case 1: _v = VERBOSITY_t::ESSENTIAL;
                        break;
                    case 2: _v = VERBOSITY_t::VERBOSE;
                        break;
                    case 3: _v = VERBOSITY_t::DEBUGGING;
                        break;
                    default:
                        throw "The verbosity level must be a value among {0, 1, 2, 3}";
                        break;
                }
            }
        }
    } else if (arg == "-c") {
        if((*_optional.find("-c")).second)
        {
            throw "Multiple -c arguments found!";
        } else {
            _optional.find("-c")->second = true;
            parseClonalCopyNumbers(value);
        }
    } else if (arg == "-o") {
        if((*_optional.find("-o")).second)
        {
            throw "Multiple -o arguments found!";
        } else {
            _optional.find("-o")->second = true;
            _o = value;
        }
    } else {
        throw "Unknown flag found!";
    }
}


void ArgParse::parseClonalCopyNumbers(const std::string &value)
{
    std::vector<std::string> tokens = split(trim(value), ',');
    if (tokens.size() == 1)
    {
        _base = 1;
        
        std::vector<std::string> parsed = split(tokens[0], ':');
        if(parsed.size() != 3 || !checkInt(parsed[1]) || !checkInt(parsed[2]))
        {
            throw "The specified clonal copy number \"" + value + "\" has wrong format!";
        } else {
            int A = std::stoi(parsed[1]);
            int B = std::stoi(parsed[2]);
            
            if (A + B != 2) {
                throw "When only one clonal copy number is given, it must be for a diploid cluster with copy number equal to 2 in the format \"ID:1:1\" or \"ID:2:0\"!";
            } else {
                assert(_cn.count(parsed[0]) == 0);
                _cn.emplace(parsed[0], CNState(A, B));
                _scalingcn.emplace(parsed[0], CNState(A, B));
            }
        }
    } else if (tokens.size() >= 2) {
        _base = 2;
        
        for(auto const& token : tokens)
        {
            std::vector<std::string> parsed = split(token, ':');
            if(parsed.size() != 3 || !checkInt(parsed[1]) || !checkInt(parsed[2]))
            {
                throw "The specified clonal copy number \"" + value + "\" has wrong format!";
            } else {
                int A = std::stoi(parsed[1]);
                int B = std::stoi(parsed[2]);
                
                if (_cn.count(parsed[0]) == 0)
                {
                    _cn.emplace(parsed[0], CNState(A, B));
                } else {
                    throw "In WGD mode, the first two clonal cluster must refer to two different clusters!";
                }
                
                if(_cn.size() == 1)
                {
                    _scalingcn.emplace(parsed[0], CNState(A, B));
                } else if (_cn.size() == 2)
                {
                    _scalingcn.emplace(parsed[0], CNState(A, B));
                }
            }
        }
        
        auto it = _scalingcn.begin();
        CNState state1 = it->second;
        std::advance(it, 1);
        CNState state2 = it->second;
        
        if((state1.first + state1.second) == (state2.first + state2.second))
        {
            throw "When two clonal copy numbers are given, they first two must be different in the two segmental clusters!";
        }
        
        for(auto const& ccluster : _cn)
        {
            if(ccluster.second.first + ccluster.second.second == 2)
            {
                warning("When a clonal cluster has copy number equal to 2, you can speficy only that!");
                break;
            }
        }
        
    } else {
        throw "Please specify only two clonal copy numbers for two clusters!";
    }
}


bool ArgParse::checkInt(const std::string &s)
{
    return s.find_first_not_of( "0123456789" ) == std::string::npos;
}

