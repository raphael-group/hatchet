#ifndef _ARGPARSE_H_
#define _ARGPARSE_H_

#include "utils.h"


class ArgParse
{
public:
    
    ArgParse(int argc, char** argv);
    
    std::string help();
    std::string toString();
    
    void setClonalCluster(const CNAssignCluster& inferred);

    std::string input() const
    {
        return _input;
    }

    std::string inputSEG() const
    {
        return _inputSEG;
    }

    std::string inputBBC() const
    {
        return _inputBBC;
    }

    int n() const
    {
        return _n;
    }
    
    int d() const
    {
        return _d;
    }
    
    int cmax() const
    {
        return _e;
    }
    
    int j() const
    {
        return _j;
    }
    
    int nrSeeds() const
    {
        return _p;
    }

    double minprop() const
    {
        return _u;
    }
    
    int m() const
    {
        return _m;
    }
    
    int s() const
    {
        return _s;
    }
    
    int maxiter() const
    {
        return _i;
    }
    
    int rnd() const
    {
        return _r;
    }
    
    CNAssignCluster cn() const
    {
        return _cn;
    }
    
    CNAssignCluster scalingCN() const
    {
        return _scalingcn;
    }
    
    SOLVE_t M() const
    {
        return _M;
    }
    
    VERBOSITY_t v() const
    {
        return _v;
    }
    
    const std::string o() const
    {
        return _o;
    }
    
    double diploidThreshold() const
    {
        return _t;
    }
    
    int base() const
    {
        return _base;
    }
    
    double forceAMPDEL() const
    {
        return _f;
    }
    
private:
    
    /// Number of arguments
    int _argc;
    /// Arguments
    char** _argv;

    /// Input .seg file
    std::string _input;
    /// Input .seg file
    std::string _inputSEG;
    /// Input .seg file
    std::string _inputBBC;
    /// Clonal segments
    CNAssignCluster _cn;
    //
    CNAssignCluster _scalingcn;
    /// Help option
    bool _h;
    /// Number of distinct clones
    int _n;
    /// Maximum number of distinct copy-number states
    int _d;
    /// Maximum copy number
    int _e;
    /// Number of parallel threads for Gurobi
    int _j;
    /// Number of seeds
    int _p;
    ///
    double _u;
    /// Maximum resident memory
    int _m;
    /// Time limit
    int _s;
    ///
    int _i;
    ///
    int _r;
    ///
    SOLVE_t _M;
    /// Verbosity of the logging
    VERBOSITY_t _v;
    ///
    std::string _o;
    ///
    double _t;
    ///
    int _base;
    ///
    bool _f;
    
    /// Set of total arguments
    std::set<std::string> _args;
    /// Map of boolean to record the presence of required arguments
    std::map<std::string, bool> _required;
    /// Map of boolean to record the present of optional arguments
    std::map<std::string, bool> _optional;
    
    void parse();
    void parseOptional(const std::string &arg, const std::string &value);
    void parseClonalCopyNumbers(const std::string &value);
    bool checkInt(const std::string &s);
};




#endif // _UTILS_GUROBI_H_
