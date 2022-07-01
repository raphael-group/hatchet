#ifndef ILPSUBSET_H
#define ILPSUBSET_H

#include "utils.h"
#include "gurobi-utils.h"
#include <gurobi_c++.h>

class IlpSubset
{

public:

    IlpSubset(const int n,
              const int m,
              const int k,
              const int maxc,
              const int d,
              const double mu,
              const int base,
              const bool ampdel,
              const CNMap& cn,
              const DoubleMatrix& FA,
              const DoubleMatrix& FB,
              const DoubleArray& w,
              const VERBOSITY_t v);


    void fixC(const IntMatrix& barCA, const IntMatrix& barCB);

    void fixU(const DoubleMatrix& barU)
    {
        _barU = barU;
        _mode = MODE_t::CARCH;
        for(int p = 0; p < _k; ++p)
        {
            double sum = 0.0;
            for(int i = 0; i < _n; ++i)
            {
                sum += barU[i][p];
                assert( (0.0 - TOL <= barU[i][p] && barU[i][p] <= 0.0 + TOL) || (barU[i][p] >= _mu - TOL) || i == 0);
            }
            assert(sum >= 1.0 - TOL && sum <= 1.0 + TOL);
        }
    }

    void init();

    void hotStart(const IntMatrix& hCA, const IntMatrix& hCB);

    int solve(const int timeLimit, const int memoryLimit, const int nrThreads);
    int solve(const int timeLimit, const int memoryLimit, const int nrThreads, const int maxIter);

    void exportModel(const std::string& filename);

    DoubleArray getObjs()
    {
        return _objs;
    }
    Int3Matrix getACNs()
    {
        return _CAs;
    }
    Int3Matrix getBCNs()
    {
        return _CBs;
    }
    Double3Matrix getUs()
    {
        return _Us;
    }

    double runtime()
    {
        return _model.get(GRB_DoubleAttr_Runtime);
    }

    double gap()
    {
        return _model.get(GRB_DoubleAttr_MIPGap);
    }

    double LB()
    {
        return _model.get(GRB_DoubleAttr_LB);
    }

    double UB()
    {
        return _model.get(GRB_DoubleAttr_UB);
    }

    void printValues();

    static std::pair<IntMatrix, IntMatrix> firstHotStart(const DoubleMatrix &FA, const DoubleMatrix &FB, const int m, const int k, const int n,
                                                         const int cmax, const int d, const int base, const bool ampdel, const CNMap& cn);

private:

    /// Build main common variables
    void buildVariables();
    /// Build the variables that depend on the chosen second objective
    void buildOptionalVariables();
    /// Build main common constraints
    void buildConstraints();
    /// Build the constraints that depend on the chosen second objective
    void buildOptionalConstraints();
    /// Build the constraints used to break the symmetry by avoid permutations of non-null columns
    void buildSymmetryBreaking();
    ///
    void fixGivenCN();
    /// Build objective
    void buildObjective();
    /// Add the constraint to avoid current solution and to find the next one
    void buildNext();
    /// Set specific values of variables for testing
    void test();
    /// Extract from current solution the values of the variables corresponding to copy numbers of clones
    IntMatrix getCA();
    /// Extract from current solution the values of the variables corresponding to B-allele copy numbers of clones
    IntMatrix getCB();
    /// Extract from current solution the values of the variables corresponding to the proportions
    DoubleMatrix getU();
    ///
    bool checkRank(const IntMatrix &hCA, const IntMatrix &hCB, const IntArray &map);
    ///
    double symmCoeff(const int index);

    /// Number of clones
    const int _n;
    /// Number of total segments
    const int _m;
    /// Number of samples
    const int _k;
    /// Maximum copy number
    const int _cmax;
    /// Maximum number of distinct copy-number states
    const int _d;
    ///
    const double _mu;
    ///
    const int _base;
    ///
    const bool _ampdel;
    ///
    const CNMap& _cn;
    /// Read depth ratio r_{t, p} for each segment t in each sample p
    const DoubleMatrix& _FA;
    /// B-allele frequency d_{t, p} for each segment t in each sample p
    const DoubleMatrix& _FB;
    ///
    const DoubleArray& _w;
    ///
    const VERBOSITY_t _v;
    ///
    IntMatrix _barCA;
    ///
    IntMatrix _barCB;
    ///
    DoubleMatrix _barU;
    ///
    MODE_t _mode;
    /// Maximum B-allele copy number
    const int _M;

    /// Gurobi environment
    GRBEnv _env;
    /// Gurobi model
    GRBModel _model;
    /// Variable yA[s][p] is the error between observed and estimated fractional A-allele copy number of segment s in sample p
    VarMatrix _yA;
    /// Variable yB[s][p] is the error between observed and estimated fractional B-allele copy number of segment s in sample p
    VarMatrix _yB;
    /// Variable fA[t][p] is the fractional A-allele copy number of chosen segment t in sample p
    VarMatrix _fA;
    /// Variable fB[t][p] is the fractional B-allele copy number of chosen segment t in sample p
    VarMatrix _fB;
    /// Variable cA[t][i] is the integer A-allele copy number of the chosen segment t in clone i
    VarMatrix _cA;
    /// Variable cB[t][i] is the integer B-allele copy number of the chosen segment t in clone i
    VarMatrix _cB;
    /// Variable cA[b][t][i] is the b-th bit in the binary representation of cA[t][i]
    Var3Matrix _bitcA;
    /// Variable cB[b][t][i] is the b-th bit in the binary representation of cB[t][i]
    Var3Matrix _bitcB;
    /// Variable u[i][p] is the proportion of clone i in sample p
    VarMatrix _u;
    ///
    VarMatrix _x;
    ///
    VarArray _adA;
    ///
    VarArray _adB;
    /// Variable z[b][t][i][p] is the product between c[b][t][i] and u[i][p]
    Var4Matrix _vA;
    /// Variable zB[b][t][i][p] is the product between cB[b][t][i] and u[i][p]
    Var4Matrix _vB;
    ///
    Var3Matrix _z;
    /// Variable barcA[b][t][i][j] = 1 iff _cA[b][t][i] and _cA[b][t][j] are different
    Var4Matrix _barcIJA;
    /// Variable barcB[b][t][i][j] = 1 iff _cB[b][t][i] and _cB[b][t][j] are different
    Var4Matrix _barcIJB;

    /// Objective value for each found solution
    DoubleArray _objs;
    /// Integer A-allele copy numbers of the clones for each solution
    Int3Matrix _CAs;
    /// Integer B-allele copy numbers of the clones for each solution
    Int3Matrix _CBs;
    /// Proportions of the clones in the samples for each solution
    Double3Matrix _Us;

};

#endif // ILPSUBSET_H
