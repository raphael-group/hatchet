#ifndef _UTILS_H_
#define _UTILS_H_


#include <string>
#include <utility>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
#include <sstream>
#include <iterator>
#include <map>
#include <set>
#include <tuple>
#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>
#include <thread>
#include <mutex>
#include <random>
#include <ctime>
#include <limits>
#include <functional>
#include <typeinfo>

#define TOL 0.001

/// Random number generator
extern std::mt19937_64 g_rng;
extern std::mutex g_mutex;
extern std::mutex g_output_mutex;

typedef enum {
    NONE = 0,
    ESSENTIAL = 1,
    VERBOSE = 2,
    DEBUGGING = 3
} VERBOSITY_t;

typedef enum {
    BOTH = 0,
    ILP = 1,
    CD = 2
} SOLVE_t;

typedef enum {
    FULL = 0,
    CARCH = 1,
    UARCH = 2
} MODE_t;

typedef std::vector<int> IntArray;
typedef std::vector<IntArray> IntMatrix;
typedef std::vector<IntMatrix> Int3Matrix;
typedef std::vector<Int3Matrix> Int4Matrix;

typedef std::vector<double> DoubleArray;
typedef std::vector<DoubleArray> DoubleMatrix;
typedef std::vector<DoubleMatrix> Double3Matrix;
typedef std::vector<Double3Matrix> Double4Matrix;

typedef std::pair<long long, long long> Bin;
typedef std::set<std::string> StringSet;
typedef std::map<std::string, size_t> StringMap;
typedef std::vector<std::map<Bin, size_t> > BinMap;
typedef std::vector<std::vector<std::map<std::string, std::pair<double, double> > > > DoubleBBRecord;
typedef std::vector<std::vector<std::map<std::string, std::tuple<int, int, int, int> > > > IntBBRecord;
typedef std::vector<std::vector<std::string> > ClusterBBRecord;

typedef std::set<int> IntSet;
typedef std::pair<IntSet, IntSet> IntSetPair;

typedef std::pair<int, int> CNState;
typedef std::map<std::string, CNState> CNAssignCluster;
typedef std::map<int, CNState> CNMap;

int countElements(const IntMatrix arg);
int countElements(const Int3Matrix arg);
int countElements(const Int4Matrix arg);

int countElements(const DoubleMatrix arg);
int countElements(const Double3Matrix arg);
int countElements(const Double4Matrix arg);

int sum_of_elements(const IntArray arg);

std::string timestamp();

void log(const std::string& msg, const VERBOSITY_t lb, const VERBOSITY_t curr);
void warning(const std::string& msg);

bool getBit(const int num, const int pos);

void split(const std::string &s, const char d, std::vector<std::string> &out);
std::vector<std::string> split(const std::string &s, const char d);

std::string rtrim(const std::string &s);
std::string ltrim(const std::string &s);
std::string trim(const std::string &s);

IntArray sort_indexes(const IntArray &v);
DoubleMatrix transpose(const DoubleMatrix &M);

std::string toString(const IntArray &M);
std::string toString(const DoubleArray &M);
std::string toString(const IntMatrix &M);
std::string toString(const DoubleMatrix &M);

int maxCeil(const DoubleMatrix &M);

std::string base_name(std::string const& path);
std::string remove_extension(std::string const& filename);

#endif // _UTILS_H_
