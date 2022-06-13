#include "utils.h"

std::mt19937_64 g_rng;
std::mutex g_mutex;
std::mutex g_output_mutex;

int countElements(const IntMatrix arg)
{
    int res = 0;
    for(IntMatrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += (*it).size();
    }
    return res;
}

int countElements(const Int3Matrix arg)
{
    int res = 0;
    for(Int3Matrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int countElements(const Int4Matrix arg)
{
    int res = 0;
    for(Int4Matrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int countElements(const DoubleMatrix arg)
{
    int res = 0;
    for(DoubleMatrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += (*it).size();
    }
    return res;
}

int countElements(const Double3Matrix arg)
{
    int res = 0;
    for(Double3Matrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int countElements(const Double4Matrix arg)
{
    int res = 0;
    for(Double4Matrix::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        res += countElements(*it);
    }
    return res;
}

int sum_of_elements(const IntArray arg)
{
    int sum = 0;
    for(IntArray::const_iterator it = arg.begin();
        it != arg.end();
        ++it)
    {
        sum += *it;
    }
    return sum;
}

std::string timestamp()
{
  time_t now = time(NULL);
  struct tm tstruct;
  char buf[40];
  tstruct = *localtime(&now);
  //format: HH:mm:ss
  strftime(buf, sizeof(buf), "[%X]", &tstruct);
  return buf;
    // std::stringstream ss;
    // std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    // std::time_t now_c = std::chrono::system_clock::to_time_t(now - std::chrono::hours(24));
    // ss << std::put_time(std::localtime(&now_c), "[%F %T]");
    // return ss.str();
    //std::time_t result = std::time(nullptr);
    //return std::asctime(std::localtime(&result));
    /**
     time_facet *facet = new time_facet("%d-%b-%Y %H:%M:%S");
     out.imbue(locale(out.getloc(), facet));
     out << "[" << second_clock::local_time() << "]";
     **/
}

void log(const std::string& msg, const VERBOSITY_t lb, const VERBOSITY_t curr)
{
    if(curr >= lb) {std::cerr << ((lb == VERBOSITY_t::ESSENTIAL) ? "\033[95m\033[1m" + timestamp() + "### " :
                                  (lb == VERBOSITY_t::VERBOSE) ? "\033[92m"+ timestamp() + "## " :
                                  "\033[96m# ") << msg << "\t\033[0m" << std::endl;}
}

void warning(const std::string& msg)
{
   std::cerr << "\033[93m\033[1m" << timestamp() << "### " << msg << "\t\033[0m" << std::endl;
}


bool getBit(const int num, const int pos)
{
    return (num >> pos) & 1;
}


void split(const std::string &s, const char d, std::vector<std::string> &out)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, d)) {
        if(!item.empty())
        {
            out.push_back(item);
        }
    }
}


std::vector<std::string> split(const std::string &s, const char d)
{
    std::vector<std::string> tokens;
    split(s, d, tokens);
    return tokens;
}


std::string rtrim(const std::string &s)
{
    return std::string(s).erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}


std::string ltrim(const std::string &s)
{
    return std::string(s).erase(0, s.find_first_not_of(" \t\n\r\f\v"));
}


std::string trim(const std::string &s)
{
    std::string res(s);
    res.erase(std::remove(res.begin(), res.end(), ' '), res.end());
    res.erase(std::remove(res.begin(), res.end(), '\t'), res.end());
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    return res;
}

IntArray sort_indexes(const IntArray &v) {

    // initialize original index locations
    IntArray idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

DoubleMatrix transpose(const DoubleMatrix &M)
{
    DoubleMatrix T(M[0].size(), DoubleArray(M.size(), 0.0));
    for(int i = 0; i < M.size(); ++i)
    {
        for(int j = 0; j < M[i].size(); ++j)
        {
            T[j][i] = M[i][j];
        }
    }
    return T;
}

std::string toString(const IntArray &M)
{
    std::stringstream s;
    for(IntArray::const_iterator it = M.begin(); it != M.end(); ++it)
    {
        s << (*it) << " ";
    }
    s << std::endl;
    return s.str();
}

std::string toString(const DoubleArray &M)
{
    std::stringstream s;
    for(DoubleArray::const_iterator it = M.begin(); it != M.end(); ++it)
    {
        s << (*it) << " ";
    }
    s << std::endl;
    return s.str();
}

std::string toString(const IntMatrix &M)
{
    std::stringstream s;
    for(IntMatrix::const_iterator it = M.begin(); it != M.end(); ++it)
    {
        for(IntArray::const_iterator itt = (*it).begin(); itt != (*it).end(); ++itt)
        {
            s << (*itt) << " ";
        }
        s << std::endl;
    }
    return s.str();
}

std::string toString(const DoubleMatrix &M)
{
    std::stringstream s;
    for(DoubleMatrix::const_iterator it = M.begin(); it != M.end(); ++it)
    {
        for(DoubleArray::const_iterator itt = (*it).begin(); itt != (*it).end(); ++itt)
        {
            s << (*itt) << " ";
        }
        s << std::endl;
    }
    return s.str();
}

int maxCeil(const DoubleMatrix &M)
{
    int max = (int)M[0][0];
    for(DoubleMatrix::const_iterator it = M.begin(); it != M.end(); ++it)
    {
        for(DoubleArray::const_iterator itt = (*it).begin(); itt != (*it).end(); ++itt)
        {
            max = std::max((int)std::ceil(*itt), max);
        }
    }
    return max;
}

std::string base_name(std::string const& path)
{
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string remove_extension(std::string const& filename)
{
    typename std::string::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != std::string::npos ? filename.substr(0, p) : filename;
}
