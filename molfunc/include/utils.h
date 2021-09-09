#ifndef MOLFUNC_EXT_UTILS_H
#define MOLFUNC_EXT_UTILS_H
#include "iostream"
#include "string"
#include "vector"
#include "atoms.h"


using namespace std;

namespace molfunc::utils{

        vector<string> split(const string &s, char delim);

        bool ends_with(const string &s, const string &ending);

        bool is_close(double a, double b);
        bool is_close(double a, double b, double tol);
        bool is_close(Coordinate a, Coordinate b);

        string to_lower(string s);
}


#endif //MOLFUNC_EXT_UTILS_H
