#include "string"
#include <sstream>
#include <fstream>
#include <cmath>
#include "utils.h"


using namespace std;

namespace molfunc{

    vector<string> utils::split(const string &s, char delim) {
        /*********************************************************
         *  Split a string by a delimiter into a vector of strings
         *
         *  Arguments:
         *      s (string):
         *
         *      delim (char):
         *
         *  Returns:
         *      (list(string)):
         *
         *  Example:
         *      items = split("A phrase", ' ');
         *      // items = ["A", "phrase"]
         ********************************************************/
        stringstream stream(s);
        string item;
        vector<string> elems;

        while (getline(stream, item, delim)) {
            if (!item.empty()){
                elems.push_back(item);
            }
        }
        return elems;
    }

    bool utils::is_close(double a, double b, double tol){
        return abs(a - b) < tol;
    }

    bool utils::is_close(double a, double b){
        // Are two numbers close with a default tolerance
        return utils::is_close(a, b, 1E-8);
    }

    bool utils::is_close(vector<double> a, vector<double> b){

        if (a.size() != b.size()){
            throw runtime_error("Size of the two vectors not identical, "
                                "therefore not the same!");
        }

        // Ensure all elements are close
        for (int i=0; i<a.size(); i++){
            if (!is_close(a[i], b[i])){
                return false;
            }
        }

        return true;
    }

}
