#include "string"
#include <sstream>
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

    bool utils::ends_with(const string &s,
                   const string &ending){
        /*********************************************************
         *  Does a string end with another string?
         *  Shamelessly stolen from:
         *  https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
         *
         *  Arguments:
         *      s (string):
         *
         *      ending (string):
         *
         *  Returns:
         *      (bool):
         *
         *  Example:
         *      ends_with("filename.xyz", ".xyz") // -> true
         ********************************************************/

        if (s.length() < ending.length()){
            return false;
        }

        return (0 == s.compare(s.length() - ending.length(), ending.length(), ending));
    }

    bool utils::is_close(double a, double b, double tol){
        return abs(a - b) < tol;
    }

    bool utils::is_close(double a, double b){
        // Are two numbers close with a default tolerance
        return utils::is_close(a, b, 1E-8);
    }

    bool utils::is_close(Coordinate a, Coordinate b) {

        for (int i=0; i<3; i++){
            if (!is_close(a[i], b[i])) return false;
        }

        return true;
    }

    string utils::to_lower(string s){
        // In-place lower case a string
        for (char &c : s){
            c = tolower(c);
        }

        return s;
    }

}
