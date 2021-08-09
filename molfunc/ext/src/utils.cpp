#include "string"
#include <sstream>
#include <fstream>
#include "utils.h"


using namespace std;

namespace molfunc{

    namespace utils{
        vector<string> split(const string &s, char delim) {
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
    }
}
