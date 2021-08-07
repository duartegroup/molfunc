#ifndef MOLFUNC_ATOMS_H
#define MOLFUNC_ATOMS_H
#include "string"
#include "vector"


using namespace std;

namespace molfunc{
    class Atom {

        public:
            string symbol;
            vector<double> coord;

            void translate(vector<double> vec);

            explicit Atom(string symbol,
                          double x,
                          double y,
                          double z);
    };
}

#endif //MOLFUNC_ATOMS_H
