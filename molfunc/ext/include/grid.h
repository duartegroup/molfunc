#ifndef MOLFUNC_EXT_GRID_H
#define MOLFUNC_EXT_GRID_H
#include "array"
#include "vector"

using namespace std;


namespace molfunc{

    class Grid3D: public vector<array<double, 3>>{

        double min_value;
        double max_value;

        public:
            Grid3D(double min_value, double max_value, unsigned int num);

    };
}



#endif //MOLFUNC_EXT_GRID_H
