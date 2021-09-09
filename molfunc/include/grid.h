#ifndef MOLFUNC_EXT_GRID_H
#define MOLFUNC_EXT_GRID_H
#include "iostream"
#include "string"
#include "array"
#include "vector"

using namespace std;


namespace molfunc{

    class GridPoint: public array<double, 3>{

        public:
            double energy = 0.0;

    };

    class Grid3D: public vector<GridPoint>{

        double min_value;
        double max_value;

        public:
            Grid3D(double min_value, double max_value, unsigned int num);

            GridPoint minimum_energy_point();
            GridPoint random_point();

    };
}



#endif //MOLFUNC_EXT_GRID_H
