#include "grid.h"


namespace molfunc{

    Grid3D::Grid3D(double min_value, double max_value, unsigned int num){
        /********************************************************
         *  Construct a 3D grid of values as a flat vector:
         *
         *  Grid3D(-1, 1, 3) --> {{-1, -1, -1}, {-1, -1, 0}, ...}
         *
         *
         *  Arguments:
         *      min_value (float):
         *
         *      max_value (float):
         *
         *      num (int): Number of points in each dimension
         ********************************************************/
         this->min_value = min_value;
         this->max_value = max_value;

         double diff = max_value - min_value;

         for (int i=0; i<num; i++){
             for (int j=0; j<num; j++){
                 for (int k=0; k<num; k++){

                     array<double, 3> arr = {min_value + diff*static_cast<double>(i)/(num-1),
                                             min_value + diff*static_cast<double>(j)/(num-1),
                                             min_value + diff*static_cast<double>(k)/(num-1)};
                     this->push_back(arr);

                 }// k
             }// j
         }// i
    }

}

