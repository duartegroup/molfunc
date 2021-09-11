#include "grid.h"
#include "iostream"


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

         for (unsigned int i=0; i<num; i++){
             for (unsigned int j=0; j<num; j++){
                 for (unsigned int k=0; k<num; k++){

                     GridPoint arr = {min_value + diff*static_cast<double>(i)/(num-1),
                                      min_value + diff*static_cast<double>(j)/(num-1),
                                      min_value + diff*static_cast<double>(k)/(num-1)};
                     this->push_back(arr);

                 }// k
             }// j
         }// i
    }

    GridPoint Grid3D::minimum_energy_point(){
        /*****************************************
         *  Get the grid point at which the
         *  energy is a minimum
         *
         *  Returns:
         *      (GridPoint):
         ***************************************/

        if (this->empty()){
            throw out_of_range("Cannot find the minimum in a "
                               "empty grid");
        }

        auto min_e_point = this->begin();


        for (auto it = this->begin(); it != this->end(); ++it) {
            if (it->energy < min_e_point->energy) {
                min_e_point = it;
            }
        }

        return *min_e_point;
    }

    GridPoint Grid3D::random_point(){
        /*********************************
         * Generate a random point within
         * the rotation grid
         ********************************/

        return this->data()[random() % this->size()];
    }
}

