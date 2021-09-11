#ifndef MOLFUNC_ATOMS_H
#define MOLFUNC_ATOMS_H
#include "string"
#include "array"
#include "memory"
#include "iostream"
#include "vector3d.h"


using namespace std;

namespace molfunc{

    class Coordinate: public array<double, 3>{

        public:

            friend std::ostream& operator<< ( std::ostream& os,
                                             const Coordinate& c );

            Vector3D operator- (const Coordinate& other);
            Coordinate operator-= (const Coordinate& other);
            Vector3D operator+ (const Coordinate& other);
            Coordinate operator+= (const Coordinate& other);

            double x();
            double y();
            double z();
    };


    class Atom {

        public:
            string symbol = "X";
            bool masked = false;

            unsigned int atomic_number() const;
            bool is_dummy() const;

            explicit Atom();
            explicit Atom(const string& symbol);

            double covalent_radius() const;

            double phi0(unsigned long n_neighbours) const;
    };

    class Atom3D: public Atom{

        private:
            Coordinate coord = {0.0, 0.0, 0.0};

        public:
            Atom3D(const string& symbol, double x, double y, double z);
            Atom3D(const string& symbol,  Coordinate &coord);

            double x();
            double y();
            double z();
    };

}

#endif //MOLFUNC_ATOMS_H
