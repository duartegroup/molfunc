#include "atoms.h"
#include <array>
#include <vector>
#include <stdexcept>


using namespace std;


namespace molfunc{

    // Default constructors are required for Cython wrapping
    Atom::Atom() = default;

    Atom::Atom(const string& symbol) {
        /*********************************************************
         * Construct an atom from an atomic symbol and its position
         * in 3D space
         *
         * Arguments:
         *      symbol (string): Atomic symbol
         ********************************************************/

        this->symbol = symbol;

        if (symbol == "R") masked = true;

    }

    unsigned int Atom::atomic_number() const{
        /*********************************************************
         * Determine the atomic number of this element by its
         * symbol.
         *
         * Returns:
         *      (int): Atomic number [1, 118]
         *
         * Raises:
         *      (runtime_error): If the symbol cannot be found in
         *                    the periodic table of known elements
         ********************************************************/

        // Defines a special dummy atom ("R") with atomic number 0
        vector<string> elements = {"R", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
                                   "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
                                   "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
                                   "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                                   "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                                   "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
                                   "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
                                   "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                   "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
                                   "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
                                   "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

        for (unsigned int i=0; i<elements.size(); i++){
            if (elements[i] == symbol){
                return i;
            }
        }

        throw runtime_error("Failed to find atomic number for "+symbol);
    }

    bool Atom::is_dummy() const{return (symbol == "R");}

    double Atom::covalent_radius() const{
        /*********************************************************
         * Determine the covalent radius of this atom
         *
         * Returns:
         *      (float): Radius in Å
         *
         * Raises:
         *      (out_of_range): If the symbol cannot be found in
         *                      the periodic table of known radii
         ********************************************************/
        //                       R     H     He    B   ...
        vector<double> radii = {0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11,
                                1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.50, 1.24,
                                1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54,
                                1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15,
                                2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.95, 1.94, 1.92, 1.92, 1.89, 1.90,
                                1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48,
                                1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69};
        //                                                                                                ^ Cm

        if (atomic_number() > radii.size()){
            throw out_of_range("Unknown element "+symbol+
                               " Not in known covalent radii");
        }

        return radii[atomic_number()];
    }

    double Atom::phi0(unsigned long n_neighbours) const {
        /******************************************
         * Determine an approximate bond angle
         * between two atoms bonded to this one
         * given a number of neighbours
         *
         * Returns:
         *      (float): Optimal angle in degrees
         *****************************************/
        if (n_neighbours == 2){
            if (symbol == "O") return 110.0;
            return 180.0;
        }

        if (n_neighbours == 3){
            if (symbol == "N" || symbol == "P") return 110.0;
            return 120.0;
        }

        if (n_neighbours == 4) return 109.5;

        if (n_neighbours == 6) return 90.0;

        return 110.0;
    }

    Atom3D::Atom3D(const string& symbol,  double x, double y, double z){
        /*********************************************************
         * Construct an atom from an atomic symbol and its position
         * in 3D space
         *
         * Arguments:
         *      symbol (string): Atomic symbol
         *      x (float):       Cartesian x coordinate (Å)
         *      y (float):                 y coordinate (Å)
         *      z (float):                 z coordinate (Å)
         ********************************************************/
        this->symbol = symbol;
        this->coord = {x, y, z};
    }

    Atom3D::Atom3D(const string& symbol,  Coordinate &coord){
        /*********************************************************
         * Construct an atom from an atomic symbol and its position
         * in 3D space
         *
         * Arguments:
         *      symbol (string): Atomic symbol
         *      coord (float):       Cartesian x, y, z coordinates (Å)
         ********************************************************/
        this->symbol = symbol;
        this->coord = coord;
    }

    double Atom3D::x() {return coord[0];}

    double Atom3D::y() {return coord[1];}

    double Atom3D::z() {return coord[2];}

    double Coordinate::x() {return this->data()[0];}

    double Coordinate::y() {return this->data()[1];}

    double Coordinate::z() {return this->data()[2];}

    std::ostream &operator<<(ostream &os, const Coordinate &coordinate) {
        os << to_string(coordinate[0]) + " \t"
              + to_string(coordinate[1]) + " \t"
              + to_string(coordinate[2]) + " \t";

        return os;
    }

    Vector3D Coordinate::operator+(const Coordinate& other) {
        return {this->data()[0] + other[0],
                this->data()[1] + other[1],
                this->data()[2] + other[2]};
    }

    Vector3D Coordinate::operator-(const Coordinate& other) {
        return {this->data()[0] - other[0],
                this->data()[1] - other[1],
                this->data()[2] - other[2]};
    }

    Coordinate Coordinate::operator-=(const Coordinate& other) {
        this->data()[0] -= other[0],
        this->data()[1] -= other[1],
        this->data()[2] -= other[2];
        return *this;
    }

    Coordinate Coordinate::operator+=(const Coordinate& other) {
        this->data()[0] += other[0],
        this->data()[1] += other[1],
        this->data()[2] += other[2];
        return *this;
    }
}
