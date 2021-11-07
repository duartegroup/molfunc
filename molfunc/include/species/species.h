#ifndef MOLFUNC_EXT_SPECIES_H
#define MOLFUNC_EXT_SPECIES_H
#include "iostream"
#include "string"
#include "vector"
#include "atoms.h"
#include "rotation.h"


namespace molfunc{

    class Species {

        public:
            Species();

            vector<Coordinate> coordinates;
            vector<Atom> atoms;

            unsigned long n_atoms() const;
            unsigned long n_masked_atoms();
            unsigned long n_unmasked_atoms();

            unsigned long no_masked_idx(unsigned long idx);

            double distance(unsigned long i, unsigned long j);

            void translate(const Vector3D &vec);

            void rotate(const RotationMatrix &R);
            void rotate(const RotationMatrix &R, unsigned long atom_idx);
            void rotate(const RotationMatrix &R, const Coordinate &origin);

            Vector3D n_vector(unsigned long i, unsigned long j);

            vector<unsigned long> masked_atom_idxs();

            void print_xyz_file(const string& filename);
            void append_xyz_file(const string& filename);

    private:
        void print_xyz_file(const string& filename, bool append);
    };
}


#endif //MOLFUNC_EXT_SPECIES_H
