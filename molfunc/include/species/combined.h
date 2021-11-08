#ifndef MOLFUNC_EXT_COMBINED_H
#define MOLFUNC_EXT_COMBINED_H
#include "molecules.h"
#include "fragments.h"
#include "angles.h"
#include "iostream"
#include "string"


namespace molfunc{

    void print_all_combined_molecules(const string& xyz_filename,
                                      const string& core_xyz_filename,
                                      const vector<unsigned int>& atom_idxs_to_del);

    void print_combined_molecule_from_names(const string& xyz_filename,
                                            const string& core_xyz_filename,
                                            const vector<unsigned int>& atom_idxs_to_del,
                                            const vector<string>& frag_names);

    void print_combined_molecule_from_xyz_filenames(const string& xyz_filename,
                                                    const string& core_xyz_filename,
                                                    const vector<unsigned int>& atom_idxs_to_del,
                                                    const vector<string>& frag_xyz_filenames);
    class CombinedMolecule {

        public:

            CoreMolecule core;
            vector<Fragment> fragments;

            vector<vector<unsigned long>> fragments_atom_idxs;
            vector<unsigned long> fragment_origin_idxs;

            AnglePotentials angle_potentials;

            CombinedMolecule();
            CombinedMolecule(CoreMolecule core,
                             vector<Fragment> fragments);

            double repulsive_energy();
            static double repulsive_energy(const vector<Coordinate> &coords);
            double repulsive_energy(const Fragment& fragment);

            double total_energy(vector<Coordinate> &coords);

            void rotate_fragment(int fragment_idx,
                                 RotationMatrix &R,
                                 vector<Coordinate> &coords);

            double dE_dw(int axis_idx,
                         vector<Coordinate> &coords,
                         int fragment_idx,
                         double curr_energy);

            void gen_angle_potentials();

            void gen_fragment_idxs();

            void build();

            vector<Coordinate> coordinates();
            void set_coordinates(vector<Coordinate> &coords);

            Molecule to_molecule();


        protected:
            void rotate_fragments_global();
            void exclude_rotational_space(Fragment &fragment,
                                          double threshold);

            void translate_fragment(Fragment &fragment,
                                    unsigned long Ra_idx);

            void minimise_total_energy();
    };

}


#endif //MOLFUNC_EXT_COMBINED_H