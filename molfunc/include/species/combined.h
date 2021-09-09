#ifndef MOLFUNC_EXT_COMBINED_H
#define MOLFUNC_EXT_COMBINED_H
#include "molecules.h"
#include "fragments.h"
#include "angles.h"
#include "iostream"
#include "string"



namespace molfunc{

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

            CombinedMolecule();
            CombinedMolecule(CoreMolecule core,
                             vector<Fragment> fragments);

            double repulsive_energy();
            double repulsive_energy(const Fragment& fragment);

            vector<AnglePotential> angle_potentials();

            void build();

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