#include <iostream>
#include <algorithm>
#include "fragments.h"
#include "utils.h"


using namespace std;


namespace molfunc{

    Fragment::Fragment() = default;

    Fragment::Fragment(const string& xyz_filename):Molecule(xyz_filename) {
        /*********************************************************
         * Construct a Fragment molecule from a standard .xyz file,
         * which may be of the form:
         *
         *      7
         *      [*]C(C)=O acetyl,come,coch3,ac    <- title
         *      R   1.56360   1.04770 -0.18810
         *      C   0.81620  -0.01900  0.15360
         *      .    .          .          .
         *
         *  where the title line contains the SMILES string and
         *  a set of name aliases of the fragment. R is the atom
         *  that will be deleted in favour of the core molecule.
         *
         *                 O
         *               //
         *          R---C -- Me
         *
         *  Here the index of the dummy (R) atom is 0 and have
         *  aliases acetyl,come,coch3,ac
         *
         * Arguments:
         *      xyz_filename (string):
         ********************************************************/

        // Populate the name aliases of this fragment
        vector<string> smiles_aliases = utils::split(xyz_title_line, ' ');

        if (smiles_aliases.size() == 2){
            // Assume aliases are second item in the space separated list
            aliases = utils::split(smiles_aliases[1], ',');
        }

        if (n_masked_atoms() != 1){
            throw runtime_error("Cannot construct a fragment molecule with "
                                "no or more than one dummy (R) atom");
        }

        this->dummy_idx = masked_atom_idxs()[0];
    }

    Fragment::Fragment(const Fragment &fragment) : Molecule(fragment) {
        // Copy constructor
        this->rot_grid_w = fragment.rot_grid_w;
        this->aliases = fragment.aliases;
        this->dummy_idx = fragment.dummy_idx;
    }

    FragmentLib::FragmentLib() {
        /*********************************************************
         * Construct the fragment library from .xyz files in the
         * data/ directory
         ********************************************************/

        string file_path = __FILE__;
        string dir_path = file_path.substr(0, file_path.rfind('/'));

        fragments = {Fragment(dir_path+"/data/Ac.xyz"),
                     Fragment(dir_path+"/data/Bn.xyz"),
                     Fragment(dir_path+"/data/Boc.xyz"),
                     Fragment(dir_path+"/data/Br.xyz"),
                     Fragment(dir_path+"/data/Bz.xyz"),
                     Fragment(dir_path+"/data/CF3.xyz"),
                     Fragment(dir_path+"/data/CH2OH.xyz"),
                     Fragment(dir_path+"/data/Cl.xyz"),
                     Fragment(dir_path+"/data/CN.xyz"),
                     Fragment(dir_path+"/data/CO2Et.xyz"),
                     Fragment(dir_path+"/data/CO2Me.xyz"),
                     Fragment(dir_path+"/data/Et.xyz"),
                     Fragment(dir_path+"/data/F.xyz"),
                     Fragment(dir_path+"/data/H.xyz"),
                     Fragment(dir_path+"/data/I.xyz"),
                     Fragment(dir_path+"/data/iPr.xyz"),
                     Fragment(dir_path+"/data/Me.xyz"),
                     Fragment(dir_path+"/data/Mes.xyz"),
                     Fragment(dir_path+"/data/Ms.xyz"),
                     Fragment(dir_path+"/data/NH2.xyz"),
                     Fragment(dir_path+"/data/NMe2.xyz"),
                     Fragment(dir_path+"/data/NO2.xyz"),
                     Fragment(dir_path+"/data/OH.xyz"),
                     Fragment(dir_path+"/data/OMe.xyz"),
                     Fragment(dir_path+"/data/OPh.xyz"),
                     Fragment(dir_path+"/data/OtBu.xyz"),
                     Fragment(dir_path+"/data/Ph.xyz"),
                     Fragment(dir_path+"/data/tBu.xyz"),
                     Fragment(dir_path+"/data/Tf.xyz"),
                     Fragment(dir_path+"/data/TMS.xyz")
        };
    }

    Fragment FragmentLib::fragment(const string& name){
        /*********************************************************
         * Get a fragment from the library given a name, which must
         * match one of the aliases of the fragment
         *
         * Arguments:
         *      name (str):
         *
         * Raises:
         *      (domain_error):
         ********************************************************/
        string l_name = utils::to_lower(name);

        for (auto &fragment : fragments){
            for (auto &alias : fragment.aliases){

                if (alias == l_name){

                    // Need to copy the fragment as it will be modified in-place
                    return Fragment(fragment);
                }
            }// aliases
        }// fragments

        throw domain_error("Failed to find a fragment with an alias "+
                           name+" in the library");
    }

}
