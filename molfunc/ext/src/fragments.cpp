#include <iostream>
#include "fragments.h"
#include "utils.h"

#include "iostream"


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
         *  that will be deleted in favour of the core molecule
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
    }

    Fragment::Fragment(const vector<Atom>& atoms,
                       const string& title):Molecule(atoms) {
        /*********************************************************
         * Construct a Fragment molecule from a set of atoms
         *
         * Arguments:
         *      atoms (list(Atom)):
         *
         *      title (string): e.g. [*]C(C)=O acetyl,come,coch3,ac
         ********************************************************/

        // Populate the name aliases of this fragment
        vector<string> smiles_aliases = utils::split(title, ' ');

        if (smiles_aliases.size() == 2){
            // Assume aliases are second item in the space separated list
            aliases = utils::split(smiles_aliases[1], ',');
        }
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

}
