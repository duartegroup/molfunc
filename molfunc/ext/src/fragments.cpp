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
         *      [*]C(C)=O acetyl,come,coch3,ac
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

    constexpr FragmentLib::FragmentLib() {
        /*********************************************************
         * Construct a fragment library from .xyz files at
         * compile time
         ********************************************************/

        const char* curr_path = __FILE__;



    }

}
