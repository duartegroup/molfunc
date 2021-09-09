#include "iostream"
#include <algorithm>
#include "species/fragments.h"
#include "utils.h"


using namespace std;


namespace molfunc{

    Fragment::Fragment() = default;

    Fragment::Fragment(const string& xyz_filename): Molecule(xyz_filename) {
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
        vector<string> smiles_and_aliases = utils::split(xyz_title_line, ' ');

        if (smiles_and_aliases.size() == 2){
            // Assume aliases are second item in the space separated list
            aliases = utils::split(smiles_and_aliases[1], ',');
        }

        if (n_masked_atoms() != 1){
            throw runtime_error("Cannot construct a fragment molecule with "
                                "no or more than one dummy (R) atom");
        }

        this->dummy_idx = masked_atom_idxs()[0];
        this->dummy_nn_idx = graph.first_neighbour(dummy_idx);
    }

    Fragment::Fragment(const Fragment &fragment): Molecule(fragment) {
        // Copy constructor
        this->rot_grid_w = fragment.rot_grid_w;
        this->aliases = fragment.aliases;
        this->dummy_idx = fragment.dummy_idx;
        this->dummy_nn_idx = fragment.dummy_nn_idx;
        this->cached_coordinates = vector<Coordinate>(fragment.coordinates);
    }

    Fragment::Fragment(vector<Atom3D> atoms, vector<string> aliases)
             : Molecule(atoms){
        /************************************************************
         * Construct a fragment from a set of atoms and
         * a list of aliases e.g. Me, CH3...
         ***********************************************************/
        this->aliases = move(aliases);

        if (n_masked_atoms() != 1){
            throw runtime_error("Cannot construct a fragment molecule with "
                                "no or more than one dummy (R) atom");
        }

        this->dummy_idx = masked_atom_idxs()[0];
        this->dummy_nn_idx = graph.first_neighbour(dummy_idx);
    }

    void Fragment::cache_coordinates(){
        /********************************************
         * Cache the coordinates so that they may
         * be reset following a translation/rotation
         *******************************************/

        this->cached_coordinates = vector<Coordinate>(coordinates);
    }

    void Fragment::reset_coordinates(){
        /************************************
         * Reset the coordinates using the
         * cached values
         ***********************************/
        if (cached_coordinates.empty()){
            throw runtime_error("Cannot reset the coordinates, no "
                                "cached coordinates found.");
        }

         for (unsigned long i=0; i<n_atoms(); i++){
             coordinates[i] = cached_coordinates[i];
         }
    }

    void Fragment::rotate(GridPoint &grid_point){
        /************************************************
         * Rotate this fragment using a point in the grid
         ***********************************************/
         rotation_matrix.update(grid_point);
         Species::rotate(rotation_matrix);
    }

    void Fragment::rotate_about_dummy_nn(GridPoint &grid_point){
        /************************************************
         * Rotate this fragment using a point in the grid
         * on a defined atom index (i.e. coordinate in
         * space)
         ***********************************************/
        rotation_matrix.update(grid_point);
        Species::rotate(rotation_matrix, dummy_nn_idx);
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
