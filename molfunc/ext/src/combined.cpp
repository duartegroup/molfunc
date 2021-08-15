#include "vector"
#include "stdexcept"
#include "combined.h"

using namespace std;


namespace molfunc{

    CombinedMolecule::CombinedMolecule() = default;

    CombinedMolecule::CombinedMolecule(CoreMolecule &core,
                                       vector<Fragment> &fragments) {
        /*********************************************
         * Generate a combined molecule from a core
         * and a number of fragments
         *
         * Arguments:
         *      core (CoreMolecule):
         *
         *      fragments (list(Fragment)): Modified
         *                                  in place
         ********************************************/

        this->core = core;
        this->fragments = fragments;

        if (!fragments.empty()){
            build();
        }
    }

    void CombinedMolecule::build() {
        /*********************************************
         * Rotate and translate the fragments such that
         * the R atoms are coincident, then translated
         * by the rough bond distance, and finally the
         * energy of the whole system minimised. e.g.
         *
         *         H
         *          \
         *      H---C---R             R---Br
         *         /
         *       H
         *
         *        ^                     ^
         *      core                 fragment
         *
         * in the most simple case. More complex
         * scenarios with polyatomic fragments, and
         * where there is more than just a single
         * fragment. For >1 fragments they are added
         * in ascending order of atom index
         *
         * NOTE: will modify the fragments in-place
         ********************************************/

        if (core.n_masked_atoms() != fragments.size()){
            throw runtime_error("Cannot add fragments to core. Number of "
                                "fragments was not equal to the number of "
                                "dummy (atoms_to_del) atoms in the core.");
        }
        vector<unsigned long> dummy_idxs = core.masked_atom_idxs();

        for (int i=0; i<fragments.size(); i++){
            translate_fragment(fragments[i],
                               dummy_idxs[i]);

            //prune_rotational_space(fragments[i]);
        }

        //rotate_fragments();

    }

    void CombinedMolecule::translate_fragment(Fragment &fragment,
                                              unsigned long core_dummy_atom_idx){
        /*********************************************
         * Translate the fragments such that the R
         * (dummy) atoms are coincident, then translated
         * by the rough bond distance
         *
         *         H                            H
         *          \                            \
         *      H---C---Ra    Rb---Br   ==>   H---C---Br
         *         /                             /
         *       H                              H
         *
         *        ^            ^
         *      core        fragment
         *
         *  Where, Ra and Rb are the core and fragment
         *  dummy atoms respectively, and here the
         *  carbon is 'x' and the bromine 'y'
         *
         * NOTE: will modify the fragment in-place
         ********************************************/

        unsigned long x_idx = core.graph.first_neighbour(core_dummy_atom_idx);
        array<double, 3> x_pos = core.coordinates[x_idx];

        unsigned long y_idx = fragment.graph.first_neighbour(fragment.dummy_idx);
        array<double, 3> y_pos = fragment.coordinates[y_idx];

        // Translate so x and y are coincident
        fragment.translate({x_pos[0] - y_pos[0],
                            x_pos[1] - y_pos[1],
                            x_pos[2] - y_pos[2]});

        double xy_dist = (core.atoms[x_idx].covalent_radius()
                          + fragment.atoms[y_idx].covalent_radius());

        // TODO this is v ugly. Eigen?
        auto n_vec = core.n_vector(x_idx, core_dummy_atom_idx);
        fragment.translate({xy_dist * n_vec[0],
                            xy_dist * n_vec[2],
                            xy_dist * n_vec[1]});
    }

    Molecule CombinedMolecule::to_molecule() {
        // Construct a standard molecule from this combined molecule
        vector<Atom> atoms;
        for (auto atom : core.atoms){
            if (!atom.masked) atoms.push_back(atom);
        }

        for (auto &fragment : fragments){
            for (auto atom : fragment.atoms){
                if (!atom.masked) atoms.push_back(atom);
            }
        }

        return Molecule(atoms);
    }

}


