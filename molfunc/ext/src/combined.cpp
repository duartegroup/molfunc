#include "vector"
#include "stdexcept"
#include "combined.h"
#include "iostream"

using namespace std;


namespace molfunc{

    CombinedMolecule::CombinedMolecule() = default;

    CombinedMolecule::CombinedMolecule(CoreMolecule core,
                                       vector<Fragment> fragments) {
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

        this->core = move(core);
        this->fragments = move(fragments);

        if (!this->fragments.empty()){
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

        rotate_fragments();

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
        vector<Atom3D> atoms;

        for (int i=0; i<core.n_atoms(); i++){

            auto atom = core.atoms[i];

            if (!atom.masked){
                atoms.emplace_back(atom.symbol, core.coordinates[i]);
            }
        }

        for (auto &fragment : fragments){
            for (int i=0; i<fragment.n_atoms(); i++){

                auto atom = fragment.atoms[i];

                if (!atom.masked){
                    atoms.emplace_back(atom.symbol, fragment.coordinates[i]);
                }
            }
        }

        return Molecule(atoms);
    }

    double CombinedMolecule::repulsive_energy() {
        /*********************************************************
         * Calculate the repulsive energy between the fragments
         * and the core
         *
         *      E = Σ_fragments Σ_ij 1/(r_ij^4)
         *
         * for all unique pairs where i enumerates atoms in the
         * core and j atoms in the specific fragment.
         *
         * Returns:
         *      E (float): Repulsive energy
         ********************************************************/

        double energy = 0.0;

        for (const auto &fragment : fragments){
            energy += repulsive_energy(fragment);
        }

        return energy;
    }

    double CombinedMolecule::repulsive_energy(const Fragment& fragment) {
        /*********************************************************
         * Calculate the repulsive energy between the core and a
         * single fragment
         *
         *      E =  Σ_ij 1/r_ij
         *
         * for all unique pairs where i enumerates atoms in the
         * core and j atoms in the fragment.
         *
         * Returns:
         *      E (float): Repulsive energy
         ********************************************************/

        double energy = 0.0;

        for (int i=0; i<core.n_atoms(); i++){

            // TODO: something that is contigious in memory
            if (core.atoms[i].masked) continue;  // Skip masked atoms

            for (int j=0; j<fragment.n_atoms(); j++){

                if (fragment.atoms[j].masked) continue;

                // inline r^2 evaluation for speed(?)
                double r_sq = 0.0;
                for (int k=0; k<3; k++){
                    double tmp = (core.coordinates[i][k]
                            - fragment.coordinates[j][k]);
                    r_sq += tmp * tmp;
                }

                energy += 1.0 / (r_sq * r_sq);

            }// fragment coordinates
        }// core coordinates

        return energy;
    }

    void CombinedMolecule::exclude_rotational_space(Fragment &fragment,
                                                    double threshold){
        /*****************************************************
         * Rotate a fragment to exclude rotational space
         * based on fragment–core repulsion only
         *
         * Arguments:
         *      fragment:
         ****************************************************/

        // Enumerate backwards through the vector, so that the indexing
        // remains valid
        auto it = fragment.rot_grid_w.rbegin();

        while (it != fragment.rot_grid_w.rend()){

            fragment.rotate()

            if (repulsive_energy(fragment) > threshold){
                fragment.rot_grid_w.erase(it);
            }

            it++;
        }

    }

    void CombinedMolecule::rotate_fragments(){
        /*****************************************************
         * Rotate the fragments to minimise the total energy
         * using the repulsion and ...
         * TODO: angle potential
         *
         *
         ****************************************************/

        for (auto &fragment : fragments){
            exclude_rotational_space(fragment);
        }


    }
}


