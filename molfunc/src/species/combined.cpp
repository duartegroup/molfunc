#include "vector"
#include "algorithm"
#include "stdexcept"
#include "species/combined.h"
#include "iostream"
#include "random"


using namespace std;


namespace molfunc{

    void print_combined_molecule_from_names(const string& xyz_filename,
                                            const string& core_xyz_filename,
                                            const vector<unsigned int>& atom_idxs_to_del,
                                            const vector<string>& frag_names){
        /*************************************************************
         * Generate and print a .xyz file for a combined molecule
         * from a core 3D structure and names of fragments that may,
         * or may not be present in the fragment library
         *
         *  xyz_filename: Name of the generated .xyz file
         *
         *  core_xyz_filename:
         *
         *  frag_names: Names of the fragments e.g. {'Me', ...}
         ************************************************************/

        auto core = CoreMolecule(core_xyz_filename,
                                 atom_idxs_to_del);

        vector<Fragment> fragments;
        fragments.reserve(frag_names.size());
        for (auto &name : frag_names){
            fragments.push_back(FragmentLib::instance().fragment(name));
        }

        CombinedMolecule(core, fragments).to_molecule().print_xyz_file(xyz_filename);
    }

    void print_combined_molecule_from_xyz_filenames(const string& xyz_filename,
                                                    const string& core_xyz_filename,
                                                    const vector<unsigned int>& atom_idxs_to_del,
                                                    const vector<string>& frag_xyz_filenames){
        /******************************************************************
         * Generate and print a .xyz file for a combined molecule
         * from a core 3D structure and filenames of .xys files
         * which contain an R atom
         *
         *  xyz_filename: Name of the generated .xyz file
         *
         *  core_xyz_filename:
         *
         *  frag_xyz_filenames: Names of the fragments e.g. {'Me.xyz', ...}
         *****************************************************************/

        auto core = CoreMolecule(core_xyz_filename,
                                 atom_idxs_to_del);

        vector<Fragment> fragments;
        fragments.reserve(frag_xyz_filenames.size());

        for (auto &filename : frag_xyz_filenames){
            fragments.emplace_back(filename);
        }

        CombinedMolecule(core, fragments).to_molecule().print_xyz_file(xyz_filename);
    }

    CombinedMolecule::CombinedMolecule() = default;

    CombinedMolecule::CombinedMolecule(CoreMolecule core,
                                       vector<Fragment> fragments) {
        /************************************************************
         * Generate a combined molecule from a core and a number of
         * fragments e.g.
         *
         *
         *            H      H                 H       F_1
         *             C -- C                   C -- C
         *           //     \\                //     \\
         *        H C        C H     -->   H C        C  F_1
         *           \      /                 \      /
         *            C == C                   C == C
         *           H      H                 H      H
         *
         *
         *  where F_1 and F_2 are two fragments that replace hydrogens
         *
         * Arguments:
         *      core (CoreMolecule):
         *
         *      fragments (list(Fragment)): Modified in place
         ***********************************************************/

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

        for (unsigned long i=0; i<fragments.size(); i++){
            translate_fragment(fragments[i], dummy_idxs[i]);
            exclude_rotational_space(fragments[i], 2.0);
        }

        rotate_fragments_global();
    }

    void CombinedMolecule::translate_fragment(Fragment &fragment,
                                              unsigned long Ra_idx){
        /***********************************************************
         * Translate the fragments such that the R (dummy) atoms are
         * coincident, then translated by the rough bond distance
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
         *  Where, Ra and Rb are the core and fragment dummy atoms
         *  respectively, and here the carbon is 'x' and the
         *  bromine 'y'.
         *
         * NOTE: will modify the fragment in-place
         ********************************************/

        unsigned long x_idx = core.graph.first_neighbour(Ra_idx);
        //                                                              Rb
        unsigned long y_idx = fragment.graph.first_neighbour(fragment.dummy_idx);

        // Translate so x and y atoms are coincident
        fragment.translate(core.coordinates[x_idx] - fragment.coordinates[y_idx]);

        double xy_dist = (core.atoms[x_idx].covalent_radius()
                          + fragment.atoms[y_idx].covalent_radius());

        // Then translate so the x-y distance is reasonable
        fragment.translate(core.n_vector(x_idx, Ra_idx) * xy_dist);
    }

    Molecule CombinedMolecule::to_molecule() {
        /***********************************************************
         * Construct a standard molecule from this combined molecule
        **********************************************************/
        vector<Atom3D> atoms;

        for (unsigned long i=0; i<core.n_atoms(); i++){

            auto atom = core.atoms[i];  // Use a copy of the atom

            if (!atom.masked){
                atoms.emplace_back(atom.symbol, core.coordinates[i]);
            }
        }

        for (auto &fragment : fragments){
            for (unsigned long i=0; i<fragment.n_atoms(); i++){

                auto atom = fragment.atoms[i];  // Use a copy of the atom

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
         *      E =  Σ_ij 1/(r_ij^4)
         *
         * for all unique pairs. NOTE: bonded repulsions are
         * included, but just provide a constant shift to the
         * rigid-body energy. Also there are relatively few bonded
         * pairs making the wasted evaluations not too severe an
         * overhead.
         *
         * Returns:
         *      E (float): Repulsive energy
         ********************************************************/

        vector<Coordinate> coords;

        // Reserve at least enough space to store all the coordinates
        unsigned long total_n_atoms = core.n_atoms();
        for (auto &frag : fragments) total_n_atoms += frag.n_atoms();
        coords.reserve(total_n_atoms);

        // Populate the coordinates with the unmasked core atoms
        for (unsigned long i=0; i<core.n_atoms(); i++){
            if (!core.atoms[i].masked) coords.push_back(core.coordinates[i]);
        }
        // and from each of the fragments
        for (auto &frag : fragments){
            for (unsigned long i=0; i<frag.n_atoms(); i++){
                if (!frag.atoms[i].masked) coords.push_back(frag.coordinates[i]);
            }
        }

        double energy = 0.0;

        // Enumerate all pairwise
        for (unsigned long i=0; i<coords.size(); i++){
            for (unsigned long j=i+1; j<coords.size(); j++){

                // inline r^2 evaluation for speed(?)
                double r_sq = 0.0;

                for (auto k=0; k<3; k++){
                    double tmp = (coords[i][k] - coords[j][k]);
                    r_sq += tmp * tmp;
                }

                energy += 1.0 / (r_sq * r_sq);
            }
        }

        return energy;
    }

    double CombinedMolecule::repulsive_energy(const Fragment &fragment) {
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

        for (unsigned long i=0; i<core.n_atoms(); i++){

            // TODO: something that is contiguous in memory
            if (core.atoms[i].masked) continue;  // Skip masked atoms

            for (unsigned long j=0; j<fragment.n_atoms(); j++){

                if (fragment.atoms[j].masked) continue;

                // inline r^2 evaluation for speed(?)
                double r_sq = 0.0;
                for (auto k=0; k<3; k++){
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
        auto y_idx = fragment.graph.first_neighbour(fragment.dummy_idx);
        auto y_coord = fragment.coordinates[y_idx];

        // Shift both the core and the fragment to the new origin
        for (auto &coord: core.coordinates) coord -= y_coord;
        for (auto &coord: fragment.coordinates) coord -= y_coord;

        // After each rotation we'll want to reset the coordinates, so cache
        // their current state
        fragment.cache_coordinates();

        // Enumerate backwards through the vector, so that the indexing
        // remains valid while deleting elements
        int end_idx = static_cast<int>(fragment.rot_grid_w.size() - 1);

        for (int i=end_idx; i>=0; i--){

            fragment.rotate(fragment.rot_grid_w[i]);
            fragment.rot_grid_w[i].energy = repulsive_energy(fragment);

            if (fragment.rot_grid_w[i].energy > threshold){
                fragment.rot_grid_w.erase(fragment.rot_grid_w.begin() + i);
            }

            fragment.reset_coordinates();
        }

        if (fragment.rot_grid_w.empty()) throw runtime_error("Deleted all points!");

        // Shift back, such that the core remains in the same position
        for (auto &coord: core.coordinates) coord += y_coord;
        for (auto &coord: fragment.coordinates) coord += y_coord;
    }

    void CombinedMolecule::rotate_fragments_global(){
        /*****************************************************
         * Rotate the fragments to minimise the total energy
         * using the repulsion and ...
         * TODO: angle potential
         *
         *
         ****************************************************/

        if (fragments.size() == 1){
            // No global optimisation needs to be done - simply
            // use the minimum energy rotation of the fragment
            auto point = fragments[0].rot_grid_w.minimum_energy_point();
            fragments[0].rotate_about_dummy_nn(point);
            return;
        }

        int max_iters = 1000;
        double min_energy = INFINITY;

        vector<GridPoint> min_points, points;
        points.reserve(fragments.size());

        for (auto &frag : fragments) frag.cache_coordinates();

        for (int iter=0; iter<max_iters; iter++){

            points.clear();   // The points used to rotate

            for (auto &frag : fragments){

                auto point = frag.rot_grid_w.random_point();
                points.push_back(point);
                frag.rotate_about_dummy_nn(point);
            }

            if (repulsive_energy() < min_energy){
                min_points = points;
                min_energy = repulsive_energy();
            }

            for (auto &frag : fragments){
                frag.reset_coordinates();
            }
        }


        // Finally, apply the minimum energy rotation
        for (unsigned long i=0; i<fragments.size(); i++){
            fragments[i].rotate_about_dummy_nn(min_points[i]);
        }

    }
}


