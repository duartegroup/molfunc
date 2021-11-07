# distutils: language = c++
# distutils: sources = [molfunc/src/atoms.cpp, molfunc/src/graph.cpp, molfunc/src/grid.cpp, molfunc/src/rotation.cpp, molfunc/src/utils.cpp, molfunc/src/vector3d.cpp, molfunc/src/species/combined.cpp, molfunc/src/species/fragments.cpp, molfunc/src/species/molecules.cpp, molfunc/src/species/species.cpp, molfunc/src/angles.cpp]
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from "include/species/combined.h" namespace "molfunc":

    void print_all_combined_molecules(string& xyz_filename,
                                      string& core_xyz_filename,
                                      vector[unsigned int]& atom_idxs_to_del) except +

    void print_combined_molecule_from_names(string& xyz_filename,
                                            string& core_xyz_filename,
                                            vector[unsigned int]& atom_idxs_to_del,
                                            vector[string]& frag_names) except +


    void print_combined_molecule_from_xyz_filenames(string& xyz_filename,
                                                    string& core_xyz_filename,
                                                    vector[unsigned int]& atom_idxs_to_del,
                                                    vector[string]& frag_xyz_filenames) except +


def c_print_all_combined_molecules(string& xyz_filename, string& core_xyz_filename, vector[unsigned int]& atom_idxs_to_del):
    print_all_combined_molecules(xyz_filename, core_xyz_filename, atom_idxs_to_del)

def c_print_combined_from_names(string& xyz_filename, string& core_xyz_filename, vector[unsigned int]& atom_idxs_to_del, vector[string]& frag_names):
    print_combined_molecule_from_names(xyz_filename, core_xyz_filename, atom_idxs_to_del, frag_names)

def c_print_combined_from_xyz_filenames(string& xyz_filename, string& core_xyz_filename, vector[unsigned int]& atom_idxs_to_del, vector[string]& frag_xyz_filenames):
    print_combined_molecule_from_xyz_filenames(xyz_filename, core_xyz_filename, atom_idxs_to_del, frag_xyz_filenames)
