# distutils: language = c++
# distutils: sources = [molfunc/src/atoms.cpp, molfunc/src/graph.cpp, molfunc/src/grid.cpp, molfunc/src/rotation.cpp, molfunc/src/utils.cpp, molfunc/src/vector3d.cpp, molfunc/src/species/combined.cpp, molfunc/src/species/fragments.cpp, molfunc/src/species/molecules.cpp, molfunc/src/species/species.cpp]
from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from "include/species/combined.h" namespace "molfunc":

    void print_combined_molecule_from_names(string& xyz_filename,
                                            string& core_xyz_filename,
                                            vector[unsigned int]& atom_idxs_to_del,
                                            vector[string]& frag_names) except +


    void print_combined_molecule_from_xyz_filenames(string& xyz_filename,
                                                    string& core_xyz_filename,
                                                    vector[unsigned int]& atom_idxs_to_del,
                                                    vector[string]& frag_xyz_filenames) except +

