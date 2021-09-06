# distutils: language = c++
from libcpp.vector cimport vector
from libcpp cimport bool as bool_t


cdef extern from "include/molecule.h" namespace "autode":
    cdef cppclass Molecule:
        # Constructors
        Molecule()
        Molecule(vector[double])

        # Attributes exposed to Python
        int n_atoms
