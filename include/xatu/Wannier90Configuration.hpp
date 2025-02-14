#pragma once
#include "xatu/SystemConfiguration.hpp"

typedef std::vector<std::vector<std::vector<double>>> cube_vector;

namespace xatu {

class Wannier90Configuration : public SystemConfiguration {

    public:
        int mSize, nFock, ndim, filling;
        arma::mat Rn, lattice, bravaisVectors;             // Real-space lattice vectors (3x3)
        arma::imat iRn;                                      // index of neighbors
        arma::cx_dcube fockMatrices;               // Hamiltonian cube (nFock, ndim, ndim)
        arma::field<arma::cx_cube> Rhop;         // Orbital localization cube (3, nFock, ndim, ndim)
        arma::mat motif, bravaisLattice;        // Diagonal elements for motif localization
        arma::irowvec Degen;                     // Degeneracies for the states

    public:
        Wannier90Configuration(std::string, int electronNum = 1);
        ~Wannier90Configuration(){};

        void parseContent(int);

    protected:
        void mapContent(bool debug = false);

};

} 
