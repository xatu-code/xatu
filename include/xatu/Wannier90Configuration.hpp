#pragma once
#include "xatu/SystemConfiguration.hpp"

typedef std::vector<std::vector<std::vector<double>>> cube_vector;

namespace xatu {

class Wannier90Configuration : public SystemConfiguration {

    public:
        int mSize, nFock, ndim, electronNum;
        arma::mat Rn;                           // Real-space lattice vectors (3x3)
        arma::mat iRn;                       // index of neighbors
        arma::cx_cube H;                        // Hamiltonian cube (nFock, ndim, ndim)
        arma::cx_cube Rhop;                     // Orbital localization cube (3, nFock, ndim, ndim)
        arma::rowvec Degen;                     // Degeneracies for the states
        arma::mat motif, bravaisLattice;        // Diagonal elements for motif localization

    public:
        Wannier90Configuration(std::string, int electronNum = 26);
        ~Wannier90Configuration(){};

        void parseContent();

    /* to be added: organize parsers into functions */
    // private:
    //     void parseBravaisLattice();
    //     void parseHamiltonian();
    //     void parseDegeneracy();
    //     void parseMotif();

    private:
       void parseBravaisLattice();
       void extractDimension();

    protected:
        void mapContent(bool debug = false, int electronNum = 26);


};

} 
