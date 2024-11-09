#include "xatu/Wannier90Configuration.hpp"
// #include <armadillo>
// #include <fstream>
// #include <sstream>
// #include <complex>
// #include <stdexcept>
// #include <vector>
// #include <iostream>

namespace xatu {

    /**
     * File constructor for Wannier90Configuration. It extracts the relevant information
     * from the Wannier90 output file and stores it in the appropriate format for further calculations.
     * @details This class is specifically designed to parse Wannier90 output files with Hamiltonian
     * and orbital localization data (motif) for exciton calculations. It initializes with an optional 
     * file path and parses the content according to the structure defined in the Wannier90 output.
     * @param file Name of the output file from Wannier90.
     */
    Wannier90Configuration::Wannier90Configuration(std::string file, int electronNum) : ConfigurationBase{file} {
        parseContent();
        mapContent(false, electronNum);
    }

    /**
     * Method to parse the content of the Wannier90 output file.
     * @details This function reads the Hamiltonian matrix and motif localization data (Rhop)
     * from the Wannier90 output file. The structure assumes that Rn matrix, ndim, nFock, Degen,
     * Hamiltonian elements, and orbital localization data are in specific order within the file.
     * @return void
     */
    void Wannier90Configuration::parseContent() {
        std::string line;
        // if (!fileStream.is_open()) {
        //     throw std::runtime_error("Failed to open file: " + file);
        // }

        // Parse Rn matrix (3x3)
        Rn = arma::mat(3, 3);
        std::getline(m_file, line);  // Skip header line (assuming there's a header to skip)
        for (int i = 0; i < 3; ++i) {
            std::getline(m_file, line);
            std::istringstream iss(line);
            iss >> Rn(i, 0) >> Rn(i, 1) >> Rn(i, 2);
            // Ensure we have exactly three values in each line
            // if (!(iss >> Rn(i, 0) >> Rn(i, 1) >> Rn(i, 2))) {
            //     throw std::runtime_error("Error reading Rn matrix values from file.");
            // }
        }

        // Parse ndim and nFock values
        std::istringstream iss(line);
        iss >> ndim;
        
        std::getline(m_file, line);  // Move to the next line
        
        // std::istringstream iss(line);
        iss >> nFock;
        
        std::getline(m_file, line);  // Move to the next line


        // Allocate matrices and vectors based on parsed sizes
        H = arma::cx_cube(nFock, ndim, ndim, arma::fill::zeros);
        Degen = arma::colvec(nFock, arma::fill::zeros); // armadillo is faster with columns
        iRn = arma::colvec(nFock, 3, arma::fill::zeros);

        //Rhop = arma::cx_cube(3, nFock, ndim, ndim, arma::fill::zeros); // don't work on armadillo
        Rhop = arma::cx_cube(3, nFock, ndim * ndim, arma::fill::zeros); // option 1: flatten one axis
        // option 2: nested vectors
        // std::vector<std::vector<arma::cx_mat>> Rhop(3, std::vector<arma::cx_mat>(nFock, arma::cx_mat(ndim, ndim, arma::fill::zeros)));

        // parsing DEGENERACIES
        int degenPerLine = 15;
        int idx = 0;

        // Read full lines of 15 elements
        for (int i = 0; i < nFock / degenPerLine; ++i) {
            for (int j = 0; j < degenPerLine; ++j) {
                std::istringstream iss(line);
                iss >> Degen[idx];
                idx++;
            }
        }

        // Read the last partial line if needed
        int remaining = nFock % degenPerLine;
        for (int j = 0; j < remaining; ++j) {
            std::istringstream iss(line);
            iss >> Degen[idx];
            idx++;
        }

        std::getline(m_file, line);  // Move to the next line

        // Parse Hamiltonian (H) data
        for (int i = 0; i < nFock; ++i) {
            std::istringstream iss(line);
            iss >> iRn(i, 0) >> iRn(i, 1) >> iRn(i, 2);
            
            // Loop through rows and columns for ndim x ndim block
            for (int row = 0; row < ndim; ++row) {
                for (int col = 0; col < ndim; ++col) {
                    double R, Im;
                    std::istringstream iss(line);
                    iss >> R >> Im;
                    H(i, row, col) = std::complex<double>(R, Im);
                }
            }
            if (i < nFock - 1) std::getline(m_file, line); // Skip blank line if not the last
        }


        // Parse orbital localization data (Rhop) - diagonal is motif
        for (int i = 0; i < nFock; ++i) {
            std::getline(m_file, line);  // Skip iRn line for Rhop

            for (int ii = 0; ii < ndim; ++ii) {
                for (int jj = 0; jj < ndim; ++jj) {
                    std::getline(m_file, line);
                    std::istringstream iss(line);

                    double a1, a1j, a2, a2j, a3, a3j;
                    iss >> a1 >> a1j >> a2 >> a2j >> a3 >> a3j;

                    // Flatten last two indices ii, jj to a single (ii * ndim + jj) index
                    Rhop(0, i, ii * ndim + jj) = std::complex<double>(a1, a1j);
                    Rhop(1, i, ii * ndim + jj) = std::complex<double>(a2, a2j);
                    Rhop(2, i, ii * ndim + jj) = std::complex<double>(a3, a3j);
                }
            }

            if (i < nFock - 1) std::getline(m_file, line);  // Skip blank line if not last
        }

        // Ensure motif is complex to handle the complex entries from Rhop
        motif = arma::mat(3, nFock, arma::fill::zeros);

        for (int i = 0; i < nFock; ++i) {
            for (int k = 0; k < ndim; ++k) {
                motif(0, i) = std::real(Rhop(0, i, k * ndim + k));  // Diagonal element for Rhop(0)
                motif(1, i) = std::real(Rhop(1, i, k * ndim + k));  // Diagonal element for Rhop(1)
                motif(2, i) = std::real(Rhop(2, i, k * ndim + k));  // Diagonal element for Rhop(2)
            }
        }


        m_file.close();
        std::cout << "Parsed Wannier90 file content successfully." << std::endl;
    }

        /**
         * Method to map parsed content into internal data structures.
         * @details This function processes the parsed data, mapping them to appropriate
         * representations within the Wannier90Configuration class for downstream usage.
         * @return void
         */
        /**
     * Method to write all the extracted information into a struct.
     * @return void 
     */
    void Wannier90Configuration::mapContent(bool debug, int electronNum) {
        
        // Fill in the system info structure with relevant data
    /**
    *    Necessary data:
    *    systemInfo.ndim           = ndim;
    *    systemInfo.bravaisLattice = bravaisLattice;
    *    systemInfo.motif          = motif;
    *    systemInfo.filling        = totalElectrons/2;
    *    systemInfo.bravaisVectors = bravaisVectors;
    *    systemInfo.overlap        = overlapMatrices;
    *    systemInfo.hamiltonian    = fockMatrices;
    **/

        systemInfo.ndim = ndim;                  // system/matrix size
        systemInfo.bravaisVectors = Rn;          // Store bravais vectors
        systemInfo.motif = motif;                // Store motif localization data
        systemInfo.hamiltonian = H;              // Store Hamiltonian
        // systemInfo.overlap = Rhop;               // ? check this
        systemInfo.filling = electronNum;        /* missing --> totalElectrons on w90 file */

        // Debug output
        if (debug) {
            // std::cout << "Number of Fock matrices (nFock): " << systemInfo.nFock << "\n";
            std::cout << "Matrix size (ndim): " << systemInfo.ndim << "\n";
            std::cout << "Motif: " << "\n";
            std::cout << systemInfo.motif << "\n";
            std::cout << "Overlap (Rhop): " << std::endl;
            std::cout << systemInfo.overlap << "\n";
            std::cout << "Hamiltonian matrix: " << "\n";
            std::cout << systemInfo.hamiltonian << "\n";
        }
    }
}