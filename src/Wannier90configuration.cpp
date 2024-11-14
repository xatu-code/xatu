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
        mapContent();
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
        std::getline(m_file, line);  // Skip header line (assuming there's a header to skip)
        // if (!fileStream.is_open()) {
        //     throw std::runtime_error("Failed to open file: " + file);
        // }

        Rn = arma::mat(3, 3, arma::fill::zeros); // Initialize as a 3x3 zero matrix
        for (int i = 0; i < 3; i++) {
            if (!std::getline(m_file, line)) {
                throw std::runtime_error("Error reading Rn matrix: insufficient lines in file.");
            }
            std::istringstream iss(line);
            double val1, val2, val3;
            iss >> val1 >> val2 >> val3;
            Rn(i, 0) = val1;
            Rn(i, 1) = val2;
            Rn(i, 2) = val3;

            std::cout << "Line " << i+1 << ": " << line << std::endl;
        }
        std::getline(m_file, line);  // Move to the next line

        // Parse ndim and nFock values 
        std::istringstream iss(line);
        iss >> mSize; // this is not ndim!!
        std::cout << "mSize: " << mSize << std::endl;
        iss.clear();  // Clear any flags (like EOF)


        std::getline(m_file, line);  // Move to the next line
        iss.str(line);  // Assign new content to the stream

        std::istringstream iss2(line);
        iss >> nFock;
        std::cout << "nFock: " << nFock << std::endl;
        iss.clear();  // Clear any flags (like EOF)

        std::getline(m_file, line);  // Move to the next line

        // Allocate matrices and vectors based on parsed sizes
        H = arma::cx_cube(nFock, mSize, mSize, arma::fill::zeros);
        std::cout << "Declared H matrix... "<< std::endl;

        Degen = arma::rowvec(nFock, arma::fill::zeros); // armadillo is faster with columns
        std::cout << "Declared degeneracies... "<< std::endl;

        iRn = arma::mat(nFock, 3, arma::fill::zeros);
        std::cout << "Declared iRn... "<< std::endl;

        //Rhop = arma::cx_cube(3, nFock, ndim, ndim, arma::fill::zeros); // don't work on armadillo
        Rhop = arma::cx_cube(3, nFock, mSize * mSize, arma::fill::zeros); // option 1: flatten one axis
        // option 2: nested vectors
        // std::vector<std::vector<arma::cx_mat>> Rhop(3, std::vector<arma::cx_mat>(nFock, arma::cx_mat(mSize, mSize, arma::fill::zeros)));

        // parsing DEGENERACIES
        int degenPerLine = 15;
        int idx = 0;

        // Read full lines of 15 elements
        for (int i = 0; i < nFock / degenPerLine; i++) {
            for (int j = 0; j < degenPerLine; j++) {
                // std::istringstream iss(line);
                iss.str(line);
                iss >> Degen[idx];
                idx++;
            }
            std::getline(m_file, line);
        }

        // Read the last partial line if needed
        int remaining = nFock % degenPerLine;
        for (int j = 0; j < remaining; j++) {
            // std::istringstream iss(line);
            iss.str(line);
            iss >> Degen[idx];
            idx++;
        }
        std::getline(m_file, line);  // Move to the next line

        //-------------------HAMILTONIAN-------------------------//
        // Parse Hamiltonian (H) data
        for (int i = 0; i < nFock; i++) {
            std::istringstream iss(line);
            iss >> iRn(i, 0) >> iRn(i, 1) >> iRn(i, 2);
            
            // Loop through rows and columns for mSize x mSize block
            for (int row = 0; row < mSize; row++) {
                for (int col = 0; col < mSize; col++) {
                    double R, Im;
                    std::istringstream iss(line);
                    iss >> R >> Im;
                    H(i, row, col) = std::complex<double>(R, Im);
                }
            }
            if (i < nFock - 1) std::getline(m_file, line); // Skip blank line if not the last
        }


        //----------------------PosMatrixElements-----------------------------//
        // Parse orbital localization data (Rhop) - diagonal is motif
        for (int i = 0; i < nFock; i++) {
            std::getline(m_file, line);  // Skip iRn line for Rhop

            for (int ii = 0; ii < mSize; ii++) {
                for (int jj = 0; jj < mSize; jj++) {
                    std::getline(m_file, line);
                    std::istringstream iss(line);

                    double a1, a1j, a2, a2j, a3, a3j;
                    iss >> a1 >> a1j >> a2 >> a2j >> a3 >> a3j;

                    // Flatten last two indices ii, jj to a single (ii * mSize + jj) index
                    Rhop(0, i, ii * mSize + jj) = std::complex<double>(a1, a1j);
                    Rhop(1, i, ii * mSize + jj) = std::complex<double>(a2, a2j);
                    Rhop(2, i, ii * mSize + jj) = std::complex<double>(a3, a3j);
                }
            }

            if (i < nFock - 1) std::getline(m_file, line);  // Skip blank line if not last
        }

        //----------------------MOTIF-----------------------------//
        // Ensure motif is complex to handle the complex entries from Rhop
        motif = arma::mat(3, mSize, arma::fill::zeros);
        int diag = 0;
        for (int i = 0; i < nFock; i++) {   
            if (iRn(i, 0) == 0 && iRn(i, 1) == 0 && iRn(i, 2) == 0) {
                diag = i;
                break;
            }
        }

        for (int k = 0; k < mSize; k++) {
            motif(0, k) = std::real(Rhop(0, diag, k * mSize + k));  // Diagonal element for Rhop(x)
            motif(1, k) = std::real(Rhop(1, diag, k * mSize + k));  // Diagonal element for Rhop(y)
            motif(2, k) = std::real(Rhop(2, diag, k * mSize + k));  // Diagonal element for Rhop(z)
        }

        m_file.close();
        std::cout << "Parsed Wannier90 file content successfully." << std::endl;
    }


/**
 * Method to parse and format the Bravais basis vectors from the file. 
 * @return void
 */
void Wannier90Configuration::parseBravaisLattice(){
    std::string line;
    std::vector<std::string> vectors;
    for(int i = 0; i < 3; i++){
        std::getline(m_file, line);
        vectors.push_back(line);     
    } 
    this->bravaisLattice = parseVectors(vectors);
    extractDimension();
}

/**
 * Method to obtain the dimension (1D,2D,3D) of the system.
 * @return void 
 */
void Wannier90Configuration::extractDimension(){
    arma::rowvec R0 = arma::zeros<arma::rowvec>(3);
    arma::rowvec R2 = {0.0, 500.0, 0.0};
    arma::rowvec R3 = {0.0, 0.0, 500.0}; 
                 
    if ( approx_equal(bravaisLattice.row(2),R3,"absdiff",0.1) || 
         approx_equal(bravaisLattice.row(2),R0,"absdiff",0.1) ){
        bravaisLattice.shed_row(2);
        if ( approx_equal(bravaisLattice.row(1),R2,"absdiff",0.1) || 
             approx_equal(bravaisLattice.row(1),R0,"absdiff",0.1) ){
            bravaisLattice.shed_row(1); 
        }
    }
    ndim = bravaisLattice.n_rows;
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
        if (true) {
            // std::cout << "Number of Fock matrices (nFock): " << systemInfo.nFock << "\n";
            std::cout << "Matrix size (mSize): " << mSize << "\n";
            std::cout << "nFock: " << nFock << "\n";
            // std::cout << "bravaisVectors: " << "\n";
            // std::cout << systemInfo.bravaisVectors << "\n";
            std::cout << "System Dimension: " << systemInfo.ndim << "\n";
            std::cout << "Filling: " << systemInfo.filling << "\n";
            std::cout << "===================" << "\n";
            std::cout << "Motif: " << "\n";
            std::cout << systemInfo.motif << "\n";
            std::cout << "===================" << "\n";
            // std::cout << "Rhop: " << std::endl;
            // // std::cout << Rhop << "\n";
            // std::cout << "===================" << "\n";
            // std::cout << "Hamiltonian matrix: " << "\n";
            // std::cout << systemInfo.hamiltonian << "\n";
        }
    }
}