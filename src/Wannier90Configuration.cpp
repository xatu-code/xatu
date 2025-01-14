#include "xatu/Wannier90Configuration.hpp"
#include <iomanip> // Required for std::setprecision


namespace xatu {

    /**
     * File constructor for Wannier90Configuration. It extracts the relevant information
     * from the Wannier90 output file and stores it in the appropriate format for further calculations.
     * @details This class is specifically designed to parse Wannier90 output files with Hamiltonian
     * and orbital localization data (motif) for exciton calculations. It initializes with an optional 
     * file path and parses the content according to the structure defined in the Wannier90 output.
     * @param file Name of the output file from Wannier90.
     */
    Wannier90Configuration::Wannier90Configuration(std::string file, int electronNum): ConfigurationBase{file} { // Initialize member
        parseContent(electronNum);
        mapContent(); // Now passes the correct value
    }


    /**
     * Method to parse the content of the Wannier90 output file.
     * @details This function reads the Hamiltonian matrix and motif localization data (Rhop)
     * from the Wannier90 output file. The structure assumes that Rn matrix, ndim, nFock, Degen,
     * Hamiltonian elements, and orbital localization data are in specific order within the file.
     * @return void
     */
    void Wannier90Configuration::parseContent(int electronNum) {
        filling = electronNum;
        std::string line;
        std::getline(m_file, line);  // Skip header line (assuming there's a header to skip)
        // if (!fileStream.is_open()) {
        //     throw std::runtime_error("Failed to open file: " + file);
        // }
        std::cout << std::fixed << std::showpoint;
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
            iss.clear();
            // std::cout << "Line " << i+1 << ": " << line << std::endl;
        }

        // change to check a threshold between highest
        // 2D: motif(z) and Rn(2,2)
        // 1D: motif(x or y) and Rn(1,1) or Rn(0,0)
        if (Rn(2,2) >= 100) {
            ndim = 2;
            if (Rn(1,1) >= 100 || Rn(0,0) >= 100) {
                ndim = 1;
            }
        }

        lattice = arma::mat(ndim, 3, arma::fill::zeros);
        for (int i = 0; i < ndim; i++) {
            lattice(i, 0) = Rn(i,0);
            lattice(i, 1) = Rn(i,1);
            lattice(i, 2) = Rn(i,2);
        }



        // Parse ndim and nFock values 
        std::getline(m_file, line);  // Move to the next line
        std::istringstream iss(line);
        iss >> mSize; // this is not ndim!!
        iss.clear();  // Clear any flags (like EOF)


        std::getline(m_file, line);  // Move to the next line
        iss.str(line);  // Assign new content to the stream
        iss >> nFock;
        iss.clear();  // Clear any flags (like EOF)

        std::getline(m_file, line);  // Move to the next line

        Degen = arma::irowvec(nFock, arma::fill::zeros); // armadillo is faster with columns

        // ---------------- DEGENERACIES ----------------//
        int degenPerLine = 15;
        int numLines = nFock / degenPerLine; // Calculate how many full chunks of 15
        int remainder = nFock % degenPerLine; // Remaining elements

        // Read full chunks of 15 elements
        for (int i = 0; i < numLines; i++) {
            iss.str(line);

            for (int j = 0; j < degenPerLine; j++) {
                iss >> Degen(i * degenPerLine + j);
            }
            std::getline(m_file, line);  // Move to the next line
            iss.clear();
        }

        iss.str(line);
        // Read the remaining elements for the last chunk
        for (int j = 0; j < remainder; ++j) {
            iss >> Degen(numLines * degenPerLine + j);
        }
        std::getline(m_file, line);  // Move to the next line
        iss.clear(); 



        iRn = arma::imat(nFock, 3, arma::fill::zeros);                           // indexes of the neighbors
        fockMatrices = arma::cx_dcube(mSize, mSize, nFock, arma::fill::zeros);   // Hamiltonian


        //-------------------HAMILTONIAN-------------------------//
        // Parse Hamiltonian (fockMatrices) data
        for (int i = 0; i < nFock; i++) {
            std::getline(m_file, line); // Skip blank line if not the last
            iss.str(line);
            iss >> iRn(i, 0) >> iRn(i, 1) >> iRn(i, 2);
            iss.clear(); 
            
            // Loop through rows and columns for mSize x mSize block
            for (int j = 0; j < mSize; j++) {
                for (int k = 0; k < mSize; k++) {
                    std::getline(m_file, line);
                    iss.str(line);

                    int ii, jj;
                    double R, Im;
                    iss >> ii >> jj >> R >> Im;
                    fockMatrices(ii-1, jj-1, i) = std::complex<double>(R, Im);
                    iss.clear(); 
                }
            }
            // fockMatrices.slice(i) /= Degen(i);  // 
            std::getline(m_file, line); // Skip blank line if not the last
        }

        bravaisVectors = arma::dmat(nFock, 3, arma::fill::zeros);

        for (int i = 0; i < nFock; i++) {
            double a1 = iRn(i, 0) * Rn(0, 0) + iRn(i, 1) * Rn(1, 0) + iRn(i, 2) * Rn(2, 0);
            double a2 = iRn(i, 0) * Rn(0, 1) + iRn(i, 1) * Rn(1, 1) + iRn(i, 2) * Rn(2, 1);
            double a3 = iRn(i, 0) * Rn(0, 2) + iRn(i, 1) * Rn(1, 2) + iRn(i, 2) * Rn(2, 2);

            bravaisVectors(i, 0) = a1;
            bravaisVectors(i, 1) = a2;
            bravaisVectors(i, 2) = a3;
        }


        Rhop = arma::field<arma::cx_cube>(nFock); // Initialize field with nFock elements
        // Initialize each field element (cx_cube) individually
        for (int i = 0; i < nFock; i++) {
            Rhop(i) = arma::cx_cube(mSize, mSize, 3, arma::fill::zeros);
        }


        //----------------------PositionMatrixElements-----------------------------//
        // Parse orbital localization data (Rhop) - diagonal is motif
        for (int i = 0; i < nFock; i++) {

            std::getline(m_file, line);  // Skip iRn line for Rhop
            iss.clear();

            for (int j = 0; j < mSize; j++) {
                for (int k = 0; k < mSize; k++) {
                    std::getline(m_file, line); 
                    iss.str(line);

                    int ii, jj;
                    double a1, a1j, a2, a2j, a3, a3j;
                    iss >> ii >> jj >> a1 >> a1j >> a2 >> a2j >> a3 >> a3j;

                    Rhop(i)(ii-1, jj-1, 0) = std::complex<double>(a1, a1j);
                    Rhop(i)(ii-1, jj-1, 1) = std::complex<double>(a2, a2j);
                    Rhop(i)(ii-1, jj-1, 2) = std::complex<double>(a3, a3j);
                    iss.clear(); 
                }
            }

            std::getline(m_file, line);  // Skip blank line if not last
        }


        //----------------------MOTIF-----------------------------//
        // Ensure motif is complex to handle the complex entries from Rhop
        motif = arma::mat(mSize, 4, arma::fill::zeros);
        int diag = 0;
        for (int i = 0; i < nFock; i++) {   
            if (iRn(i, 0) == 0 && iRn(i, 1) == 0 && iRn(i, 2) == 0) {
                diag = i;
                break;
            }
        }

        for (int k = 0; k < mSize; k++) {
            motif(k, 0) = std::real(Rhop(diag)(k, k, 0));  // Diagonal element for Rhop(x)
            motif(k, 1) = std::real(Rhop(diag)(k, k, 1));  // Diagonal element for Rhop(y)
            motif(k, 2) = std::real(Rhop(diag)(k, k, 2));  // Diagonal element for Rhop(z)
            motif(k, 3) = k; // species?
        }

        m_file.close();
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
    void Wannier90Configuration::mapContent(bool debug) {
        
        // Fill in the system info structure with relevant data

        systemInfo.ndim = ndim;                        // system dimension
        systemInfo.bravaisLattice = lattice;           // Store bravais lattice
        systemInfo.bravaisVectors = bravaisVectors;    // Store bravais nieghbors vectors
        systemInfo.motif = motif;                      // Store motif localization data
        systemInfo.hamiltonian = fockMatrices;         // Store Hamiltonian
        systemInfo.filling = filling;                  /* missing on w90 file */
        // systemInfo.extendedMotif = Rhop;               // plans ?

        // int norbitals_size = motif.n_rows;
        arma::urowvec norbitals = arma::zeros<arma::urowvec>(mSize);
        for (int i = 0; i < mSize; i++){
            // norbitals(i) = (std::floor(std::sqrt(mSize)))/4;
            norbitals(i) = 1;
        }
        systemInfo.norbitals      = norbitals;

        // Debug output
        if (false) {
            std::cout << "========================" << std::endl;
            std::cout << "||     PARSED VALUES  ||" << std::endl;
            std::cout << "========================" << std::endl;
            // std::cout << "Number of Fock matrices (nFock): " << systemInfo.nFock << "\n";
            std::cout << "Matrix size (mSize): " << mSize               << std::endl;

            std::cout << "No. of Fock matrices: "               << nFock    << std::endl;

            // std::cout << "-------------------"                          << std::endl;
            // std::cout << "bravaisLattice: "                             << std::endl;
            // std::cout << systemInfo.bravaisLattice                      << std::endl;
            // std::cout << "-------------------"                          << std::endl;

            std::cout << "Degen: "                             << std::endl;
            std::cout <<  Degen                                 << std::endl;
            std::cout << "-------------------"                          << std::endl;

            // std::cout << "-------------------"                          << std::endl;
            // std::cout << "bravaisVectors: "                             << std::endl;
            // std::cout << systemInfo.bravaisVectors                      << std::endl;
            // std::cout << "-------------------"                          << std::endl;

            std::cout << "System Dimension: "   << systemInfo.ndim      << std::endl;

            std::cout << "Filling: "            << systemInfo.filling   << std::endl;

            // std::cout << "-------------------"                          << std::endl;
            // std::cout << "Motif: "                                      << std::endl;
            // std::cout << systemInfo.motif                           << std::endl;
            // std::cout << "-------------------"                          << std::endl;

            // std::cout << "Rhop: "                                       << std::endl;
            // std::cout << Rhop                                           << std::endl;;
            // std::cout << "==================="                          << std::endl;;

            // std::cout << "iRn index: "                         << std::endl;
            // std::cout <<  iRn                         << std::endl;
            // std::cout << "========================" << std::endl;

            // std::cout << "Hamiltonian matrix: "                         << std::endl;
            // std::cout << systemInfo.hamiltonian                         << std::endl;
            // std::cout << "========================" << std::endl;
            std::cout << std::setprecision(15) << systemInfo.hamiltonian.slice(0)                         << std::endl;
        }
    }
}