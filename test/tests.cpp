#define CATCH_CONFIG_MAIN

#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include <string>
#include <vector>

#include <xatu.hpp>
#include <catch.hpp>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


TEST_CASE("Modelfile parsing", "[model_parsing]"){

    std::cout << std::setw(40) << std::left << "Testing model file parsing... ";
    std::cout.setstate(std::ios_base::failbit);
    
    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    int expectedFilling = 1;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {1, 1};

    arma::vec expectedBravaisLatticeHash = {4.55422, 2.05422};
    arma::vec expectedMotifHash = {0, 3.90211};
    arma::vec expectedBravaisVectorsHash = {0, -0.110843, 2.38916, 2.05422, 4.55422};
    arma::vec expectedHamiltionianHash = {-1.34375, -0.9, -0.9, -0.0375, -0.0375};
    
    REQUIRE(config.systemInfo.ndim == expectedDim);
    REQUIRE(config.systemInfo.filling == expectedFilling);
    for(uint i = 0; i < config.systemInfo.norbitals.n_cols; i++){
        REQUIRE(config.systemInfo.norbitals(i) == expectedOrbitals(i));
    }

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisLatticeHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedMotifHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisVectorsHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedHamiltionianHash(i), 1E-4));
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("CRYSTAL file parsing", "[CRYSTAL_parsing]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing CRYSTAL file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config = xatu::CRYSTALConfiguration(modelfile);

    int expectedFilling = 6;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {33, 36};
    uint expectedFockMatricesNumber = 20;

    arma::vec expectedBravaisLatticeHash = {2.05977, 4.51667};
    arma::vec expectedMotifHash = {0.068875, 2.84563};
    arma::vec expectedBravaisVectorsHash = {0, -0.113953, 4.56977, 2.39605, 2.05977,
                                            -0.503333, 4.51667, -0.39124, 3.9562, -0.95062,
                                            8.7531, 6.57938, 1.2231, -0.894573, 8.47287,
                                            4.12543, 3.45287, -1.34, 8.7, -0.838527};
    arma::vec expectedHamiltionianHash = {61.2007, 2.87756, -0.679479, 0.246344, -0.981652,
                                          0.447539, -1.41873, 0.32476, 0.408354, 0.596467,
                                          0.309267, 0.191303, 0.107262, 0.244018, 0.220732,
                                          0.165454, 0.248117, 0.32759, 0.169674, 0.0856007};
    arma::vec expectedOverlapHash = {46.302, -1.47274, -0.128825, 0.444089, 1.43755,
                                    -0.87285, 1.00188, 0.0256365, 0.0421825, -0.0733795,
                                    0.0628083, 0.113745, 0.340853, 0.0733232, 0.0869901,
                                    0.0867042, 0.151335, 0.0713888, 0.0795013, 0.0507182};
    
    REQUIRE(config.systemInfo.ndim == expectedDim);
    REQUIRE(config.systemInfo.filling == expectedFilling);
    REQUIRE(config.systemInfo.hamiltonian.n_slices == expectedFockMatricesNumber);
    REQUIRE(config.systemInfo.overlap.n_slices == expectedFockMatricesNumber);
    for(uint i = 0; i < config.systemInfo.norbitals.n_cols; i++){
        REQUIRE(config.systemInfo.norbitals(i) == expectedOrbitals(i));
    }

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisLatticeHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedMotifHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisVectorsHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedHamiltionianHash(i), 1E-2));
    } 
    for(uint i = 0; i < config.systemInfo.overlap.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.overlap.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedOverlapHash(i), 1E-2));
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("Wannier90 file parsing", "[wannier_parsing]") {
    // Suppress output during test
    std::cout << std::setw(40) << std::left << "Testing Wannier90 file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    // Load Wannier90 file
    std::string modelfile = "../examples/material_models/wannier/hBN_tb.dat";
    std::cout << "got file: " << modelfile << std::endl;
    xatu::Wannier90Configuration config = xatu::Wannier90Configuration(modelfile, 4);

    int expectedDim = 2;                         // Dimension of the system
    int expectedFilling = 4;                     // Number of filled states
    int expectedmSize = 6;                       // Number of unit cells
    int expectednFock = 83;                       // Number of Fock matrices
    
    // Hashes for real-space lattice vectors (flattened 3x3 matrix)
    arma::vec expectedBravaisLatticeHash = {2.425222, 4.080734};
    
    arma::vec expectedBravaisVectorsHash = {-2.971610, -2.456147, -1.940685, -3.584157, -3.068694, -2.553232, -2.037769, -1.522307,
    -1.340178,  2.407223, -4.014574, -3.165779, -2.650316, -2.134854, -1.619392, -1.103929,
    -0.921800,  2.825601,  6.239668, -1.910645, -2.565272, -1.716476, -1.201014, -0.685551,
    -0.503422,  3.243979,  6.658046, 10.072114,  0.181244, -0.140049, -0.461343, -1.115969,
    -0.267174, -0.085044,  3.662356,  7.076424, 10.490492,  1.951840,  1.630546,  1.309253,
     0.987960,  0.000000,  4.080734,  7.494802, 10.908869, 14.322937,  3.722436,  3.401142,
     3.079849,  2.425222,  5.335868,  7.579846, 11.327247, 14.741315, 18.155382,  5.814325,
     5.493031,  5.171738,  4.517112,  7.427757, 10.005069, 12.582381, 14.826359, 18.573760,
     7.584920,  7.263627,  6.609001,  9.519646, 12.096958, 14.674270, 17.251582, 19.828894,
    22.072872,  9.355516,  8.700890, 11.611535, 14.188847, 16.766159, 19.343471, 21.920783,
    13.703424, 16.280736, 18.858048
};

    // Hashes for motif positions (flattened matrix)
    arma::vec expectedMotifHash = {2.566478, 5.266827, 5.661396, 10.201080, 13.255683, 14.551132};

    // arma::vec expectedDegenHash = {1, 1, 2, 1, 1};

    // Hashes for Fock matrices (real/imag parts handled separately)
    arma::vec expectedHamiltonianHash = {0.978634, 0.990920, 1.046875, 1.054984, 0.996068, 0.968801, 1.101170, 0.961007,
        0.968796, 1.043775, 1.032634, 1.017193, 0.940727, 1.172638, 0.890277, 1.110180,
        1.088959, 1.027729, 0.962491, 0.970735, 1.036763, 1.089453, 1.079106, -0.531977,
        0.946280, 0.965448, 1.031659, 1.011340, 1.016294, 1.009594, 0.999775, 0.951863,
        0.458796, 0.579032, 1.125541, 1.059755, 0.985202, 0.986659, 1.042266, 1.042844,
        2.030694, -11.551247, 2.313774, 1.014736, 1.049511, 0.994099, 0.997162, 1.026182,
        1.131210, 2.914045, 6.239713, 0.937042, 1.013736, 0.995272, 1.019086, 0.998125,
        1.032483, 0.973986, 1.026194, 1.357262, 1.224249, 1.032389, 1.027929, 0.987408,
        0.982629, 1.003563, 1.071731, 1.044472, 0.880948, 1.118926, 0.963565, 0.998474,
        1.029212, 1.046875, 0.982783, 0.964047, 1.057069, 0.999934, 0.973145, 1.045742,
        1.043775, 1.005556, 0.966527
    };
    
    // Hashes for Rhop matrices (if used)
    arma::vec expectedRhopHash = {0.997114, 0.990306, 0.997694, 0.991232, 1.005100, 1.004297, 0.989814, 0.977408,
    1.005649, 0.983461, 1.006747, 1.007933, 1.012504, 0.994354, 0.991596, 1.013520,
    1.007218, 1.005119, 0.996115, 0.995221, 1.012035, 1.011999, 1.006950, 0.997909,
    1.000703, 0.989825, 1.002780, 1.002820, 0.985856, 1.007615, 0.986269, 0.993846,
    1.007167, 1.002167, 0.989714, 0.993715, 1.014779, 1.004536, 0.999528, 1.020440,
    1.012040, 1.008106, 0.994788, 1.035806, 0.989133, 0.991763, 0.972419, 0.991202,
    1.022495, 1.038578, 1.006076, 1.005136, 1.000803, 0.992482, 0.990968, 1.014487,
    1.003521, 1.004732, 0.985285, 1.003408, 0.988708, 0.953025, 1.004636, 1.021256,
    1.001702, 0.998779, 1.022129, 0.971661, 1.021564, 1.088059, 1.024697, 1.026738,
    1.006775, 0.963720, 0.993811, 1.026752, 1.029950, 1.005488, 0.957017, 1.021292,
    1.008024, 0.997984, 1.019589, 0.994519, 1.053920, 1.020668, 1.007176, 0.985443,
    1.008066, 0.993171, 1.007998, 0.976088, 1.000136, 0.992240, 0.992790, 1.035847,
    1.024869, 1.066353, 1.095552, 1.019990, 1.021736, 1.087952, 0.970075, 0.995309,
    1.004839, 1.015675, 1.005783, 0.995738, 0.980129, 0.997494, 1.004768, 0.983212,
    0.996325, 1.002814, 1.005364, 1.001616, 1.002397, 0.954128, 1.003217, 1.018674,
    1.085609, 1.078255, 1.012354, 4.548316, 2.873857, 1.197828, 1.115675, 1.006616,
    0.964120, 0.984689, 0.978937, 1.024818, 1.025584, 1.003606, 1.005235, 1.017286,
    0.988046, 0.995220, 0.996773, 0.999681, 0.996953, 1.013851, 0.991732, 1.009270,
    0.957047, 1.006585, 0.995265, 0.992805, 0.997985, 1.040597, 0.972494, 0.999376,
    1.115256, 1.006251, 0.996555, 1.018850, 1.021810, 0.972551, 0.998658, 0.984656,
    1.009052, 1.005229, 1.040141, 1.028630, 1.004017, 1.007448, 1.002345, 1.001194,
    0.991527, 1.005624, 1.007665, 1.025050, 1.006177, 1.011064, 0.970850, 0.992055,
    0.981971, 1.049555, 1.007129, 1.029757, 1.030265, 0.961084, 0.987074, 1.007111,
    0.998259, 1.021026, 0.994860, 1.013045, 1.006976, 0.999025, 1.012405, 0.992590,
    1.009446, 0.996579, 0.992565, 1.007522, 0.998227, 1.001843, 1.009624, 1.028827,
    1.022804, 1.022722, 1.002362, 0.977078, 0.983739, 1.012481, 0.999343, 1.005116,
    1.003181, 1.022401, 0.997381, 0.990140, 0.993093, 1.004442, 1.016081, 0.998645,
    0.993441, 1.001894, 1.008756, 0.989814, 0.977408, 1.005649, 1.006771, 1.007816,
    0.996179, 0.981866, 1.020622, 1.012108, 1.056906, 1.033585, 1.020818, 1.011852,
    0.997104, 1.001710, 1.015061, 0.997085, 1.002824, 0.989569, 1.023114, 1.005145,
    1.002820, 0.985856, 1.007615, 1.010556, 1.008685, 0.990653, 1.017227, 1.020576,
    1.002143
    };

    // --- Validation ---
    REQUIRE(config.ndim == expectedDim);
    REQUIRE(config.filling == expectedFilling);
    REQUIRE(config.mSize == expectedmSize);
    REQUIRE(config.nFock == expectednFock);

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisLatticeHash(i), 1E-4));
    }

    // double degenHash = xatu::array2hash(arma::vectorise(config.degeneracies));
    // REQUIRE_THAT(degenHash, Catch::Matchers::WithinAbs(expectedDegenHash), 1E-4);

    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedMotifHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisVectorsHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedHamiltonianHash(i), 1E-2));
    } 

    uint hashIndex = 0;
    for(uint i = 0; i < config.Rhop.n_elem; i++) {
        for(uint s = 0; s < config.Rhop(i).n_slices; s++) {
            arma::cx_mat slice = config.Rhop(i).slice(s);
            double sliceHash = xatu::array2hash(slice);
            REQUIRE_THAT(sliceHash, Catch::Matchers::WithinAbs(expectedRhopHash(hashIndex), 1E-2));
            hashIndex++;
        }
    }

    // Restore output
    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("Exciton file parsing", "[excitonfile_parsing]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing exciton file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string excitonconfig = "./data/hBN_spinless.txt";
    xatu::ExcitonConfiguration config = xatu::ExcitonConfiguration(excitonconfig);

    std::string expectedName = "hBN_N30";
    int expectedNcell = 30;
    int expectedNbands = 1;
    arma::vec expectedDielectric = {1, 1, 10};

    REQUIRE(config.excitonInfo.ncell == expectedNcell);
    REQUIRE(config.excitonInfo.nbands == expectedNbands);
    for (uint i = 0; i < expectedDielectric.n_elem; i++){
        REQUIRE(config.excitonInfo.eps(i) == expectedDielectric(i));
    }
    REQUIRE(xatu::array2hash(config.excitonInfo.Q) == 0);

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("Exciton file parsing - anisotropic", "[excitonfile_parsing_anisotropic]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing anisotropic exciton file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string excitonconfig = "./data/hBN_anisotropic.txt";
    xatu::ExcitonConfiguration config = xatu::ExcitonConfiguration(excitonconfig);

    std::string expectedName = "hBN_ani";
    int expectedNcell = 30;
    int expectedNbands = 1;
    arma::vec expectedDielectric = {1, 1, 10, 25, 35};

    REQUIRE(config.excitonInfo.ncell == expectedNcell);
    REQUIRE(config.excitonInfo.nbands == expectedNbands);
    for (uint i = 0; i < expectedDielectric.n_elem; i++){
        REQUIRE(config.excitonInfo.eps(i) == expectedDielectric(i));
    }
    REQUIRE(xatu::array2hash(config.excitonInfo.Q) == 0);

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}
    
TEST_CASE("TB hBN energies (full diagonalization)", "[tb-hBN-fulldiag]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (fulldiag)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (davidson)", "[tb-hBN-davidson]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (Davidson)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("davidson", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (Lanczos)", "[hBN-lanczos]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (Lanczos)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("sparse", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (reciprocal)", "[tb-hBN-reciprocal]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (reciprocal)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.setMode("reciprocalspace");
    exciton.setReciprocalVectors(5);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{6.234291, 1}, 
                                                         {6.236636, 1},
                                                         {6.731819, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN reciprocal w.f.", "[tb-hBN-kwf]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN reciprocal w.f... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }

    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.086145105;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN real-space w.f.", "[tb-hBN-rswf]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN real-space w.f... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::rowvec holePosition = exciton.system->motif.row(holeIndex).subvec(0, 2) + holeCell;

    double radius = arma::norm(exciton.system->bravaisLattice.row(0)) * exciton.ncell;
    arma::mat cellCombinations = exciton.system->truncateSupercell(exciton.ncell, radius);
    arma::vec rswf = arma::zeros(cellCombinations.n_rows*exciton.system->motif.n_rows);

    // Compute probabilities
    for(int n = 0; n < nstates; n++){
        int it = 0;
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.system->motif.n_rows; atomIndex++){
            rswf(it) += results->realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
        }
    }

    double rswfHash = xatu::array2hash(rswf);
    double expectedRSwfHash = 3.4185866144;
    REQUIRE_THAT(rswfHash, Catch::Matchers::WithinAbs(expectedRSwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN absorption", "[tb-hBN-kubo]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN conductivity... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 47.4140063784;
    REQUIRE_THAT(cum_norm_vme_ex, Catch::Matchers::WithinAbs(expectedTotalOscillator, 1E-7));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (spinful)", "[tb-hBN-spinful]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN spinful... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 12;

    std::string modelfile = "../examples/material_models/hBN_spinful.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 8}, 
                                                         {6.074062, 4}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}
    
TEST_CASE("TB hBN spin", "[tb-hBN-spin]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN spin... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 16;
    int nstates = 4;

    std::string modelfile = "../examples/material_models/hBN_spinful.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1, 1, 10});
    
    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::mat expectedSpin = arma::vec{-1, 0, 0, 1};  

    for(uint i = 0; i < nstates; i++){
        arma::cx_vec spin = results->spinX(i);
        REQUIRE(spin.n_elem == 3);
        REQUIRE_THAT(std::real(spin(0)), Catch::Matchers::WithinAbs(expectedSpin(i), 1E-2));
        REQUIRE_THAT(std::abs(spin(1)), Catch::Matchers::WithinAbs(0.5, 1E-2));
        REQUIRE_THAT(std::abs(spin(2)), Catch::Matchers::WithinAbs(0.5, 1E-2));
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("DFT hBN", "[dft-hBN]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing DFT hBN... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;
    bool triangular = true;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config = xatu::CRYSTALConfiguration(modelfile, 100);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.system->setAU(true);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{4.442317, 1}, 
                                                         {4.442427, 1}};
                                                         
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    // Check reciprocal w.f.
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }

    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.0864895707;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    // Check realspace w.f.
    arma::rowvec holePosition = exciton.system->motif.row(holeIndex).subvec(0, 2) + holeCell;

    double radius = arma::norm(exciton.system->bravaisLattice.row(0)) * exciton.ncell;
    arma::mat cellCombinations = exciton.system->truncateSupercell(exciton.ncell, 2);
    arma::vec rswf = arma::zeros(cellCombinations.n_rows*exciton.system->motif.n_rows);

    // Compute probabilities
    for(int n = 0; n < nstates; n++){
        int it = 0;
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.system->motif.n_rows; atomIndex++){
            rswf(it) += results->realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
        }
    }

    double rswfHash = xatu::array2hash(rswf);
    double expectedRSwfHash = 83.2242560463;
    REQUIRE_THAT(rswfHash, Catch::Matchers::WithinAbs(expectedRSwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("Wannier hBN", "[w90-hBN]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing Wannierized hBN... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 30;
    int nstates = 2;
    bool triangular = true;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/wannier/hBN_tb.dat";
    xatu::Wannier90Configuration config = xatu::Wannier90Configuration(modelfile, 4);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.system->setAU(true);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{2.833347, 1}, 
                                                         {2.834596, 1}};
                                                         
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    // check oscillator strength for absorption
    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 103.0312591808;
    REQUIRE_THAT(cum_norm_vme_ex, Catch::Matchers::WithinAbs(expectedTotalOscillator, 1E-7));

    // Check reciprocal w.f.
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }


    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.057666;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    // Check realspace w.f.
    arma::rowvec holePosition = exciton.system->motif.row(holeIndex).subvec(0, 2) + holeCell;

    double radius = arma::norm(exciton.system->bravaisLattice.row(0)) * exciton.ncell;
    arma::mat cellCombinations = exciton.system->truncateSupercell(exciton.ncell, 2);
    arma::vec rswf = arma::zeros(cellCombinations.n_rows*exciton.system->motif.n_rows);

    // Compute probabilities
    for(int n = 0; n < nstates; n++){
        int it = 0;
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.system->motif.n_rows; atomIndex++){
            rswf(it) += results->realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
        }
    }

    double rswfHash = xatu::array2hash(rswf);
    double expectedRSwfHash = 1.003540;
    REQUIRE_THAT(rswfHash, Catch::Matchers::WithinAbs(expectedRSwfHash, 1E-5));

    std::cout.clear();
    std::cout << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN-ani energies (full diagonalization)", "[tb-hBN-ani-fulldiag]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN-ani energies (fulldiag)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10, 25, 35});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{6.428535, 1}, 
                                                         {6.456359, 1},
                                                         {6.731320, 1}};

    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN absorption - anisotropic", "[tb-hBN-kubo-ani]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN anisotropic conductivity... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10, 25, 35});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 47.4140063784;
    REQUIRE_THAT(cum_norm_vme_ex, Catch::Matchers::WithinAbs(expectedTotalOscillator, 1E-7));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 energies", "[MoS2-energies]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 energies... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.768783, 2}, 
                                                         {1.780562, 2}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 reciprocal w.f.", "[MoS2-kwf]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 reciprocal w.f... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }

    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.1814790902;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 spin", "[MoS2-spin]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 spin... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    int factor = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.system->reducedBrillouinZoneMesh(ncell, factor);
    exciton.system->shiftBZ({0.6628, -1.1480, 0});

    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::mat expectedSpin = arma::mat{{-9.915e-01, -5.000e-01, -4.915e-01},
                                       {-6.966e-04, -5.000e-01,  4.993e-01},
                                       { 8.297e-03,  4.998e-01, -4.915e-01},
                                       { 9.991e-01,  4.998e-01,  4.993e-01}}; 
                              
    for(uint i = 0; i < nstates; i++){
        arma::cx_vec spin = results->spinX(i);
        for(uint j = 0; j < 3; j++){
            double spinValue = real(spin(j));
            REQUIRE_THAT(spinValue, Catch::Matchers::WithinAbs(expectedSpin(i, j), 1E-4));
        }
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 exciton bands", "[MoS2-Q]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 exciton bands... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);
    
    arma::vec Q_values = {-0.1, -0.05, 0.0, 0.05, 0.1};
    arma::rowvec Q = {0., 0., 0.};
    std::vector<std::vector<double>> energies;

    for (const auto& q : Q_values){
        Q(1) = q;
        xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55}, Q);

        exciton.brillouinZoneMesh(ncell);
        exciton.initializeHamiltonian();
        exciton.BShamiltonian();

        auto results = exciton.diagonalize("diag", nstates);
        auto resultingEnergies = xatu::detectDegeneracies(results->eigval, nstates, 6);

        energies.push_back(resultingEnergies[0]);
    }
    
    std::vector<std::vector<double>> expectedEnergies = {{1.810464, 2}, 
                                                         {1.780106, 2},
                                                         {1.768783, 2},
                                                         {1.780106, 2},
                                                         {1.810464, 2}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-5));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 exchange", "[MoS2-exchange]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 with exchange... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    arma::rowvec Q = {0., 0.1, 0.};

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55}, Q);
    exciton.setExchange(true);
    
    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();

    auto results = exciton.diagonalize("diag", nstates);
    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.810464, 2}, 
                                                         {1.825740, 1},
                                                         {1.855544, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-5));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 reduced BZ", "[MoS2-reducedBZ]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 with reduced BZ... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    int factor = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.system->reducedBrillouinZoneMesh(ncell, factor);
    exciton.system->shiftBZ({0.6628, -1.1480, 0});

    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.795378, 1}, 
                                                         {1.809810, 1},
                                                         {1.935880, 1},
                                                         {1.950567, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 absorption", "[MoS2-kubo]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 conductivity... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 37.2792093358;
    REQUIRE_THAT(cum_norm_vme_ex, Catch::Matchers::WithinAbs(expectedTotalOscillator, 1E-7));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

