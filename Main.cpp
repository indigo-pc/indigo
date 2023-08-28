#include "Spectrum.cpp"
#include <iostream>
#include <vector>

int main() {

// uWatt per cm^2: 210714
// uWatt: 421429
// Joules: 1.26429
// eV: 7.89193e+18
// Lux: 2.61276
// Candella: 52255.1
// Photons per cm^2 per s: 6.11848e+17
// Total photons: 3.67109e+18

    std::vector<double> photopicEfficiencyData = { 0.027, 0.082, 0.082, 0.27, 0.826, 2.732, 7.923, 15.709, 25.954, 40.98, 62.139, 94.951, 142.078, 220.609, 303.464, 343.549, 484.93, 588.746, 651.582, 679.551, 683, 679.585, 650.216, 594.21, 517.031, 430.973, 343.549, 260.223, 180.995, 119.525, 73.081, 41.663, 21.856, 11.611, 5.607, 2.802, 1.428, 0.715, 0.355, 0.17, 0.082, 0.041, 0.02 };
    double peakPhotopicLuminosity = 683.0;

    std::vector<double> wavelengths = { 380, 390, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 507, 510, 520, 530, 540, 550, 555, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770 };
    std::vector<double> data = { 79, 363, 50, 821, 984, 948, 694, 855, 214, 197, 695, 6, 402, 639, 101, 242, 672, 590, 808, 459, 853, 971, 264, 293, 443, 608, 453, 17, 923, 108, 533, 874, 633, 773, 51, 723, 747, 623, 723, 545, 999, 116, 481 };
    Spectrum s { data, wavelengths };
    Spectrum photopic { photopicEfficiencyData, wavelengths };

    double collectionArea = 2;
    double solidAngle = 1;
    double integrationTime = 3;

    std::cout << "uWatt per cm^2: " << s.computeMicrowattsPerCmSquared( s.lowerBound, s.upperBound ) << "\n";
    std::cout << "uWatt: " << s.computeMicrowatts( s.lowerBound, s.upperBound, collectionArea ) << "\n";
    std::cout << "Joules: " << s.computeJoules( s.lowerBound, s.upperBound, collectionArea, integrationTime ) << "\n";
    std::cout << "eV: " << s.computeElectronVolts( s.lowerBound, s.upperBound, collectionArea, integrationTime ) << "\n";
    std::cout << "Lux: " << s.computeLuxPhotopic( ) << "\n";
    std::cout << "Candella: " << s.computeCandella( photopic, peakPhotopicLuminosity, collectionArea, solidAngle ) << "\n";
    std::cout << "Photons per cm^2 per s: " << s.computePhotonsPerCmSquaredPerSecond( s.lowerBound, s.upperBound ) << "\n";
    std::cout << "Total photons: " << s.computeTotalPhotons( s.lowerBound, s.upperBound, collectionArea, integrationTime ) << "\n";

    return 0;
}



