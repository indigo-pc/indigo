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

    std::vector<double> data = { 79, 363, 50, 821, 984, 948, 694, 855, 214, 197, 695, 6, 402, 639, 101, 242, 672, 590, 808, 459, 853, 971, 264, 293, 443, 608, 453, 17, 923, 108, 533, 874, 633, 773, 51, 723, 747, 623, 723, 545, 999, 116, 481 };
    Spectrum s { data, Spectrum::efficiencyFunctionWavelengths };
    Spectrum photopic { Spectrum::photopicEfficiencyData, Spectrum::efficiencyFunctionWavelengths };

    double collectionArea = 2;
    double solidAngle = 1;
    double integrationTime = 3;

    std::cout << "uWatt per cm^2: " << s.computeMicrowattsPerCmSquared( Spectrum::lowerBound, Spectrum::upperBound ) << "\n";
    std::cout << "uWatt: " << s.computeMicrowatts( Spectrum::lowerBound, Spectrum::upperBound, collectionArea ) << "\n";
    std::cout << "Joules: " << s.computeJoules( Spectrum::lowerBound, Spectrum::upperBound, collectionArea, integrationTime ) << "\n";
    std::cout << "eV: " << s.computeElectronVolts( Spectrum::lowerBound, Spectrum::upperBound, collectionArea, integrationTime ) << "\n";
    std::cout << "Lux: " << s.computeLuxPhotopic( ) << "\n";
    std::cout << "Candella: " << s.computeCandella( photopic, Spectrum::peakPhotopicLuminosity, collectionArea, solidAngle ) << "\n";
    std::cout << "Photons per cm^2 per s: " << s.computePhotonsPerCmSquaredPerSecond( Spectrum::lowerBound, Spectrum::upperBound ) << "\n";
    std::cout << "Total photons: " << s.computeTotalPhotons( Spectrum::lowerBound, Spectrum::upperBound, collectionArea, integrationTime ) << "\n";

    return 0;
}



