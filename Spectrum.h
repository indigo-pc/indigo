#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include <tuple>

class Spectrum {

    public:

        // constructors

        explicit Spectrum( std::vector<double> data, std::vector<double> wavelengths );

        // getters & setters

        std::vector<double> getSpectralData();
        void setSpectralData( const int & i, const double & d );
        std::vector<double> getSpectralWavelengths();
        std::tuple<int, int> getWavelengthIndices( const double & lowerWavelength, const double & upperWavelength );

        // peaks

        std::vector<int> getPeakIndices( const double & tolerance );
        std::vector<double> getPeakAmplitudes( const double & tolerance );
        double peakFullWidthHalfMax( const int & index );
        double peakFullWidthHalfMax( const double & wavelength );
        bool isSaturated();

        // radiometry and photometry

        double computeMicrowattsPerCmSquared( int & lowerBoundWavelength, int & upperBoundWavelength );
        double computeMicrowatts( int & lowerBoundWavelength, int & upperBoundWavelength, double & collectionArea );
        Spectrum computeMicrowattsPerNanometer( double & collectionArea );

        double computeJoules( int & lowerBoundWavelength, int & upperBoundWavelength, double & collectionArea, double & integrationTime );
        double computeElectronVolts( int & lowerBoundWavelength, int & upperBoundWavelength, double & collectionArea, double & integrationTime );

        double computeLux( Spectrum & luminousEfficiencyFunction, double & maxLuminousEfficiencyCoefficient );
        double computeLuxPhotopic( );
        double computeLuxScotopic( );
        double computeLumens( Spectrum & luminousEfficiencyFunction, double & maxLuminousEfficiencyCoefficient, double & collectionArea );

        double computeCandella( Spectrum & luminousEfficiencyFunction, double & maxLuminousEfficiencyCoefficient, double & collectionArea, double & solidAngle );

        double computePhotonsPerCmSquaredPerSecond( int & lowerBoundWavelength, int & upperBoundWavelength );
        double computeTotalPhotons( int & lowerBoundWavelength, int & upperBoundWavelength, double & collectionArea, double & integrationTime );

        // Spectral Objects

        Spectrum splice( Spectrum & s );
        void crop( double & lowerWavelength, double & upperWavelength );
        void crop( int & lowerWavelengthIndex, int & upperWavelengthIndex );
        std::vector<double> unitlessSpectralCalibration( std::vector<double> & idealReferenceData );
        Spectrum boxCarAveraged( const int & pixels );

        // load vectors from file

        static std::vector<double> loadFromFile( const std::string & filePath, const int & headerRowCount );

        // maths

        double computeIntegral( double & lowerWavelength, double & upperWavelength );
        double computeIntegral( int & lowerWavelengthIndex, int & upperWavelengthIndex );
        double getMaxValue( const std::vector<double> & v );
        double getMinValue( const std::vector<double> & v );
        double percentDifference( const double & v1, const double & v2 );

        // constants

        int lowerBound = 0;
        int upperBound = -1;

    private:

        std::vector<double> spectralData;
        std::vector<double> spectralWavelengths;

        // Photometry constants courtesy of Williamson & Cummins, Light and Color in Nature and Art, Wiley, 1983. Accessed through http://hyperphysics.phy-astr.gsu.edu/hbase/vision/efficacy.html
        std::vector<double> photopicEfficiencyData = { 0.027, 0.082, 0.082, 0.27, 0.826, 2.732, 7.923, 15.709, 25.954, 40.98, 62.139, 94.951, 142.078, 220.609, 303.464, 343.549, 484.93, 588.746, 651.582, 679.551, 683, 679.585, 650.216, 594.21, 517.031, 430.973, 343.549, 260.223, 180.995, 119.525, 73.081, 41.663, 21.856, 11.611, 5.607, 2.802, 1.428, 0.715, 0.355, 0.17, 0.082, 0.041, 0.02 };
        std::vector<double> scotopicEfficiencyData = { 1.001, 3.755, 3.755, 15.793, 59.228, 164.22, 339.66, 557.77, 773.5, 963.9, 1149.2, 1348.1, 1536.8, 1669.4, 1700, 1694.9, 1589.5, 1378.7, 1105, 817.7, 683, 558.96, 352.92, 206.04, 111.35, 56.355, 27.081, 12.529, 5.67, 2.545, 1.151, 0.532, 0.252, 0.122, 0.06, 0.03, 0.016, 0.008, 0.004, 0.002, 0.001, 0, 0 };
        std::vector<double> efficiencyFunctionWavelengths = { 380, 390, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 507, 510, 520, 530, 540, 550, 555, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770 };
        double peakPhotopicLuminosity = 683.0;
        double peakScotopicLuminosity = 1700.0;

        double electricChargeElectron = 1.602e-19;
};

// whole-spectrum vector arithmetic

Spectrum operator*=( Spectrum s, std::vector<double> & v );
Spectrum operator*=( std::vector<double> & v, Spectrum s );
Spectrum operator*=( Spectrum s1, Spectrum s2 );
Spectrum operator/=( Spectrum s, std::vector<double> & v );

// spectral constant vector arithmetic

Spectrum operator*( Spectrum s, double & t );
Spectrum operator*( double & t, Spectrum s );
std::vector<double> operator*( std::vector<double> & v, double & t );
std::vector<double> operator*( double & t, std::vector<double> & v );

// other operators

std::ostream & operator<<( std::ostream & out, Spectrum & s );

#endif
