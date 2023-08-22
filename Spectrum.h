#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include <tuple>

class Spectrum {

    public:

        // constructors

        explicit Spectrum( std::vector<double> data, std::vector<double> wavelengths );
        explicit Spectrum( );

        // getters & setters

        std::vector<double> getSpectralData();
        void setSpectralData( int i, double d );
        std::vector<double> getSpectralWavelengths();
        std::tuple<int, int> getWavelengthIndices( double lowerWavelength, double upperWavelength );

        // peaks

        std::vector<int> getPeakIndices( double tolerance );
        std::vector<double> getPeakAmplitudes( double tolerance );
        double peakFullWidthHalfMax( int index );
        double peakFullWidthHalfMax( double wavelength );
        bool isSaturated();

        // maths

        double computeIntegral( double lowerWavelength, double upperWavelength );
        double computeIntegral( int lowerWavelengthIndex, int upperWavelengthIndex );
        double getMaxValue( std::vector<double> v );
        double getMinValue( std::vector<double> v );
        double percentDifference( double v1, double v2 );

        // advanced photometrics otherwise available in Java from Ocean Insight
        // https://www.oceaninsight.com/globalassets/catalog-blocks-and-images/software-downloads-installers/javadocs-api/spam/com/oceanoptics/spam/advancedprocessing/advancedphotometrics.html
        // https://www.oceaninsight.com/globalassets/catalog-blocks-and-images/software-downloads-installers/javadocs-api/spam/index.html
        // TODO implement 'Advanced Color' from https://www.oceaninsight.com/globalassets/catalog-blocks-and-images/software-downloads-installers/javadocs-api/spam/index.html
        // TODO implement 'XYZ Color' from https://www.oceaninsight.com/globalassets/catalog-blocks-and-images/software-downloads-installers/javadocs-api/spam/index.html
        // also see: https://wp.optics.arizona.edu/jpalmer/radiometry/radiometry-and-photometry-faq/

        // radiometry and photometry

        double computePower( int lowerBoundWavelength, int upperBoundWavelength );
        double computePowerPerSquareArea( int lowerBoundWavelength, int upperBoundWavelength, double collectionArea );
        Spectrum computePowerPerSquareAreaPerWavelength( double collectionArea );
        
        double computeEnergy( int lowerBoundWavelength, int upperBoundWavelength, double integrationTime );
        double computeEnergyPerSquareArea( int lowerBoundWavelength, int upperBoundWavelength, double collectionArea, double integrationTime );
        Spectrum computeEnergyPerSquareAreaPerWavelength( double collectionArea, double integrationTime );

        Spectrum computeIlluminancePerSquareAreaPerWavelength( double collectionArea, Spectrum luminousEfficiencyFunction, double maxLuminousEfficiencyCoefficient );
        Spectrum computeIlluminancePerSquareAreaPerWavelength( double collectionArea );

        double computeLuminanceCandelaPerSquareMeter( std::vector<double> wavelengths, std::vector<double> energyWattsPerNanometer, std::vector<double> V_wavelengths, std::vector<double> V, double K_m, double steradians, double areaSquareMeters );
        double computeLuminousIntensityCandela( std::vector<double> wavelengths, std::vector<double> energyWattsPerNanometer, std::vector<double> V_wavelengths, std::vector<double> V, double K_m, double steradians );

        double computeTotalPhotons( std::vector<double> wavelengths, std::vector<double> uWPerCmSquaredPerNm, double startingWavelength, double endingWavelength, int integrationMethod, double integrationTimeSeconds, double surfaceAreaCmSquared );
        double computePhotonsPerCmSquaredPerSecond( std::vector<double> wavelengths, std::vector<double> uWPerCmSquaredPerNm, double startingWavelength, double endingWavelength, int integrationMethod );

        // Spectral Objects

        Spectrum splice( Spectrum s );
        void crop( double lowerWavelength, double upperWavelength );
        void crop( int lowerWavelengthIndex, int upperWavelengthIndex );
        std::vector<double> unitlessSpectralCalibration( std::vector<double> idealReferenceData );
        Spectrum boxCarAveraged( int pixels );

        // load vectors from file

        static std::vector<double> loadFromFile( std::string const& filePath, int headerRowCount );

    private:

        std::vector<double> spectralData;
        std::vector<double> spectralWavelengths;

        // data fields for photometry calculations in photopic regime
        std::vector<double> photopicEfficiencyData = { 0.027, 0.082, 0.082, 0.27, 0.826, 2.732, 7.923, 15.709, 25.954, 40.98, 62.139, 94.951, 142.078, 220.609, 303.464, 343.549, 484.93, 588.746, 651.582, 679.551, 683, 679.585, 650.216, 594.21, 517.031, 430.973, 343.549, 260.223, 180.995, 119.525, 73.081, 41.663, 21.856, 11.611, 5.607, 2.802, 1.428, 0.715, 0.355, 0.17, 0.082, 0.041, 0.02 };
        std:: vector<double> photopicEfficiencyWavelengths = { 380, 390, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 507, 510, 520, 530, 540, 550, 555, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770 };
        double peakLuminosityNormalization = 683.002;       // units: [ lm / [power] ]

};

// whole-spectrum vector arithmetic

Spectrum operator*=( Spectrum s, std::vector<double> v );
Spectrum operator*=( std::vector<double> v, Spectrum s );
Spectrum operator*=( Spectrum s1, Spectrum s2 );
Spectrum operator/=( Spectrum s, std::vector<double> v );

// spectral constant vector arithmetic

Spectrum operator*( Spectrum s, double t );
Spectrum operator*( double t, Spectrum s );

// other operators

std::ostream & operator<<( std::ostream & out, Spectrum & s );

#endif
