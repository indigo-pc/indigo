#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include <tuple>

class Spectrum {

    public:

        // constructors

        explicit Spectrum( std::vector<double> data, std::vector<double> wavelengths, double integrationTime );
        explicit Spectrum( );

        // getters & setters

        std::vector<double> getSpectralData();
        void setSpectralData( int i, double d );
        std::vector<double> getSpectralWavelengths();
        double getIntegrationTime();
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

        // radiometry

        double computePower( int lowerBoundWavelength, int upperBoundWavelength );
        double computePowerPerSquareArea( int lowerBoundWavelength, int upperBoundWavelength, double collectionArea );
        Spectrum computePowerPerSquareAreaPerWavelength( double collectionArea );
        double computeEnergy( int lowerBoundWavelength, int upperBoundWavelength );
        double computeEnergyPerSquareArea( int lowerBoundWavelength, int upperBoundWavelength, double collectionArea );
        Spectrum computeEnergyPerSquareAreaPerWavelength( double collectionArea );

        double computeIlluminanceLux( std::vector<double> wavelengths, std::vector<double> energyWattsPerNanometer, std::vector<double> V_wavelengths, std::vector<double> V, double K_m, double areaSquareMeters );
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
        double integrationTime;

};

// whole-spectrum vector arithmetic

Spectrum operator*=( Spectrum s, std::vector<double> v );
Spectrum operator*=( std::vector<double> v, Spectrum s );
Spectrum operator/=( Spectrum s, std::vector<double> v );

// spectral index-wise vector arithmetic

Spectrum operator*( Spectrum s, double t );
Spectrum operator*( double t, Spectrum s );

// other operators

std::ostream & operator<<( std::ostream & out, Spectrum & s );

#endif
