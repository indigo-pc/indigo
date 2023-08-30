#include "Spectrum.h"

#include <vector>
#include <system_error>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>
#include <tuple>
#include <ostream>
#include <cmath>

//
// constructors
//

/**
 * Radiometrc and photometric calculations assume #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
*/
Spectrum::Spectrum( std::vector<double> data,
                    std::vector<double> wavelengths)
    : spectralData(std::move(data)), spectralWavelengths(std::move(wavelengths)) {
    if ( spectralData.size() != spectralWavelengths.size() )  {
        throw std::runtime_error( "Data and wavelength vectors must have matching length." );
    }
    for ( int i = 0; i < getSpectralWavelengths().size()-1; i++ )  {
        if ( ! ( getSpectralWavelengths().at(i) <= getSpectralWavelengths().at(i+1) ) )  { 
            throw std::runtime_error( "Spectral wavelengths must be a vector of unique values in increasing order" );
        }
        if ( getSpectralWavelengths().at(i) < 0 ) {
            throw std::runtime_error( "Spectral wavelength values must be positive." );
        }
    }
}

//
// getters & setters
//

std::vector<double> Spectrum::getSpectralData() { return this -> spectralData; }

void Spectrum::setSpectralData( const int & i, const double & d ) { ( this -> spectralData[i] ) = d; }

std::vector<double> Spectrum::getSpectralWavelengths() { return this -> spectralWavelengths; }

std::tuple<int, int> Spectrum::getWavelengthIndices( const double & lowerWavelength, const double & upperWavelength ) {
    if ( lowerWavelength >= upperWavelength ) {
        throw std::runtime_error( "Lower wavelength argument must not be equal or larger than upper." );
    }
    int lowerWavelengthIndex;
    int upperWavelengthIndex;
    for ( int i = 0; i < getSpectralWavelengths().size(); i++ )  {
        if ( getSpectralWavelengths().at(i) == lowerWavelength )  {
            lowerWavelengthIndex = i; continue;
        }
        if ( getSpectralWavelengths().at(i) == upperWavelength ) {
            upperWavelengthIndex = i;
        }
    }
    if ( lowerWavelengthIndex >= upperWavelengthIndex )  {
        throw std::runtime_error( "Lower wavelength argument is equal or larger than upper." );
    }
    return  { lowerWavelengthIndex, upperWavelengthIndex };
}

//
// peaks
//

/**
 * Return vector of peaks (amplitudes) within some normalized percent difference tolerance.
 * E.g., percentTolerance=0 returns only the largest peak(s) and percentTolerance=1 returns _all_ peaks.
*/
std::vector<double> Spectrum::getPeakAmplitudes( const double & percentTolerance ) {
    if ( percentTolerance < 0 || percentTolerance > 1 ) {
        throw std::runtime_error( "percentTolerance argument to Spectrum#getPeakAmplitudes must be [0 1] inclusive." );
    }
    std::vector<double> peaks = {};
    bool peakFound;
    int index;
    double smallestPeak = DBL_MAX;
    double largestPeak = -DBL_MAX; 
    for ( int i = 1; i < getSpectralData().size()-1; i++ ) {
        peakFound = false;
        if ( i == 1 && getSpectralData().at(i-1) > getSpectralData().at(i) ) {
            peakFound = true; index = i-1;
        } else if ( getSpectralData().at(i-1) < getSpectralData().at(i) && getSpectralData().at(i) > getSpectralData().at(i+1) ) {
            peakFound = true; index = i;
        } else if ( i == getSpectralData().size()-2 && getSpectralData().at(i+1) > getSpectralData().at(i) ) {
            peakFound = true; index = i+1;
        }
        if ( peakFound ) {
            if ( getSpectralData().at( index ) < smallestPeak ) { smallestPeak = getSpectralData().at( index ); }
            if ( getSpectralData().at( index ) > largestPeak ) { largestPeak = getSpectralData().at( index ); }
            peaks.push_back( getSpectralData().at( index ) );
        }
    } 
    for ( int k = peaks.size()-1; k >= 0 ; k-- ) {
        if ( (  percentDifference( largestPeak, peaks.at(k) ) /
                percentDifference( largestPeak, smallestPeak ) ) > 
                percentTolerance ) {
            peaks.erase( peaks.begin() + k );
        }
    }
    return peaks;
}

/**
 * Return vector of peaks as vector indices. See #getPeakAmplitudes.
*/
std::vector<int> Spectrum::getPeakIndices( const double & tolerance ) {
    std::vector<double> peakAmplitudes = getPeakAmplitudes( tolerance );
    std::vector<int> indices = {};
    std::vector<double> data = getSpectralData();
    for ( int i = 0; i < peakAmplitudes.size(); i++ ) {
        for ( int j = 0; j < data.size(); j++ ) {
            if ( data.at(j) == peakAmplitudes.at(i) ) {
                indices.push_back( j );
                data.at(j) = 1.0 / 0;
                break;
            }
        }
    }
    return indices;
}

/**
 * Return full width half max for given peak.
 * If the peak occurs at the first or last index, that peak is assumed symmetric.
 * Spectra with large DC offset discontinuities at peak locations will return less accurate FWHM estimates.
*/
double Spectrum::peakFullWidthHalfMax( const int & index ) {
    std::vector<int> peaks = getPeakIndices( 1 );
    bool notFound = true;
    for ( int i = 0; i < peaks.size(); i++ ) {
        if ( peaks.at(i) == index ) { notFound = false; }
    }
    if ( notFound ) { throw std::runtime_error( "No peak found at specificed index." ); }
    int leftIndex = index;
    bool leftFringePeak = true;
    while ( leftIndex > 0 ) {
        leftFringePeak = false;
        if ( ! (getSpectralData().at( leftIndex - 1 ) < getSpectralData().at( leftIndex )) ) { break; }
        leftIndex--;
    }
    int rightIndex = index;
    bool rightFringePeak = true;
    while ( rightIndex < (getSpectralData().size() - 2) ) {
        rightFringePeak = false;
        if ( ! (getSpectralData().at( rightIndex + 1 ) < getSpectralData().at( rightIndex )) ) { break; }
        rightIndex++;
    }
    if ( leftFringePeak ) { 
        return getSpectralWavelengths().at( rightIndex ) / 2; 
    }
    if ( rightFringePeak ) { 
        return getSpectralWavelengths().at( getSpectralWavelengths().size()-1 ) - getSpectralWavelengths().at( leftFringePeak ); 
    }
    return ( getSpectralWavelengths().at( rightIndex ) - getSpectralWavelengths().at( leftIndex ) ) / 2;
}

/**
 * Return full width half max for a peak at specified wavelength. See overloaded method.
*/
double Spectrum::peakFullWidthHalfMax( const double & wavelength ) {
    for ( int i = 0; i < getSpectralWavelengths().size(); i++ ) {
        if ( wavelength == getSpectralWavelengths().at(i) ) {
            return peakFullWidthHalfMax( i );
        }
    }
    throw std::runtime_error( "Wavelength not found." );
}

/**
 * Determine if spectrum is saturated, defined as three contiguous maximum values located anywhere in the spectrum.
*/
bool Spectrum::isSaturated() {
    double maxValue = getMaxValue( getSpectralData() );
    for ( int i = 0; i < getSpectralData().size() - 2; i++ )  {
        if (    ( getSpectralData().at(i)   == getSpectralData().at(i+1) ) & 
                ( getSpectralData().at(i+1) == getSpectralData().at(i+2) ) & 
                ( getSpectralData().at(i) == maxValue ) )  {
            return true;
        }
    }
    return false;
}

//
// radiometry and photometry
//

/**
 * Compute irradiance in [uW/cm^2]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#compute_uWattPerCmSquared.
*/
double Spectrum::computeMicrowattsPerCmSquared( int & lowerBoundWavelength, int & upperBoundWavelength ) {
    return computeIntegral( lowerBoundWavelength, upperBoundWavelength );
}

/**
 * Compute absolute power in [uW]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#compute_uWatt.
 * collectionArea in [cm^2].
*/
double Spectrum::computeMicrowatts( int & lowerBoundWavelength, int & upperBoundWavelength, double & collectionArea ) {
    return computeMicrowattsPerCmSquared( lowerBoundWavelength, upperBoundWavelength ) * collectionArea;
}

/**
 * Compute spectral power in [uW/nm]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#compute_uWattPerNm.
 * collectionArea in [cm^2].
*/
Spectrum Spectrum::computeMicrowattsPerNanometer( double & collectionArea ) {
    return *this * collectionArea;
}

/**
 * Compute energy in Joules. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#compute_Joules.
 * collectionArea in [cm^2], integrationTime in [seconds].
*/
double Spectrum::computeJoules( int & lowerBoundWavelength, 
                                int & upperBoundWavelength, 
                                double & collectionArea, 
                                double & integrationTime ) {
    return computeMicrowatts( lowerBoundWavelength, upperBoundWavelength, collectionArea ) * integrationTime * (1.0/1000000.0);
}

/**
 * Compute energy in electron volts. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#computeElectronVolts_eV.
 * collectionArea in [cm^2], integrationTime in [seconds].
*/
double Spectrum::computeElectronVolts(  int & lowerBoundWavelength, 
                                        int & upperBoundWavelength, 
                                        double & collectionArea, 
                                        double & integrationTime ) {
    return computeJoules( lowerBoundWavelength, upperBoundWavelength, collectionArea, integrationTime ) * ( 1 / electricChargeElectron );
}

/**
 * Compute lux in [lumen/m^2] for an arbitrary luminousEfficiencyFunction. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#computeIlluminanceLux.
 * Note: user need not specify an incident area in m^2, as method returns lux in lumens/m^2 based upon the collection area in cm^2 for original data.
 * Note: #getSpectralWavelengths must match efficiencyFunctionWavelengths. Else, #crop.
*/
double Spectrum::computeLux(    Spectrum & luminousEfficiencyFunction, 
                                double & maxLuminousEfficiencyCoefficient ) {
    double convertMicrowattToWatt = 1.0 / 1000000.0;
    double convertPerCmSquaredToPerMeterSquared = 1.0 / 10000.0;
    Spectrum wattsPerMeterSquaredPerNanometer = *this * convertMicrowattToWatt * convertPerCmSquaredToPerMeterSquared;
    Spectrum lumensPerMeterSquaredPerNanometer = wattsPerMeterSquaredPerNanometer *= luminousEfficiencyFunction;
    return maxLuminousEfficiencyCoefficient * lumensPerMeterSquaredPerNanometer.computeIntegral( lowerBound, upperBound );
}

/**
 * Compute photopic lux in [lumen/m^2]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * See base method #computeLux.
*/
double Spectrum::computeLuxPhotopic( ) {
    Spectrum photopic { photopicEfficiencyData, efficiencyFunctionWavelengths };
    return computeLux( photopic, peakPhotopicLuminosity );
}

/**
 * Compute scotopic lux in [lumen/m^2]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * See base method #computeLux.
*/
double Spectrum::computeLuxScotopic( ) {
    Spectrum scotopic { scotopicEfficiencyData, efficiencyFunctionWavelengths };
    return computeLux( scotopic, peakScotopicLuminosity );
}

/**
 * Compute lumens. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * collectionArea in [cm^2].
*/
double Spectrum::computeLumens( Spectrum & luminousEfficiencyFunction, 
                                double & maxLuminousEfficiencyCoefficient,
                                double & collectionArea ) {
    double convertPerMeterSquaredToPerCentimeterSquared = 10000.0;
    return collectionArea * convertPerMeterSquaredToPerCentimeterSquared * computeLux( luminousEfficiencyFunction, maxLuminousEfficiencyCoefficient );
}

/**
 * Compute candella in [lumen/steradians]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#computeLuminousIntensityCandela.
 * solidAngle in [steradians]. collectionArea in [cm^2].
*/
double Spectrum::computeCandella(   Spectrum & luminousEfficiencyFunction, 
                                    double & maxLuminousEfficiencyCoefficient,
                                    double & collectionArea,
                                    double & solidAngle ) {
    return computeLumens( luminousEfficiencyFunction, maxLuminousEfficiencyCoefficient, collectionArea ) / solidAngle;
}

/**
 * Compute luminous flux in [photons/cm^2/s]. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#computePhotonsPerCmSquaredPerSecond.
*/
double Spectrum::computePhotonsPerCmSquaredPerSecond( int & lowerBoundWavelength, int & upperBoundWavelength ) {
    std::vector<double> photonsPerCmSquaredPerSecondPerNanometer = {};
    double plancksConstant = 6.62607015e-34;
    double velocityOfLight = 2.998e8;
    double metersToNanometer = 1e-9;
    double wattsToMicrowatts = 1e-6;
    for ( int i = 0; i < getSpectralData().size(); i++ ) {
        photonsPerCmSquaredPerSecondPerNanometer.push_back( (getSpectralData().at(i) * getSpectralWavelengths().at(i) * metersToNanometer * wattsToMicrowatts) / ( plancksConstant * velocityOfLight ) );
    }
    return Spectrum{ photonsPerCmSquaredPerSecondPerNanometer, getSpectralWavelengths() }.computeIntegral( lowerBound, upperBound );
}

/**
 * Compute total photons. Calculation assumes #getSpectralData in units of [uW/cm^2/nm] (absolute spectral irradiance).
 * Calculation replicates Ocean Optics Spectral Processing and Math AdvancedPhotometrics#computeTotalPhotons.
 * integrationTime in [seconds]. collectionArea in [cm^2].
*/
double Spectrum::computeTotalPhotons( int & lowerBoundWavelength, int & upperBoundWavelength, double & collectionArea, double & integrationTime ) {
    return computePhotonsPerCmSquaredPerSecond( lowerBoundWavelength, upperBoundWavelength ) * collectionArea * integrationTime;
}

//
// Spectral Objects
//

/**
 * Splice two spectral objects. 
 * Spectra can be merged, such that an overlapping indices are necessarily averaged. Overlapping regions without matching wavelength values are simply interleaved.
 * Two spectra that have no overlap, e.g., an abscissa gap exists between the two, will be spliced as though there is no gap at all. E.g., { 0, 1 } and { 5, 6, 7 } merge as { 0, 1, 5, 6, 7 }.
*/
Spectrum Spectrum::splice( Spectrum & s ) {
    std::vector<double> data = getSpectralData();
    std::vector<double> wavelengths = getSpectralWavelengths();
    int parameterIndexer = 0;
    for ( int j = parameterIndexer; j < s.getSpectralWavelengths().size(); j++ ) {
        for ( int k = 0; k < wavelengths.size(); k++ ) {
            if ( s.getSpectralWavelengths().at(j) > wavelengths.at(k) & k == wavelengths.size() - 1 ) {
                wavelengths.push_back( s.getSpectralWavelengths().at(j) );
                data.push_back( s.getSpectralData().at(j) );
                parameterIndexer++;
                break;
            } else if ( s.getSpectralWavelengths().at(j) < wavelengths.at(k) ) {
                wavelengths.insert( wavelengths.begin() + k, s.getSpectralWavelengths().at(j) );
                data.insert( data.begin() + k, s.getSpectralData().at(j) );
                parameterIndexer++;
                break;
            } else if ( s.getSpectralWavelengths().at(j) == wavelengths.at(k) ) {
                data.at(k) = ( data.at(k) + s.getSpectralData().at(j) ) / 2;
                parameterIndexer++;
                break;
            }
        }
    }
    return Spectrum{ data, wavelengths };
}

/**
 * Crop a spectrum using wavelength values rather than vector indices. See overloaded method.
*/
void Spectrum::crop( double & lowerWavelength, double & upperWavelength ) {
    auto [ lowerWavelengthIndex, upperWavelengthIndex ] = getWavelengthIndices( lowerWavelength, upperWavelength );
    Spectrum::crop( lowerWavelengthIndex, upperWavelengthIndex );
}

/**
 * Crop a spectrum from an initial index lowerWavelengthIndex to a final index upperWavelengthIndex.
*/
void Spectrum::crop( int & lowerWavelengthIndex, int & upperWavelengthIndex ) {
    if ( upperWavelengthIndex == -1 ) { upperWavelengthIndex = getSpectralData().size() - 1; }
    if ( lowerWavelengthIndex < 0 || upperWavelengthIndex > getSpectralData().size() - 1 ) {
        throw std::runtime_error( "Index argument to #crop is out-of-bounds." );
    }
    if ( lowerWavelengthIndex >= upperWavelengthIndex ) {
        throw std::runtime_error( "Index argument to #crop is out-of-bounds." );
    }
    // erase upper indices
    spectralData.erase(         spectralData.begin() + upperWavelengthIndex+1, 
                                spectralData.begin() + spectralData.size() );
    spectralWavelengths.erase(  spectralWavelengths.begin() + upperWavelengthIndex+1, 
                                spectralWavelengths.begin() + spectralWavelengths.size() );
    // erase lower indices
    spectralData.erase(         spectralData.begin() + 0, 
                                spectralData.begin() + lowerWavelengthIndex );
    spectralWavelengths.erase(  spectralWavelengths.begin() + 0, 
                                spectralWavelengths.begin() + lowerWavelengthIndex );
}

/**
 * Derive a vector that, when multiplied with a spectrum, results in some idealReferenceData vector.
*/
std::vector<double> Spectrum::unitlessSpectralCalibration( std::vector<double> & idealReferenceData ) {
    return ( *this /= idealReferenceData ).getSpectralData();
}

/**
 * Conduct spectral smoothing with a boxcar. Data resulting from this method call matches MatLab function movmean().
*/
Spectrum Spectrum::boxCarAveraged( const int & pixelWindow ) {
    if ( pixelWindow <= 0 ) { throw std::runtime_error( "Box car pixel width must be >= 1." ); }
    std::vector<double> data = {}; 
    double sum = 0; 
    int sumCount = 0; 
    int centeringOffset;
    if ( pixelWindow % 2 == 0 ) {
        centeringOffset = ( pixelWindow / 2 ) - 1;
    } else {
        centeringOffset = ( (pixelWindow + 1) / 2 ) - 1;
    }
    for ( int i = centeringOffset; i < getSpectralData().size() + centeringOffset; i++ ) {
        for ( int j = 0; j < pixelWindow; j++ ) {
            if ( i - j < 0 ) { break; }
            try {
                sum += getSpectralData().at( i - j );
                sumCount += 1;
            } catch ( const std::out_of_range& ) {
                continue;
            }
        }
        data.push_back( sum / sumCount );
        sum = 0;
        sumCount = 0;
    }
    return Spectrum{ data, getSpectralWavelengths() };
}

//
// load vectors from file
//

/**
 * Load carriage-return-delimited data from file. Specify headerRowCount>0 if that file contains a header which must be omitted.
*/
std::vector<double> Spectrum::loadFromFile( const std::string & filePath, const int & headerRowCount = 0 )
{
    std::vector<double> v;
    std::ifstream file( filePath );
    std::string line;
    if ( file.is_open() ) {
        int i = 1;
        while ( file ) {
            std::getline( file, line );
            if ( i <= headerRowCount ) { i++; continue; }            
            v.push_back( stod(line) );
        }
    } else {
        std::cerr << "Error: Cannot locate or open file.\n";
    }
    file.close();
    return v;
}

//
// maths
//

/**
 * Create a Gaussian plotted along s.getSpectralWavelengths() with some maxValue centered on centerValue with rmsWidth.
*/
Spectrum Spectrum::createGaussian( std::vector<double> & abscissa, const double & maxValue, const double & centerValue, const double & rmsWidth ) {
    std::vector<double> gaussWuzHere = {};
    double e = 2.71828;
    for ( int i = 0; i < abscissa.size(); i++ ) {
        gaussWuzHere.push_back( maxValue * std::pow( e, -( std::pow( (abscissa.at(i)-centerValue), 2 ) / ( 2.0 * std::pow( rmsWidth, 2 )) ) ) );
    }
    return Spectrum{ gaussWuzHere, abscissa };
}

/**
 * Perform a trapezoidal integration using wavelengths as boundaries.
*/
double Spectrum::computeIntegral( double & lowerWavelength, double & upperWavelength ) {
    auto [ lowerWavelengthIndex, upperWavelengthIndex ] = getWavelengthIndices( lowerWavelength, upperWavelength );
    return Spectrum::computeIntegral( lowerWavelengthIndex, upperWavelengthIndex );
}

/**
 * Perform a trapezoidal integration using abscissa vector indices as boudnaries. 
*/
double Spectrum::computeIntegral( int & lowerWavelengthIndex, int & upperWavelengthIndex ) {
    if ( upperWavelengthIndex == -1 ) {
        upperWavelengthIndex = getSpectralData().size()-1;
    }
    if ( lowerWavelengthIndex < 0 || upperWavelengthIndex > getSpectralData().size()-1 ) {
        throw std::runtime_error( "Index argument to #computeIntegral is out-of-bounds." );
    }
    double dy = 0.0;
    for ( int i = 1; i < getSpectralWavelengths().size(); ++i ) {
        double dx = getSpectralWavelengths().at(i) - getSpectralWavelengths().at( i-1 );
        dy += 0.5 * (getSpectralData().at(i-1) + getSpectralData().at(i) ) * dx;
    }
    return dy;
}

/**
 * Identify maximum value for any std::vector<double>.
*/
double Spectrum::getMaxValue( const std::vector<double> & v ) {
    double maxValue = -DBL_MAX; 
    for ( int i = 0; i < v.size(); i ++ ) {
        if ( v.at(i) > maxValue ) { maxValue = v.at(i); }
    }
    return maxValue;
}

/**
 * Identify minimum value for any std::vector<double>.
*/
double Spectrum::getMinValue( const std::vector<double> & v ) {
    double minValue = DBL_MAX; 
    for ( int i = 0; i < v.size(); i ++ ) {
        if ( v.at(i) < minValue ) { minValue = v.at(i); }
    }
    return minValue;
}

/**
 * Compute percent difference of two values in arbitrary order.
*/
double Spectrum::percentDifference( const double & v1, const double & v2 ) {
    return ( std::abs( v1 - v2 ) / ((v1 + v2) / 2.0) ) * 100;
}

//
// whole-spectrum vector arithmetic
//

/**
 * Perform an index-wise multiplication of spectral data with a vector of matching length.
*/
Spectrum operator*=( Spectrum s, std::vector<double> & v )  {
    for (int i = 0; i < s.getSpectralData().size(); i++) {
        s.setSpectralData( i, s.getSpectralData()[i] * v.at(i) );
    }
    return Spectrum( s.getSpectralData(), s.getSpectralWavelengths() );
}

/**
 * See overloaded method.
*/
Spectrum operator*=( std::vector<double> & v, Spectrum s ) { return s *= v; }

/**
 * Perform an index-wise multiplication of data within respective spectral objects. Abscissa must match and will persist in the return object.
*/
Spectrum operator*=( Spectrum s1, Spectrum s2 )  {
    for ( int k = 0; k < s1.getSpectralWavelengths().size(); k++ ) {
        if ( s1.getSpectralWavelengths().at(k) != s2.getSpectralWavelengths().at(k) ) {
            throw std::runtime_error( "Spectral-wise multiplication requires matching wavelength vectors." );
        }
    }
    if ( s1.getSpectralData().size() != s2.getSpectralData().size() ) {
        throw std::runtime_error( "Spectral-wise multiplication requires matching data dimensionality." );
    }
    std::vector<double> v = {};
    for ( int i = 0; i < s1.getSpectralData().size(); i++ ) {
        v.push_back( s1.getSpectralData().at(i) * s2.getSpectralData().at(i) );
    }
    return Spectrum( v, s1.getSpectralWavelengths() );
}

/**
 * Perform an index-wise division of data within a spectrum with a vector of matching length.
*/
Spectrum operator/=( Spectrum s, std::vector<double> & v )  {
    for (int i = 0; i < s.getSpectralData().size(); i++) {
        s.setSpectralData( i, s.getSpectralData()[i] / v.at(i) );
    }
    return Spectrum( s.getSpectralData(), s.getSpectralWavelengths() );
}

//
// spectral constant vector arithmetic
//

/**
 * Multiply some constant t with each index of spectral data s.
*/
Spectrum operator*( Spectrum s, double & t ) {
    for (int i = 0; i < s.getSpectralData().size(); i++) {
        s.setSpectralData( i, s.getSpectralData()[i] * t );
    }
    // Spectrum s { s.getSpectralData(), s.getSpectralWavelengths() };
    return Spectrum { s.getSpectralData(), s.getSpectralWavelengths() };
}

/**
 * See overloaded method.
*/
Spectrum operator*( double & t, Spectrum s ) { return s * t; }

/**
 * Multiply some constant t with each index of std::vector v.
*/
std::vector<double> operator*( std::vector<double> & v, double & t ) {
    for (int i = 0; i < v.size(); i++) {
        v.at(i) *= t;
    }
    return v;
}

/**
 * See overloaded method.
*/
std::vector<double> operator*( double & t, std::vector<double> & v ) { return v * t; }

//
// other operators
//

/**
 * Print a spectrum object in a tidy format.
*/
std::ostream & operator<<( std::ostream & out, Spectrum & s ) {
    out << "\nAmplitudes : [";
    for ( size_t i = 0; i < s.getSpectralData().size(); ++i ) {
        out << s.getSpectralData()[i];
        if ( i != s.getSpectralData().size()-1 ) 
            out << ", ";
    }
    out << "]";
    out << "\nWavelengths: [";
    for ( size_t i = 0; i < s.getSpectralWavelengths().size(); ++i ) {
        out << s.getSpectralWavelengths()[i];
        if ( i != s.getSpectralWavelengths().size()-1 ) 
            out << ", ";
    }
    out << "]";
    return out;
}
