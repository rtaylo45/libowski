//*****************************************************************************
// Author: Zack Taylor
//
// List of unit conversions
//*****************************************************************************
#ifndef UTILCONVERSIONS_H
#define UTILCONVERSIONS_H

//*****************************************************************************
// General SI conversions
//*****************************************************************************
double baseToGiga = 1e9;
double baseToMega = 1e6;
double baseToKilo = 1e3;
double baseToHecto = 1e2;
double baseToDeka = 1e1;
double baseToDeci = 1e-1;
double baseToCenti = 1e-2;
double basetoMilli = 1e-3;
double baseToMicro = 1e-6;
double basetoNano = 1e-9;

//*****************************************************************************
// Length conversions
//*****************************************************************************
double mTocm = baseToCenti;
double cmTom = 1./mTocm;
double ftTom = 0.3048;
double mToft = 1./ftTom;

//*****************************************************************************
// Mass conversions
//*****************************************************************************
double gTokg = baseToKilo;
double kgTog = 1./baseToKilo;
double gTomilg = baseToMilli;
double milgTog = 1./baseToMilli;
double lbTokg = 0.4536;
double kgTolb = 1./lbToKg;

//*****************************************************************************
// Pressure conversions
//*****************************************************************************
double atmTokpa = 101.325;
double kpaToatm = 1./atmTokpa;
double paTokpa = baseToKilo;
double kpaTopa = 1./paTokpa;

#endif
