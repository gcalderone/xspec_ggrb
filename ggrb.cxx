/*
  Copyright 2014-2015 Giorgio Calderone (giorgio.calderone@gmail.com)

  This file is part of XSPEC_GGRB.

  XSPEC_GGRB is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  XSPEC_GGRB is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <xsTypes.h>
#include <XSUtil/FunctionUtils/funcType.h>
#include <XSFunctions/functionMap.h>
#include <XSUtil/Numerics/IncGamma.h>
#include <XSUtil/Numerics/ExpInt.h>



//Uncomment this line to print debug informations:
//#define GGRB_DEBUG


#ifdef GGRB_DEBUG
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#endif






//--------------------------------------------------------------------
Real ggrb_int_cpl(Real a,
                  Real Ec,
                  Real Emin,
                  Real Emax) {
  Numerics::IncGamma incGamma;
  //incGamma(2+PhoIndex, Energy) is the primitive of the cutoff power law.
  Real i1 = incGamma(2+a, Emin / Ec);
  Real i2 = incGamma(2+a, Emax / Ec);

#ifdef GGRB_DEBUG
  cout << "GGRB: i1-i2 = " <<  i1-i2 << endl;
#endif

  return -Ec*Ec * (i2-i1);
}



//--------------------------------------------------------------------
Real ggrb_int_pl(Real a,
                 Real b,
                 Real Ec,
                 Real Emin,
                 Real Emax) {
  Real pre = pow(a-b, a-b) * exp(b-a) / pow(Ec, b);

  if (b != -2)
    return pre / (2+b) * (pow(Emax, 2+b) - pow(Emin, 2+b));
  else
    return pre * log(Emax / Emin);
}







//--------------------------------------------------------------------
extern "C" void ggrb (const RealArray& energy,
		      const RealArray& parameter,
		      int spectrum,
		      RealArray& flux,
		      RealArray& fluxError,
		      const string& init)
{

  /*
    The GGRB model implements either the Band or the cutoff power law
    spectral models.

    The value of the "opt" parameter (to be freezed during fit) select
    the appropriate model:
    - opt=0: Cutoff power law ("b" parameter must be freezed)
    - opt=1: Band model

    This model has the integrated flux as a free parameter (log_Flux).
    The nuFnu spectrum has a peak (at E=10^logEp) if alpha is greater
    than beta.  To ensure this condition set up appropriate limits in
    XSPEC.
  */


  //Model parameters
  const Real a        = parameter[0]; //First photon index
  const Real b        = parameter[1]; //Second photon index
  const Real log_Ep   = parameter[2]; //log10 of peak energy [keV]
  const Real log_Flux = parameter[3]; //log10 of integrated flux in energy range Emin-Emax [erg s^-1 cm^-2]
  const Real Emin     = parameter[4]; //Flux integration limits [keV] (to be freezed during fit)
  const Real Emax     = parameter[5];
  const int  opt      = parameter[6]; //Option parameter (either 0 or 1, to be freezed during fit)




  //Cutoff energy
  Real Ec;
  if (a == -2) Ec = pow(10, log_Ep) / 0.0001; //TRICK: avoid a=-2
  else         Ec = pow(10, log_Ep) / (2+a);


  //Split energy
  Real Esplit = (a-b) * Ec;

  

  //Evaluate model integrated flux and normalization
  Real intflux;

  if (opt == 0) { //Cutoff power law
    intflux = ggrb_int_cpl(a, Ec, Emin, Emax);
  }
  else { //Band model
    if (Emax < Esplit) {
      intflux = ggrb_int_cpl(a, Ec, Emin, Emax);
    }
    else {
      if ( (Emin   <=  Esplit)   &&
	   (Esplit <=  Emax  ) ) {
	intflux =
	  ggrb_int_cpl(a   , Ec, Emin  , Esplit) +
	  ggrb_int_pl (a, b, Ec, Esplit, Emax  );
      }
      else {
	if (Esplit < Emin) {
	  intflux = ggrb_int_pl(a, b, Ec, Emin, Emax);
	}
	else {
	  intflux = 1/0;  //?????
	}
      }
    }
  }
  Real norm = pow(10, log_Flux) * 6.24151e8 / intflux;





  //Loop over energy bins
  size_t N(energy.size()-1);
  flux.resize(N);
  for (size_t k=0; k<N; ++k) { 
    Real de = energy[k+1] - energy[k]; //width of energy bin [keV]
    Real E  = pow(10, (log10(energy[k]) + log10(energy[k+1])) / 2); //(logarithmic) average energy in bin [keV]

    flux[k] = norm * de;

    if (opt == 0) { //Cutoff power law
      flux[k] *= pow(E / Ec, a) * exp(-E / Ec);
    }
    else { //Band model
      Real Esplit = (a-b) * Ec;
      if (E < Esplit) flux[k] *= pow(E / Ec, a) * exp(-E / Ec);
      else            flux[k] *= pow(a-b, a-b) * exp(b-a) * pow(E / Ec, b);
    }
#ifdef GGRB_DEBUG
    cout << "GGRB: photon flux [ph / s / cm2]"  << "\t" 
	 << E       << "\t"
	 << flux[k] << endl;
#endif
  }


#ifdef GGRB_DEBUG
  cout << "GGRB: " 
       << a             << "\t"
       << b             << "\t"
       << log_Ep        << "\t"
       << log_Flux      << "\t"
       << Emin          << "\t"
       << Emax          << "\t"
       << opt           << "\t"
       << intflux       << "\t"
       << norm          << endl << endl;
#endif


  fluxError.resize(0);
  return;
}
