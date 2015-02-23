// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Mon 23 Feb 2015 13:31:13

#ifndef NUTNMSSM_INPUT_PARAMETERS_H
#define NUTNMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NUTNMSSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   double Azero;
   double LambdaInput;
   double KappaInput;
   double ALambdaInput;
   double AKappaInput;
   double MuEff;

   NUTNMSSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), Azero(0), LambdaInput(0), KappaInput(0),
   ALambdaInput(0), AKappaInput(0), MuEff(0)

   {}
};

std::ostream& operator<<(std::ostream&, const NUTNMSSM_input_parameters&);

} // namespace flexiblesusy

#endif
