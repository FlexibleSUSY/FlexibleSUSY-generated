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

// File generated at Tue 8 Sep 2015 13:06:04

#ifndef NUTSMSSM_INPUT_PARAMETERS_H
#define NUTSMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NUTSMSSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   double Azero;
   double LambdaInput;
   double KappaInput;
   double LambdaSInput;
   double L1Input;
   double MSInput;
   double BInput;
   double MuInput;
   double LInput;

   NUTSMSSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), Azero(0), LambdaInput(0), KappaInput(0),
   LambdaSInput(0), L1Input(0), MSInput(0), BInput(0), MuInput(0), LInput(0)

   {}
};

std::ostream& operator<<(std::ostream&, const NUTSMSSM_input_parameters&);

} // namespace flexiblesusy

#endif
