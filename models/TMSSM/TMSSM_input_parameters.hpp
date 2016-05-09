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

// File generated at Mon 9 May 2016 12:11:42

#ifndef TMSSM_INPUT_PARAMETERS_H
#define TMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct TMSSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   int SignMu;
   double Azero;
   double MTinput;
   double Qin;
   double LambdaInput;
   double vTInput;

   TMSSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), SignMu(1), Azero(0), MTinput(0), Qin(0),
   LambdaInput(0), vTInput(0)

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const TMSSM_input_parameters&);

} // namespace flexiblesusy

#endif
