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

// File generated at Sun 10 Jan 2016 15:35:22

#ifndef UMSSM_INPUT_PARAMETERS_H
#define UMSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct UMSSM_input_parameters {
   double m0;
   double m12;
   double TanBeta;
   double Azero;
   double LambdaInput;
   double ALambdaInput;
   double vSInput;
   double Qq;
   double Ql;
   double QHd;
   double QHu;
   double Qd;
   double Qu;
   double Qe;
   double Qs;
   double Qv;

   UMSSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), Azero(0), LambdaInput(0), ALambdaInput(0),
   vSInput(0), Qq(0), Ql(0), QHd(0), QHu(0), Qd(0), Qu(0), Qe(0), Qs(0), Qv(0)

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const UMSSM_input_parameters&);

} // namespace flexiblesusy

#endif
