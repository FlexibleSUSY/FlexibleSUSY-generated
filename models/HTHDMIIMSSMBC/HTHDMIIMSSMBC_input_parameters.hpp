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

// File generated at Thu 15 Dec 2016 12:40:51

#ifndef HTHDMIIMSSMBC_INPUT_PARAMETERS_H
#define HTHDMIIMSSMBC_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct HTHDMIIMSSMBC_input_parameters {
   double TanBeta;
   double MSUSY;
   double MEWSB;
   double MuInput;
   double MAInput;
   double AtInput;
   double AbInput;
   double AtauInput;
   double LambdaLoopOrder;

   HTHDMIIMSSMBC_input_parameters()
      : TanBeta(0), MSUSY(0), MEWSB(0), MuInput(0), MAInput(0), AtInput(0),
   AbInput(0), AtauInput(0), LambdaLoopOrder(0)

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const HTHDMIIMSSMBC_input_parameters&);

} // namespace flexiblesusy

#endif
