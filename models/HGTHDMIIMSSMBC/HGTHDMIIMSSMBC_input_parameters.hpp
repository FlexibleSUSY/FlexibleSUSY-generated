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

// File generated at Sun 4 Aug 2019 19:38:50

#ifndef HGTHDMIIMSSMBC_INPUT_PARAMETERS_H
#define HGTHDMIIMSSMBC_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct HGTHDMIIMSSMBC_input_parameters {
   double TanBeta{};
   double MSUSY{};
   double MEWSB{};
   double MuInput{};
   double M1Input{};
   double M2Input{};
   double M3Input{};
   double MAInput{};
   double AtInput{};
   double AbInput{};
   double AtauInput{};
   double LambdaLoopOrder{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const HGTHDMIIMSSMBC_input_parameters&);

} // namespace flexiblesusy

#endif
