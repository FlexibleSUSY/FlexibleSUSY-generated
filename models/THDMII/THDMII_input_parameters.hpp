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

// File generated at Sun 26 Aug 2018 14:09:24

#ifndef THDMII_INPUT_PARAMETERS_H
#define THDMII_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct THDMII_input_parameters {
   double Lambda1IN{};
   double Lambda2IN{};
   double Lambda3IN{};
   double Lambda4IN{};
   double Lambda5IN{};
   double Lambda6IN{};
   double Lambda7IN{};
   double M122IN{};
   double TanBeta{};
   double Qin{};


   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const THDMII_input_parameters&);

} // namespace flexiblesusy

#endif
