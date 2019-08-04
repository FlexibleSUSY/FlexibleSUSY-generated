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

// File generated at Sun 4 Aug 2019 18:14:08

/**
 * @file MRSSMEFTHiggs_ewsb_solver_interface.hpp
 *
 * @brief contains interface for EWSB solver
 *
 * This file was generated at Sun 4 Aug 2019 18:14:08 with FlexibleSUSY
 * 2.4.0 (git commit: 544c83a2e6b5f23da8d0b6ccdb06f1c91f75d6eb) and SARAH 4.14.2 .
 */

#ifndef MRSSMEFTHiggs_EWSB_SOLVER_INTERFACE_H
#define MRSSMEFTHiggs_EWSB_SOLVER_INTERFACE_H

namespace flexiblesusy {

class MRSSMEFTHiggs_mass_eigenstates;

/**
 * @class MRSSMEFTHiggs_ewsb_solver_interface
 * @brief interface for EWSB solvers to be used to solve the EWSB equations
 */
class MRSSMEFTHiggs_ewsb_solver_interface {
public:
   virtual ~MRSSMEFTHiggs_ewsb_solver_interface() {}

   virtual void set_loop_order(int) = 0;
   virtual void set_number_of_iterations(int) = 0;
   virtual void set_precision(double) = 0;

   virtual int get_loop_order() const = 0;
   virtual int get_number_of_iterations() const = 0;
   virtual double get_precision() const = 0;

   virtual int solve(MRSSMEFTHiggs_mass_eigenstates&) = 0;
};

} // namespace flexiblesusy

#endif
