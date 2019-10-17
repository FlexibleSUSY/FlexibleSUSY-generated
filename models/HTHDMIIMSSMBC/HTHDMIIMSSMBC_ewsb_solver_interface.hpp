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

// File generated at Wed 16 Oct 2019 21:26:47

/**
 * @file HTHDMIIMSSMBC_ewsb_solver_interface.hpp
 *
 * @brief contains interface for EWSB solver
 *
 * This file was generated at Wed 16 Oct 2019 21:26:47 with FlexibleSUSY
 * 2.4.1 (git commit: 3e3a10f4fde301d99a4732f29d14f4ac1c646945) and SARAH 4.14.3 .
 */

#ifndef HTHDMIIMSSMBC_EWSB_SOLVER_INTERFACE_H
#define HTHDMIIMSSMBC_EWSB_SOLVER_INTERFACE_H

namespace flexiblesusy {

class HTHDMIIMSSMBC_mass_eigenstates;

/**
 * @class HTHDMIIMSSMBC_ewsb_solver_interface
 * @brief interface for EWSB solvers to be used to solve the EWSB equations
 */
class HTHDMIIMSSMBC_ewsb_solver_interface {
public:
   virtual ~HTHDMIIMSSMBC_ewsb_solver_interface() {}

   virtual void set_loop_order(int) = 0;
   virtual void set_number_of_iterations(int) = 0;
   virtual void set_precision(double) = 0;

   virtual int get_loop_order() const = 0;
   virtual int get_number_of_iterations() const = 0;
   virtual double get_precision() const = 0;

   virtual int solve(HTHDMIIMSSMBC_mass_eigenstates&) = 0;
};

} // namespace flexiblesusy

#endif
