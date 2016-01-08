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

// File generated at Fri 8 Jan 2016 13:28:26

#ifndef MSSMatMGUT_OBSERVABLES_H
#define MSSMatMGUT_OBSERVABLES_H

#include "observables.hpp"

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class MSSMatMGUT_mass_eigenstates;

Observables calculate_observables(const MSSMatMGUT_mass_eigenstates&, const softsusy::QedQcd&);

} // namespace flexiblesusy

#endif
