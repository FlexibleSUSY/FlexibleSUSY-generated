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

// File generated at Wed 16 Oct 2019 21:36:18

/**
 * @file cxx_qft/SM_context_base.hpp
 *
 * This file was generated at Wed 16 Oct 2019 21:36:18 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.3 .
 */

#ifndef SM_CXXQFT_CONTEXT_BASE_H
#define SM_CXXQFT_CONTEXT_BASE_H

#include "SM_mass_eigenstates.hpp"

#include "SM_fields.hpp"
#include "SM_mass_eigenstates.hpp"

namespace flexiblesusy {
namespace SM_cxx_diagrams {

   struct context_base {
      SM_mass_eigenstates model; ///< The model object.

      template <class Field>
      double mass(const typename field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename fields::remove_lorentz_conjugation<Field>::type;
         return mass_impl<CleanField>(indices);
      }

      context_base(const SM_mass_eigenstates& m) : model(m) {}
      context_base(const context_base&) = default;
      context_base(context_base&&) = default;

      context_base& operator=(const context_base&) = default;
      context_base& operator=(context_base&&) = default;

      virtual ~context_base(void) = default;

   private:
      template <class Field>
      double
      mass_impl(const typename field_indices<Field>::type& indices) const;
   };

   template<> inline
double context_base::mass_impl<fields::VG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::gG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::Hp>(const std::array<int, 0>& indices) const
{ return model.get_MHp(); }

template<> inline
double context_base::mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_MFv(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Ah>(const std::array<int, 0>& indices) const
{ return model.get_MAh(); }

template<> inline
double context_base::mass_impl<fields::hh>(const std::array<int, 0>& indices) const
{ return model.get_Mhh(); }

template<> inline
double context_base::mass_impl<fields::VP>(const std::array<int, 0>& indices) const
{ return model.get_MVP(); }

template<> inline
double context_base::mass_impl<fields::VZ>(const std::array<int, 0>& indices) const
{ return model.get_MVZ(); }

template<> inline
double context_base::mass_impl<fields::gP>(const std::array<int, 0>& indices) const
{ return model.get_MVP(); }

template<> inline
double context_base::mass_impl<fields::gZ>(const std::array<int, 0>& indices) const
{ return model.get_MVZ(); }

template<> inline
double context_base::mass_impl<fields::gWp>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

template<> inline
double context_base::mass_impl<fields::gWpC>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

template<> inline
double context_base::mass_impl<fields::Fd>(const std::array<int, 1>& indices) const
{ return model.get_MFd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fu>(const std::array<int, 1>& indices) const
{ return model.get_MFu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_MFe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::VWp>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

} // namespace SM_cxx_diagrams
} // namespace flexiblesusy

#include "SM_vertices.hpp"

#endif
