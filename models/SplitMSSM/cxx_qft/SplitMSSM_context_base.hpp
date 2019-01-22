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

// File generated at Tue 22 Jan 2019 16:49:16

/**
 * @file cxx_qft/SplitMSSM_context_base.hpp
 *
 * This file was generated at Tue 22 Jan 2019 16:49:16 with FlexibleSUSY
 * 2.3.0 and SARAH 4.14.1 .
 */

#ifndef SplitMSSM_CXXQFT_CONTEXT_BASE_H
#define SplitMSSM_CXXQFT_CONTEXT_BASE_H

#include "SplitMSSM_fields.hpp"

namespace flexiblesusy
{
namespace SplitMSSM_cxx_diagrams
{
   struct context_base {
      SplitMSSM_mass_eigenstates& model; ///< The model object.

      template <class Field>
      double mass(const typename field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename fields::remove_lorentz_conjugation<Field>::type;
         return mass_impl<CleanField>(indices);
      }

      context_base(SplitMSSM_mass_eigenstates& m) : model(m) {}
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
double context_base::mass_impl<fields::VG>(const std::array<int, 1>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::gG>(const std::array<int, 1>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::Hp>(const std::array<int, 0>& indices) const
{ return model.get_MHp(); }

template<> inline
double context_base::mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_MFv(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Glu>(const std::array<int, 1>& indices) const
{ return model.get_MGlu(); }

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
double context_base::mass_impl<fields::Fd>(const std::array<int, 2>& indices) const
{ return model.get_MFd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fu>(const std::array<int, 2>& indices) const
{ return model.get_MFu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_MFe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Chi>(const std::array<int, 1>& indices) const
{ return model.get_MChi(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Cha>(const std::array<int, 1>& indices) const
{ return model.get_MCha(indices[0]); }

template<> inline
double context_base::mass_impl<fields::VWp>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

} // namespace SplitMSSM_cxx_diagrams
} // namespace flexiblesusy

#include "SplitMSSM_vertices.hpp"

#endif
