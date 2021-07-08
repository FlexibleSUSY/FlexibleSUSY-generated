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


/**
 * @file cxx_qft/TMSSM_context_base.hpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#ifndef TMSSM_CXXQFT_CONTEXT_BASE_H
#define TMSSM_CXXQFT_CONTEXT_BASE_H

#include "TMSSM_fields.hpp"
#include "TMSSM_mass_eigenstates_interface.hpp"

namespace flexiblesusy {
namespace TMSSM_cxx_diagrams {

   struct context_base {
      TMSSM_mass_eigenstates_interface const& model; ///< The model object.

      template <class Field>
      double mass(const typename field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename fields::remove_lorentz_conjugation<Field>::type;
         return mass_impl<CleanField>(indices);
      }
      template<class Field>
      double physical_mass(const typename field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename fields::remove_lorentz_conjugation<Field>::type;
         return physical_mass_impl<CleanField>(indices);
      }

      context_base(TMSSM_mass_eigenstates_interface const& m) : model(m) {}
      context_base(TMSSM_mass_eigenstates_interface const* const m) : model(*m) {}

      context_base(const context_base&) = default;
      context_base(context_base&&) = default;

      virtual ~context_base(void) = default;

   private:
      template <class Field>
      double
      mass_impl(const typename field_indices<Field>::type& indices) const;
      template<class Field>
      double physical_mass_impl(const typename field_indices<Field>::type& indices) const;
   };

template<> inline
double context_base::mass_impl<fields::VG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::gG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::Glu>(const std::array<int, 0>& indices) const
{ return model.get_MGlu(); }

template<> inline
double context_base::mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_MFv(indices[0]); }

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
double context_base::mass_impl<fields::gWm>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

template<> inline
double context_base::mass_impl<fields::gWmC>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

template<> inline
double context_base::mass_impl<fields::Sd>(const std::array<int, 1>& indices) const
{ return model.get_MSd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Sv>(const std::array<int, 1>& indices) const
{ return model.get_MSv(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Su>(const std::array<int, 1>& indices) const
{ return model.get_MSu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Se>(const std::array<int, 1>& indices) const
{ return model.get_MSe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::hh>(const std::array<int, 1>& indices) const
{ return model.get_Mhh(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Ah>(const std::array<int, 1>& indices) const
{ return model.get_MAh(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Hpm>(const std::array<int, 1>& indices) const
{ return model.get_MHpm(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Chi>(const std::array<int, 1>& indices) const
{ return model.get_MChi(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Cha>(const std::array<int, 1>& indices) const
{ return model.get_MCha(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_MFe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fd>(const std::array<int, 1>& indices) const
{ return model.get_MFd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fu>(const std::array<int, 1>& indices) const
{ return model.get_MFu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::VWm>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

template<> inline
double context_base::physical_mass_impl<fields::VG>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVG; }

template<> inline
double context_base::physical_mass_impl<fields::gG>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVG; }

template<> inline
double context_base::physical_mass_impl<fields::Glu>(const std::array<int, 0>& indices) const
{ return model.get_physical().MGlu; }

template<> inline
double context_base::physical_mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFv[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::VP>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVP; }

template<> inline
double context_base::physical_mass_impl<fields::VZ>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVZ; }

template<> inline
double context_base::physical_mass_impl<fields::gP>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVP; }

template<> inline
double context_base::physical_mass_impl<fields::gZ>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVZ; }

template<> inline
double context_base::physical_mass_impl<fields::gWm>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWm; }

template<> inline
double context_base::physical_mass_impl<fields::gWmC>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWm; }

template<> inline
double context_base::physical_mass_impl<fields::Sd>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSd[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Sv>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSv[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Su>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSu[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Se>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSe[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::hh>(const std::array<int, 1>& indices) const
{ return model.get_physical().Mhh[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Ah>(const std::array<int, 1>& indices) const
{ return model.get_physical().MAh[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Hpm>(const std::array<int, 1>& indices) const
{ return model.get_physical().MHpm[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Chi>(const std::array<int, 1>& indices) const
{ return model.get_physical().MChi[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Cha>(const std::array<int, 1>& indices) const
{ return model.get_physical().MCha[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFe[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Fd>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFd[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Fu>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFu[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::VWm>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWm; }

} // namespace TMSSM_cxx_diagrams
} // namespace flexiblesusy

#include "TMSSM_vertices.hpp"

#endif
