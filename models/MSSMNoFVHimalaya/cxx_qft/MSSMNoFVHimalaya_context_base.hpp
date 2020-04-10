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

// File generated at Fri 10 Apr 2020 17:35:25

/**
 * @file cxx_qft/MSSMNoFVHimalaya_context_base.hpp
 *
 * This file was generated at Fri 10 Apr 2020 17:35:25 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef MSSMNoFVHimalaya_CXXQFT_CONTEXT_BASE_H
#define MSSMNoFVHimalaya_CXXQFT_CONTEXT_BASE_H

#include "MSSMNoFVHimalaya_mass_eigenstates.hpp"

#include "MSSMNoFVHimalaya_fields.hpp"
#include "MSSMNoFVHimalaya_mass_eigenstates.hpp"

namespace flexiblesusy {
namespace MSSMNoFVHimalaya_cxx_diagrams {

   struct context_base {
      MSSMNoFVHimalaya_mass_eigenstates model; ///< The model object.

      template <class Field>
      double mass(const typename field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename fields::remove_lorentz_conjugation<Field>::type;
         return mass_impl<CleanField>(indices);
      }

      context_base(const MSSMNoFVHimalaya_mass_eigenstates& m) : model(m) {}
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
double context_base::mass_impl<fields::Glu>(const std::array<int, 0>& indices) const
{ return model.get_MGlu(); }

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
double context_base::mass_impl<fields::Fd>(const std::array<int, 0>& indices) const
{ return model.get_MFd(); }

template<> inline
double context_base::mass_impl<fields::Fs>(const std::array<int, 0>& indices) const
{ return model.get_MFs(); }

template<> inline
double context_base::mass_impl<fields::Fb>(const std::array<int, 0>& indices) const
{ return model.get_MFb(); }

template<> inline
double context_base::mass_impl<fields::Fu>(const std::array<int, 0>& indices) const
{ return model.get_MFu(); }

template<> inline
double context_base::mass_impl<fields::Fc>(const std::array<int, 0>& indices) const
{ return model.get_MFc(); }

template<> inline
double context_base::mass_impl<fields::Ft>(const std::array<int, 0>& indices) const
{ return model.get_MFt(); }

template<> inline
double context_base::mass_impl<fields::Fve>(const std::array<int, 0>& indices) const
{ return model.get_MFve(); }

template<> inline
double context_base::mass_impl<fields::Fvm>(const std::array<int, 0>& indices) const
{ return model.get_MFvm(); }

template<> inline
double context_base::mass_impl<fields::Fvt>(const std::array<int, 0>& indices) const
{ return model.get_MFvt(); }

template<> inline
double context_base::mass_impl<fields::Fe>(const std::array<int, 0>& indices) const
{ return model.get_MFe(); }

template<> inline
double context_base::mass_impl<fields::Fm>(const std::array<int, 0>& indices) const
{ return model.get_MFm(); }

template<> inline
double context_base::mass_impl<fields::Ftau>(const std::array<int, 0>& indices) const
{ return model.get_MFtau(); }

template<> inline
double context_base::mass_impl<fields::SveL>(const std::array<int, 0>& indices) const
{ return model.get_MSveL(); }

template<> inline
double context_base::mass_impl<fields::SvmL>(const std::array<int, 0>& indices) const
{ return model.get_MSvmL(); }

template<> inline
double context_base::mass_impl<fields::SvtL>(const std::array<int, 0>& indices) const
{ return model.get_MSvtL(); }

template<> inline
double context_base::mass_impl<fields::Sd>(const std::array<int, 1>& indices) const
{ return model.get_MSd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Su>(const std::array<int, 1>& indices) const
{ return model.get_MSu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Se>(const std::array<int, 1>& indices) const
{ return model.get_MSe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Sm>(const std::array<int, 1>& indices) const
{ return model.get_MSm(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Stau>(const std::array<int, 1>& indices) const
{ return model.get_MStau(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Ss>(const std::array<int, 1>& indices) const
{ return model.get_MSs(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Sc>(const std::array<int, 1>& indices) const
{ return model.get_MSc(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Sb>(const std::array<int, 1>& indices) const
{ return model.get_MSb(indices[0]); }

template<> inline
double context_base::mass_impl<fields::St>(const std::array<int, 1>& indices) const
{ return model.get_MSt(indices[0]); }

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
double context_base::mass_impl<fields::VWm>(const std::array<int, 0>& indices) const
{ return model.get_MVWm(); }

} // namespace MSSMNoFVHimalaya_cxx_diagrams
} // namespace flexiblesusy

#include "MSSMNoFVHimalaya_vertices.hpp"

#endif
