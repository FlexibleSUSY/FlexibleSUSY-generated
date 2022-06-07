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
 * @file cxx_qft/NUHMSSMNoFVHimalaya_context_base.hpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef NUHMSSMNoFVHimalaya_CXXQFT_CONTEXT_BASE_H
#define NUHMSSMNoFVHimalaya_CXXQFT_CONTEXT_BASE_H

#include "NUHMSSMNoFVHimalaya_fields.hpp"
#include "NUHMSSMNoFVHimalaya_mass_eigenstates_interface.hpp"

namespace flexiblesusy {
namespace NUHMSSMNoFVHimalaya_cxx_diagrams {

   struct context_base {
      NUHMSSMNoFVHimalaya_mass_eigenstates_interface const& model; ///< The model object.

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

      context_base(NUHMSSMNoFVHimalaya_mass_eigenstates_interface const& m) : model(m) {}
      context_base(NUHMSSMNoFVHimalaya_mass_eigenstates_interface const* const m) : model(*m) {}

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
double context_base::physical_mass_impl<fields::Fd>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFd; }

template<> inline
double context_base::physical_mass_impl<fields::Fs>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFs; }

template<> inline
double context_base::physical_mass_impl<fields::Fb>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFb; }

template<> inline
double context_base::physical_mass_impl<fields::Fu>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFu; }

template<> inline
double context_base::physical_mass_impl<fields::Fc>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFc; }

template<> inline
double context_base::physical_mass_impl<fields::Ft>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFt; }

template<> inline
double context_base::physical_mass_impl<fields::Fve>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFve; }

template<> inline
double context_base::physical_mass_impl<fields::Fvm>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFvm; }

template<> inline
double context_base::physical_mass_impl<fields::Fvt>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFvt; }

template<> inline
double context_base::physical_mass_impl<fields::Fe>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFe; }

template<> inline
double context_base::physical_mass_impl<fields::Fm>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFm; }

template<> inline
double context_base::physical_mass_impl<fields::Ftau>(const std::array<int, 0>& indices) const
{ return model.get_physical().MFtau; }

template<> inline
double context_base::physical_mass_impl<fields::SveL>(const std::array<int, 0>& indices) const
{ return model.get_physical().MSveL; }

template<> inline
double context_base::physical_mass_impl<fields::SvmL>(const std::array<int, 0>& indices) const
{ return model.get_physical().MSvmL; }

template<> inline
double context_base::physical_mass_impl<fields::SvtL>(const std::array<int, 0>& indices) const
{ return model.get_physical().MSvtL; }

template<> inline
double context_base::physical_mass_impl<fields::Sd>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSd[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Su>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSu[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Se>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSe[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Sm>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSm[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Stau>(const std::array<int, 1>& indices) const
{ return model.get_physical().MStau[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Ss>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSs[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Sc>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSc[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Sb>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSb[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::St>(const std::array<int, 1>& indices) const
{ return model.get_physical().MSt[indices[0]]; }

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
double context_base::physical_mass_impl<fields::VWm>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWm; }

} // namespace NUHMSSMNoFVHimalaya_cxx_diagrams
} // namespace flexiblesusy

#include "NUHMSSMNoFVHimalaya_vertices.hpp"

#endif
