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
 * @file CE6SSM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef CE6SSM_DECAYS_H
#define CE6SSM_DECAYS_H

#include "CE6SSM_decay_table.hpp"
#include "CE6SSM_mass_eigenstates.hpp"
#include "CE6SSM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/CE6SSM_qft.hpp"
#include "CE6SSM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename CE6SSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename CE6SSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename CE6SSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename CE6SSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class CE6SSM_decays {
public:
   CE6SSM_decays() = default;
   CE6SSM_decays(CE6SSM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   CE6SSM_decays(const CE6SSM_decays&) = default;
   CE6SSM_decays(CE6SSM_decays&&) = default;
   ~CE6SSM_decays() = default;
   CE6SSM_decays& operator=(const CE6SSM_decays&) = default;
   CE6SSM_decays& operator=(CE6SSM_decays&&) = default;

   const CE6SSM_decay_table& get_decay_table() const;
   const FlexibleDecay_problems& get_problems() const;

   void clear();
   void clear_problems();
   void calculate_decays();

   const Decays_list& get_hh_decays(int i) const { return decay_table.
      get_hh_decays(i); }
   const Decays_list& get_Hpm_decays(int i) const { return decay_table.
      get_Hpm_decays(i); }
   const Decays_list& get_Ah_decays(int i) const { return decay_table.
      get_Ah_decays(i); }
   void calculate_hh_decays();
   void calculate_Hpm_decays();
   void calculate_Ah_decays();

double partial_width_hh_to_SdconjSd(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SDXconjSDX(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjHpm(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHI0SHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHI0conjSHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHIpconjSHIp(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SSI0SSI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SSI0conjSSI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHp0SHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHp0conjSHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHppconjSHpp(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSHI0conjSHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSSI0conjSSI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSHp0conjSHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhVP(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZp(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZp(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjHpmVWm(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjVWm(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_VGVG(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZp(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZp(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZpVZp(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaPChaP(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_ChiChi(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaCha(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFDXFDX(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaIChaI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiIChiI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_FSIFSI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiPChiP(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SdconjSu(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHI0SHIp(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHIpconjSHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHp0SHpp(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHppconjSHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhVWm(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZpVWm(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChaPChiP(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_ChiCha(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_ChaIChiI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SdconjSd(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SDXconjSDX(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjHpm(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHI0SHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHI0conjSHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHIpconjSHIp(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SSI0SSI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SSI0conjSSI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHp0SHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHp0conjSHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHppconjSHpp(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSHI0conjSHI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSSI0conjSSI0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSHp0conjSHp0(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhVP(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZp(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjHpmVWm(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjVWm(CE6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_VGVG(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZp(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZp(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZpVZp(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaPChaP(CE6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_ChiChi(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaCha(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFDXFDX(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaIChaI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiIChiI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_FSIFSI(CE6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiPChiP(CE6SSM_mass_eigenstates_interface*, int, int, int) const;

private:
   CE6SSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   CE6SSM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const CE6SSM_cxx_diagrams::context_base&,
      const typename CE6SSM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename CE6SSM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename CE6SSM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(CE6SSM_cxx_diagrams::context_base const& context,
                  typename CE6SSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename CE6SSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename CE6SSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const CE6SSM_cxx_diagrams::context_base&,
      typename CE6SSM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename CE6SSM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename CE6SSM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Sd, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Sd >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Sv, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Sv >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Su, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Su >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Se, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Se >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SDX, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SDX >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Ah>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Ah>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SHI0, CE6SSM_cxx_diagrams::fields::SHI0>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SHI0, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SHIp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHIp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SSI0, CE6SSM_cxx_diagrams::fields::SSI0>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SSI0 >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SSI0, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SHp0, CE6SSM_cxx_diagrams::fields::SHp0>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SHp0, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::SHpp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHpp>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHpp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHpp>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VG, CE6SSM_cxx_diagrams::fields::VG>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZp, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Glu, CE6SSM_cxx_diagrams::fields::Glu>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Glu >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type, CE6SSM_cxx_diagrams::fields::Fv>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaP>::type, CE6SSM_cxx_diagrams::fields::ChaP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaP>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChaP >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::Chi, CE6SSM_cxx_diagrams::fields::Chi>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Chi >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type, CE6SSM_cxx_diagrams::fields::Cha>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type, CE6SSM_cxx_diagrams::fields::Fd>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type, CE6SSM_cxx_diagrams::fields::Fu>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type, CE6SSM_cxx_diagrams::fields::FDX>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::FDX >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type, CE6SSM_cxx_diagrams::fields::ChaI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChaI >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::ChiI, CE6SSM_cxx_diagrams::fields::ChiI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiI >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiI >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::FSI, CE6SSM_cxx_diagrams::fields::FSI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::FSI >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::FSI >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::ChiP, CE6SSM_cxx_diagrams::fields::ChiP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiP >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::Sd, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Sd >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::Se, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Se >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::SHI0, CE6SSM_cxx_diagrams::fields::SHIp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHIp >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::SHIp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHIp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::SHp0, CE6SSM_cxx_diagrams::fields::SHpp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHpp >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::SHpp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHpp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::VZp, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::ChaP, CE6SSM_cxx_diagrams::fields::ChiP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChaP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiP >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::Chi, CE6SSM_cxx_diagrams::fields::Cha>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Chi >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type, CE6SSM_cxx_diagrams::fields::Fe>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type, CE6SSM_cxx_diagrams::fields::Fd>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Hpm, CE6SSM_cxx_diagrams::fields::ChaI, CE6SSM_cxx_diagrams::fields::ChiI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChaI >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiI >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Sd, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Sd >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Sv, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Sv >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Su, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Su >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Se, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Se >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SDX, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SDX >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SDX>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::hh>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SHI0, CE6SSM_cxx_diagrams::fields::SHI0>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SHI0, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SHIp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHIp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHIp>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SSI0, CE6SSM_cxx_diagrams::fields::SSI0>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SSI0 >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SSI0, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SHp0, CE6SSM_cxx_diagrams::fields::SHp0>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SHp0, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::SHpp, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHpp>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::SHpp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHpp>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::hh, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Hpm, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Hpm >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VG, CE6SSM_cxx_diagrams::fields::VG>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VP, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::VZ>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VZ, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::VZp, CE6SSM_cxx_diagrams::fields::VZp>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type, CE6SSM_cxx_diagrams::fields::VWm>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Glu, CE6SSM_cxx_diagrams::fields::Glu>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Glu >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type, CE6SSM_cxx_diagrams::fields::Fv>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fv>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaP>::type, CE6SSM_cxx_diagrams::fields::ChaP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaP>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChaP >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::Chi, CE6SSM_cxx_diagrams::fields::Chi>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Chi >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type, CE6SSM_cxx_diagrams::fields::Cha>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Cha>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type, CE6SSM_cxx_diagrams::fields::Fe>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type, CE6SSM_cxx_diagrams::fields::Fd>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type, CE6SSM_cxx_diagrams::fields::Fu>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type, CE6SSM_cxx_diagrams::fields::FDX>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::FDX>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::FDX >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type, CE6SSM_cxx_diagrams::fields::ChaI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::ChaI>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChaI >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::ChiI, CE6SSM_cxx_diagrams::fields::ChiI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiI >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiI >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::FSI, CE6SSM_cxx_diagrams::fields::FSI>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::FSI >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::FSI >::type&) const;

template<>
Decay_amplitude_SFF CE6SSM_decays::calculate_amplitude<CE6SSM_cxx_diagrams::fields::Ah, CE6SSM_cxx_diagrams::fields::ChiP, CE6SSM_cxx_diagrams::fields::ChiP>(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::ChiP >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
CE6SSM_decays::amplitude_squared(CE6SSM_cxx_diagrams::context_base const& context,
                  typename CE6SSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename CE6SSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename CE6SSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
{

   const auto mat_elem = calculate_amplitude<FieldIn, FieldOut1, FieldOut2>(
      context, indexIn, indexOut1, indexOut2);
   return mat_elem.square();
}

// returns a squared color generator for a 3 point amplitude with FieldIn, FieldOut1 and FieldOut2
// averaged over inital state colors
// the generator is guessed from color representations of FieldIn, FieldOut1 and FieldOut2
// This is not a bulletproof solution and might fail in general but is enough for
// decays of color singlets

// 1 -> 1, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(CE6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   CE6SSM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (CE6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   CE6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// 8 -> 8, 8 with identical particles in the final state
// because of symmetry of the final state it must be proportional to d^2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
// average: 1/8
squared_color_generator() {return 40/24.;}

// 8 -> 8, 8 with differnt particles in the final state
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
CE6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
!std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   f^2 = 2 Nc (-1 + Nc^2) TR = 24
// average: 1/8
squared_color_generator() {return 3.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double CE6SSM_decays::get_partial_width(
   const CE6SSM_cxx_diagrams::context_base& context,
   typename CE6SSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename CE6SSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename CE6SSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
   ) const
{

   // on-shell masses
   const double mIn = context.physical_mass<FieldIn>(indexIn);
   const double mOut1 = context.physical_mass<FieldOut1>(indexOut1);
   const double mOut2 = context.physical_mass<FieldOut2>(indexOut2);

   // is decay kinematically allowed?
   if(mIn < mOut1 + mOut2) {
      WARNING("Called kinematically forbidden decay");
      return 0.;
   }

   // phase space without symmetry factor
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mOut1/mIn), Sqr(mOut2/mIn)));

   // phase space symmetry factor
   const double ps_symmetry =
      final_state_symmetry_factor<FieldOut1, FieldOut2>(indexOut1, indexOut2);

   // color factor
   constexpr double color_factor = squared_color_generator<FieldIn, FieldOut1, FieldOut2>();

   // matrix element squared
   const auto mat_elem_sq = amplitude_squared<FieldIn, FieldOut1, FieldOut2>(
      context, indexIn, indexOut1, indexOut2);

   // flux * phase space factor * symmetry factor * color factor * |matrix element|^2
   const auto result = 0.5/mIn * ps * ps_symmetry * color_factor * mat_elem_sq;

   return result;
}

template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,CE6SSM_cxx_diagrams::fields::VZ,CE6SSM_cxx_diagrams::fields::VZ >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type,CE6SSM_cxx_diagrams::fields::VWm >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::conj<CE6SSM_cxx_diagrams::fields::VWm>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VWm >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,CE6SSM_cxx_diagrams::fields::VG,CE6SSM_cxx_diagrams::fields::VG >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,CE6SSM_cxx_diagrams::fields::VP,CE6SSM_cxx_diagrams::fields::VP >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,CE6SSM_cxx_diagrams::fields::VP,CE6SSM_cxx_diagrams::fields::VZ >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type,CE6SSM_cxx_diagrams::fields::Fu >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type,CE6SSM_cxx_diagrams::fields::Fd >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::hh,typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type,CE6SSM_cxx_diagrams::fields::Fe >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::hh >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fe >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::Ah,typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type,CE6SSM_cxx_diagrams::fields::Fd >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::Ah,typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type,CE6SSM_cxx_diagrams::fields::Fu >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::Ah,CE6SSM_cxx_diagrams::fields::VG,CE6SSM_cxx_diagrams::fields::VG >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::Ah,CE6SSM_cxx_diagrams::fields::VP,CE6SSM_cxx_diagrams::fields::VP >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::Ah,CE6SSM_cxx_diagrams::fields::VP,CE6SSM_cxx_diagrams::fields::VZ >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VP >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double CE6SSM_decays::get_partial_width<CE6SSM_cxx_diagrams::fields::Ah,typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type,CE6SSM_cxx_diagrams::fields::Fe >(const CE6SSM_cxx_diagrams::context_base&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Ah >::type&, const typename CE6SSM_cxx_diagrams::field_indices<typename CE6SSM_cxx_diagrams::fields::bar<CE6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename CE6SSM_cxx_diagrams::field_indices<CE6SSM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
