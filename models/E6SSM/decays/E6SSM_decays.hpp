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
 * @file E6SSM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.5 .
 */

#ifndef E6SSM_DECAYS_H
#define E6SSM_DECAYS_H

#include "E6SSM_decay_table.hpp"
#include "E6SSM_mass_eigenstates.hpp"
#include "E6SSM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/E6SSM_qft.hpp"
#include "E6SSM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename E6SSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename E6SSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename E6SSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename E6SSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class E6SSM_decays {
public:
   E6SSM_decays() = default;
   E6SSM_decays(E6SSM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   E6SSM_decays(const E6SSM_decays&) = default;
   E6SSM_decays(E6SSM_decays&&) = default;
   ~E6SSM_decays() = default;
   E6SSM_decays& operator=(const E6SSM_decays&) = default;
   E6SSM_decays& operator=(E6SSM_decays&&) = default;

   const E6SSM_decay_table& get_decay_table() const;
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

double partial_width_hh_to_SdconjSd(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SDXconjSDX(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjHpm(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHI0SHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHI0conjSHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHIpconjSHIp(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SSI0SSI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SSI0conjSSI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHp0SHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHp0conjSHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SHppconjSHpp(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSHI0conjSHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSSI0conjSSI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSHp0conjSHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhVP(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZp(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZp(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjHpmVWm(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjVWm(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_VGVG(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZp(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZp(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZpVZp(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaPChaP(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_ChiChi(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaCha(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFDXFDX(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaIChaI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiIChiI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_FSIFSI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiPChiP(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SdconjSu(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHI0SHIp(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHIpconjSHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHp0SHpp(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SHppconjSHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhVWm(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZpVWm(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChaPChiP(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_ChiCha(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_ChaIChiI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SdconjSd(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SDXconjSDX(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjHpm(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHI0SHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHI0conjSHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHIpconjSHIp(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SSI0SSI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SSI0conjSSI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHp0SHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHp0conjSHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SHppconjSHpp(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSHI0conjSHI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSSI0conjSSI0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSHp0conjSHp0(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhVP(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZp(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjHpmVWm(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjVWm(E6SSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_VGVG(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZp(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZp(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZpVZp(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaPChaP(E6SSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_ChiChi(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaCha(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFDXFDX(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaIChaI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiIChiI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_FSIFSI(E6SSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiPChiP(E6SSM_mass_eigenstates_interface*, int, int, int) const;

private:
   E6SSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   E6SSM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const E6SSM_cxx_diagrams::context_base&,
      const typename E6SSM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename E6SSM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename E6SSM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(E6SSM_cxx_diagrams::context_base const& context,
                  typename E6SSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename E6SSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename E6SSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const E6SSM_cxx_diagrams::context_base&,
      typename E6SSM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename E6SSM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename E6SSM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Sd, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Sd >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Sv, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Sv >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Su, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Su >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Se, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Se >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SDX, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SDX >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Ah>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Ah>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SHI0, E6SSM_cxx_diagrams::fields::SHI0>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SHI0, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SHIp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHIp >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SSI0, E6SSM_cxx_diagrams::fields::SSI0>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SSI0 >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SSI0, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SHp0, E6SSM_cxx_diagrams::fields::SHp0>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SHp0, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::SHpp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHpp>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHpp >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHpp>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VG, E6SSM_cxx_diagrams::fields::VG>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZp, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Glu, E6SSM_cxx_diagrams::fields::Glu>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Glu >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type, E6SSM_cxx_diagrams::fields::Fv>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaP>::type, E6SSM_cxx_diagrams::fields::ChaP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaP>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChaP >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::Chi, E6SSM_cxx_diagrams::fields::Chi>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Chi >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type, E6SSM_cxx_diagrams::fields::Cha>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type, E6SSM_cxx_diagrams::fields::Fd>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type, E6SSM_cxx_diagrams::fields::Fu>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type, E6SSM_cxx_diagrams::fields::FDX>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::FDX >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type, E6SSM_cxx_diagrams::fields::ChaI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChaI >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::ChiI, E6SSM_cxx_diagrams::fields::ChiI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiI >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiI >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::FSI, E6SSM_cxx_diagrams::fields::FSI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::FSI >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::FSI >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::ChiP, E6SSM_cxx_diagrams::fields::ChiP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiP >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::Sd, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Sd >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::Se, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Se >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::SHI0, E6SSM_cxx_diagrams::fields::SHIp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHIp >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::SHIp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHIp >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::SHp0, E6SSM_cxx_diagrams::fields::SHpp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHpp >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::SHpp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHpp >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::VZp, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::ChaP, E6SSM_cxx_diagrams::fields::ChiP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChaP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiP >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::Chi, E6SSM_cxx_diagrams::fields::Cha>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Chi >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type, E6SSM_cxx_diagrams::fields::Fe>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type, E6SSM_cxx_diagrams::fields::Fd>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Hpm, E6SSM_cxx_diagrams::fields::ChaI, E6SSM_cxx_diagrams::fields::ChiI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChaI >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiI >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Sd, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Sd >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Sv, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Sv >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Su, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Su >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Se, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Se >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SDX, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SDX >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SDX>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::hh>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SHI0, E6SSM_cxx_diagrams::fields::SHI0>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SHI0, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SHIp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHIp >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHIp>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SSI0, E6SSM_cxx_diagrams::fields::SSI0>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SSI0 >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SSI0, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SSI0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SHp0, E6SSM_cxx_diagrams::fields::SHp0>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SHp0, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHp0 >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::SHpp, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHpp>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::SHpp >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHpp>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SSI0>::type >::type&) const;

template<>
Decay_amplitude_SSS E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::SHp0>::type >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::hh, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Hpm, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Hpm >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VG, E6SSM_cxx_diagrams::fields::VG>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VP, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::VZ>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VZ, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::VZp, E6SSM_cxx_diagrams::fields::VZp>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZp >::type&) const;

template<>
Decay_amplitude_SVV E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type, E6SSM_cxx_diagrams::fields::VWm>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Glu, E6SSM_cxx_diagrams::fields::Glu>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Glu >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type, E6SSM_cxx_diagrams::fields::Fv>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fv>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaP>::type, E6SSM_cxx_diagrams::fields::ChaP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaP>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChaP >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::Chi, E6SSM_cxx_diagrams::fields::Chi>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Chi >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type, E6SSM_cxx_diagrams::fields::Cha>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Cha>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type, E6SSM_cxx_diagrams::fields::Fe>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type, E6SSM_cxx_diagrams::fields::Fd>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type, E6SSM_cxx_diagrams::fields::Fu>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type, E6SSM_cxx_diagrams::fields::FDX>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::FDX>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::FDX >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type, E6SSM_cxx_diagrams::fields::ChaI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::ChaI>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChaI >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::ChiI, E6SSM_cxx_diagrams::fields::ChiI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiI >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiI >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::FSI, E6SSM_cxx_diagrams::fields::FSI>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::FSI >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::FSI >::type&) const;

template<>
Decay_amplitude_SFF E6SSM_decays::calculate_amplitude<E6SSM_cxx_diagrams::fields::Ah, E6SSM_cxx_diagrams::fields::ChiP, E6SSM_cxx_diagrams::fields::ChiP>(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::ChiP >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
E6SSM_decays::amplitude_squared(E6SSM_cxx_diagrams::context_base const& context,
                  typename E6SSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename E6SSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename E6SSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   E6SSM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
E6SSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
E6SSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
E6SSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(E6SSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
E6SSM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   E6SSM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (E6SSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   E6SSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double E6SSM_decays::get_partial_width(
   const E6SSM_cxx_diagrams::context_base& context,
   typename E6SSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename E6SSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename E6SSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,E6SSM_cxx_diagrams::fields::VZ,E6SSM_cxx_diagrams::fields::VZ >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type,E6SSM_cxx_diagrams::fields::VWm >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::conj<E6SSM_cxx_diagrams::fields::VWm>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VWm >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,E6SSM_cxx_diagrams::fields::VG,E6SSM_cxx_diagrams::fields::VG >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,E6SSM_cxx_diagrams::fields::VP,E6SSM_cxx_diagrams::fields::VP >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,E6SSM_cxx_diagrams::fields::VP,E6SSM_cxx_diagrams::fields::VZ >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type,E6SSM_cxx_diagrams::fields::Fu >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type,E6SSM_cxx_diagrams::fields::Fd >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::hh,typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type,E6SSM_cxx_diagrams::fields::Fe >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::hh >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fe >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::Ah,typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type,E6SSM_cxx_diagrams::fields::Fd >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fd>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::Ah,typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type,E6SSM_cxx_diagrams::fields::Fu >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fu>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::Ah,E6SSM_cxx_diagrams::fields::VG,E6SSM_cxx_diagrams::fields::VG >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::Ah,E6SSM_cxx_diagrams::fields::VP,E6SSM_cxx_diagrams::fields::VP >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::Ah,E6SSM_cxx_diagrams::fields::VP,E6SSM_cxx_diagrams::fields::VZ >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VP >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double E6SSM_decays::get_partial_width<E6SSM_cxx_diagrams::fields::Ah,typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type,E6SSM_cxx_diagrams::fields::Fe >(const E6SSM_cxx_diagrams::context_base&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Ah >::type&, const typename E6SSM_cxx_diagrams::field_indices<typename E6SSM_cxx_diagrams::fields::bar<E6SSM_cxx_diagrams::fields::Fe>::type >::type&, const typename E6SSM_cxx_diagrams::field_indices<E6SSM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
