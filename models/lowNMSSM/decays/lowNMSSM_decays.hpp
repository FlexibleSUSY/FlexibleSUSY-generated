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
 * @file lowNMSSM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef lowNMSSM_DECAYS_H
#define lowNMSSM_DECAYS_H

#include "lowNMSSM_decay_table.hpp"
#include "lowNMSSM_mass_eigenstates.hpp"
#include "lowNMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/lowNMSSM_qft.hpp"
#include "lowNMSSM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename lowNMSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename lowNMSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename lowNMSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename lowNMSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class lowNMSSM_decays {
public:
   lowNMSSM_decays() = default;
   lowNMSSM_decays(lowNMSSM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   lowNMSSM_decays(const lowNMSSM_decays&) = default;
   lowNMSSM_decays(lowNMSSM_decays&&) = default;
   ~lowNMSSM_decays() = default;
   lowNMSSM_decays& operator=(const lowNMSSM_decays&) = default;
   lowNMSSM_decays& operator=(lowNMSSM_decays&&) = default;

   const lowNMSSM_decay_table& get_decay_table() const;
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

double partial_width_hh_to_SdconjSd(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjHpm(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhVP(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjHpmVWm(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjVWm(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_VGVG(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiChi(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaCha(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SdconjSu(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhVWm(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChiCha(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SdconjSd(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhAh(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhAh(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjHpm(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhVP(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVP(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVZ(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjHpmVWm(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjVWm(lowNMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_VGVG(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(lowNMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiChi(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaCha(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(lowNMSSM_mass_eigenstates_interface*, int, int, int) const;

private:
   lowNMSSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   lowNMSSM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const lowNMSSM_cxx_diagrams::context_base&,
      const typename lowNMSSM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename lowNMSSM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename lowNMSSM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(lowNMSSM_cxx_diagrams::context_base const& context,
                  typename lowNMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename lowNMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename lowNMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const lowNMSSM_cxx_diagrams::context_base&,
      typename lowNMSSM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename lowNMSSM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename lowNMSSM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Sd, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sd>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Sd >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Sv, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sv>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Sv >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Su, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Su>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Su >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Se, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Se>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Se >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::hh>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Ah>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Ah>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Hpm, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VP>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VP>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Hpm, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VG, lowNMSSM_cxx_diagrams::fields::VG>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VP, lowNMSSM_cxx_diagrams::fields::VP>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VP, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VZ, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Glu, lowNMSSM_cxx_diagrams::fields::Glu>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Glu >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fv>::type, lowNMSSM_cxx_diagrams::fields::Fv>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Chi, lowNMSSM_cxx_diagrams::fields::Chi>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Chi >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Cha>::type, lowNMSSM_cxx_diagrams::fields::Cha>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type, lowNMSSM_cxx_diagrams::fields::Fe>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type, lowNMSSM_cxx_diagrams::fields::Fd>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::hh, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type, lowNMSSM_cxx_diagrams::fields::Fu>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::Sd, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Su>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Sd >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::Se, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sv>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Se >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::VP, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::VZ, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, lowNMSSM_cxx_diagrams::fields::Chi, lowNMSSM_cxx_diagrams::fields::Cha>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Chi >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fv>::type, lowNMSSM_cxx_diagrams::fields::Fe>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Hpm, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type, lowNMSSM_cxx_diagrams::fields::Fd>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Sd, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sd>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Sd >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Sv, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sv>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Sv >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Su, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Su>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Su >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Se, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Se>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Se >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::hh>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::Ah>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Ah>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Hpm, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VP>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VP>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::hh, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Hpm, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Hpm >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VG, lowNMSSM_cxx_diagrams::fields::VG>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VP, lowNMSSM_cxx_diagrams::fields::VP>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VP, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::VZ, lowNMSSM_cxx_diagrams::fields::VZ>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type, lowNMSSM_cxx_diagrams::fields::VWm>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Glu, lowNMSSM_cxx_diagrams::fields::Glu>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Glu >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fv>::type, lowNMSSM_cxx_diagrams::fields::Fv>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, lowNMSSM_cxx_diagrams::fields::Chi, lowNMSSM_cxx_diagrams::fields::Chi>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Chi >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Cha>::type, lowNMSSM_cxx_diagrams::fields::Cha>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type, lowNMSSM_cxx_diagrams::fields::Fe>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type, lowNMSSM_cxx_diagrams::fields::Fd>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF lowNMSSM_decays::calculate_amplitude<lowNMSSM_cxx_diagrams::fields::Ah, typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type, lowNMSSM_cxx_diagrams::fields::Fu>(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fu >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
lowNMSSM_decays::amplitude_squared(lowNMSSM_cxx_diagrams::context_base const& context,
                  typename lowNMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename lowNMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename lowNMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(lowNMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (lowNMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   lowNMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// 8 -> 8, 8 with identical particles in the final state
// because of symmetry of the final state it must be proportional to d^2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
// average: 1/8
squared_color_generator() {return 40/24.;}

// 8 -> 8, 8 with differnt particles in the final state
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
lowNMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
!std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   f^2 = 2 Nc (-1 + Nc^2) TR = 24
// average: 1/8
squared_color_generator() {return 3.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double lowNMSSM_decays::get_partial_width(
   const lowNMSSM_cxx_diagrams::context_base& context,
   typename lowNMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename lowNMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename lowNMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,lowNMSSM_cxx_diagrams::fields::VZ,lowNMSSM_cxx_diagrams::fields::VZ >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type,lowNMSSM_cxx_diagrams::fields::VWm >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::conj<lowNMSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VWm >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,lowNMSSM_cxx_diagrams::fields::VG,lowNMSSM_cxx_diagrams::fields::VG >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,lowNMSSM_cxx_diagrams::fields::VP,lowNMSSM_cxx_diagrams::fields::VP >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,lowNMSSM_cxx_diagrams::fields::VP,lowNMSSM_cxx_diagrams::fields::VZ >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type,lowNMSSM_cxx_diagrams::fields::Fu >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type,lowNMSSM_cxx_diagrams::fields::Fd >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::hh,typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type,lowNMSSM_cxx_diagrams::fields::Fe >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::hh >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fe >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::Ah,typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type,lowNMSSM_cxx_diagrams::fields::Fd >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::Ah,typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type,lowNMSSM_cxx_diagrams::fields::Fu >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::Ah,lowNMSSM_cxx_diagrams::fields::VG,lowNMSSM_cxx_diagrams::fields::VG >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::Ah,lowNMSSM_cxx_diagrams::fields::VP,lowNMSSM_cxx_diagrams::fields::VP >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::Ah,lowNMSSM_cxx_diagrams::fields::VP,lowNMSSM_cxx_diagrams::fields::VZ >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VP >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double lowNMSSM_decays::get_partial_width<lowNMSSM_cxx_diagrams::fields::Ah,typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type,lowNMSSM_cxx_diagrams::fields::Fe >(const lowNMSSM_cxx_diagrams::context_base&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Ah >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<typename lowNMSSM_cxx_diagrams::fields::bar<lowNMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename lowNMSSM_cxx_diagrams::field_indices<lowNMSSM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
