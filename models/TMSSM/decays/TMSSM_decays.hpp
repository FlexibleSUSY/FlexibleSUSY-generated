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
 * @file TMSSM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#ifndef TMSSM_DECAYS_H
#define TMSSM_DECAYS_H

#include "TMSSM_decay_table.hpp"
#include "TMSSM_mass_eigenstates.hpp"
#include "TMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/TMSSM_qft.hpp"
#include "TMSSM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename TMSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename TMSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename TMSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename TMSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class TMSSM_decays {
public:
   TMSSM_decays() = default;
   TMSSM_decays(TMSSM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   TMSSM_decays(const TMSSM_decays&) = default;
   TMSSM_decays(TMSSM_decays&&) = default;
   ~TMSSM_decays() = default;
   TMSSM_decays& operator=(const TMSSM_decays&) = default;
   TMSSM_decays& operator=(TMSSM_decays&&) = default;

   const TMSSM_decay_table& get_decay_table() const;
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

double partial_width_hh_to_SdconjSd(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjHpm(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhVP(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjHpmVWm(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjVWm(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_VGVG(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiChi(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaCha(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SdconjSu(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhHpm(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_AhHpm(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_HpmVP(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_HpmVZ(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_hhVWm(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChiCha(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SdconjSd(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhAh(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhAh(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjHpm(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhVP(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVP(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVZ(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjHpmVWm(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjVWm(TMSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_VGVG(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(TMSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiChi(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaCha(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(TMSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(TMSSM_mass_eigenstates_interface*, int, int, int) const;

private:
   TMSSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   TMSSM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const TMSSM_cxx_diagrams::context_base&,
      const typename TMSSM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename TMSSM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename TMSSM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(TMSSM_cxx_diagrams::context_base const& context,
                  typename TMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename TMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename TMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const TMSSM_cxx_diagrams::context_base&,
      typename TMSSM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename TMSSM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename TMSSM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Sd, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sd>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Sd >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Sv, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sv>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Sv >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Su, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Su>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Su >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Se, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Se>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Se >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::hh>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Ah>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Ah>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Hpm, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Hpm, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VG, TMSSM_cxx_diagrams::fields::VG>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VP, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VP, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VZ, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Glu, TMSSM_cxx_diagrams::fields::Glu>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Glu >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fv>::type, TMSSM_cxx_diagrams::fields::Fv>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Chi, TMSSM_cxx_diagrams::fields::Chi>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Chi >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Cha>::type, TMSSM_cxx_diagrams::fields::Cha>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type, TMSSM_cxx_diagrams::fields::Fe>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type, TMSSM_cxx_diagrams::fields::Fd>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::hh, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type, TMSSM_cxx_diagrams::fields::Fu>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Sd, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Su>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Sd >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Se, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sv>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Se >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Hpm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Hpm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::VP, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::VZ, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, TMSSM_cxx_diagrams::fields::Chi, TMSSM_cxx_diagrams::fields::Cha>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Chi >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fv>::type, TMSSM_cxx_diagrams::fields::Fe>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Hpm, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type, TMSSM_cxx_diagrams::fields::Fd>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Sd, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sd>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Sd >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Sv, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sv>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Sv >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Su, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Su>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Su >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Se, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Se>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Se >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::hh>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::Ah>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Ah>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Hpm, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::hh, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Hpm, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Hpm >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VG, TMSSM_cxx_diagrams::fields::VG>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VP, TMSSM_cxx_diagrams::fields::VP>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VP, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::VZ, TMSSM_cxx_diagrams::fields::VZ>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type, TMSSM_cxx_diagrams::fields::VWm>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Glu, TMSSM_cxx_diagrams::fields::Glu>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Glu >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fv>::type, TMSSM_cxx_diagrams::fields::Fv>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, TMSSM_cxx_diagrams::fields::Chi, TMSSM_cxx_diagrams::fields::Chi>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Chi >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Cha>::type, TMSSM_cxx_diagrams::fields::Cha>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type, TMSSM_cxx_diagrams::fields::Fe>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type, TMSSM_cxx_diagrams::fields::Fd>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF TMSSM_decays::calculate_amplitude<TMSSM_cxx_diagrams::fields::Ah, typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type, TMSSM_cxx_diagrams::fields::Fu>(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fu >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
TMSSM_decays::amplitude_squared(TMSSM_cxx_diagrams::context_base const& context,
                  typename TMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename TMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename TMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   TMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
TMSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
TMSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
TMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(TMSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   TMSSM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (TMSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   TMSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// 8 -> 8, 8 with identical particles in the final state
// because of symmetry of the final state it must be proportional to d^2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
TMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
// average: 1/8
squared_color_generator() {return 40/24.;}

// 8 -> 8, 8 with differnt particles in the final state
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
TMSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
TMSSM_cxx_diagrams::fields::is_octet_v<FieldOut2> &&
!std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   f^2 = 2 Nc (-1 + Nc^2) TR = 24
// average: 1/8
squared_color_generator() {return 3.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double TMSSM_decays::get_partial_width(
   const TMSSM_cxx_diagrams::context_base& context,
   typename TMSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename TMSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename TMSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,TMSSM_cxx_diagrams::fields::VZ,TMSSM_cxx_diagrams::fields::VZ >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type,TMSSM_cxx_diagrams::fields::VWm >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::conj<TMSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VWm >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,TMSSM_cxx_diagrams::fields::VG,TMSSM_cxx_diagrams::fields::VG >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,TMSSM_cxx_diagrams::fields::VP,TMSSM_cxx_diagrams::fields::VP >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,TMSSM_cxx_diagrams::fields::VP,TMSSM_cxx_diagrams::fields::VZ >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type,TMSSM_cxx_diagrams::fields::Fu >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type,TMSSM_cxx_diagrams::fields::Fd >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::hh,typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type,TMSSM_cxx_diagrams::fields::Fe >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::hh >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fe >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::Ah,typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type,TMSSM_cxx_diagrams::fields::Fd >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::Ah,typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type,TMSSM_cxx_diagrams::fields::Fu >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::Ah,TMSSM_cxx_diagrams::fields::VG,TMSSM_cxx_diagrams::fields::VG >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::Ah,TMSSM_cxx_diagrams::fields::VP,TMSSM_cxx_diagrams::fields::VP >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::Ah,TMSSM_cxx_diagrams::fields::VP,TMSSM_cxx_diagrams::fields::VZ >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VP >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double TMSSM_decays::get_partial_width<TMSSM_cxx_diagrams::fields::Ah,typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type,TMSSM_cxx_diagrams::fields::Fe >(const TMSSM_cxx_diagrams::context_base&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Ah >::type&, const typename TMSSM_cxx_diagrams::field_indices<typename TMSSM_cxx_diagrams::fields::bar<TMSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename TMSSM_cxx_diagrams::field_indices<TMSSM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
