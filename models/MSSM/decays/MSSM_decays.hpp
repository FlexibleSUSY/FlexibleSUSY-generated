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
 * @file MSSM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#ifndef MSSM_DECAYS_H
#define MSSM_DECAYS_H

#include "MSSM_decay_table.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "MSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/MSSM_qft.hpp"
#include "MSSM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename MSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename MSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename MSSM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename MSSM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class MSSM_decays {
public:
   MSSM_decays() = default;
   MSSM_decays(MSSM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   MSSM_decays(const MSSM_decays&) = default;
   MSSM_decays(MSSM_decays&&) = default;
   ~MSSM_decays() = default;
   MSSM_decays& operator=(const MSSM_decays&) = default;
   MSSM_decays& operator=(MSSM_decays&&) = default;

   const MSSM_decay_table& get_decay_table() const;
   const FlexibleDecay_problems& get_problems() const;

   void clear();
   void clear_problems();
   void calculate_decays();

   const Decays_list& get_hh_decays(int i) const { return decay_table.
      get_hh_decays(i); }
   const Decays_list& get_Ah_decays(int i) const { return decay_table.
      get_Ah_decays(i); }
   const Decays_list& get_Hpm_decays(int i) const { return decay_table.
      get_Hpm_decays(i); }
   const Decays_list& get_Su_decays(int i) const { return decay_table.
      get_Su_decays(i); }
   const Decays_list& get_Sd_decays(int i) const { return decay_table.
      get_Sd_decays(i); }
   const Decays_list& get_Se_decays(int i) const { return decay_table.
      get_Se_decays(i); }
   const Decays_list& get_Sv_decays(int i) const { return decay_table.
      get_Sv_decays(i); }
   void calculate_hh_decays();
   void calculate_Ah_decays();
   void calculate_Hpm_decays();
   void calculate_Su_decays();
   void calculate_Sd_decays();
   void calculate_Se_decays();
   void calculate_Sv_decays();

double partial_width_hh_to_SdconjSd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjHpm(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjHpmVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_VGVG(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiChi(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaCha(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SdconjSd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjHpm(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjHpmVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_VGVG(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiChi(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaCha(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SdconjSu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(MSSM_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChiCha(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Su_to_SdconjHpm(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Su_to_Suhh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Su_to_SuAh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Su_to_SuVG(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Su_to_SuVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Su_to_SuVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Su_to_SdconjVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Su_to_GluFu(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Su_to_ChiFu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Su_to_barChaFd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sd_to_Sdhh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sd_to_SdAh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sd_to_SuHpm(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sd_to_SdVG(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sd_to_SdVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sd_to_SdVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sd_to_SuVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sd_to_GluFd(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sd_to_ChiFd(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sd_to_ChaFu(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Se_to_SvHpm(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Se_to_Sehh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Se_to_SeAh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Se_to_SeVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Se_to_SeVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Se_to_SvVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Se_to_FvCha(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Se_to_ChiFe(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sv_to_Svhh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sv_to_SvAh(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sv_to_SeconjHpm(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sv_to_SvVP(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sv_to_SvVZ(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sv_to_SeconjVWm(MSSM_mass_eigenstates_interface*, int, int) const;
double partial_width_Sv_to_FvChi(MSSM_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Sv_to_barChaFe(MSSM_mass_eigenstates_interface*, int, int, int) const;

private:
   MSSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   MSSM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const MSSM_cxx_diagrams::context_base&,
      const typename MSSM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename MSSM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename MSSM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(MSSM_cxx_diagrams::context_base const& context,
                  typename MSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename MSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename MSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const MSSM_cxx_diagrams::context_base&,
      typename MSSM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename MSSM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename MSSM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Sd, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sd>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Sv, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sv>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Su, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Su>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Se, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Se>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::hh>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Ah>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Ah>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Hpm, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Hpm, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VG, MSSM_cxx_diagrams::fields::VG>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VP, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VP, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VZ, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Glu, MSSM_cxx_diagrams::fields::Glu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Glu >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fv>::type, MSSM_cxx_diagrams::fields::Fv>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::Chi, MSSM_cxx_diagrams::fields::Chi>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type, MSSM_cxx_diagrams::fields::Cha>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type, MSSM_cxx_diagrams::fields::Fe>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type, MSSM_cxx_diagrams::fields::Fd>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::hh, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type, MSSM_cxx_diagrams::fields::Fu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Sd, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sd>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Sv, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sv>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Su, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Su>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Se, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Se>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::hh>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Hpm, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Hpm, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VG, MSSM_cxx_diagrams::fields::VG>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VP, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VP, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VZ, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Glu, MSSM_cxx_diagrams::fields::Glu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Glu >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fv>::type, MSSM_cxx_diagrams::fields::Fv>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::Chi, MSSM_cxx_diagrams::fields::Chi>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type, MSSM_cxx_diagrams::fields::Cha>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type, MSSM_cxx_diagrams::fields::Fe>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type, MSSM_cxx_diagrams::fields::Fd>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Ah, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type, MSSM_cxx_diagrams::fields::Fu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::Sd, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Su>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::Se, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sv>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::hh, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::Ah, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::VP, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::VZ, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, MSSM_cxx_diagrams::fields::Chi, MSSM_cxx_diagrams::fields::Cha>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fv>::type, MSSM_cxx_diagrams::fields::Fe>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fv>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Hpm, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type, MSSM_cxx_diagrams::fields::Fd>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Sd, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::hh>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Ah>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::VG>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Sd, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Glu, MSSM_cxx_diagrams::fields::Fu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Glu >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Chi, MSSM_cxx_diagrams::fields::Fu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Su, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type, MSSM_cxx_diagrams::fields::Fd>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::hh>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Ah>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::Hpm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::VG>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Su, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Su >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Glu, MSSM_cxx_diagrams::fields::Fd>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Glu >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Chi, MSSM_cxx_diagrams::fields::Fd>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sd, MSSM_cxx_diagrams::fields::Cha, MSSM_cxx_diagrams::fields::Fu>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sd >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Cha >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Hpm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::hh>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Ah>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::VWm>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Fv, MSSM_cxx_diagrams::fields::Cha>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Se, MSSM_cxx_diagrams::fields::Chi, MSSM_cxx_diagrams::fields::Fe>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::hh>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Ah>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Se, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::VP>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::VZ>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Se, typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Se >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, MSSM_cxx_diagrams::fields::Fv, MSSM_cxx_diagrams::fields::Chi>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fv >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MSSM_decays::calculate_amplitude<MSSM_cxx_diagrams::fields::Sv, typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type, MSSM_cxx_diagrams::fields::Fe>(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Sv >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Cha>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
MSSM_decays::amplitude_squared(MSSM_cxx_diagrams::context_base const& context,
                  typename MSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename MSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename MSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
MSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   MSSM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MSSM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
MSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
MSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
MSSM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(MSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(MSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MSSM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((MSSM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MSSM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MSSM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   MSSM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (MSSM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   MSSM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double MSSM_decays::get_partial_width(
   const MSSM_cxx_diagrams::context_base& context,
   typename MSSM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename MSSM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename MSSM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,MSSM_cxx_diagrams::fields::VZ,MSSM_cxx_diagrams::fields::VZ >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type,MSSM_cxx_diagrams::fields::VWm >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::conj<MSSM_cxx_diagrams::fields::VWm>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VWm >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,MSSM_cxx_diagrams::fields::VG,MSSM_cxx_diagrams::fields::VG >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,MSSM_cxx_diagrams::fields::VP,MSSM_cxx_diagrams::fields::VP >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,MSSM_cxx_diagrams::fields::VP,MSSM_cxx_diagrams::fields::VZ >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type,MSSM_cxx_diagrams::fields::Fu >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type,MSSM_cxx_diagrams::fields::Fd >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::hh,typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type,MSSM_cxx_diagrams::fields::Fe >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::hh >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::Ah,typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type,MSSM_cxx_diagrams::fields::Fd >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fd>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::Ah,typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type,MSSM_cxx_diagrams::fields::Fu >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fu>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::Ah,MSSM_cxx_diagrams::fields::VG,MSSM_cxx_diagrams::fields::VG >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VG >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::Ah,MSSM_cxx_diagrams::fields::VP,MSSM_cxx_diagrams::fields::VP >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::Ah,MSSM_cxx_diagrams::fields::VP,MSSM_cxx_diagrams::fields::VZ >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VP >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MSSM_decays::get_partial_width<MSSM_cxx_diagrams::fields::Ah,typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type,MSSM_cxx_diagrams::fields::Fe >(const MSSM_cxx_diagrams::context_base&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Ah >::type&, const typename MSSM_cxx_diagrams::field_indices<typename MSSM_cxx_diagrams::fields::bar<MSSM_cxx_diagrams::fields::Fe>::type >::type&, const typename MSSM_cxx_diagrams::field_indices<MSSM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
