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
 * @file NMSSMEFTHiggs_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#ifndef NMSSMEFTHiggs_DECAYS_H
#define NMSSMEFTHiggs_DECAYS_H

#include "NMSSMEFTHiggs_decay_table.hpp"
#include "NMSSMEFTHiggs_mass_eigenstates.hpp"
#include "NMSSMEFTHiggs_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/NMSSMEFTHiggs_qft.hpp"
#include "NMSSMEFTHiggs_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename NMSSMEFTHiggs_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename NMSSMEFTHiggs_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename NMSSMEFTHiggs_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename NMSSMEFTHiggs_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class NMSSMEFTHiggs_decays {
public:
   NMSSMEFTHiggs_decays() = default;
   NMSSMEFTHiggs_decays(NMSSMEFTHiggs_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   NMSSMEFTHiggs_decays(const NMSSMEFTHiggs_decays&) = default;
   NMSSMEFTHiggs_decays(NMSSMEFTHiggs_decays&&) = default;
   ~NMSSMEFTHiggs_decays() = default;
   NMSSMEFTHiggs_decays& operator=(const NMSSMEFTHiggs_decays&) = default;
   NMSSMEFTHiggs_decays& operator=(NMSSMEFTHiggs_decays&&) = default;

   const NMSSMEFTHiggs_decay_table& get_decay_table() const;
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

double partial_width_hh_to_SdconjSd(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjHpm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhVP(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjHpmVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_VGVG(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiChi(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChaCha(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SdconjSu(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChiCha(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SdconjSd(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhAh(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhAh(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjHpm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhVP(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVP(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjHpmVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_VGVG(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(NMSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiChi(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChaCha(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(NMSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;

private:
   NMSSMEFTHiggs_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   NMSSMEFTHiggs_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const NMSSMEFTHiggs_cxx_diagrams::context_base&,
      const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(NMSSMEFTHiggs_cxx_diagrams::context_base const& context,
                  typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const NMSSMEFTHiggs_cxx_diagrams::context_base&,
      typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Sd, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sd>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Sd >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Sv, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sv>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Sv >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Su, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Su>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Su >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Se, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Se>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Se >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::hh>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Ah>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Ah>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VP>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VP>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VG, NMSSMEFTHiggs_cxx_diagrams::fields::VG>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VP, NMSSMEFTHiggs_cxx_diagrams::fields::VP>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VP, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VZ, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Glu, NMSSMEFTHiggs_cxx_diagrams::fields::Glu>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fv>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fv>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Chi, NMSSMEFTHiggs_cxx_diagrams::fields::Chi>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Cha>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Cha>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Cha>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fe>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fd>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::hh, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fu>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::Sd, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Su>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Sd >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::Se, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sv>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Se >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::VP, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::VZ, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, NMSSMEFTHiggs_cxx_diagrams::fields::Chi, NMSSMEFTHiggs_cxx_diagrams::fields::Cha>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fe>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fv>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fd>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Sd, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sd>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Sd >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Sv, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sv>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Sv >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Su, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Su>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Su >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Se, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Se>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Se >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::hh>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::Ah>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Ah>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VP>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VP>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::hh, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VG, NMSSMEFTHiggs_cxx_diagrams::fields::VG>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VP, NMSSMEFTHiggs_cxx_diagrams::fields::VP>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VP, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::VZ, NMSSMEFTHiggs_cxx_diagrams::fields::VZ>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, NMSSMEFTHiggs_cxx_diagrams::fields::VWm>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Glu, NMSSMEFTHiggs_cxx_diagrams::fields::Glu>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fv>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fv>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, NMSSMEFTHiggs_cxx_diagrams::fields::Chi, NMSSMEFTHiggs_cxx_diagrams::fields::Chi>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Cha>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Cha>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Cha>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Cha >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fe>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fd>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF NMSSMEFTHiggs_decays::calculate_amplitude<NMSSMEFTHiggs_cxx_diagrams::fields::Ah, typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, NMSSMEFTHiggs_cxx_diagrams::fields::Fu>(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
NMSSMEFTHiggs_decays::amplitude_squared(NMSSMEFTHiggs_cxx_diagrams::context_base const& context,
                  typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(NMSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   NMSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (NMSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   NMSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double NMSSMEFTHiggs_decays::get_partial_width(
   const NMSSMEFTHiggs_cxx_diagrams::context_base& context,
   typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename NMSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,NMSSMEFTHiggs_cxx_diagrams::fields::VZ,NMSSMEFTHiggs_cxx_diagrams::fields::VZ >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type,NMSSMEFTHiggs_cxx_diagrams::fields::VWm >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::conj<NMSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,NMSSMEFTHiggs_cxx_diagrams::fields::VG,NMSSMEFTHiggs_cxx_diagrams::fields::VG >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,NMSSMEFTHiggs_cxx_diagrams::fields::VP,NMSSMEFTHiggs_cxx_diagrams::fields::VP >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,NMSSMEFTHiggs_cxx_diagrams::fields::VP,NMSSMEFTHiggs_cxx_diagrams::fields::VZ >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type,NMSSMEFTHiggs_cxx_diagrams::fields::Fu >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type,NMSSMEFTHiggs_cxx_diagrams::fields::Fd >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::hh,typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type,NMSSMEFTHiggs_cxx_diagrams::fields::Fe >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::Ah,typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type,NMSSMEFTHiggs_cxx_diagrams::fields::Fd >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::Ah,typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type,NMSSMEFTHiggs_cxx_diagrams::fields::Fu >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::Ah,NMSSMEFTHiggs_cxx_diagrams::fields::VG,NMSSMEFTHiggs_cxx_diagrams::fields::VG >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::Ah,NMSSMEFTHiggs_cxx_diagrams::fields::VP,NMSSMEFTHiggs_cxx_diagrams::fields::VP >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::Ah,NMSSMEFTHiggs_cxx_diagrams::fields::VP,NMSSMEFTHiggs_cxx_diagrams::fields::VZ >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;
template <>
double NMSSMEFTHiggs_decays::get_partial_width<NMSSMEFTHiggs_cxx_diagrams::fields::Ah,typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type,NMSSMEFTHiggs_cxx_diagrams::fields::Fe >(const NMSSMEFTHiggs_cxx_diagrams::context_base&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<typename NMSSMEFTHiggs_cxx_diagrams::fields::bar<NMSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename NMSSMEFTHiggs_cxx_diagrams::field_indices<NMSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
