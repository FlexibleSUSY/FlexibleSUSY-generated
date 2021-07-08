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
 * @file MRSSMEFTHiggs_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#ifndef MRSSMEFTHiggs_DECAYS_H
#define MRSSMEFTHiggs_DECAYS_H

#include "MRSSMEFTHiggs_decay_table.hpp"
#include "MRSSMEFTHiggs_mass_eigenstates.hpp"
#include "MRSSMEFTHiggs_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/MRSSMEFTHiggs_qft.hpp"
#include "MRSSMEFTHiggs_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename MRSSMEFTHiggs_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename MRSSMEFTHiggs_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename MRSSMEFTHiggs_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename MRSSMEFTHiggs_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class MRSSMEFTHiggs_decays {
public:
   MRSSMEFTHiggs_decays() = default;
   MRSSMEFTHiggs_decays(MRSSMEFTHiggs_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   MRSSMEFTHiggs_decays(const MRSSMEFTHiggs_decays&) = default;
   MRSSMEFTHiggs_decays(MRSSMEFTHiggs_decays&&) = default;
   ~MRSSMEFTHiggs_decays() = default;
   MRSSMEFTHiggs_decays& operator=(const MRSSMEFTHiggs_decays&) = default;
   MRSSMEFTHiggs_decays& operator=(MRSSMEFTHiggs_decays&&) = default;

   const MRSSMEFTHiggs_decay_table& get_decay_table() const;
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

double partial_width_hh_to_SRdpSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SRdpHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_SRdpconjSRdp(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SRumconjSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SRumconjHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_sigmaOsigmaO(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_sigmaOphiO(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_phiOphiO(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_SdconjSd(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SvconjSv(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SuconjSu(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_SeconjSe(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhhh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhAh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_hhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhAh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_AhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_RhRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_RhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_HpmconjSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_HpmconjHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_conjSRdpconjSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjSRdpconjHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjRhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_sigmaOVG(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_phiOVG(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_hhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_RhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjRhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_hhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_AhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_RhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjRhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_SRdpVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjSRumVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjHpmVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_SRumconjVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_HpmconjVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_conjSRdpconjVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VGVG(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VPVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_VZVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_conjVWmVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_GluGlu(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barGluGlu(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barFvFv(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_ChiChi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barChiChi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_Cha1Cha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barCha1Cha1(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barCha2Cha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFeFe(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFdFd(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barFuFu(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barGlubarGlu(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_hh_to_barChibarChi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_hh_to_barCha1barCha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SRumhh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumAh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SdconjSu(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_SeconjSv(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_hhconjSRdp(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_AhconjSRdp(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_RhHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_RhconjSRdp(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_HpmconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_conjSRdpconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_SRumVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_HpmVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_conjSRdpVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_SRumVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_HpmVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_conjSRdpVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_hhVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_AhVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_RhVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_conjRhVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Hpm_to_VPVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_VZVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Hpm_to_ChiCha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barCha1Chi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barChiCha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFvFe(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barFuFd(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Hpm_to_barChibarCha1(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SRdpSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SRdpHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_SRdpconjSRdp(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SRumconjSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SRumconjHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_sigmaOsigmaO(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_sigmaOphiO(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_phiOphiO(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_SdconjSd(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SvconjSv(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SuconjSu(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_SeconjSe(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhhh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhAh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_hhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhAh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_AhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_RhRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_RhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_HpmconjSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_HpmconjHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_conjSRdpconjSRum(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjSRdpconjHpm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjRhconjRh(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_sigmaOVG(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_phiOVG(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_hhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_RhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjRhVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_hhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_AhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_RhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjRhVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_SRdpVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjSRumVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjHpmVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_SRumconjVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_HpmconjVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int) const;
double partial_width_Ah_to_conjSRdpconjVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VGVG(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVP(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VPVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_VZVZ(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_conjVWmVWm(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_GluGlu(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barGluGlu(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barFvFv(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_ChiChi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barChiChi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_Cha1Cha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barCha1Cha1(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barCha2Cha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFeFe(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFdFd(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barFuFu(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barGlubarGlu(MRSSMEFTHiggs_mass_eigenstates_interface*, int) const;
double partial_width_Ah_to_barChibarChi(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;
double partial_width_Ah_to_barCha1barCha2(MRSSMEFTHiggs_mass_eigenstates_interface*, int, int, int) const;

private:
   MRSSMEFTHiggs_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   MRSSMEFTHiggs_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const MRSSMEFTHiggs_cxx_diagrams::context_base&,
      const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(MRSSMEFTHiggs_cxx_diagrams::context_base const& context,
                  typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const MRSSMEFTHiggs_cxx_diagrams::context_base&,
      typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, MRSSMEFTHiggs_cxx_diagrams::fields::SRum>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO, MRSSMEFTHiggs_cxx_diagrams::fields::phiO>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::phiO, MRSSMEFTHiggs_cxx_diagrams::fields::phiO>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Sd, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Sd >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Sv, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Sv >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Su, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Su >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Se, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Se >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO, MRSSMEFTHiggs_cxx_diagrams::fields::VG>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::phiO, MRSSMEFTHiggs_cxx_diagrams::fields::VG>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VG, MRSSMEFTHiggs_cxx_diagrams::fields::VG>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Glu, MRSSMEFTHiggs_cxx_diagrams::fields::Glu>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Glu>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fv>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Chi, MRSSMEFTHiggs_cxx_diagrams::fields::Chi>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Chi>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1 >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fd>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fu>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, MRSSMEFTHiggs_cxx_diagrams::fields::hh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, MRSSMEFTHiggs_cxx_diagrams::fields::Ah>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Sd, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Sd >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Se, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Se >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, MRSSMEFTHiggs_cxx_diagrams::fields::Chi, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Chi>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fd>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, MRSSMEFTHiggs_cxx_diagrams::fields::SRum>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO, MRSSMEFTHiggs_cxx_diagrams::fields::phiO>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::phiO, MRSSMEFTHiggs_cxx_diagrams::fields::phiO>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Sd, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Sd >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sd>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Sv, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Sv >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Sv>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Su, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Su >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Su>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Se, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Se >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Se>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::hh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Ah>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::hh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::Rh>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&) const;

template<>
Decay_amplitude_SSS MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO, MRSSMEFTHiggs_cxx_diagrams::fields::VG>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::sigmaO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::phiO, MRSSMEFTHiggs_cxx_diagrams::fields::VG>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::phiO >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::hh, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Rh, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Rh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Rh>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRdp, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRum>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::SRum, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::SRum >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Hpm, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Hpm >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SSV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::SRdp>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VG, MRSSMEFTHiggs_cxx_diagrams::fields::VG>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::VP>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VP, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::VZ, MRSSMEFTHiggs_cxx_diagrams::fields::VZ>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type, MRSSMEFTHiggs_cxx_diagrams::fields::VWm>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Glu, MRSSMEFTHiggs_cxx_diagrams::fields::Glu>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Glu>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Glu >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fv>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fv>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Chi, MRSSMEFTHiggs_cxx_diagrams::fields::Chi>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Chi>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Chi >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1 >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2 >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fe>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fd>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type, MRSSMEFTHiggs_cxx_diagrams::fields::Fu>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Glu>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Chi>::type >::type&) const;

template<>
Decay_amplitude_SFF MRSSMEFTHiggs_decays::calculate_amplitude<MRSSMEFTHiggs_cxx_diagrams::fields::Ah, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type, typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type>(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha1>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Cha2>::type >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
MRSSMEFTHiggs_decays::amplitude_squared(MRSSMEFTHiggs_cxx_diagrams::context_base const& context,
                  typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
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
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(MRSSMEFTHiggs_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   MRSSMEFTHiggs_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (MRSSMEFTHiggs_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   MRSSMEFTHiggs_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double MRSSMEFTHiggs_decays::get_partial_width(
   const MRSSMEFTHiggs_cxx_diagrams::context_base& context,
   typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename MRSSMEFTHiggs_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
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
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,MRSSMEFTHiggs_cxx_diagrams::fields::VZ,MRSSMEFTHiggs_cxx_diagrams::fields::VZ >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type,MRSSMEFTHiggs_cxx_diagrams::fields::VWm >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::conj<MRSSMEFTHiggs_cxx_diagrams::fields::VWm>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VWm >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,MRSSMEFTHiggs_cxx_diagrams::fields::VG,MRSSMEFTHiggs_cxx_diagrams::fields::VG >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,MRSSMEFTHiggs_cxx_diagrams::fields::VP,MRSSMEFTHiggs_cxx_diagrams::fields::VP >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,MRSSMEFTHiggs_cxx_diagrams::fields::VP,MRSSMEFTHiggs_cxx_diagrams::fields::VZ >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type,MRSSMEFTHiggs_cxx_diagrams::fields::Fu >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type,MRSSMEFTHiggs_cxx_diagrams::fields::Fd >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::hh,typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type,MRSSMEFTHiggs_cxx_diagrams::fields::Fe >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::hh >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::Ah,typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type,MRSSMEFTHiggs_cxx_diagrams::fields::Fd >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fd>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fd >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::Ah,typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type,MRSSMEFTHiggs_cxx_diagrams::fields::Fu >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fu>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fu >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::Ah,MRSSMEFTHiggs_cxx_diagrams::fields::VG,MRSSMEFTHiggs_cxx_diagrams::fields::VG >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VG >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::Ah,MRSSMEFTHiggs_cxx_diagrams::fields::VP,MRSSMEFTHiggs_cxx_diagrams::fields::VP >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::Ah,MRSSMEFTHiggs_cxx_diagrams::fields::VP,MRSSMEFTHiggs_cxx_diagrams::fields::VZ >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VP >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::VZ >::type&) const;
template <>
double MRSSMEFTHiggs_decays::get_partial_width<MRSSMEFTHiggs_cxx_diagrams::fields::Ah,typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type,MRSSMEFTHiggs_cxx_diagrams::fields::Fe >(const MRSSMEFTHiggs_cxx_diagrams::context_base&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Ah >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<typename MRSSMEFTHiggs_cxx_diagrams::fields::bar<MRSSMEFTHiggs_cxx_diagrams::fields::Fe>::type >::type&, const typename MRSSMEFTHiggs_cxx_diagrams::field_indices<MRSSMEFTHiggs_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
