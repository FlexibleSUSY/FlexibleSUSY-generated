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


#ifndef lowNMSSMTanBetaAtMZ_WEINBERG_ANGLE_H
#define lowNMSSMTanBetaAtMZ_WEINBERG_ANGLE_H

#include "lowNMSSMTanBetaAtMZ_mass_eigenstates.hpp"
#include "ew_input.hpp"
#include <utility>

namespace flexiblesusy {

/**
 * @class lowNMSSMTanBetaAtMZ_weinberg_angle
 * @brief Class to calculate the DR-bar weak mixing angle from muon decay
 */
class lowNMSSMTanBetaAtMZ_weinberg_angle {
public:
   /**
    * @class Sm_parameters
    * @brief SM parameters necessary for calculating the weak mixing angle
    */
   struct Sm_parameters {
      double fermi_constant{0.}; ///< Fermi constant
      double mw_pole{0.};        ///< W pole mass
      double mz_pole{0.};        ///< Z pole mass
      double mt_pole{0.};        ///< top quark pole mass
      double mh_pole{Electroweak_constants::MH}; ///< Higgs pole mass
      double alpha_s{0.};        ///< strong coupling at Q = mt_pole
      double alpha_s_mz{0.};     ///< strong coupling at Q = mz_pole
      double dalpha_s_5_had{Electroweak_constants::delta_alpha_s_5_had}; ///< 5-flavour hadronic contributions
      int higgs_index{0};        ///< index of SM-like Higgs
   };

   lowNMSSMTanBetaAtMZ_weinberg_angle(const lowNMSSMTanBetaAtMZ_mass_eigenstates*, const Sm_parameters&);

   void set_number_of_iterations(int);       ///< set maximum number of iterations
   void set_number_of_loops(int);            ///< set number of loops
   void set_precision_goal(double);          ///< set precision goal
   void enable_dvb_bsm();                    ///< enable bsm wave, vertex and box corrections
   void disable_dvb_bsm();                   ///< disable bsm wave, vertex and box corrections
   void set_model(const lowNMSSMTanBetaAtMZ_mass_eigenstates*);  ///< set pointer to investigated model
   void set_sm_parameters(const Sm_parameters&);  ///< set sm_parameters member variable

   /// calculates and returns the sine of the Weinberg angle and the W pole mass
   std::pair<double,double> calculate(double sinThetaW_start = 0.48);
   double calculate_mw_pole() const;

private:
   int number_of_iterations{20};      ///< maximum number of iterations
   int number_of_loops{2};            ///< number of loops
   double precision_goal{1e-8};       ///< precision goal
   bool include_dvb_bsm{true};        ///< bsm wave, vertex and box corrections are included or not
   const lowNMSSMTanBetaAtMZ_mass_eigenstates* model{nullptr}; ///< pointer to investigated model
   Sm_parameters sm_parameters{};     ///< SM parameters
   double pizzt_MZ{0.};               ///< transverse Z self-energy at p^2 = MZ^2
   double piwwt_MW{0.};               ///< transverse W self-energy at p^2 = MW^2
   double piwwt_0{0.};                ///< transverse W self-energy at p^2 = 0

   double calculate_rho_hat_tree() const;
   double calculate_delta_rho_hat(double) const;
   double calculate_delta_r_hat(double, double) const;
   double calculate_delta_vb(double, double) const;
   double calculate_delta_vb_sm(double) const;
   int get_neutrino_index(int) const;
   double calculate_delta_vb_bsm(double) const;
   double calculate_mw_pole(double) const;
   double calculate_delta_alpha_hat_bsm(double) const;

   std::complex<double> CpAhHpmconjVWm(int gI2, int gI1) const;
   std::complex<double> CphhHpmconjVWm(int gI2, int gI1) const;
   std::complex<double> CpbarFvFeconjVWmPL(int gI1, int gI2) const;
   std::complex<double> CpSeconjSvconjVWm(int gI2, int gI1) const;
   std::complex<double> CpChiChaconjVWmPL(int gI1, int gI2) const;
   std::complex<double> CpChiChaconjVWmPR(int gI1, int gI2) const;
   double CpbarFvFeconjHpmPL(int , int , int ) const;
   std::complex<double> CpbarFvFeconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChabarFvSePR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFvChiSvPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const;
   double CpbarFeFvHpmPR(int , int , int ) const;
   std::complex<double> CpbarFeChaSvPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeFehhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFeChiSePR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpChaFvconjSePL(int gI1, int gO1, int gI2) const;
   std::complex<double> CpChiFvconjSvPL(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarChaFeconjSvPL(int gI1, int gO1, int gI2) const;
   std::complex<double> CpChiFeconjSePL(int gI1, int gO1, int gI2) const;
   std::complex<double> delta_vb_wave_Fv(int gO1) const;
   std::complex<double> delta_vb_wave_Fe(int gO1) const;
   std::complex<double> delta_vb_vertex(int gO1, int gO2) const;
   std::complex<double> delta_vb_box(int gO1, int gO2, int gO3, int gO4) const;

   // Passarino-Veltman loop functions
   double B0(double, double, double) const noexcept;
   double B1(double, double, double) const noexcept;
   double C0(double, double, double) const noexcept;
   double D0(double, double, double, double) const noexcept;
   double D27(double, double, double, double) const noexcept;

   static double rho_2(double);

   double calculate_self_energy_VZ(double) const;
   double calculate_self_energy_VWm(double) const;
   double calculate_self_energy_VZ_top(double, double) const;
   double calculate_self_energy_VWm_top(double, double) const;
};

} // namespace flexiblesusy

#endif
