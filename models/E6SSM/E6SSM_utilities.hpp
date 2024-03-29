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


#ifndef E6SSM_UTILITIES_H
#define E6SSM_UTILITIES_H

#include "E6SSM_info.hpp"

#include <Eigen/Core>
#include <array>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class E6SSM_mass_eigenstates;
struct E6SSM_observables;
class Physical_input;

class E6SSM_parameter_getter {
private:
   static std::vector<std::string> get_mass_names(const std::string& head = "");

   static std::array<std::string, E6SSM_info::NUMBER_OF_MIXINGS> get_mixing_names();

public:
   /// returns DR-bar parameters
   static Eigen::ArrayXd get_parameters(const E6SSM_mass_eigenstates& model);

   /// returns names of DR-bar parameters
   static std::array<std::string, E6SSM_info::NUMBER_OF_PARAMETERS> get_parameter_names();

   /// returns names of particles
   static std::array<std::string, E6SSM_info::NUMBER_OF_PARTICLES> get_particle_names();

   /// returns names of DR-bar masses
   static std::vector<std::string> get_DRbar_mass_names();

   /// returns names of pole masses
   static std::vector<std::string> get_pole_mass_names();

   /// returns names of DR-bar mixing matrices
   static std::array<std::string, E6SSM_info::NUMBER_OF_MIXINGS> get_DRbar_mixing_names();

   /// returns names of pole mixing matrices
   static std::array<std::string, E6SSM_info::NUMBER_OF_MIXINGS> get_pole_mixing_names();

   /// returns names of input parameters
   static std::array<std::string, E6SSM_info::NUMBER_OF_INPUT_PARAMETERS> get_input_parameter_names();

   /// returns names of input parameters
   static std::array<std::string, E6SSM_info::NUMBER_OF_EXTRA_PARAMETERS> get_extra_parameter_names();

   /// returns number of particle masses
   static decltype(E6SSM_info::NUMBER_OF_MASSES) get_number_of_masses();
};

class E6SSM_spectrum_plotter {
public:
   E6SSM_spectrum_plotter() = default;
   explicit E6SSM_spectrum_plotter(const E6SSM_mass_eigenstates&);
   void extract_spectrum(const E6SSM_mass_eigenstates&);
   void write_to_file(const std::string&) const;

private:
   struct TParticle {
      std::string name;
      std::string latex_name;
      std::valarray<double> masses;
      TParticle(const std::string& name_, const std::string& latex_name_,
                const std::valarray<double>& masses_)
         : name(name_)
         , latex_name(latex_name_)
         , masses(masses_)
         {}
   };
   using TSpectrum = std::vector<TParticle>;
   TSpectrum spectrum{};
   double scale{0.};
   int width{16};

   void write_spectrum(const TSpectrum&, std::ofstream&) const;
};

namespace E6SSM_database {

/// append parameter point to database
void to_database(
   const std::string&,
   const E6SSM_mass_eigenstates&,
   const softsusy::QedQcd* qedqcd = nullptr,
   const Physical_input* physical_input = nullptr,
   const E6SSM_observables* observables = nullptr);

/// fill model from an entry of the database
E6SSM_mass_eigenstates from_database(
   const std::string&,
   long long,
   softsusy::QedQcd* qedqcd = nullptr,
   Physical_input* physical_input = nullptr,
   E6SSM_observables* observables = nullptr);

} // namespace E6SSM_database

} // namespace flexiblesusy

#endif
