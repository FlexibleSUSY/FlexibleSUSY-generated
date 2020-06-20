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

#include "E6SSMEFTHiggs_model_slha.hpp"
#include "ckm.hpp"
#include "linalg2.hpp"
#include "mixings.hpp"
#include "pmns.hpp"
#include "wrappers.hpp"

#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) this->p
#define PHYSICAL_SLHA(p) physical_slha.p
#define PHYSICAL_SLHA_REAL(p) Re(physical_slha.p)

namespace flexiblesusy {

E6SSMEFTHiggs_slha::E6SSMEFTHiggs_slha(const E6SSMEFTHiggs_input_parameters& input_,
                                   bool do_convert_masses_to_slha)
   : E6SSMEFTHiggs_mass_eigenstates(input_)
   , convert_masses_to_slha(do_convert_masses_to_slha)
{
}

/**
 * Copy constructor.  Copies from base class (model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 * @param do_convert_masses_to_slha whether to convert majorana
 *    fermion masses to SLHA convention (allow them to be negative)
 */
E6SSMEFTHiggs_slha::E6SSMEFTHiggs_slha(const E6SSMEFTHiggs_mass_eigenstates& model_, bool do_convert_masses_to_slha)
   : E6SSMEFTHiggs_mass_eigenstates(model_)
   , convert_masses_to_slha(do_convert_masses_to_slha)
{
   convert_to_slha();
}

void E6SSMEFTHiggs_slha::clear()
{
   E6SSMEFTHiggs_mass_eigenstates::clear();
   physical_slha.clear();
}

void E6SSMEFTHiggs_slha::calculate_spectrum()
{
   E6SSMEFTHiggs_mass_eigenstates::calculate_spectrum();
   convert_to_slha();
}

void E6SSMEFTHiggs_slha::convert_to_slha()
{
   physical_slha = this->get_physical();

   if (convert_masses_to_slha)
      physical_slha.convert_to_slha();

   convert_yukawa_couplings_to_slha();
   calculate_ckm_matrix();
   calculate_pmns_matrix();
   convert_trilinear_couplings_to_slha();
   convert_soft_squared_masses_to_slha();
}

void E6SSMEFTHiggs_slha::calculate_ckm_matrix()
{
   ckm = ZUL_slha * ZDL_slha.adjoint();
   CKM_parameters::to_pdg_convention(ckm, ZUL_slha, ZDL_slha, ZUR_slha, ZDR_slha);

}

void E6SSMEFTHiggs_slha::calculate_pmns_matrix()
{
   pmns << 1, 0, 0, 0, 1, 0, 0, 0, 1;

}

/**
 * Convert Yukawa couplings to SLHA convention
 */
void E6SSMEFTHiggs_slha::convert_yukawa_couplings_to_slha()
{
   fs_svd(MODELPARAMETER(Yu), Yu_slha, ZUR_slha, ZUL_slha);
   fs_svd(MODELPARAMETER(Yd), Yd_slha, ZDR_slha, ZDL_slha);
   fs_svd(MODELPARAMETER(Ye), Ye_slha, ZER_slha, ZEL_slha);

}

/**
 * Convert trilinear couplings to SLHA convention
 */
void E6SSMEFTHiggs_slha::convert_trilinear_couplings_to_slha()
{
   TYu_slha = (ZUR_slha.conjugate() * MODELPARAMETER(TYu) * ZUL_slha.adjoint()).real();
   TYd_slha = (ZDR_slha.conjugate() * MODELPARAMETER(TYd) * ZDL_slha.adjoint()).real();
   TYe_slha = (ZER_slha.conjugate() * MODELPARAMETER(TYe) * ZEL_slha.adjoint()).real();

}

/**
 * Convert soft-breaking squared mass parameters to SLHA convention
 */
void E6SSMEFTHiggs_slha::convert_soft_squared_masses_to_slha()
{
   mq2_slha = (ZDL_slha * MODELPARAMETER(mq2) * ZDL_slha.adjoint()).real();
   mu2_slha = (ZUR_slha.conjugate() * MODELPARAMETER(mu2) * ZUR_slha.transpose()).real();
   md2_slha = (ZDR_slha.conjugate() * MODELPARAMETER(md2) * ZDR_slha.transpose()).real();
   ml2_slha = (ZEL_slha * MODELPARAMETER(ml2) * ZEL_slha.adjoint()).real();
   me2_slha = (ZER_slha.conjugate() * MODELPARAMETER(me2) * ZER_slha.transpose()).real();

}

const E6SSMEFTHiggs_physical& E6SSMEFTHiggs_slha::get_physical_slha() const
{
   return physical_slha;
}

E6SSMEFTHiggs_physical& E6SSMEFTHiggs_slha::get_physical_slha()
{
   return physical_slha;
}

void E6SSMEFTHiggs_slha::print(std::ostream& ostr) const
{
   E6SSMEFTHiggs_mass_eigenstates::print(ostr);

   ostr << "----------------------------------------\n"
           "SLHA convention:\n"
           "----------------------------------------\n";
   physical_slha.print(ostr);
}

void E6SSMEFTHiggs_slha::set_convert_masses_to_slha(bool flag)
{
   convert_masses_to_slha = flag;
}

std::ostream& operator<<(std::ostream& ostr, const E6SSMEFTHiggs_slha& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
