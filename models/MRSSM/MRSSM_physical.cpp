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

// File generated at Sun 31 May 2015 12:28:19

#include "MRSSM_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

MRSSM_physical::MRSSM_physical()
   :
    MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MsigmaO(0), MphiO(
       0), MVP(0), MVZ(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(Eigen::Array
       <double,3,1>::Zero()), MSu(Eigen::Array<double,6,1>::Zero()), MSe(
       Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,4,1>::Zero()),
       MAh(Eigen::Array<double,4,1>::Zero()), MRh(Eigen::Array<double,2,1>::Zero()
       ), MHpm(Eigen::Array<double,4,1>::Zero()), MRpm(Eigen::Array<double,2,1>
       ::Zero()), MChi(Eigen::Array<double,4,1>::Zero()), MCha1(Eigen::Array<
       double,2,1>::Zero()), MCha2(Eigen::Array<double,2,1>::Zero()), MFe(
       Eigen::Array<double,3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()),
       MFu(Eigen::Array<double,3,1>::Zero()), MVWm(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,4,4>::Zero()), ZA(Eigen::Matrix<double,4,
      4>::Zero()), ZHR(Eigen::Matrix<double,2,2>::Zero()), ZP(Eigen::Matrix<double
      ,4,4>::Zero()), ZRP(Eigen::Matrix<double,2,2>::Zero()), ZN1(Eigen::Matrix<
      std::complex<double>,4,4>::Zero()), ZN2(Eigen::Matrix<std::complex<double>,4
      ,4>::Zero()), UM1(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP1(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), UM2(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), UP2(Eigen::Matrix<std::complex<double>,2
      ,2>::Zero()), ZEL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZUL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(
      Eigen::Matrix<std::complex<double>,3,3>::Zero())

{
}

void MRSSM_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MsigmaO = 0.;
   MphiO = 0.;
   MVP = 0.;
   MVZ = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,4,1>::Zero();
   ZH = Eigen::Matrix<double,4,4>::Zero();
   MAh = Eigen::Matrix<double,4,1>::Zero();
   ZA = Eigen::Matrix<double,4,4>::Zero();
   MRh = Eigen::Matrix<double,2,1>::Zero();
   ZHR = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,4,1>::Zero();
   ZP = Eigen::Matrix<double,4,4>::Zero();
   MRpm = Eigen::Matrix<double,2,1>::Zero();
   ZRP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN1 = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   ZN2 = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha1 = Eigen::Matrix<double,2,1>::Zero();
   UM1 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP1 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MCha2 = Eigen::Matrix<double,2,1>::Zero();
   UM2 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP2 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWm = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void MRSSM_physical::convert_to_hk()
{

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void MRSSM_physical::convert_to_slha()
{

}

void MRSSM_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MsigmaO = " << MsigmaO << '\n';
   ostr << "MphiO = " << MphiO << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MRh = " << MRh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MRpm = " << MRpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha1 = " << MCha1.transpose() << '\n';
   ostr << "MCha2 = " << MCha2.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZHR = " << ZHR << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZRP = " << ZRP << '\n';
   ostr << "ZN1 = " << ZN1 << '\n';
   ostr << "ZN2 = " << ZN2 << '\n';
   ostr << "UM1 = " << UM1 << '\n';
   ostr << "UP1 = " << UP1 << '\n';
   ostr << "UM2 = " << UM2 << '\n';
   ostr << "UP2 = " << UP2 << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const MRSSM_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
