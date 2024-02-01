! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/basis/qvszp.f90
!>
!> Implements the charge dependent single-zeta q-vSZP basis set based on: 
!> Marcel Müller, Andreas Hansen, and Stefan Grimme, 
!> "An atom-in-molecule adaptive polarized valence single-ζ 
!> atomic orbital basis for electronic structure calculations"
!> J. Chem. Phys. 159, 164108 (2023). DOI: 10.1063/5.0172373
module tblite_basis_qvszp
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type, integral_cutoff
   use tblite_basis_type, only : cgto_type
   implicit none
   private

   public :: new_basis, add_qvszp_cgtos
   !> Maximum contraction length of basis functions.
   !> The limit is chosen as twice the maximum size returned by the STO-NG expansion
   integer, parameter :: maxg = 12

   !> Maximal number of primitives per CGTO
   integer, parameter :: max_prim = 7
   !> Maximal number of elements
   integer, parameter :: max_elem = 86
   !> Maximal number of shells
   integer, parameter :: max_shell = 3

   !> Number of functions
   integer, parameter :: nf = 15
   !> Two over pi
   real(wp), parameter :: top = 2.0_wp / pi
   !> Double factorial, see OEIS A001147
   real(wp), parameter :: dfactorial(8) = &
      & [1.0_wp,1.0_wp,3.0_wp,15.0_wp,105.0_wp,945.0_wp,10395.0_wp,135135.0_wp]

   !> Contracted Gaussian type basis function of qvszp type
   type, public, extends(cgto_type) :: qvszp_cgto_type
      !> Derivative of the Contraction coefficient w.r.t. R
      real(wp) :: dcoeff(maxg) = 0.0_wp
      !> Contraction coefficients of the primitive Gaussian functions,
      real(wp) :: coeff0(maxg) = 0.0_wp
      !> Primitive-specific strength of the environment specific scaling
      real(wp) :: coeff1(maxg) = 0.0_wp
      !> Atom-specific parameter scaling the charge-dependence of the coefficient
      real(wp) :: k1 = 0.0_wp
      !> Atom-specific parameter scaling the CN-dependence of the coefficient
      real(wp) :: k2 = 0.0_wp
      !> Global parameter scaling the mixed charge- and CN-dependence of the coefficient
      real(wp) :: k3 = 0.0_wp
   end type qvszp_cgto_type

   !> Collection of information regarding the basis set of a system
   type, public, extends(basis_type) :: qvszp_basis_type 
   end type qvszp_basis_type


   integer, parameter :: nshell(max_elem) = [ &
   & 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, &
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, &
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 1, 1, 1, &
   & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
   & 3, 3, 3, 3, 3, 3]

   integer, parameter :: n_prim(max_elem, max_shell) = reshape([&
   & 8, 3, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 5, 4, 0, 0, 0, 0, 0, & ! up to element: 3
   & 6, 4, 0, 0, 0, 0, 0, 6, 5, 0, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 0, & ! up to element: 6
   & 6, 6, 0, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 0, & ! up to element: 9
   & 6, 6, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 4, 3, 2, 0, 0, 0, 0, & ! up to element: 12
   & 5, 4, 2, 0, 0, 0, 0, 6, 4, 2, 0, 0, 0, 0, 5, 5, 2, 0, 0, 0, 0, & ! up to element: 15
   & 5, 5, 2, 0, 0, 0, 0, 5, 5, 2, 0, 0, 0, 0, 5, 5, 2, 0, 0, 0, 0, & ! up to element: 18
   & 4, 3, 0, 0, 0, 0, 0, 4, 3, 3, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, & ! up to element: 21
   & 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, & ! up to element: 24
   & 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, & ! up to element: 27
   & 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, 5, 3, 0, 0, 0, 0, 0, & ! up to element: 30
   & 5, 4, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, & ! up to element: 33
   & 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, & ! up to element: 36
   & 4, 3, 0, 0, 0, 0, 0, 4, 3, 3, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, & ! up to element: 39
   & 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, & ! up to element: 42
   & 5, 2, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, & ! up to element: 45
   & 5, 3, 6, 0, 0, 0, 0, 5, 3, 6, 0, 0, 0, 0, 5, 3, 2, 0, 0, 0, 0, & ! up to element: 48
   & 5, 4, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, & ! up to element: 51
   & 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, & ! up to element: 54
   & 4, 3, 0, 0, 0, 0, 0, 4, 3, 3, 0, 0, 0, 0, 4, 2, 5, 0, 0, 0, 0, & ! up to element: 57
   & 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, & ! up to element: 60
   & 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, & ! up to element: 63
   & 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, & ! up to element: 66
   & 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, & ! up to element: 69
   & 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 4, 3, 5, 0, 0, 0, 0, & ! up to element: 72
   & 5, 3, 5, 0, 0, 0, 0, 5, 3, 5, 0, 0, 0, 0, 5, 3, 5, 0, 0, 0, 0, & ! up to element: 75
   & 5, 3, 5, 0, 0, 0, 0, 5, 3, 5, 0, 0, 0, 0, 5, 3, 5, 0, 0, 0, 0, & ! up to element: 78
   & 5, 3, 5, 0, 0, 0, 0, 5, 3, 2, 0, 0, 0, 0, 5, 4, 2, 0, 0, 0, 0, & ! up to element: 81
   & 5, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0, & ! up to element: 84
   & 6, 5, 2, 0, 0, 0, 0, 6, 5, 2, 0, 0, 0, 0], shape(n_prim)) 

   integer, parameter :: ang_shell(max_elem, max_shell) = reshape([&
   & 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, & ! up to element: 3
   & 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, & ! up to element: 6
   & 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, & ! up to element: 9
   & 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 12
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 15
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 18
   & 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 21
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 24
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 27
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, & ! up to element: 30
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 33
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 36
   & 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 39
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 42
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 45
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 48
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 51
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 54
   & 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 57
   & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! up to element: 60
   & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! up to element: 63
   & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! up to element: 66
   & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! up to element: 69
   & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 72
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 75
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 78
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 81
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, & ! up to element: 84
   & 0, 1, 2, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0], shape(ang_shell))

   real(wp), parameter :: p_k1(86) = [&
   &    0.4617920882_wp,    0.8057551376_wp,    0.8345391366_wp,    0.9352502005_wp, & ! -4
   &    0.7374526520_wp,    0.9201468994_wp,    0.1998619300_wp,    0.3873692821_wp, & ! -8
   &    0.2182421148_wp,    0.0546184554_wp,    0.6729111637_wp,    2.8289503248_wp, & ! -12
   &    1.8715641031_wp,    1.3825557914_wp,    0.2784619343_wp,    0.1677047996_wp, & ! -16
   &    0.3419763327_wp,    0.3551862388_wp,    0.9054671925_wp,    3.7657222920_wp, & ! -20
   &    6.3138005556_wp,    2.3435161030_wp,    1.8406529083_wp,    0.9003587849_wp, & ! -24
   &    2.6962491168_wp,    0.6521680519_wp,    0.8224560018_wp,    0.0363835510_wp, & ! -28
   &    2.0113709201_wp,    1.3348285691_wp,    3.4835976689_wp,    1.8461172060_wp, & ! -32
   &    0.2325456263_wp,    0.2766785064_wp,    0.3788137064_wp,    0.4102186578_wp, & ! -36
   &    1.5190059967_wp,    5.9675064723_wp,   10.5283425352_wp,    1.2818353604_wp, & ! -40
   &    2.6822755248_wp,    0.9124266051_wp,   11.7334043524_wp,    1.4552418701_wp, & ! -44
   &    4.5786060844_wp,   -0.5025607912_wp,    3.0992007331_wp,    2.3129482626_wp, & ! -48
   &    1.7994519869_wp,    1.4293697892_wp,    0.4163692776_wp,    0.5259871447_wp, & ! -52
   &    0.6225453181_wp,    0.3340518516_wp,    1.7058377174_wp,    7.2594862271_wp, & ! -56
   &   16.9619460644_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -60
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -64
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -68
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    1.3023989549_wp, & ! -72
   &    2.2999117078_wp,    1.9875510656_wp,    4.8555129660_wp,    0.7289354285_wp, & ! -76
   &    1.4025407706_wp,    1.0915811660_wp,    0.5758670225_wp,    2.6169472024_wp, & ! -80
   &    2.2071772428_wp,    1.2614862277_wp,    0.9806462683_wp,    0.5731434555_wp, & ! -84
   &    0.6474666807_wp,    2.1015118685_wp]

   real(wp), parameter :: p_k2(86) = [&
   &    0.1606087330_wp,   -0.0363837083_wp,   -0.7587353255_wp,    0.2639711545_wp, & ! -4
   &    0.1472401862_wp,    0.0467082703_wp,    0.2551447257_wp,   -0.0320674298_wp, & ! -8
   &   -0.0127660402_wp,    0.0220132140_wp,   -0.3310736225_wp,   -0.4490921824_wp, & ! -12
   &    0.2615371148_wp,    0.0315036442_wp,    0.2918038834_wp,    0.1991241881_wp, & ! -16
   &   -0.0253627196_wp,    0.0098630949_wp,   -0.2956604942_wp,    0.0845845631_wp, & ! -20
   &   -0.8051851949_wp,   -0.3615693527_wp,    0.2057919905_wp,    0.6286786148_wp, & ! -24
   &    0.5252884020_wp,    0.4304116762_wp,    0.3763488499_wp,    0.0613884486_wp, & ! -28
   &    0.0193253403_wp,    0.3758993048_wp,    0.1140403624_wp,    0.0784027913_wp, & ! -32
   &    0.5206927335_wp,    0.2545917256_wp,    0.0195973022_wp,    0.0047452656_wp, & ! -36
   &   -0.5584320885_wp,    2.8567145885_wp,   -0.6903215172_wp,   -0.0377926748_wp, & ! -40
   &    0.1193001538_wp,    0.6390784511_wp,    0.8367669986_wp,    0.2761479903_wp, & ! -44
   &    0.2378114224_wp,   -0.1911775305_wp,    0.1147030576_wp,    0.5451792336_wp, & ! -48
   &    0.0728390152_wp,    0.1212865540_wp,    0.5803216588_wp,    0.2882676174_wp, & ! -52
   &    0.0382709261_wp,    0.0948718823_wp,   -0.7750240287_wp,    0.3568416480_wp, & ! -56
   &   -0.0174939003_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -60
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -64
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -68
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,   -0.3775706870_wp, & ! -72
   &    0.0277595072_wp,   -0.0374034897_wp,   -0.0688439928_wp,   -0.0683856481_wp, & ! -76
   &   -0.0130270022_wp,    0.8396809160_wp,   -0.1477660485_wp,    0.0466701560_wp, & ! -80
   &    0.0456239770_wp,    0.1072620782_wp,    0.1946437724_wp,    0.3596951321_wp, & ! -84
   &    0.0569460627_wp,    0.0019532518_wp]

   real(wp), parameter :: p_k3(86) = [&
   &   -0.3554218004_wp,    0.0621617162_wp,    0.0356870097_wp,    0.4265235744_wp, & ! -4
   &   -0.0997140767_wp,   -0.0388421664_wp,   -0.0206605361_wp,   -0.0210898295_wp, & ! -8
   &    0.0333257503_wp,   -0.0235939405_wp,   -0.1896374365_wp,    0.0715400075_wp, & ! -12
   &    0.0205039808_wp,    0.0101036400_wp,    0.0487293183_wp,   -0.0701037722_wp, & ! -16
   &   -0.0276715176_wp,   -0.0638928630_wp,    0.0550363489_wp,    0.0231904439_wp, & ! -20
   &    0.2637721778_wp,    0.0350118510_wp,   -0.0535516981_wp,   -0.2512403937_wp, & ! -24
   &   -0.3675145209_wp,   -0.2126762024_wp,   -0.1957737999_wp,    1.4069644813_wp, & ! -28
   &   -0.5414717014_wp,   -0.1650303356_wp,    0.1523555425_wp,    0.0454393729_wp, & ! -32
   &    0.1420448115_wp,   -0.0301446752_wp,    0.0497931804_wp,    0.1403022143_wp, & ! -36
   &    0.1378014984_wp,    0.1226099841_wp,    0.7241025569_wp,    0.0137677251_wp, & ! -40
   &    0.0220997343_wp,    0.0150237080_wp,   -0.4439942620_wp,   -0.0168728417_wp, & ! -44
   &    0.0325055850_wp,    0.3477903062_wp,    0.0185097447_wp,    0.0575171291_wp, & ! -48
   &    0.0217778288_wp,    0.0947835658_wp,   -0.0295816729_wp,   -0.0712062447_wp, & ! -52
   &    0.0183619481_wp,    0.1608275620_wp,    0.2765529721_wp,    0.0101489920_wp, & ! -56
   &    2.4324136854_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -60
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -64
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! -68
   &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0721523731_wp, & ! -72
   &   -0.0350627611_wp,   -0.0560672732_wp,   -0.0686776453_wp,    0.0143641773_wp, & ! -76
   &   -0.1359450351_wp,    0.8420828052_wp,   -0.0110374429_wp,    0.0244227610_wp, & ! -80
   &    0.0200322157_wp,    0.0992535713_wp,    0.0124519511_wp,    0.0603826428_wp, & ! -84
   &   -0.0797804279_wp,   -0.7043039945_wp]

   !> CGTO exponents (Exponent of the primitive Gaussian functions)
   real(wp), protected :: exponents(max_prim, max_shell, max_elem) = 0.0_wp

   !> CGTO coefficients (Contraction coefficients of the primitive Gaussian functions,
   !> might contain normalization)
   real(wp), protected :: coefficients(max_prim, max_shell, max_elem) = 0.0_wp

   !> CGTO coefficients (Contraction coefficients of the primitive Gaussian functions,
   !> might contain normalization)
   real(wp), protected :: coefficients_env(max_prim, max_shell, max_elem) = 0.0_wp

contains


   !> Create a new basis set
   subroutine new_basis(self, mol, nsh_id, cgto, acc)
      !> Instance of the basis set data
      type(qvszp_basis_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Number of shells per species
      integer, intent(in) :: nsh_id(:)
      !> Contracted Gaussian basis functions for each shell and species
      class(qvszp_cgto_type), intent(in) :: cgto(:, :)
      !> Calculation accuracy
      real(wp), intent(in) :: acc
   
      integer :: iat, isp, ish, iao, ii
      real(wp) :: min_alpha
   
      self%nsh_id = nsh_id
      self%cgto = cgto
      self%intcut = integral_cutoff(acc)
   
      ! Make count of shells for each atom
      self%nsh_at = nsh_id(mol%id)
   
      ! Create mapping between atoms and shells
      self%nsh = sum(self%nsh_at)
      allocate(self%ish_at(mol%nat), self%sh2at(self%nsh))
      ii = 0
      do iat = 1, mol%nat
         self%ish_at(iat) = ii
         do ish = 1, self%nsh_at(iat)
            self%sh2at(ii+ish) = iat
         end do
         ii = ii + self%nsh_at(iat)
      end do
   
      ! Make count of spherical orbitals for each shell
      allocate(self%nao_sh(self%nsh))
      do iat = 1, mol%nat
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            self%nao_sh(ii+ish) = 2*cgto(ish, iat)%ang + 1
         end do
      end do
   
      ! Create mapping between shells and spherical orbitals, also map directly back to atoms
      self%nao = sum(self%nao_sh)
      allocate(self%iao_sh(self%nsh), self%ao2sh(self%nao), self%ao2at(self%nao))
      ii = 0
      do ish = 1, self%nsh
         self%iao_sh(ish) = ii
         do iao = 1, self%nao_sh(ish)
            self%ao2sh(ii+iao) = ish
            self%ao2at(ii+iao) = self%sh2at(ish)
         end do
         ii = ii + self%nao_sh(ish)
      end do
   
      ii = 0
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nsh_id(isp)
            self%iao_sh(ish+self%ish_at(iat)) = ii
            ii = ii + 2*cgto(ish, iat)%ang + 1
         end do
      end do
   
      min_alpha = huge(acc)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nsh_id(isp)
            self%maxl = max(self%maxl, cgto(ish, iat)%ang)
            min_alpha = min(min_alpha, minval(cgto(ish, iat)%alpha(:cgto(ish, iat)%nprim)))
         end do
      end do
   
      self%min_alpha = min_alpha
   
   end subroutine new_basis
   
   subroutine scale_basis(self, mol, qat, cn, expscal)
      !> Instance of the basis set data
      class(qvszp_basis_type), intent(inout) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Atomic charges for the charge scaling of the basis set 
      real(wp), intent(in) :: qat(:)
      !> Coordination number
      real(wp), intent(in) :: cn(:)
      !> Exponent scaling factor
      real(wp), intent(in), optional :: expscal

      integer :: iat, ish
   
      do iat = 1, mol%nat
         do ish = 1, self%nsh_at(iat)
            call self%cgto(ish, iat)%scale_cgto(qat(iat), cn(iat), expscal)
         end do
      end do
   
   end subroutine scale_basis

   subroutine add_qvszp_cgtos(self, mol, nsh_id)
      !> Instance of the cgto data
      type(qvszp_cgto_type), allocatable, intent(inout) :: self(:,:)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Number of shells at each atom
      integer, allocatable, intent(inout) :: nsh_id(:)
   
      integer :: iat, isp, izp, ish, il, nprim
   
      !> Initialize full parameter set
      !> set up the array of CGTOs for the molecule of interest
      call setCGTOexponents()
      call setCGTOcoefficients()
   
      nsh_id = nshell(mol%num)
      allocate(self(maxval(nsh_id), mol%nat))
      do iat = 1, mol%nat
         isp = mol%id(iat)
         izp = mol%num(isp)
         do ish = 1, nsh_id(isp)
            il = ang_shell(ish, izp)
            nprim = n_prim(ish, izp)
            
            ! Gather all necessary parameters
            self(ish, iat)%ang = il
            self(ish, iat)%nprim = nprim
            self(ish, iat)%alpha(1:nprim) = exponents(1:nprim, ish, izp)
            self(ish, iat)%coeff(1:nprim) = coefficients(1:nprim, ish, izp)
            self(ish, iat)%coeff0(1:nprim) = coefficients(1:nprim, ish, izp)
            self(ish, iat)%coeff1(1:nprim) = coefficients_env(1:nprim, ish, izp)
            self(ish, iat)%k1 = p_k1(izp)
            self(ish, iat)%k2 = p_k2(izp)
            self(ish, iat)%k3 = p_k3(izp)
   
         end do
      end do
   
   end subroutine add_qvszp_cgtos
   
   !> Procedure to scale the contraction coefficients of the CGTOs
   subroutine scale_cgto(self, qat, cn, expscal)
      !> Instance of the cgto data
      class(qvszp_cgto_type), intent(inout) :: self
      !> Atomic charge for the charge scaling of the basis set 
      real(wp), intent(in) :: qat
      !> Coordination number
      real(wp), intent(in) :: cn
      !> Exponent scaling factor
      real(wp), intent(in), optional :: expscal
   
      integer :: ipr
      real(wp) :: qeff
      
      ! Scale the coeffients by the charge-dependence
      do ipr = 1, maxg
         qeff = qat + self%k1 * qat**2 + self%k2 * sqrt(cn) + self%k3 * cn * qat
         self%coeff(ipr) = self%coeff0(ipr) + self%coeff1(ipr) * qeff
      end do    
      
      ! Optionally scale the exponent 
      if(present(expscal)) then
         do ipr = 1, maxg
            self%alpha(ipr) = self%alpha(ipr) * expscal
         end do        
      end if
   
   end subroutine scale_cgto
   
   !> Procedure to scale the contraction coefficients of the CGTOs
   subroutine scale_cgto_grad(self, q, cn)
      !> Instance of the basis set data
      class(qvszp_cgto_type), intent(inout) :: self
      !> Atomic charge for the charge scaling of the basis set 
      real(wp), intent(in) :: q
      !> Coordination number
      real(wp), intent(in) :: cn
      
      integer :: ipr
      real(wp) :: qeff
   
      do ipr = 1, maxg
         qeff = q + self%k1 * q**2 + self%k2 * sqrt(cn) + self%k3 * cn * q
         self%coeff(ipr) = self%coeff0(ipr) + self%coeff1(ipr) * qeff
      end do    
   
   end subroutine scale_cgto_grad

   subroutine setCGTOexponents()
      ! set up the array of CGTO exponents during initialization of the basis set data
      
      exponents(:, :,  1) = reshape([&
      &  296.0629766090_wp,   47.9780935580_wp,   11.9388524106_wp,    3.4589164031_wp, &
      &    1.1269931341_wp,    0.4677535907_wp,    0.2098893356_wp,    0.0741956361_wp, & ! shell 1
      &    1.9212296822_wp,    0.6145584683_wp,    0.2273061668_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  2) = reshape([&
      &  325.3193889837_wp,   53.8529443290_wp,   13.0508309735_wp,    3.8713602188_wp, &
      &    1.2108391704_wp,    0.4533424114_wp,    0.2475298319_wp,    0.2002843134_wp, & ! shell 1
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  3) = reshape([&
      &    3.7043076907_wp,    0.5637486123_wp,    0.1455350857_wp,    0.0686608449_wp, &
      &    0.0295010199_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.7437867540_wp,    0.4730880728_wp,    0.1570436642_wp,    0.0599773554_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  4) = reshape([&
      &    7.4077694199_wp,    1.2355438040_wp,    0.4433048670_wp,    0.1788375974_wp, &
      &    0.0626826282_wp,    0.0194239695_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.7655075121_wp,    0.6565159830_wp,    0.2142912316_wp,    0.0605181540_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  5) = reshape([&
      &    1.5439131356_wp,    1.2089824249_wp,    0.3058442147_wp,    0.1160249703_wp, &
      &    0.0682235936_wp,    0.0502268935_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    8.9769010628_wp,    2.3657980483_wp,    0.6943666432_wp,    0.2328484561_wp, &
      &    0.0761862667_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  6) = reshape([&
      &    2.1969435617_wp,    1.9249701500_wp,    0.5018527654_wp,    0.1966308943_wp, &
      &    0.1377690883_wp,    0.0762566102_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   22.5294566414_wp,    6.1533099644_wp,    2.1678631339_wp,    0.7723850897_wp, &
      &    0.2806512667_wp,    0.0954367884_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  7) = reshape([&
      &    2.6986954458_wp,    2.2504031153_wp,    0.6171177991_wp,    0.3116793580_wp, &
      &    0.1638031516_wp,    0.0675314174_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   27.2287637240_wp,    7.0606760739_wp,    2.3483638135_wp,    0.8485780333_wp, &
      &    0.2986799649_wp,    0.0909296536_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  8) = reshape([&
      &    3.4962893786_wp,    2.8288974669_wp,    0.7249120926_wp,    0.3361163612_wp, &
      &    0.2387426762_wp,    0.0985625339_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   37.3050123505_wp,    9.6142969083_wp,    3.0827510750_wp,    1.0931055447_wp, &
      &    0.3840580235_wp,    0.1212507025_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :,  9) = reshape([&
      &    6.8229606760_wp,    5.1904412429_wp,    1.3692346791_wp,    0.7016091148_wp, &
      &    0.4259574125_wp,    0.1605751395_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   53.6967237849_wp,   14.8893967124_wp,    4.7495843078_wp,    1.6898244647_wp, &
      &    0.5688111793_wp,    0.1780790714_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 10) = reshape([&
      &    8.8362115995_wp,    7.4896279798_wp,    1.7523581081_wp,    0.6931614114_wp, &
      &    0.3108356666_wp,    0.1748493239_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   71.4752440352_wp,   18.5157491263_wp,    6.0351158600_wp,    2.1607109774_wp, &
      &    0.7339973868_wp,    0.2355636322_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 11) = reshape([&
      &    2.2049566801_wp,    0.5378058890_wp,    0.0715844411_wp,    0.0275628967_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4248136101_wp,    0.2581316954_wp,    0.0800762228_wp,    0.0236633030_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 12) = reshape([&
      &    2.6334836883_wp,    0.7908890874_wp,    0.1303989301_wp,    0.0517718329_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8515995867_wp,    0.1724280102_wp,    0.0576679218_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4186096077_wp,    0.1009160378_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 13) = reshape([&
      &    3.2922192739_wp,    1.0408800818_wp,    0.3712308404_wp,    0.1676012256_wp, &
      &    0.0649014274_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8835859483_wp,    0.4073805077_wp,    0.1424168358_wp,    0.0541120697_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.5058995677_wp,    0.1330757688_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 14) = reshape([&
      &    4.2484659124_wp,    1.3254502628_wp,    0.3071780296_wp,    0.1292848447_wp, &
      &    0.0965719453_wp,    0.0524151279_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2065748118_wp,    0.4996563298_wp,    0.1889116493_wp,    0.0701085381_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.0589057689_wp,    0.2871608270_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 15) = reshape([&
      &    7.9611244782_wp,    1.4589712233_wp,    0.4345871479_wp,    0.1865322669_wp, &
      &    0.0663466643_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.1713337396_wp,    0.7711286095_wp,    0.3376093165_wp,    0.1602097044_wp, &
      &    0.0671957686_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.6580618305_wp,    0.3649059874_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 16) = reshape([&
      &    7.4651982988_wp,    1.9329553350_wp,    0.5205100977_wp,    0.2302754276_wp, &
      &    0.0860210463_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.6375502697_wp,    1.1362521699_wp,    0.4223961389_wp,    0.1860968380_wp, &
      &    0.0794613763_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.6181982268_wp,    0.3735700452_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 17) = reshape([&
      &   16.4379690803_wp,    2.1310561497_wp,    0.7018504516_wp,    0.2950510846_wp, &
      &    0.1113181312_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.6852262094_wp,    1.6334658243_wp,    0.6471540154_wp,    0.2632798726_wp, &
      &    0.1040962011_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    2.0965021978_wp,    0.4851267355_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 18) = reshape([&
      &   17.2077820733_wp,    2.2811018859_wp,    1.7479478805_wp,    0.4654674147_wp, &
      &    0.1682728781_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    3.2323664724_wp,    1.1261631764_wp,    0.4398790097_wp,    0.2334854311_wp, &
      &    0.1266670515_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.7350463941_wp,    0.3634120679_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 19) = reshape([&
      &    0.9105152366_wp,    0.3037322089_wp,    0.0443108826_wp,    0.0183243318_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2069767309_wp,    0.0741990636_wp,    0.0245000321_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 20) = reshape([&
      &    0.8594832837_wp,    0.4727901462_wp,    0.0508507361_wp,    0.0398733003_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3586477973_wp,    0.0770419912_wp,    0.0242772945_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    1.5264809172_wp,    0.4506518991_wp,    0.1196551892_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 21) = reshape([&
      &    0.8306707451_wp,    0.4654634186_wp,    0.1571842701_wp,    0.1095324675_wp, &
      &    0.0403752058_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4520553793_wp,    0.1576499805_wp,    0.1159015610_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   35.3993939379_wp,    9.2237130122_wp,    2.9802097002_wp,    1.0603213126_wp, &
      &    0.3620479701_wp,    0.1250409202_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 22) = reshape([&
      &    1.0787305631_wp,    0.5227318658_wp,    0.2659156930_wp,    0.1049459311_wp, &
      &    0.0422872604_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5796056355_wp,    0.1165172934_wp,    0.0804870873_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   40.5510394549_wp,   10.9301936751_wp,    3.6094276078_wp,    1.3098640939_wp, &
      &    0.4493494621_wp,    0.1458767654_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 23) = reshape([&
      &    0.9722923426_wp,    0.5540477155_wp,    0.2614418621_wp,    0.1161041521_wp, &
      &    0.0386015787_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.7502260397_wp,    0.1332704136_wp,    0.1121052883_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   39.0451550260_wp,   11.0367333189_wp,    3.7466021315_wp,    1.3652031999_wp, &
      &    0.4659298449_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 24) = reshape([&
      &    0.8613981186_wp,    0.6083737405_wp,    0.2430342508_wp,    0.1025460689_wp, &
      &    0.0277793094_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8097904722_wp,    0.1174447647_wp,    0.0294785765_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   43.5624715979_wp,   12.0038600987_wp,    4.0518662943_wp,    1.4659017591_wp, &
      &    0.4965798804_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 25) = reshape([&
      &    1.1674696963_wp,    0.7522820330_wp,    0.1923117018_wp,    0.1480221484_wp, &
      &    0.0588205302_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.8293026407_wp,    0.1526550794_wp,    0.1320923271_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   43.2774686387_wp,   12.2091517240_wp,    4.2127545096_wp,    1.5234390169_wp, &
      &    0.5068887384_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 26) = reshape([&
      &    1.4345008149_wp,    0.8384222876_wp,    0.2536490631_wp,    0.1100629629_wp, &
      &    0.0262617693_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.9127217962_wp,    0.1534939292_wp,    0.0429652564_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   47.5163871319_wp,   13.4829424357_wp,    4.6789918295_wp,    1.6857476940_wp, &
      &    0.5518382943_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 27) = reshape([&
      &    1.3655362096_wp,    0.8333687882_wp,    0.3101103317_wp,    0.1304100197_wp, &
      &    0.0369947013_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.9769874338_wp,    0.1768560243_wp,    0.0894661613_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   49.2955004674_wp,   14.2172829758_wp,    4.9571417489_wp,    1.7755767421_wp, &
      &    0.5763062135_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 28) = reshape([&
      &    1.4470300713_wp,    0.9412032436_wp,    0.2887395603_wp,    0.1147674161_wp, &
      &    0.0282883490_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.1286187088_wp,    0.1709788701_wp,    0.0527416617_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   56.6626460481_wp,   15.8963422320_wp,    5.4516937144_wp,    1.9237061681_wp, &
      &    0.6107858432_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 29) = reshape([&
      &    1.5113423344_wp,    0.9558846550_wp,    0.2445102708_wp,    0.1579736623_wp, &
      &    0.0555112704_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.1122194202_wp,    0.1890200227_wp,    0.0970534280_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   56.3912849171_wp,   16.0819697857_wp,    5.5684785009_wp,    1.9524273096_wp, &
      &    0.6101381651_wp,    0.1521373750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 30) = reshape([&
      &    1.8083224514_wp,    1.1341539195_wp,    0.3691137315_wp,    0.1495968016_wp, &
      &    0.0519122231_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2178742063_wp,    0.1832738968_wp,    0.0479838131_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 31) = reshape([&
      &    3.2781463741_wp,    1.7353677348_wp,    0.2217785045_wp,    0.1795022133_wp, &
      &    0.0671222115_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.0599015561_wp,    0.7776843644_wp,    0.1909339923_wp,    0.0654922435_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3445113790_wp,    0.1059643396_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 32) = reshape([&
      &    3.1784020104_wp,    1.9172406836_wp,    0.2879063648_wp,    0.1916456611_wp, &
      &    0.0716382866_wp,    0.0594866492_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2110886910_wp,    0.8143183661_wp,    0.2411958564_wp,    0.0842675992_wp, &
      &    0.0793262087_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3742430139_wp,    0.1199288234_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 33) = reshape([&
      &    3.3959834521_wp,    1.8774718178_wp,    0.3773073232_wp,    0.2216454940_wp, &
      &    0.1374100487_wp,    0.0728534130_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.3412056183_wp,    0.9368162421_wp,    0.3278044176_wp,    0.1685573872_wp, &
      &    0.0644280149_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2871650230_wp,    0.0838070691_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 34) = reshape([&
      &    3.6933360117_wp,    2.1914868281_wp,    0.3575158144_wp,    0.1912816782_wp, &
      &    0.1287895455_wp,    0.0677083089_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.4583881774_wp,    1.0124457661_wp,    0.3678204131_wp,    0.1859870696_wp, &
      &    0.0771928990_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3235036177_wp,    0.1119420008_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 35) = reshape([&
      &    4.5291369073_wp,    2.2964861296_wp,    0.4258542828_wp,    0.2450543837_wp, &
      &    0.1359057290_wp,    0.0764368580_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.7721665463_wp,    1.3086216630_wp,    0.4642724437_wp,    0.2046820847_wp, &
      &    0.0876297039_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4597465232_wp,    0.1838304254_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 36) = reshape([&
      &    4.1551034747_wp,    2.8214253386_wp,    0.4617262010_wp,    0.2232537521_wp, &
      &    0.1621058406_wp,    0.0698534641_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0891432410_wp,    1.4555990158_wp,    0.5559447836_wp,    0.2514666827_wp, &
      &    0.1041667197_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4897230858_wp,    0.1737069992_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 37) = reshape([&
      &    0.5749546742_wp,    0.3420127048_wp,    0.0259494074_wp,    0.0185252242_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1568178078_wp,    0.0849536651_wp,    0.0221134395_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 38) = reshape([&
      &    0.5954071331_wp,    0.4152975457_wp,    0.0450108930_wp,    0.0331722748_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2046367275_wp,    0.0864704826_wp,    0.0336153577_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.6180913968_wp,    0.2004636179_wp,    0.0707531984_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 39) = reshape([&
      &    0.6764765958_wp,    0.3800665845_wp,    0.1035044711_wp,    0.0516749760_wp, &
      &    0.0322119947_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2890110460_wp,    0.1008735199_wp,    0.0904433617_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.2384546244_wp,    2.2832789230_wp,    0.9010403870_wp,    0.3650599630_wp, &
      &    0.1579709951_wp,    0.0713026669_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 40) = reshape([&
      &    0.5392761395_wp,    0.3444598486_wp,    0.2827234759_wp,    0.1371494096_wp, &
      &    0.0478289077_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3489791992_wp,    0.1142484465_wp,    0.0537780877_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   10.0810460009_wp,    2.2487707552_wp,    0.9188527433_wp,    0.3640823834_wp, &
      &    0.1477812228_wp,    0.0643919646_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 41) = reshape([&
      &    0.7785562920_wp,    0.4374177807_wp,    0.1773701335_wp,    0.1052866004_wp, &
      &    0.0490395758_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3948551018_wp,    0.1043303055_wp,    0.0895188884_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   10.1908225383_wp,    1.8681908504_wp,    0.8854314957_wp,    0.3812857837_wp, &
      &    0.1528881480_wp,    0.0755313966_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 42) = reshape([&
      &    0.7129465040_wp,    0.4619138528_wp,    0.1223419482_wp,    0.0534987192_wp, &
      &    0.0329945617_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4374831659_wp,    0.1919532059_wp,    0.1138171731_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   11.0794690319_wp,    3.5990576140_wp,    1.3542860564_wp,    0.5770253070_wp, &
      &    0.2493898096_wp,    0.1096952739_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 43) = reshape([&
      &    0.9632757947_wp,    0.4747292576_wp,    0.1554463392_wp,    0.1130245747_wp, &
      &    0.0540159936_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4559170006_wp,    0.1381856089_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   14.5007348902_wp,    4.3376558008_wp,    1.7609699217_wp,    0.8787802854_wp, &
      &    0.3878093526_wp,    0.1506201387_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 44) = reshape([&
      &    1.0732283615_wp,    0.5846425374_wp,    0.1234713112_wp,    0.0781899278_wp, &
      &    0.0382670200_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5253287557_wp,    0.1433580153_wp,    0.0486047832_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   14.1908594409_wp,    4.0523812745_wp,    1.7820655178_wp,    0.8420248803_wp, &
      &    0.3657897622_wp,    0.1388794428_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 45) = reshape([&
      &    1.2534322956_wp,    0.6353695688_wp,    0.1520223523_wp,    0.0686786828_wp, &
      &    0.0507321743_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4980038655_wp,    0.1563880086_wp,    0.0700952839_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   16.8067481665_wp,    5.7277788794_wp,    1.8830331580_wp,    0.8498356288_wp, &
      &    0.3662734564_wp,    0.1454309021_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 46) = reshape([&
      &    1.3719793790_wp,    0.6441996411_wp,    0.1073426140_wp,    0.0490090004_wp, &
      &    0.0185761967_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5880612185_wp,    0.1463459073_wp,    0.0670539791_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   18.2560629442_wp,    3.7819041632_wp,    1.9082188860_wp,    0.8834523111_wp, &
      &    0.3743221295_wp,    0.1358460286_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 47) = reshape([&
      &    1.6883392445_wp,    0.6556694926_wp,    0.1366027982_wp,    0.0666102898_wp, &
      &    0.0440540239_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5870492284_wp,    0.1544383746_wp,    0.0736802768_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   17.6539941832_wp,    4.5274784839_wp,    2.3322407702_wp,    1.0889648174_wp, &
      &    0.4485553389_wp,    0.1686430070_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 48) = reshape([&
      &    1.3503008671_wp,    0.8614574884_wp,    0.4043122590_wp,    0.1280005499_wp, &
      &    0.0473981278_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.7660012761_wp,    0.1547342249_wp,    0.0387159717_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1058020500_wp,    0.0351860655_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 49) = reshape([&
      &    1.3228663135_wp,    0.9442427410_wp,    0.3271063888_wp,    0.1394146133_wp, &
      &    0.0566012650_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.4305338328_wp,    1.2150155205_wp,    0.1866977828_wp,    0.0585700111_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1186291731_wp,    0.0757732953_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 50) = reshape([&
      &    2.2509720035_wp,    1.1830513117_wp,    0.4603611431_wp,    0.1833289016_wp, &
      &    0.1277913368_wp,    0.0507002081_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.2944208688_wp,    1.5983237423_wp,    0.2217739778_wp,    0.0773321582_wp, &
      &    0.0655381886_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1963477682_wp,    0.0898761566_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 51) = reshape([&
      &    1.9546443939_wp,    1.3972682527_wp,    0.3396924855_wp,    0.2227489251_wp, &
      &    0.1456212820_wp,    0.0526070607_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.0926403119_wp,    1.3927912485_wp,    0.2714431291_wp,    0.1275834042_wp, &
      &    0.0533836071_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1879166901_wp,    0.1551913077_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 52) = reshape([&
      &    2.2674705626_wp,    1.5934378114_wp,    0.3356684935_wp,    0.2210544927_wp, &
      &    0.1794207556_wp,    0.0802802007_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    2.3190783635_wp,    1.5500085119_wp,    0.3581443938_wp,    0.1651351817_wp, &
      &    0.0654337713_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1882412330_wp,    0.0412880147_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 53) = reshape([&
      &    2.5927987167_wp,    1.6200376522_wp,    0.3419234720_wp,    0.2747073514_wp, &
      &    0.2059248675_wp,    0.1012654473_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.9468683265_wp,    0.6564654372_wp,    0.2189842016_wp,    0.1853603803_wp, &
      &    0.0835096198_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2458963824_wp,    0.1229705157_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 54) = reshape([&
      &    2.6742885822_wp,    1.7567645576_wp,    0.3883299888_wp,    0.3315361303_wp, &
      &    0.1787231702_wp,    0.0718749422_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.0399175475_wp,    0.7659928838_wp,    0.2722487515_wp,    0.1862661077_wp, &
      &    0.0761656134_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2044347483_wp,    0.0671786274_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 55) = reshape([&
      &    0.3458918276_wp,    0.2906279979_wp,    0.0229171732_wp,    0.0172463943_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1079351513_wp,    0.0889086049_wp,    0.0188609591_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 56) = reshape([&
      &    0.3851681594_wp,    0.2975636012_wp,    0.0370011696_wp,    0.0346330487_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1329954135_wp,    0.1153725672_wp,    0.0304309682_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2822785789_wp,    0.0882479945_wp,    0.0261349312_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 57) = reshape([&
      &    0.3652067118_wp,    0.2135532475_wp,    0.1199987183_wp,    0.0354386658_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1668554422_wp,    0.1091657654_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    3.1193600378_wp,    0.5777229524_wp,    0.3291375182_wp,    0.1397870406_wp, &
      &    0.0527952624_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 72) = reshape([&
      &    0.6837191225_wp,    0.4676908224_wp,    0.1032754319_wp,    0.0343387390_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3784784425_wp,    0.1244468385_wp,    0.0505132821_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.6149683328_wp,    0.9812300925_wp,    0.3908750564_wp,    0.1582219844_wp, &
      &    0.0622857724_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 73) = reshape([&
      &    0.6960787665_wp,    0.4496551889_wp,    0.2186334602_wp,    0.0829100122_wp, &
      &    0.0390265498_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3854210047_wp,    0.1478781816_wp,    0.0795033400_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.0396335955_wp,    1.0790867168_wp,    0.4463540459_wp,    0.1836856248_wp, &
      &    0.0729214443_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 74) = reshape([&
      &    1.1230157557_wp,    0.5896280540_wp,    0.1362600773_wp,    0.0616371055_wp, &
      &    0.0377400581_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4409845972_wp,    0.1550869705_wp,    0.0671865105_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    7.4489011662_wp,    1.1050607468_wp,    0.4526042530_wp,    0.1883153444_wp, &
      &    0.0858781381_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 75) = reshape([&
      &    1.1256618658_wp,    0.6090517530_wp,    0.1414823711_wp,    0.0544499783_wp, &
      &    0.0478103973_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4345105058_wp,    0.2042208712_wp,    0.1193935478_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    7.0592225356_wp,    1.5514075726_wp,    1.1137551097_wp,    0.4458446018_wp, &
      &    0.1598958872_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 76) = reshape([&
      &    1.0331612402_wp,    0.6588862216_wp,    0.1361895010_wp,    0.0556929979_wp, &
      &    0.0307569661_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5044244667_wp,    0.1601752493_wp,    0.0643382805_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    6.8371580008_wp,    1.3038841730_wp,    0.5753299068_wp,    0.2658882909_wp, &
      &    0.1123772298_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 77) = reshape([&
      &    1.0794609092_wp,    0.6492423459_wp,    0.2199426177_wp,    0.0834922927_wp, &
      &    0.0388519894_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4487050433_wp,    0.2418587057_wp,    0.1081757627_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.4024255233_wp,    1.3298289575_wp,    0.5777613431_wp,    0.2734387652_wp, &
      &    0.1210613626_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 78) = reshape([&
      &    1.3698023838_wp,    0.7257453910_wp,    0.1319615787_wp,    0.0777821895_wp, &
      &    0.0392611160_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5479325568_wp,    0.1481247655_wp,    0.0550461136_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.4760176647_wp,    1.4257610916_wp,    0.6203299774_wp,    0.2696753513_wp, &
      &    0.1060426336_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 79) = reshape([&
      &    1.3101235185_wp,    0.7424205017_wp,    0.1492400978_wp,    0.0569819664_wp, &
      &    0.0462275609_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.5551190544_wp,    0.1821547642_wp,    0.0679459671_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    8.8724027378_wp,    1.4968892843_wp,    0.6628731563_wp,    0.3044308364_wp, &
      &    0.1324251318_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 80) = reshape([&
      &    1.4414111867_wp,    0.8859054655_wp,    0.3072833627_wp,    0.1373321497_wp, &
      &    0.0587810800_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.6016138345_wp,    0.1455530159_wp,    0.0277509758_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0777461389_wp,    0.0385457088_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 81) = reshape([&
      &    1.4274324972_wp,    0.9671436880_wp,    0.2598401553_wp,    0.1679694816_wp, &
      &    0.0606732368_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.1705842936_wp,    1.0687150265_wp,    0.1688729812_wp,    0.0518197726_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.6097949719_wp,    0.0858730236_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 82) = reshape([&
      &    1.4072150736_wp,    1.0343124310_wp,    0.5114274272_wp,    0.1876895475_wp, &
      &    0.0639014487_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.2707421905_wp,    1.0192238533_wp,    0.2215246292_wp,    0.0740688503_wp, &
      &    0.0643999270_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2034715847_wp,    0.0840830597_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 83) = reshape([&
      &    1.6286675405_wp,    1.0125231394_wp,    0.3914877717_wp,    0.1926675612_wp, &
      &    0.1651017294_wp,    0.0610936554_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.3835100086_wp,    1.1915513148_wp,    0.2761659482_wp,    0.1405150570_wp, &
      &    0.0620946082_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1092332919_wp,    0.0291157869_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 84) = reshape([&
      &    1.7976784830_wp,    1.3509248224_wp,    0.6956891781_wp,    0.2543747879_wp, &
      &    0.1466239493_wp,    0.0722609656_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.6596013907_wp,    1.2455175542_wp,    0.2518870683_wp,    0.1466558932_wp, &
      &    0.0772995342_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1564139292_wp,    0.0424878804_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 85) = reshape([&
      &    1.9638838709_wp,    1.3060666881_wp,    0.5462006822_wp,    0.3792792878_wp, &
      &    0.2338279502_wp,    0.1052257635_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.9202808257_wp,    1.2321759927_wp,    0.3262403959_wp,    0.1415497933_wp, &
      &    0.0644222185_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1908503703_wp,    0.1643201294_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      exponents(:, :, 86) = reshape([&
      &    1.9628316800_wp,    1.4642104527_wp,    0.3651899239_wp,    0.2156657792_wp, &
      &    0.1558792626_wp,    0.0951988921_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    1.7555810689_wp,    1.4562661914_wp,    0.3760332653_wp,    0.1609065878_wp, &
      &    0.0579682124_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1604814020_wp,    0.0243495180_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
   end subroutine setCGTOexponents

   subroutine setCGTOcoefficients()
      ! set up the array of CGTO coefficients during initialization of the basis set data
      coefficients(:, :,  1) = reshape([&
      &    0.0010949386_wp,    0.0072326738_wp,    0.0339004967_wp,    0.1391347923_wp, &
      &    0.3830065360_wp,    0.5467684408_wp,    0.5700479997_wp,    0.4569839155_wp, & ! shell 1
      &    0.2255341694_wp,    0.5291370409_wp,    0.8180148717_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  2) = reshape([&
      &    0.0028546294_wp,    0.0194448716_wp,    0.0906466552_wp,    0.2906836491_wp, &
      &    0.6144265753_wp,    0.6137808702_wp,    0.0485348306_wp,    0.3876820157_wp, & ! shell 1
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  3) = reshape([&
      &    0.0090146316_wp,   -0.1743650828_wp,    0.0115312555_wp,    0.4778794895_wp, &
      &    0.8608215715_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0391917981_wp,    0.0813299551_wp,    0.4010973462_wp,    0.9115757567_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  4) = reshape([&
      &    0.0107093433_wp,   -0.2456008301_wp,    0.0672621012_wp,    0.6377135888_wp, &
      &    0.7267324053_wp,    0.0149245351_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0538129598_wp,    0.2275585136_wp,    0.5279956606_wp,    0.8164201557_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  5) = reshape([&
      &   -0.5582059698_wp,    0.4311748325_wp,    0.5531335167_wp,    0.4158594078_wp, &
      &    0.0984253549_wp,    0.1179452170_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0217527399_wp,    0.1232949416_wp,    0.4835439401_wp,    0.6086038533_wp, &
      &    0.6165320619_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  6) = reshape([&
      &   -0.6444768885_wp,    0.5541287618_wp,    0.3662930406_wp,    0.3531528795_wp, &
      &    0.0358811984_wp,    0.1319691146_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0090734598_wp,    0.0694438794_wp,    0.2218132871_wp,    0.5601762237_wp, &
      &    0.5540968223_wp,    0.5701520809_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  7) = reshape([&
      &   -0.6302370972_wp,    0.5950833423_wp,    0.3876892769_wp,    0.2689377543_wp, &
      &    0.1419833316_wp,    0.0767284779_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0168702357_wp,    0.1162342364_wp,    0.3440184346_wp,    0.6562285180_wp, &
      &    0.6128335426_wp,    0.2483052462_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  8) = reshape([&
      &   -0.5986461430_wp,    0.5928903063_wp,    0.4720830846_wp,    0.2100701609_wp, &
      &    0.1274271965_wp,    0.0829112688_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0179453646_wp,    0.1189908463_wp,    0.3583415329_wp,    0.6242383100_wp, &
      &    0.6227928229_wp,    0.2820746713_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :,  9) = reshape([&
      &   -0.4801928273_wp,    0.2998085670_wp,    0.5261051184_wp,    0.3635220361_wp, &
      &    0.4963475815_wp,    0.1556723602_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0174699483_wp,    0.1011033895_wp,    0.3070724710_wp,    0.5829007661_wp, &
      &    0.6756117661_wp,    0.3145708848_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 10) = reshape([&
      &   -0.4357984747_wp,    0.2816992209_wp,    0.5153844768_wp,    0.6568024876_wp, &
      &    0.1738731947_wp,    0.0590142686_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0146484990_wp,    0.1034556349_wp,    0.3101927587_wp,    0.5830039270_wp, &
      &    0.6679197262_wp,    0.3268829544_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 11) = reshape([&
      &    0.0169546496_wp,   -0.2498540276_wp,    0.3893771186_wp,    0.8863808235_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1031491023_wp,    0.0862151969_wp,    0.3464982866_wp,    0.9283674595_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 12) = reshape([&
      &    0.0265680049_wp,   -0.2787776385_wp,    0.5186899267_wp,    0.8077981984_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0433051828_wp,    0.2572312999_wp,    0.9653790548_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0364466611_wp,    0.9993355997_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 13) = reshape([&
      &    0.0392356539_wp,   -0.3536663664_wp,    0.0936987159_wp,    0.6851483729_wp, &
      &    0.6286278092_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1365566327_wp,    0.2815261723_wp,    0.5904177966_wp,    0.7439772348_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2699053117_wp,    0.9628868691_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 14) = reshape([&
      &    0.0396899986_wp,   -0.3866611485_wp,    0.5460550687_wp,    0.6946830491_wp, &
      &    0.1627280911_wp,    0.2041488476_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1023118647_wp,    0.3129231213_wp,    0.6729852350_wp,    0.6623460394_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1533169451_wp,    0.9881770663_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 15) = reshape([&
      &    0.0290071695_wp,   -0.3842474695_wp,    0.5010544647_wp,    0.7394217341_wp, &
      &    0.2317593338_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1809957020_wp,    0.3601760449_wp,    0.6264075297_wp,    0.5746008026_wp, &
      &    0.3390594297_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2214433162_wp,    0.9751732450_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 16) = reshape([&
      &    0.0413415802_wp,   -0.4036431828_wp,    0.5216147592_wp,    0.7274812960_wp, &
      &    0.1845320071_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2511105598_wp,    0.3844791597_wp,    0.6444155690_wp,    0.5645473314_wp, &
      &    0.2348066174_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2661713796_wp,    0.9639257215_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 17) = reshape([&
      &    0.0191446759_wp,   -0.3941868394_wp,    0.4836750278_wp,    0.7486157013_wp, &
      &    0.2233455088_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1427346169_wp,    0.1733874265_wp,    0.6045076752_wp,    0.6852131917_wp, &
      &    0.3385512990_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2579774104_wp,    0.9661509487_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 18) = reshape([&
      &    0.0089657613_wp,   -0.5991609353_wp,    0.4455629873_wp,    0.6264908245_wp, &
      &    0.2234024618_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0822890990_wp,    0.4932019909_wp,    0.7649634175_wp,    0.2489536824_wp, &
      &    0.3206763700_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3140884921_wp,    0.9493937113_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 19) = reshape([&
      &    0.0594660805_wp,   -0.3778600312_wp,    0.5231862752_wp,    0.7615521674_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0651633551_wp,    0.2364090998_wp,    0.9694660771_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 20) = reshape([&
      &    0.1144183909_wp,   -0.3348129832_wp,    0.9379249930_wp,   -0.0255729982_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0336217144_wp,    0.5957686636_wp,    0.7631134251_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.2368221088_wp,    0.6209354579_wp,    0.7546571008_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 21) = reshape([&
      &    0.0802221906_wp,   -0.3790459155_wp,    0.0076002856_wp,    0.3584936744_wp, &
      &    0.8493015455_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1283153170_wp,    0.0091689160_wp,    0.9916910358_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0111931883_wp,    0.0872875924_wp,    0.2897256927_wp,    0.5204897295_wp, &
      &    0.5521281189_wp,    0.5766798016_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 22) = reshape([&
      &    0.0493099025_wp,   -0.3520965937_wp,    0.0123295592_wp,    0.4309900976_wp, &
      &    0.8292719940_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0895317589_wp,    0.2994150737_wp,    0.9499129843_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0128503084_wp,    0.0932278826_wp,    0.3101320440_wp,    0.5514879848_wp, &
      &    0.5907057862_wp,    0.4918223496_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 23) = reshape([&
      &    0.0653342813_wp,   -0.3948137155_wp,    0.0219945447_wp,    0.5756992371_wp, &
      &    0.7126992283_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1222573548_wp,    0.4225979516_wp,    0.8980334685_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0199115222_wp,    0.1187580820_wp,    0.3711833996_wp,    0.6072513991_wp, &
      &    0.5924999683_wp,    0.3576485131_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 24) = reshape([&
      &    0.1142793077_wp,   -0.4169289949_wp,    0.0153871508_wp,    0.7088011921_wp, &
      &    0.5572024396_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0805021354_wp,    0.5731977485_wp,    0.8154530932_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0216361681_wp,    0.1326529336_wp,    0.3978895195_wp,    0.6213308443_wp, &
      &    0.6123366245_wp,    0.2502215936_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 25) = reshape([&
      &    0.0990261636_wp,   -0.3774104315_wp,    0.0106093368_wp,    0.3499717733_wp, &
      &    0.8515646687_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0808022283_wp,    0.6155674395_wp,    0.7839309455_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0278068305_wp,    0.1563483089_wp,    0.4408334284_wp,    0.6509247923_wp, &
      &    0.5447197230_wp,    0.2450004350_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 26) = reshape([&
      &    0.0827431662_wp,   -0.3524443664_wp,    0.0093420971_wp,    0.8005601741_wp, &
      &    0.4774438919_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1245191264_wp,    0.6247220042_wp,    0.7708549829_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0281336102_wp,    0.1549882631_wp,    0.4329530098_wp,    0.6323111928_wp, &
      &    0.5785127309_wp,    0.2307474923_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 27) = reshape([&
      &    0.0776852630_wp,   -0.3519416744_wp,    0.0102041608_wp,    0.6911710256_wp, &
      &    0.6263230366_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1208394510_wp,    0.4504793528_wp,    0.8845711841_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0319583704_wp,    0.1675807640_wp,    0.4579670521_wp,    0.6433995101_wp, &
      &    0.5605766346_wp,    0.1815280594_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 28) = reshape([&
      &    0.1036977741_wp,   -0.3610108116_wp,    0.0118034234_wp,    0.8450923553_wp, &
      &    0.3802598529_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1341309403_wp,    0.7358567755_wp,    0.6637195920_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0296661714_wp,    0.1713974168_wp,    0.4604744883_wp,    0.6393890079_wp, &
      &    0.5486391152_wp,    0.2188216332_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 29) = reshape([&
      &    0.0745059292_wp,   -0.3236654341_wp,    0.0087324338_wp,    0.4874282539_wp, &
      &    0.8074818854_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1028890524_wp,    0.5531004383_wp,    0.8267368070_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0355607850_wp,    0.1913098794_wp,    0.4885052132_wp,    0.6445374314_wp, &
      &    0.5186744098_wp,    0.1976030702_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 30) = reshape([&
      &    0.1625382444_wp,   -0.4190411332_wp,    0.0161818715_wp,    0.7147361211_wp, &
      &    0.5356083196_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0550731667_wp,    0.5206125385_wp,    0.8520149829_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 31) = reshape([&
      &    0.1746236041_wp,   -0.4839721399_wp,    0.2214694468_wp,    0.5947171150_wp, &
      &    0.5766631616_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3082970578_wp,    0.3021951550_wp,    0.4443481081_wp,    0.7849750132_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0141417333_wp,    0.9999000007_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 32) = reshape([&
      &    0.2471153357_wp,   -0.5696353526_wp,    0.2373790503_wp,    0.6951664876_wp, &
      &    0.1363957227_wp,    0.2371508459_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3467769135_wp,    0.3142435575_wp,    0.6372033137_wp,    0.5046972824_wp, &
      &    0.3467698789_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3041010877_wp,    0.9526397685_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 33) = reshape([&
      &    0.1815335231_wp,   -0.5466552363_wp,    0.3352718201_wp,    0.6895458373_wp, &
      &    0.0965685936_wp,    0.2664722956_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.4324935168_wp,    0.4009239554_wp,    0.4844422206_wp,    0.5587009066_wp, &
      &    0.3246203501_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.8985815217_wp,    0.4388066190_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 34) = reshape([&
      &    0.2315684536_wp,   -0.5873551144_wp,    0.6780832924_wp,    0.3316244354_wp, &
      &    0.1486020098_wp,    0.0976511436_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3641011797_wp,    0.3722755657_wp,    0.5706760184_wp,    0.5595497526_wp, &
      &    0.3001236252_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9502754357_wp,    0.3114106554_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 35) = reshape([&
      &    0.1472601719_wp,   -0.5329068161_wp,    0.6931258049_wp,    0.3854770770_wp, &
      &    0.2454433556_wp,    0.0711784236_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3894522717_wp,    0.3446405463_wp,    0.5608725576_wp,    0.5737772689_wp, &
      &    0.2928334709_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.7273130959_wp,    0.6863058069_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 36) = reshape([&
      &    0.2920011133_wp,   -0.6126666701_wp,    0.6657736865_wp,    0.2875541894_wp, &
      &    0.1122161680_wp,    0.0289899865_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3400394681_wp,    0.2900250677_wp,    0.5642654314_wp,    0.6307256713_wp, &
      &    0.2899107978_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.5986358375_wp,    0.8010213068_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 37) = reshape([&
      &    0.2127722928_wp,   -0.4708945571_wp,    0.8066570194_wp,    0.2868635923_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1637526039_wp,    0.1522054130_wp,    0.9746889745_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 38) = reshape([&
      &    0.2669464427_wp,   -0.5280977879_wp,    0.7639112316_wp,    0.2509225965_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1296966669_wp,    0.5190667011_wp,    0.8549244134_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1426417156_wp,    0.6106047130_wp,    0.7991293970_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 39) = reshape([&
      &    0.1390853191_wp,   -0.4827265709_wp,    0.3808371062_wp,    0.0339543325_wp, &
      &    0.7755259722_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3008754571_wp,    0.7112811800_wp,    0.6352582485_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0119144970_wp,    0.0375378247_wp,    0.3146570996_wp,    0.5401045695_wp, &
      &    0.4507228641_wp,    0.6360627485_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 40) = reshape([&
      &    0.1082431966_wp,   -0.5047129057_wp,    0.0445274839_wp,    0.2623330475_wp, &
      &    0.8140927273_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1511037208_wp,    0.3273452886_wp,    0.9327447280_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0117581692_wp,    0.0706419371_wp,    0.4502680733_wp,    0.6883882746_wp, &
      &    0.4459400058_wp,    0.3455274504_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 41) = reshape([&
      &    0.0885707747_wp,   -0.4468474517_wp,    0.0583006427_wp,    0.3556505216_wp, &
      &    0.8140001930_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2057254347_wp,    0.8122997237_wp,    0.5457528785_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0153295379_wp,    0.1750810591_wp,    0.5174086443_wp,    0.7129592834_wp, &
      &    0.3903769847_wp,    0.2017295042_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 42) = reshape([&
      &    0.1556652881_wp,   -0.5691378361_wp,    0.6665235373_wp,    0.2220242611_wp, &
      &    0.3978718931_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.4097476002_wp,    0.2321791150_wp,    0.8821563142_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0195140967_wp,    0.0489100378_wp,    0.4716043721_wp,    0.6814125133_wp, &
      &    0.4930247980_wp,    0.2596533451_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 43) = reshape([&
      &    0.0587793758_wp,   -0.4418564963_wp,    0.3891153906_wp,    0.0394440530_wp, &
      &    0.8051963743_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1570864172_wp,    0.9875848609_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0155256107_wp,    0.0306554740_wp,    0.3565599778_wp,    0.6257365560_wp, &
      &    0.6070665876_wp,    0.3340779852_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 44) = reshape([&
      &    0.1012945357_wp,   -0.4417810154_wp,    0.5708458165_wp,    0.0508827192_wp, &
      &    0.6827261194_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2194533991_wp,    0.5442443537_wp,    0.8097149431_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0185045975_wp,    0.0436915804_wp,    0.4264744668_wp,    0.6298245591_wp, &
      &    0.5760776044_wp,    0.2955059615_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 45) = reshape([&
      &    0.0958047071_wp,   -0.4671792466_wp,    0.5228765411_wp,    0.1257231962_wp, &
      &    0.6952401098_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2809328186_wp,    0.6126331745_wp,    0.7387539137_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0193020501_wp,    0.0437373438_wp,    0.5091770490_wp,    0.6709956193_wp, &
      &    0.4935172259_wp,    0.2113263707_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 46) = reshape([&
      &    0.0671079668_wp,   -0.4156987391_wp,    0.8246669695_wp,    0.3615116480_wp, &
      &    0.1092007181_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2790018623_wp,    0.7443059066_wp,    0.6067674005_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0156999562_wp,    0.0861759796_wp,    0.5030694772_wp,    0.6471015661_wp, &
      &    0.5266115176_wp,    0.2078176741_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 47) = reshape([&
      &    0.0423667102_wp,   -0.3728619236_wp,    0.5760906574_wp,    0.0553825864_wp, &
      &    0.7240382389_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2284179488_wp,    0.6628922630_wp,    0.7130210995_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0193535384_wp,    0.0475264099_wp,    0.4469645769_wp,    0.6765233025_wp, &
      &    0.5452166514_wp,    0.2065051380_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 48) = reshape([&
      &    0.2236991485_wp,   -0.5305743184_wp,    0.0299911504_wp,    0.6755556928_wp, &
      &    0.4595373983_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1144365509_wp,    0.5510156550_wp,    0.8266111684_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1249899763_wp,    0.9921580045_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 49) = reshape([&
      &    0.3661466726_wp,   -0.6974197331_wp,    0.1706557928_wp,    0.4704979227_wp, &
      &    0.3592361830_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2545400066_wp,   -0.3622422289_wp,    0.3899138385_wp,    0.8074386361_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1223247684_wp,    0.9924901264_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 50) = reshape([&
      &    0.2124540959_wp,   -0.6382014727_wp,    0.1323435953_wp,    0.5608008094_wp, &
      &    0.3400141442_wp,    0.3161331114_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1221698649_wp,   -0.2663948920_wp,    0.6816389227_wp,    0.5825717885_wp, &
      &    0.3317932729_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.4613945456_wp,    0.8871950593_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 51) = reshape([&
      &    0.4149286722_wp,   -0.7356616832_wp,    0.1597841612_wp,    0.2596645726_wp, &
      &    0.4183845653_wp,    0.1365055750_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1632463566_wp,   -0.3847696333_wp,    0.6368772771_wp,    0.5730845125_wp, &
      &    0.3021000360_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9469045087_wp,    0.3215149317_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 52) = reshape([&
      &    0.3324994583_wp,   -0.6347709373_wp,    0.0988561656_wp,    0.6523694592_wp, &
      &    0.0343140667_wp,    0.2235487854_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1944949458_wp,   -0.4360152590_wp,    0.4745169614_wp,    0.6645373154_wp, &
      &    0.3244783810_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9177859614_wp,    0.3970754703_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 53) = reshape([&
      &    0.3239071021_wp,   -0.7313251814_wp,    0.3661865576_wp,    0.2595904319_wp, &
      &    0.2785432606_wp,    0.2849237319_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.5288357293_wp,    0.5939775709_wp,    0.4938378671_wp,    0.1123448268_wp, &
      &    0.3332059689_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.8863035179_wp,    0.4631048199_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 54) = reshape([&
      &    0.3550594150_wp,   -0.7222432520_wp,    0.4480632958_wp,    0.0823462363_wp, &
      &    0.3733805293_wp,    0.0730948521_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.5765268770_wp,    0.6396608643_wp,    0.3233428102_wp,    0.3428821884_wp, &
      &    0.1906094718_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.6667619147_wp,    0.7452707891_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 55) = reshape([&
      &    0.5458629366_wp,   -0.7071105469_wp,    0.3797788016_wp,    0.2404087994_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2285090239_wp,    0.2061888350_wp,    0.9514566676_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 56) = reshape([&
      &    0.3007680965_wp,   -0.5554577476_wp,    0.7699842603_wp,    0.0919544982_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2651457936_wp,    0.6888574357_wp,    0.6750735939_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0359640602_wp,    0.9541239476_wp,    0.2899117146_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 57) = reshape([&
      &    0.1261279229_wp,   -0.6356538911_wp,    0.4435362803_wp,    0.6191215114_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.4666972802_wp,    0.8844171237_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0207772483_wp,    0.2164081058_wp,    0.3570173007_wp,    0.6778890354_wp, &
      &    0.6047651944_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 72) = reshape([&
      &    0.1799635903_wp,   -0.5188869273_wp,    0.6194645037_wp,    0.5609217338_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1997303709_wp,    0.4459737827_wp,    0.8724764547_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0141912136_wp,    0.3697863287_wp,    0.6360347355_wp,    0.4868541790_wp, &
      &    0.4706267142_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 73) = reshape([&
      &    0.1864297890_wp,   -0.6849874901_wp,    0.4515324195_wp,    0.3145449208_wp, &
      &    0.4395634643_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2311606652_wp,    0.3966132273_wp,    0.8884045783_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0219270334_wp,    0.3869901974_wp,    0.6476551830_wp,    0.4981149452_wp, &
      &    0.4268279015_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 74) = reshape([&
      &    0.1311946769_wp,   -0.5275683655_wp,    0.6101002551_wp,    0.2362105489_wp, &
      &    0.5257773595_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2123828159_wp,    0.3774504031_wp,    0.9013460671_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0173771030_wp,    0.4833993694_wp,    0.7014563488_wp,    0.4297499141_wp, &
      &    0.2988261838_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 75) = reshape([&
      &    0.1210272768_wp,   -0.4984452243_wp,    0.6077044030_wp,    0.0717470410_wp, &
      &    0.6020402622_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2330391211_wp,    0.1907272031_wp,    0.9535805693_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0234581866_wp,    0.0438566263_wp,    0.5210786479_wp,    0.7262045523_wp, &
      &    0.4456795942_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 76) = reshape([&
      &    0.1788482313_wp,   -0.5438886728_wp,    0.7237314593_wp,    0.2120313036_wp, &
      &    0.3216425390_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3360651002_wp,    0.7442276874_wp,    0.5772221390_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0306058204_wp,    0.5235924669_wp,    0.6558373050_wp,    0.4532357205_wp, &
      &    0.2989465227_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 77) = reshape([&
      &    0.1700613814_wp,   -0.6106138562_wp,    0.4850608334_wp,    0.4708040985_wp, &
      &    0.3758847347_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.4627575159_wp,    0.5264526807_wp,    0.7132342227_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0220562773_wp,    0.5908016683_wp,    0.6490580904_wp,    0.4035321275_wp, &
      &    0.2575894539_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 78) = reshape([&
      &    0.1164964235_wp,   -0.4851170031_wp,    0.7153708614_wp,    0.0997292349_wp, &
      &    0.4789453904_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.2712554175_wp,    0.6564649229_wp,    0.7038993561_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0261291136_wp,    0.5880215462_wp,    0.6558212214_wp,    0.4315548391_wp, &
      &    0.1928908422_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 79) = reshape([&
      &    0.1188357950_wp,   -0.4866826873_wp,    0.7314457030_wp,    0.0624610098_wp, &
      &    0.4583708340_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.3254799136_wp,    0.6939358345_wp,    0.6422739941_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0272684249_wp,    0.6327056671_wp,    0.6350441380_wp,    0.4016046712_wp, &
      &    0.1853984970_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 80) = reshape([&
      &    0.2416087904_wp,   -0.5771931434_wp,    0.0362222406_wp,    0.6542021044_wp, &
      &    0.4232975590_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.1494683605_wp,    0.6517415570_wp,    0.7435671807_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.5674809204_wp,    0.8233865465_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 81) = reshape([&
      &    0.3353276294_wp,   -0.6593936371_wp,    0.0376303224_wp,    0.5941995455_wp, &
      &    0.3134745146_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4239149653_wp,   -0.5131223337_wp,    0.3758390761_wp,    0.6447841203_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0007490924_wp,    0.9999997194_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 82) = reshape([&
      &    0.3193338362_wp,   -0.6524070942_wp,    0.0703401448_wp,    0.6339057741_wp, &
      &    0.2561378889_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3583901752_wp,   -0.5651433966_wp,    0.5092833625_wp,    0.5101018224_wp, &
      &    0.1805436544_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0461899809_wp,    0.9989326732_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 83) = reshape([&
      &    0.0747268679_wp,   -0.4691036489_wp,    0.1382502251_wp,    0.8261395758_wp, &
      &    0.1003229100_wp,    0.2258496723_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.4624764084_wp,   -0.6404816003_wp,    0.3432852845_wp,    0.3161388225_wp, &
      &    0.3964081749_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.9688837855_wp,    0.0135411388_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 84) = reshape([&
      &    0.3973908977_wp,   -0.7278893798_wp,    0.0729982415_wp,    0.5203849698_wp, &
      &    0.1517219196_wp,    0.1144933373_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.2299199503_wp,   -0.4452963825_wp,    0.6976918468_wp,    0.2371877382_wp, &
      &    0.4536694963_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.7451325995_wp,    0.6669163435_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 85) = reshape([&
      &    0.3309498638_wp,   -0.7519624346_wp,    0.0629128287_wp,    0.3033617172_wp, &
      &    0.4282202060_wp,    0.2136955411_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.1106898491_wp,   -0.3310392217_wp,    0.6236023488_wp,    0.6375332371_wp, &
      &    0.2878059644_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3703323028_wp,    0.9288993409_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

      coefficients(:, :, 86) = reshape([&
      &    0.3900659815_wp,   -0.7295660473_wp,    0.5075609283_wp,    0.0470260410_wp, &
      &    0.2276273819_wp,    0.0627546266_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.3050248300_wp,   -0.4836445058_wp,    0.5033509493_wp,    0.6296240849_wp, &
      &    0.1525096019_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0020336485_wp,    0.9999979321_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

   end subroutine setCGTOcoefficients

   subroutine setCGTOcoefficients_environment()
      ! set up the array of CGTO coeffients of environment dependence during initialization of the basis set data
      coefficients_env(:, :,  1) = reshape([&
      &    0.0003808289_wp,    0.0020833625_wp,    0.0116213375_wp,    0.0363847150_wp, &
      &    0.1630499574_wp,    0.0002196400_wp,    0.0001421331_wp,   -1.1038095111_wp, & ! shell 1
      &    0.0272895011_wp,    0.0112868737_wp,   -0.9359727329_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  2) = reshape([&
      &    0.0001630339_wp,    0.0003648991_wp,    0.0003367473_wp,    0.0028730204_wp, &
      &    0.0411062036_wp,    0.0010532063_wp,   -0.4859020691_wp,    0.0010004954_wp, & ! shell 1
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  3) = reshape([&
      &   -0.0000385956_wp,   -0.1362445659_wp,    0.1489403105_wp,    0.0024206777_wp, &
      &   -0.5050188536_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0009941299_wp,    0.0748439926_wp,   -0.0018340368_wp,   -0.4217706489_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  4) = reshape([&
      &    0.0002822327_wp,   -0.0006782595_wp,    0.0006391637_wp,    0.0004324185_wp, &
      &   -0.3614388051_wp,    0.0005167150_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0103615143_wp,    0.0721093603_wp,    0.0112351816_wp,   -0.3736729367_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  5) = reshape([&
      &    0.0011591943_wp,    0.0033175575_wp,    0.0053646298_wp,    0.0124937362_wp, &
      &   -0.3083420490_wp,    0.0145328556_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0011558020_wp,    0.0132990955_wp,    0.0018143704_wp,    0.0024575772_wp, &
      &   -0.4683962549_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  6) = reshape([&
      &    0.0005085813_wp,    0.0003488988_wp,    0.0004320064_wp,    0.0001993897_wp, &
      &   -0.1893046851_wp,    0.0001362983_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0009685398_wp,    0.0005561082_wp,    0.0017895920_wp,    0.0008545133_wp, &
      &    0.0006097436_wp,   -0.3120800745_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  7) = reshape([&
      &    0.0003632028_wp,    0.0003875061_wp,    0.0003747232_wp,    0.0001952960_wp, &
      &    0.0003614710_wp,   -0.1044325009_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0014125939_wp,    0.0020890789_wp,    0.0281829275_wp,    0.0138423801_wp, &
      &    0.0008452096_wp,   -0.3399406586_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  8) = reshape([&
      &   -0.0010995753_wp,    0.0038535712_wp,    0.0086922512_wp,    0.0011954736_wp, &
      &    0.0011134570_wp,   -0.0677639799_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0012087972_wp,    0.0060242778_wp,    0.0229813921_wp,    0.0359198618_wp, &
      &    0.0013254987_wp,   -0.2847389810_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :,  9) = reshape([&
      &   -0.0037565216_wp,    0.0067597866_wp,    0.0067527898_wp,    0.0030122921_wp, &
      &    0.0032295133_wp,   -0.1073873519_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0000642194_wp,    0.0005653914_wp,    0.0021046523_wp,    0.0027943546_wp, &
      &    0.0019670367_wp,   -0.2440352710_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 10) = reshape([&
      &   -0.0007256706_wp,    0.0011302804_wp,    0.0010797578_wp,    0.0013778589_wp, &
      &    0.0013875720_wp,   -0.0871650143_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0004560226_wp,    0.0012005386_wp,    0.0010009089_wp,    0.0012710018_wp, &
      &    0.0011923267_wp,   -0.2001746332_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 11) = reshape([&
      &    0.0000580887_wp,    0.0000335097_wp,    0.0012275034_wp,   -0.6982007616_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0013249877_wp,    0.0019300598_wp,    0.0048114942_wp,   -0.7097998587_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 12) = reshape([&
      &    0.0008479145_wp,   -0.0011083746_wp,    0.0050506852_wp,   -0.1550731970_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0023707942_wp,    0.0361173347_wp,   -0.1833273166_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0937214821_wp,   -0.1348429580_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 13) = reshape([&
      &    0.0006884157_wp,   -0.0013113628_wp,    0.0006110640_wp,    0.0013895142_wp, &
      &   -0.1866390057_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0005608897_wp,    0.0117899875_wp,    0.0040932572_wp,   -0.2841357337_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0218111855_wp,   -0.2387468664_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 14) = reshape([&
      &    0.0028629210_wp,   -0.0001503051_wp,    0.0021009223_wp,    0.0015511007_wp, &
      &   -0.2578865994_wp,    0.0018320798_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0107438354_wp,    0.0417972130_wp,    0.0005666162_wp,   -0.2438936577_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0846840315_wp,   -0.0769693800_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 15) = reshape([&
      &    0.0004558321_wp,   -0.0004677458_wp,    0.0139789655_wp,    0.0018562011_wp, &
      &   -0.1002418646_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0019421566_wp,    0.0102792294_wp,    0.0015805240_wp,    0.0012925660_wp, &
      &   -0.3636186679_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0301630882_wp,   -0.0419876704_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 16) = reshape([&
      &    0.0006346874_wp,   -0.0006261349_wp,    0.0257111759_wp,    0.0020110327_wp, &
      &   -0.0875112585_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0030055451_wp,    0.0041955356_wp,    0.0029938858_wp,    0.0034840493_wp, &
      &   -0.2763737435_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0238038568_wp,   -0.0297046219_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 17) = reshape([&
      &    0.0005931335_wp,   -0.0001221984_wp,    0.0230716315_wp,    0.0014845030_wp, &
      &   -0.0881611865_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0009042565_wp,    0.0028339841_wp,    0.0007827647_wp,    0.0008306077_wp, &
      &   -0.2306632947_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0682582701_wp,   -0.0230052285_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 18) = reshape([&
      &    0.0004047334_wp,   -0.0005502849_wp,    0.0099346506_wp,    0.0010944114_wp, &
      &   -0.0438203489_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0008250201_wp,    0.0012005067_wp,    0.0010562672_wp,    0.0012675075_wp, &
      &   -0.1723985236_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0115891534_wp,   -0.0129877889_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 19) = reshape([&
      &   -0.0010723692_wp,   -0.0027400352_wp,   -0.0042587738_wp,   -0.4632300283_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0685858067_wp,    0.0233686207_wp,   -0.5083136484_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 20) = reshape([&
      &   -0.0000216751_wp,   -0.0007520554_wp,    0.0011167357_wp,   -0.0821372075_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0399902075_wp,    0.0237609502_wp,   -0.1477992118_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0218441359_wp,    0.0010957731_wp,    0.0027319660_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 21) = reshape([&
      &    0.0001513813_wp,   -0.0043525817_wp,    0.0002898642_wp,    0.0006866745_wp, &
      &   -0.0810516852_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0006298886_wp,    0.0005769238_wp,   -0.0678341151_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001346888_wp,    0.0001124630_wp,    0.0003376227_wp,    0.0006123464_wp, &
      &    0.0198279015_wp,   -0.0223928380_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 22) = reshape([&
      &    0.0010768508_wp,   -0.0008937884_wp,    0.0003478071_wp,    0.0006278550_wp, &
      &   -0.1733619862_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0008727934_wp,    0.0131008200_wp,   -0.1955240791_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0002789258_wp,    0.0003908913_wp,    0.0005945402_wp,    0.0005414832_wp, &
      &    0.0246394262_wp,   -0.0212695770_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 23) = reshape([&
      &    0.0005912209_wp,   -0.0010085180_wp,    0.0004183270_wp,    0.0008043609_wp, &
      &   -0.2249974308_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0013524684_wp,    0.0164055533_wp,   -0.2667873633_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001788220_wp,    0.0007876598_wp,    0.0006722045_wp,    0.0003694450_wp, &
      &    0.0237432112_wp,   -0.0274577146_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 24) = reshape([&
      &    0.0003328899_wp,   -0.0010296609_wp,    0.0007066219_wp,    0.0015738078_wp, &
      &   -0.4348462068_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0221533717_wp,    0.0277912524_wp,   -0.5740954283_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001787774_wp,    0.0011279011_wp,    0.0005545761_wp,    0.0002379271_wp, &
      &    0.0017866173_wp,   -0.0192102376_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 25) = reshape([&
      &    0.0018305132_wp,   -0.0003116346_wp,    0.0007064445_wp,    0.0053159860_wp, &
      &   -0.2178704553_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0250539570_wp,    0.0138614059_wp,   -0.2124367029_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001223261_wp,    0.0007714939_wp,    0.0007532423_wp,    0.0002618335_wp, &
      &    0.0171804168_wp,   -0.0188093933_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 26) = reshape([&
      &    0.0005065770_wp,   -0.0003007638_wp,    0.0004134233_wp,    0.0016157303_wp, &
      &   -0.4561613178_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0033346588_wp,    0.0101004111_wp,   -0.7344926757_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001864304_wp,    0.0027359587_wp,    0.0010359357_wp,    0.0002085049_wp, &
      &    0.0006294295_wp,   -0.0257097834_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 27) = reshape([&
      &    0.0004749423_wp,   -0.0007282872_wp,    0.0003910907_wp,    0.0010556999_wp, &
      &   -0.4427980223_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0023370590_wp,    0.0093828881_wp,   -0.5680074752_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001780441_wp,    0.0040033865_wp,    0.0006962059_wp,    0.0001951214_wp, &
      &    0.0004117783_wp,   -0.0139364355_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 28) = reshape([&
      &    0.0009726596_wp,   -0.0002178162_wp,    0.0002655605_wp,    0.0014455686_wp, &
      &   -0.3284347063_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0036569166_wp,    0.0171319017_wp,   -0.4757284063_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001798191_wp,    0.0023521630_wp,    0.0100036496_wp,    0.0002514949_wp, &
      &    0.0001780043_wp,   -0.0167995705_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 29) = reshape([&
      &    0.0004831043_wp,   -0.0008326102_wp,    0.0003399436_wp,    0.0004045033_wp, &
      &   -0.3231385465_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0023314323_wp,    0.0107939489_wp,   -0.3363085062_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0001415001_wp,    0.0005829268_wp,    0.0002980210_wp,    0.0002171128_wp, &
      &    0.0004782984_wp,   -0.0151705605_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 30) = reshape([&
      &    0.0013555077_wp,   -0.0015226887_wp,    0.0030392991_wp,    0.0031521951_wp, &
      &   -0.2480249972_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0042373970_wp,    0.0028590051_wp,   -0.4719353448_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 31) = reshape([&
      &    0.0024336683_wp,   -0.0035436911_wp,    0.0029914258_wp,    0.0020216060_wp, &
      &   -0.0976316600_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0042471149_wp,    0.0028914609_wp,    0.0037430350_wp,   -0.1473080795_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0685773508_wp,   -0.1361609662_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 32) = reshape([&
      &    0.0034992414_wp,   -0.0034219756_wp,    0.0051878575_wp,    0.0051396995_wp, &
      &   -0.1034713120_wp,    0.0028335502_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0026224709_wp,    0.0030919590_wp,    0.0051793071_wp,    0.0050643813_wp, &
      &   -0.2347366970_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0391620513_wp,   -0.2433370382_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 33) = reshape([&
      &    0.0025075701_wp,   -0.0022661062_wp,    0.0175775576_wp,    0.0027790037_wp, &
      &    0.0022501770_wp,   -0.0838840930_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0033676455_wp,    0.0062447868_wp,    0.0035476758_wp,    0.0024528555_wp, &
      &   -0.2974758717_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1823693397_wp,   -0.3864755159_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 34) = reshape([&
      &    0.0058321169_wp,   -0.0018424305_wp,    0.0094469163_wp,    0.0019988887_wp, &
      &    0.0027091038_wp,   -0.0713508831_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0017782813_wp,    0.0030692690_wp,    0.0022745804_wp,    0.0022962952_wp, &
      &   -0.2810869091_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0323592833_wp,   -0.0694028199_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 35) = reshape([&
      &    0.0045926605_wp,   -0.0009433182_wp,    0.0026764202_wp,    0.0042044861_wp, &
      &    0.0025519350_wp,   -0.0646739461_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0036338974_wp,    0.0066364939_wp,    0.0041334253_wp,    0.0035339228_wp, &
      &   -0.1937775078_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0431438353_wp,   -0.1424263046_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 36) = reshape([&
      &    0.0035020445_wp,   -0.0008581097_wp,    0.0042767452_wp,    0.0019325250_wp, &
      &    0.0019045265_wp,    0.0036620917_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0032787614_wp,    0.0049504733_wp,    0.0037102406_wp,    0.0032176967_wp, &
      &   -0.1284601587_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0165079679_wp,   -0.0168518674_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 37) = reshape([&
      &    0.0023714283_wp,   -0.0049791284_wp,    0.0031289148_wp,   -0.2875443671_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0007810684_wp,    0.0014332687_wp,   -0.3052790494_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 38) = reshape([&
      &    0.0036160855_wp,    0.0000397704_wp,    0.0002730774_wp,   -0.0719584717_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0424813739_wp,    0.0009692495_wp,   -0.1257800319_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0360331326_wp,    0.0131924994_wp,   -0.0473478320_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 39) = reshape([&
      &    0.0003189778_wp,   -0.0060949718_wp,    0.0004051687_wp,    0.0000853855_wp, &
      &   -0.0481762316_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0020369937_wp,    0.0118045691_wp,   -0.0710368590_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000083523_wp,    0.0001011423_wp,    0.0002330005_wp,    0.0033454767_wp, &
      &    0.0008964894_wp,   -0.0345153122_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 40) = reshape([&
      &    0.0009984033_wp,   -0.0006790304_wp,    0.0401052871_wp,    0.0011296314_wp, &
      &   -0.2788036683_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0018598871_wp,    0.0076029194_wp,   -0.3304709420_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0004178000_wp,    0.0004310824_wp,    0.0017713470_wp,    0.0005718097_wp, &
      &    0.0010255796_wp,   -0.1220181000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 41) = reshape([&
      &    0.0005638263_wp,   -0.0009048162_wp,    0.0190441365_wp,    0.0007399897_wp, &
      &   -0.1677502860_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0012793573_wp,    0.0081106685_wp,   -0.1875880318_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0002388125_wp,    0.0003359272_wp,    0.0005937854_wp,    0.0002686623_wp, &
      &    0.0021249787_wp,   -0.0413295060_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 42) = reshape([&
      &    0.0001456450_wp,   -0.0001116225_wp,    0.0002189920_wp,    0.0001424016_wp, &
      &   -0.2996573522_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0078478369_wp,    0.0020122103_wp,   -0.0936949106_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0003742339_wp,    0.0001823958_wp,    0.0007659715_wp,    0.0007785715_wp, &
      &    0.0238740025_wp,   -0.0677167009_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 43) = reshape([&
      &    0.0001133880_wp,   -0.0003465599_wp,   -0.0000127632_wp,    0.0000479096_wp, &
      &   -0.0439646984_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0011559357_wp,   -0.0398968634_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000707823_wp,    0.0001860024_wp,    0.0003507711_wp,    0.0002162158_wp, &
      &    0.0057668574_wp,   -0.0005310959_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 44) = reshape([&
      &    0.0000556176_wp,   -0.0001359446_wp,   -0.0000189402_wp,    0.0000592616_wp, &
      &   -0.2552416434_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0076277945_wp,    0.0175040215_wp,   -0.3010546730_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000797443_wp,    0.0002418592_wp,    0.0009726422_wp,    0.0004669221_wp, &
      &    0.0100708608_wp,   -0.0290551432_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 45) = reshape([&
      &    0.0000732756_wp,   -0.0001191377_wp,    0.0000229564_wp,    0.0001468020_wp, &
      &   -0.1028034623_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0018701235_wp,    0.0174733261_wp,   -0.0912300895_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000584147_wp,    0.0002150019_wp,    0.0012266634_wp,    0.0002456909_wp, &
      &    0.0104336948_wp,   -0.0056412239_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 46) = reshape([&
      &    0.0000491154_wp,   -0.0000876082_wp,    0.0000098325_wp,    0.0001042717_wp, &
      &   -0.1803128807_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0021948569_wp,    0.0049370505_wp,   -0.6617956640_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001339656_wp,    0.0003323879_wp,    0.0264161467_wp,    0.0011649641_wp, &
      &    0.0005102223_wp,   -0.1259780977_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 47) = reshape([&
      &    0.0000688071_wp,   -0.0004346478_wp,    0.0000188401_wp,    0.0001469365_wp, &
      &   -0.1849445219_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0016621029_wp,    0.0122924851_wp,   -0.1829829677_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000515598_wp,    0.0002215167_wp,    0.0006845323_wp,    0.0005453090_wp, &
      &    0.0012047956_wp,   -0.0007680746_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 48) = reshape([&
      &    0.0014907253_wp,   -0.0026179620_wp,    0.0051503229_wp,    0.0036859637_wp, &
      &   -0.1406178704_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0086898468_wp,    0.0033005910_wp,   -0.3037808704_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.6525413386_wp,   -0.1841249824_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 49) = reshape([&
      &    0.0043028690_wp,   -0.0039857632_wp,    0.0062405308_wp,    0.0039735309_wp, &
      &   -0.1164928818_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0066497251_wp,   -0.0082407889_wp,    0.0034846589_wp,   -0.2679328111_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1249355809_wp,   -0.3317800236_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 50) = reshape([&
      &    0.0029667858_wp,   -0.0031941302_wp,    0.0038055426_wp,    0.0030136555_wp, &
      &    0.0027444233_wp,   -0.1085732791_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0019870063_wp,   -0.0034213016_wp,    0.0052452268_wp,    0.0027469567_wp, &
      &   -0.2775268348_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0708063348_wp,   -0.3410053453_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 51) = reshape([&
      &    0.0020765131_wp,   -0.0037233325_wp,    0.0026854350_wp,    0.0026692290_wp, &
      &    0.0027472640_wp,    0.0016052866_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0032665030_wp,   -0.0038503489_wp,    0.0042825019_wp,    0.0036883149_wp, &
      &   -0.3262836861_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0458723610_wp,   -0.0485562741_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 52) = reshape([&
      &    0.0009896423_wp,   -0.0008325018_wp,    0.0116493425_wp,    0.0009439708_wp, &
      &    0.0012726315_wp,   -0.0373841139_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0039447243_wp,   -0.0033794749_wp,    0.0046837609_wp,    0.0036140165_wp, &
      &   -0.2410559632_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0897060699_wp,   -0.3274263753_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 53) = reshape([&
      &    0.0026927854_wp,   -0.0013928706_wp,    0.0087043918_wp,    0.0021821807_wp, &
      &    0.0023089281_wp,   -0.0317686874_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0016410052_wp,    0.0080151275_wp,    0.0053256794_wp,    0.0039357649_wp, &
      &   -0.1453868629_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0884636878_wp,   -0.2282551354_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 54) = reshape([&
      &    0.0046688729_wp,   -0.0008287304_wp,    0.0022166294_wp,    0.0034455032_wp, &
      &    0.0021289710_wp,   -0.0356839382_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0027619117_wp,    0.0050304093_wp,    0.0046699789_wp,    0.0034161871_wp, &
      &   -0.1037995199_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0547575281_wp,   -0.4139770255_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 55) = reshape([&
      &    0.0043727559_wp,   -0.0044178370_wp,    0.0053668790_wp,   -0.1591063835_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0064528896_wp,    0.0020271757_wp,   -0.3061057211_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 56) = reshape([&
      &    0.0092757981_wp,    0.0010340312_wp,    0.0016738215_wp,   -0.0802605035_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0871787579_wp,    0.0002257955_wp,   -0.0856975708_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1275531766_wp,    0.0014351581_wp,   -0.0115873399_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 57) = reshape([&
      &    0.0012669116_wp,   -0.0053441555_wp,    0.0000329574_wp,   -0.0238146498_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0015660976_wp,   -0.0246615421_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000396602_wp,    0.0002599301_wp,    0.0144113568_wp,    0.0001570797_wp, &
      &   -0.0201210343_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 72) = reshape([&
      &    0.0000342688_wp,   -0.0003356282_wp,    0.0000088311_wp,   -0.1686507420_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0005867764_wp,    0.0029108816_wp,   -0.2534008643_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001241592_wp,    0.0001117522_wp,    0.0010686922_wp,    0.0008482590_wp, &
      &   -0.1475274260_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 73) = reshape([&
      &    0.0000519728_wp,   -0.0002505903_wp,    0.0000146085_wp,    0.0001750814_wp, &
      &   -0.1186792437_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0018473420_wp,    0.0001603969_wp,   -0.1919991778_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001865685_wp,    0.0073607825_wp,    0.0139514023_wp,    0.0008513136_wp, &
      &   -0.0947167876_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 74) = reshape([&
      &    0.0002181772_wp,   -0.0000919114_wp,    0.0000259654_wp,    0.0001901507_wp, &
      &   -0.1564668132_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0050503627_wp,    0.0121456570_wp,   -0.2095572244_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000886208_wp,    0.0007132972_wp,    0.0009477216_wp,    0.0040617825_wp, &
      &   -0.0715531334_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 75) = reshape([&
      &    0.0001176286_wp,   -0.0003028426_wp,    0.0000396282_wp,    0.0001820327_wp, &
      &   -0.0731536902_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0008276011_wp,    0.0019619053_wp,   -0.0986684535_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000640753_wp,    0.0004656600_wp,    0.0007885931_wp,    0.0044198570_wp, &
      &   -0.0077685082_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 76) = reshape([&
      &    0.0004089258_wp,   -0.0002401379_wp,    0.0001005822_wp,    0.0005105471_wp, &
      &   -0.2292759055_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0017023725_wp,    0.0038859824_wp,   -0.3139001357_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0001014193_wp,    0.0009695151_wp,    0.0005902329_wp,    0.0028378595_wp, &
      &   -0.1148586186_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 77) = reshape([&
      &    0.0001061123_wp,   -0.0003847915_wp,    0.0000647033_wp,    0.0001742429_wp, &
      &   -0.1930539939_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0012864796_wp,    0.0040582089_wp,   -0.1118622317_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0003856785_wp,    0.0006001519_wp,    0.0007483259_wp,    0.0015307698_wp, &
      &   -0.0653962180_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 78) = reshape([&
      &    0.0001192736_wp,   -0.0003113305_wp,    0.0000500334_wp,    0.0003059971_wp, &
      &   -0.1900983500_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0028091779_wp,    0.0121153360_wp,   -0.3007932893_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000834291_wp,    0.0061276397_wp,    0.0007877087_wp,    0.0006512936_wp, &
      &   -0.0439238695_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 79) = reshape([&
      &    0.0001150373_wp,   -0.0005079669_wp,    0.0000351473_wp,    0.0001996096_wp, &
      &   -0.2216419043_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0018213635_wp,    0.0193376087_wp,   -0.3757074584_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.0000775031_wp,    0.0006961823_wp,    0.0008897587_wp,    0.0010657364_wp, &
      &   -0.0594262162_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 80) = reshape([&
      &    0.0017003633_wp,   -0.0014467132_wp,    0.0094341304_wp,    0.0044763608_wp, &
      &   -0.1033504808_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &   -0.0101005022_wp,    0.0042246110_wp,   -0.2421040932_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &   -0.1655699561_wp,    0.9695915903_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 81) = reshape([&
      &    0.0027344257_wp,   -0.0026051603_wp,    0.0022160816_wp,    0.0059115487_wp, &
      &   -0.0717491661_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0028165512_wp,   -0.0073210322_wp,    0.0090588052_wp,   -0.1819394929_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0782481319_wp,   -0.0449672907_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 82) = reshape([&
      &    0.0035517466_wp,   -0.0030336860_wp,    0.0026093813_wp,    0.0064965670_wp, &
      &   -0.0719044947_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0026981176_wp,   -0.0039250336_wp,    0.0046231167_wp,    0.0046533989_wp, &
      &   -0.2208197506_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3304205341_wp,   -0.0623814931_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 83) = reshape([&
      &    0.0009480430_wp,   -0.0023016636_wp,    0.0025909610_wp,    0.0016384491_wp, &
      &    0.0022608932_wp,    0.0002957354_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0023245533_wp,   -0.0016419241_wp,    0.0021036384_wp,    0.0024109803_wp, &
      &   -0.1881827662_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1958753681_wp,   -0.3770965841_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 84) = reshape([&
      &    0.0037072014_wp,   -0.0010593912_wp,    0.0021573370_wp,    0.0020623739_wp, &
      &    0.0017189761_wp,   -0.0372538208_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0059386951_wp,   -0.0011952929_wp,    0.0030570801_wp,    0.0025726723_wp, &
      &   -0.2500950258_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.1370111147_wp,   -0.4883795868_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 85) = reshape([&
      &    0.0024551126_wp,   -0.0010234111_wp,    0.0012970758_wp,    0.0017478145_wp, &
      &    0.0015898649_wp,   -0.0155864357_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0029341631_wp,   -0.0018312706_wp,    0.0094849533_wp,    0.0028360669_wp, &
      &   -0.1874683304_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.0569051402_wp,   -0.0523177867_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))
      
      coefficients_env(:, :, 86) = reshape([&
      &    0.0040109203_wp,   -0.0014211529_wp,    0.0013887698_wp,    0.0013345295_wp, &
      &    0.0020485980_wp,   -0.0185253889_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 1
      &    0.0025246756_wp,   -0.0032642541_wp,    0.0004296680_wp,    0.0269286179_wp, &
      &   -0.0145864333_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 2
      &    0.3303148086_wp,   -0.4087710695_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 3
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 4
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 5
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, & ! shell 6
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp, &
      &    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp,    0.0000000000_wp], (/max_prim, max_shell/))

   end subroutine setCGTOcoefficients_environment

end module tblite_basis_qvszp
