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

!> @file tblite/xtb/gp3.f90
!> Provides the parametrization for the GP3-xTB Hamiltonian

!> Implementation of the GP3-xTB Hamiltonian to parametrize an xTB calculator.
module tblite_xtb_gp3
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_symbols, only : to_symbol
   use tblite_basis_type, only : basis_type
   use tblite_basis_qvszp, only : qvszp_basis_type, new_basis, add_qvszp_cgtos, qvszp_cgto_type
   use tblite_coulomb_charge, only : new_effective_coulomb, effective_coulomb, &
      & arithmetic_average, coulomb_kernel
   use tblite_coulomb_multipole, only : new_damped_multipole
   use tblite_coulomb_thirdorder, only : new_onsite_thirdorder
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_disp, only : d4_dispersion, new_d4_dispersion
   use tblite_ncoord, only : new_ncoord
   use tblite_param, only : param_record
   use tblite_repulsion, only : new_repulsion
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : new_hamiltonian
   use tblite_xtb_spec, only : tb_h0spec
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: new_gp3_calculator
   public :: export_gp3_param

   integer, parameter :: max_elem = 103
   integer, parameter :: max_shell = 4

   ! Global parameters: -----------------------------------------------------------

   !> Global Parameters: Dispersion 
   real(wp), parameter :: d4_s6 = 1.0000000000_wp
   real(wp), parameter :: d4_s8 = 6.6071949202_wp
   real(wp), parameter :: d4_a1 = 0.8528178637_wp
   real(wp), parameter :: d4_a2 = 0.0000000000_wp
   real(wp), parameter :: d4_s9 = 1.5000000000_wp
   real(wp), parameter :: d4_beta1 = 0.0000000000_wp
   real(wp), parameter :: d4_beta2 = 0.0000000000_wp

   !> Global Parameter: Diatomic frame scaling for different orbital pairs
   real(wp), parameter :: kdiat_ss = 4.6544537664_wp
   real(wp), parameter :: kdiat_sp = 1.3671012900_wp
   real(wp), parameter :: kdiat_pp = 0.3222196243_wp
   real(wp), parameter :: kdiat_sd = 0.6328106033_wp
   real(wp), parameter :: kdiat_pd = 0.4979124643_wp
   real(wp), parameter :: kdiat_dd = 0.9061889851_wp

   !> Global Parameter: Coordination number exponents
   real(wp), parameter :: rep_exp_cn = 1.8120000000_wp
   real(wp), parameter :: se_exp_cn  = 2.0000000000_wp

   !> Global Parameter: Exponent of the repulsion polynomial
   real(wp), parameter :: rep_rexp = 2.0000000000_wp

   !> Global Parameter: Shell dependence of H0 scaling
   real(wp), parameter :: kh0ss = 0.6752395702_wp
   real(wp), parameter :: kh0pp = 0.9434008972_wp
   real(wp), parameter :: kh0dd = 1.1000000000_wp
   real(wp), parameter :: kh0ff = 2.7000000000_wp

   !> Parameter: CN-dependence of first-order IP/EA for non-H atoms
   real(wp), parameter :: kipea_cn_nonH = 0.2000000000_wp
   !> Parameter: CN-dependence of first-order IP/EA for H atoms
   real(wp), parameter :: kipea_cn_H = 0.5000000000_wp

   !> Global Parameter: Second order Klopman-Ohno ES
   real(wp), parameter :: gexp = 2.0_wp

   !> Global Parameter: R-dependence of the 3rd order KO-gamma
   real(wp), parameter :: k3rd_split   = -0.9000000000_wp
   real(wp), parameter :: k3rd_scaleMN =  0.4500208536_wp
   !> Global Parameter: Shell dependence of hubbard derivs for 3rd order 
   real(wp), parameter :: kll_hubbard_derivs(0:3) = &
      & -0.0347117344_wp * [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp]

   !> Global Parameter: Scaling of third order hubbard derivatives for the 4th order
   real(wp), parameter :: k4th_scale = -0.0112920259_wp

   !> Global Parameter: Damping of multipole electrostatics for dipole and quadrupole
   real(wp), parameter :: kmp_ddamp =  0.190_wp
   real(wp), parameter :: kmp_qdamp =  0.320_wp
   !> Global Parameter: Z-dependent scaling of the D3 radii for multipole electrostatics
   real(wp), parameter :: kmp_rc0  =  1.7500_wp
   real(wp), parameter :: kmp_rc1  = -0.0750_wp

   !> Global Parameter: base-level Fock exchange and range-separation
   real(wp), parameter :: fock_ax    = 0.150_wp
   real(wp), parameter :: fock_omega = 0.300_wp
   !> Global Parameter: Scaling of the diagonal in the Mulliken exchange matrix
   real(wp), parameter :: kdiag_fock = 1.2788811412_wp
   !> Global Parameter: Shell-dependent scaling
   !> of the off-diagonals in the Mulliken exchange matrix
   real(wp), parameter :: kll_fock(0:3) = &
      & [0.3047852845_wp, 1.1594840133_wp, 1.1594840133_wp, 1.1594840133_wp]

   !> Global Parameter: Shell-dependent scaling of the LDA exchange
   real(wp), parameter :: kll_ldax(0:3) = &
      & [-0.0024757921_wp, -0.0408433413_wp, -0.0410000000_wp, -0.0410000000_wp]

   !> Global Parameter: Shell-dpendent scaling for the shell polynomials
   real(wp), parameter :: kll_shpoly(0:3) = &
      & 0.2_wp * [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp]

   ! Atom/shell-resolved parameter: -----------------------------------------------

   !> Parameter: Exponent in GP3-xTB repulsion
   real(wp), parameter :: rep_alpha(max_elem) = [&
   &  0.2342286016_wp,  0.8803317221_wp,  0.5146531158_wp,  0.4280997510_wp, & !1-4
   &  0.5058996423_wp,  0.4868704255_wp,  0.4608909036_wp,  0.3755559418_wp, & !5-8
   &  0.4975431361_wp,  0.5733862867_wp,  0.4461313165_wp,  0.4822649381_wp, & !9-12
   &  0.4096852133_wp,  0.4439802705_wp,  0.4668800640_wp,  0.4902942052_wp, & !13-16
   &  0.4513441213_wp,  0.5404520531_wp,  0.5123109723_wp,  0.4594484167_wp, & !17-20
   &  0.5580392072_wp,  0.5942261103_wp,  0.4728199835_wp,  0.4891777544_wp, & !21-24
   &  0.4046019277_wp,  0.3003640806_wp,  0.4041305017_wp,  0.4564941003_wp, & !25-28
   &  0.5597231752_wp,  0.4746080808_wp,  0.3588844719_wp,  0.3850462064_wp, & !29-32
   &  0.4432923047_wp,  0.4299301233_wp,  0.4588430504_wp,  0.5832759445_wp, & !33-36
   &  0.2827077076_wp,  0.4848640235_wp,  0.3713045775_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.4483902791_wp, & !45-48
   &  0.3718305452_wp,  0.3485258469_wp,  0.4178838280_wp,  0.4495766534_wp, & !49-52
   &  0.4464780816_wp,  0.4200372112_wp,  0.3623188084_wp,  0.4415169995_wp, & !53-56
   &  0.4024689782_wp,  0.4675189420_wp,  0.4675189420_wp,  0.4840295743_wp, & !57-60
   &  0.5001658174_wp,  0.5159276714_wp,  0.5313151362_wp,  0.5609668985_wp, & !61-64
   &  0.5463282119_wp,  0.5609668985_wp,  0.5752311959_wp,  0.5891211043_wp, & !65-68
   &  0.6026366234_wp,  0.6157777535_wp,  0.6529548088_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp,  0.3432570211_wp,  0.0000000000_wp, & !77-80
   &  0.4135979785_wp,  0.3701652928_wp,  0.4678490284_wp,  0.3204969833_wp, & !81-84
   &  0.3665530025_wp,  0.4280887198_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp,  0.6285444944_wp,  0.6409368462_wp,  0.6529548088_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103

   !> Parameter: Reference radii in GP3-xTB repulsion
   real(wp), parameter :: rep_rcov(max_elem) = [&
   &  0.6719079170_wp,  1.7019143884_wp,  1.2509057985_wp,  1.1064591188_wp, & !1-4
   &  1.0462436357_wp,  0.9852758482_wp,  0.8489038165_wp,  0.6513920820_wp, & !5-8
   &  0.7811154636_wp,  1.1152976802_wp,  1.1127693304_wp,  1.3352730420_wp, & !9-12
   &  1.3915878165_wp,  1.2195716461_wp,  1.1057052621_wp,  1.0988287529_wp, & !13-16
   &  0.9886349126_wp,  1.1352768544_wp,  1.6450241053_wp,  1.4869197542_wp, & !17-20
   &  1.2934217462_wp,  1.2836422943_wp,  1.0606302719_wp,  1.1233557717_wp, & !21-24
   &  0.9992247942_wp,  0.8760813065_wp,  0.9949026386_wp,  1.0242878062_wp, & !25-28
   &  1.2438309971_wp,  1.3846383929_wp,  1.0792557222_wp,  1.0936895329_wp, & !29-32
   &  1.1807319955_wp,  1.0733037192_wp,  1.0438219580_wp,  1.2765949855_wp, & !33-36
   &  1.4100637807_wp,  1.3235008819_wp,  1.3264570364_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  1.5017708792_wp, & !45-48
   &  1.2647681764_wp,  1.1204697281_wp,  1.2235735183_wp,  1.2390063744_wp, & !49-52
   &  1.2083549203_wp,  1.0427995340_wp,  1.8266360341_wp,  1.5230371773_wp, & !53-56
   &  1.3086636164_wp,  1.5612482050_wp,  1.5612482050_wp,  1.5069582615_wp, & !57-60
   &  1.4568089357_wp,  1.4108002276_wp,  1.3689321372_wp,  1.2976178094_wp, & !61-64
   &  1.3312046644_wp,  1.2976178094_wp,  1.2681715721_wp,  1.2428659524_wp, & !65-68
   &  1.2217009505_wp,  1.2046765662_wp,  1.1784471196_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp,  1.0400415214_wp,  0.0000000000_wp, & !77-80
   &  1.4182273357_wp,  1.2122677379_wp,  1.2968585653_wp,  0.9899494203_wp, & !81-84
   &  1.0204817816_wp,  0.9504266427_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp,  1.1917927996_wp,  1.1830496508_wp,  1.1784471196_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103

   !> Parameter: Environment-independent effective nuclear charge in GP3-xTB repulsion
   real(wp), parameter :: rep_zeff(max_elem) = [&
   &  0.8349280246_wp,  0.5164194548_wp,  3.3186588270_wp,  2.7483927859_wp, & !1-4
   &  3.1514033253_wp,  3.2715032000_wp,  3.0825487458_wp,  2.9582589564_wp, & !5-8
   &  2.5414641410_wp,  2.1050638756_wp,  6.2584193721_wp,  5.6489204656_wp, & !9-12
   &  4.5532494572_wp,  5.5733558960_wp,  6.3218383567_wp,  5.5237011332_wp, & !13-16
   &  4.8839050775_wp,  5.8546075966_wp,  8.9371689515_wp,  7.0245743182_wp, & !17-20
   &  9.6339891257_wp,  7.3788218558_wp,  6.2288499912_wp,  4.9518119155_wp, & !21-24
   &  4.3232727235_wp,  4.1065647232_wp,  4.5207442222_wp,  3.8482474498_wp, & !25-28
   &  4.3212666849_wp,  4.0017316652_wp,  5.1243914652_wp,  5.5542716202_wp, & !29-32
   &  5.8725116279_wp,  7.0550850475_wp,  6.5301788127_wp,  7.7948428423_wp, & !33-36
   &  5.4924054624_wp, 14.8360773518_wp,  7.3982272886_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  4.0137483913_wp, & !45-48
   &  5.1606584985_wp,  6.3804205736_wp,  7.4068746083_wp,  8.1199377759_wp, & !49-52
   &  7.3291974802_wp,  9.7836714040_wp,  7.0901052467_wp, 12.6809646309_wp, & !53-56
   &  8.8916204375_wp,  9.0758898478_wp,  9.0758898478_wp, 10.1375150241_wp, & !57-60
   & 11.0851607900_wp, 11.9188271454_wp, 12.6385140903_wp, 13.7359497489_wp, & !61-64
   & 13.2442216248_wp, 13.7359497489_wp, 14.1136984625_wp, 14.3774677657_wp, & !65-68
   & 14.5272576584_wp, 14.5630681406_wp, 13.9866231247_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp,  3.2775582896_wp,  0.0000000000_wp, & !77-80
   &  5.2746307451_wp,  6.6463052645_wp,  9.1394388158_wp,  9.0177604805_wp, & !81-84
   &  9.2445480794_wp, 10.9933073945_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp, 14.4848992124_wp, 14.2927508738_wp, 13.9866231247_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103

   !> Parameter: CN-dependence of Zeff in GP3-xTB repulsion
   real(wp), parameter :: rep_cn(max_elem) = [&
   &  0.1888951308_wp, -0.0150613309_wp, -0.1711013242_wp, -0.0049403131_wp, & !1-4
   &  0.0360342527_wp,  0.0233077978_wp,  0.0218716480_wp,  0.0859251997_wp, & !5-8
   & -0.0164621829_wp,  0.0108864040_wp, -0.0356858269_wp, -0.1139920806_wp, & !9-12
   & -0.0159103005_wp,  0.0336457887_wp,  0.0227112823_wp,  0.0673624372_wp, & !13-16
   &  0.0412699810_wp,  0.0130095413_wp, -0.1619818551_wp, -0.2030414817_wp, & !17-20
   & -0.4665626052_wp,  0.0368457085_wp,  0.0138051194_wp,  0.0022702526_wp, & !21-24
   &  0.0845629411_wp, -0.0040474363_wp, -1.3164892620_wp,  0.0891036378_wp, & !25-28
   &  0.1664338757_wp, -0.0642067810_wp, -0.8348305117_wp, -0.0879580034_wp, & !29-32
   & -0.0544175890_wp, -0.0428588849_wp, -0.0031494668_wp, -0.0276481717_wp, & !33-36
   & -0.0878684118_wp, -0.2753916589_wp, -0.0314683933_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.1246211137_wp, & !45-48
   &  0.0868343002_wp,  0.0304425141_wp, -0.0451897910_wp,  0.0279787730_wp, & !49-52
   &  0.0218755705_wp,  0.1578361412_wp, -0.2129755542_wp,  0.3173254159_wp, & !53-56
   & -0.0084749605_wp, -0.0627252130_wp, -0.0627252130_wp, -0.0583072067_wp, & !57-60
   & -0.0540641768_wp, -0.0499961233_wp, -0.0461030462_wp, -0.0388418210_wp, & !61-64
   & -0.0423849454_wp, -0.0388418210_wp, -0.0354736730_wp, -0.0322805013_wp, & !65-68
   & -0.0292623061_wp, -0.0264190872_wp, -0.0189392888_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp, -0.0005691876_wp,  0.0000000000_wp, & !77-80
   &  0.2010826437_wp,  0.1174172473_wp,  0.1253622556_wp,  0.1364259880_wp, & !81-84
   & -0.0638870618_wp,  0.0563222478_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp, -0.0237508447_wp, -0.0212575786_wp, -0.0189392888_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103

   !> Parameter: Charge-dependence of Zeff in GP3-xTB repulsion
   real(wp), parameter :: rep_q(max_elem) = [&
   &  0.0008271593_wp, -0.3522609502_wp,  0.1353380841_wp,  0.0450706434_wp, & !1-4
   &  0.0525063232_wp,  0.0007776371_wp,  0.0005586691_wp,  0.0007880504_wp, & !5-8
   &  0.0101403049_wp, -0.0852263855_wp,  0.3882216335_wp,  0.2065255534_wp, & !9-12
   &  0.1501898531_wp,  0.0116197630_wp,  0.0249841827_wp,  0.0076204783_wp, & !13-16
   &  0.0075788259_wp,  0.1054809329_wp,  0.4208321354_wp,  0.1771453922_wp, & !17-20
   &  0.1950989099_wp,  0.2305898199_wp,  0.1166685877_wp,  0.0685549853_wp, & !21-24
   &  0.1243005681_wp,  0.0541374963_wp,  0.1139522924_wp,  0.1139521451_wp, & !25-28
   &  0.1739825513_wp,  0.1240625271_wp,  0.1726382176_wp,  0.0329848041_wp, & !29-32
   &  0.0228209732_wp,  0.0060841344_wp,  0.0050488787_wp, -0.0080967279_wp, & !33-36
   &  0.4066541682_wp,  0.1468425593_wp,  0.0700000003_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.1736633295_wp, & !45-48
   &  0.2389184696_wp,  0.0387652478_wp,  0.1245409830_wp,  0.0608073820_wp, & !49-52
   &  0.0047456921_wp, -0.0592459850_wp,  0.2838070252_wp,  0.1515413690_wp, & !53-56
   &  0.0299999993_wp,  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp, & !57-60
   &  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp, & !61-64
   &  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp, & !65-68
   &  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp,  0.0399999991_wp,  0.0000000000_wp, & !77-80
   &  0.1759990858_wp,  0.0493036498_wp,  0.0574318580_wp,  0.0632894629_wp, & !81-84
   &  0.0457083121_wp,  0.0916276903_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp,  0.2399999946_wp,  0.2399999946_wp,  0.2399999946_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103

   !> Parameter: Radii in the coordination number for the GP3-xTB repulsion
   real(wp), parameter :: rep_rcov_cn(max_elem) = [&
   &  0.2806743503_wp, -0.2367094066_wp,  1.3199466509_wp,  1.0928398320_wp, & !1-4
   &  1.2163297359_wp,  1.2892723926_wp,  1.1093215479_wp,  0.9522993236_wp, & !5-8
   &  0.7601560354_wp,  0.9091963275_wp,  1.7896424369_wp,  1.5939753643_wp, & !9-12
   &  1.1783917691_wp,  1.1091054929_wp,  1.6808311996_wp,  1.6602686396_wp, & !13-16
   &  1.5106000154_wp,  1.3188860170_wp,  2.5074795871_wp,  1.7241829531_wp, & !17-20
   &  0.8382276375_wp,  1.8869086378_wp,  2.0725781810_wp,  1.7976184182_wp, & !21-24
   &  1.3711201783_wp,  1.2166611527_wp,  0.3954438357_wp,  1.3227586250_wp, & !25-28
   &  1.3969700852_wp,  1.3875299693_wp,  0.4895709548_wp,  0.4346276434_wp, & !29-32
   &  1.1973353138_wp,  1.4553451923_wp,  1.6411504403_wp,  1.9254514147_wp, & !33-36
   &  2.6671685687_wp,  1.2364298085_wp,  2.6456163286_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  1.6726585096_wp, & !45-48
   &  1.8468237307_wp,  1.2056707410_wp,  0.5585515604_wp,  1.4743880961_wp, & !49-52
   &  1.7199777178_wp,  1.4948856188_wp,  2.6175814767_wp,  1.4891100087_wp, & !53-56
   &  2.6456163286_wp,  3.0235615184_wp,  3.0235615184_wp,  2.9757992142_wp, & !57-60
   &  2.9345436587_wp,  2.8997948518_wp,  2.8715527937_wp,  2.8345889235_wp, & !61-64
   &  2.8498174843_wp,  2.8345889235_wp,  2.8258671114_wp,  2.8236520481_wp, & !65-68
   &  2.8279437334_wp,  2.8387421674_wp,  2.9101779615_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp,  2.1353903224_wp,  0.0000000000_wp, & !77-80
   &  1.9740882878_wp,  1.1052631807_wp,  1.5488265474_wp,  1.5810332912_wp, & !81-84
   &  1.5011537428_wp,  2.7393421954_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp,  2.8560473501_wp,  2.8798592814_wp,  2.9101779615_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103

   !> Parameter: Covalent radii for the GP3-xTB repulsion 
   real(wp), parameter :: cov_radii(max_elem) = [&
   &  0.9482131418_wp,  0.7407305377_wp,  1.9999396679_wp,  1.6432722774_wp, & !1-4
   &  1.5667588383_wp,  1.6038422957_wp,  1.5726980342_wp,  1.3534525238_wp, & !5-8
   &  1.1734308389_wp,  1.1154356439_wp,  2.2509996732_wp,  1.9082726166_wp, & !9-12
   &  2.2772331198_wp,  2.2975649691_wp,  2.4404242428_wp,  2.3116928364_wp, & !13-16
   &  2.0703156330_wp,  1.9486489403_wp,  3.1553530678_wp,  1.5524205976_wp, & !17-20
   &  1.7325493113_wp,  1.6634604100_wp,  1.6695156712_wp,  1.4302983432_wp, & !21-24
   &  1.6905509384_wp,  1.5587268194_wp,  1.6288892257_wp,  1.9567800332_wp, & !25-28
   &  1.9810161062_wp,  1.8728556342_wp,  1.8138599974_wp,  1.8359315936_wp, & !29-32
   &  1.4595635077_wp,  2.2630478355_wp,  1.9339004881_wp,  2.3424054535_wp, & !33-36
   &  3.1228854965_wp,  1.5409484429_wp,  2.0428639708_wp,  0.0000000000_wp, & !37-40
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !41-44
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  1.9670636903_wp, & !45-48
   &  2.0081601663_wp,  1.6434376839_wp,  1.4956837610_wp,  1.6405944999_wp, & !49-52
   &  2.2834897846_wp,  2.9505129131_wp,  3.5111280503_wp,  3.2947196497_wp, & !53-56
   &  2.4000000000_wp,  2.8477440674_wp,  2.8477440674_wp,  2.6824512091_wp, & !57-60
   &  2.5315619541_wp,  2.3950763026_wp,  2.2729942545_wp,  2.0720409684_wp, & !61-64
   &  2.1653158097_wp,  2.0720409684_wp,  1.9931697305_wp,  1.9287020960_wp, & !65-68
   &  1.8786380648_wp,  1.8429776371_wp,  1.8224179744_wp,  0.0000000000_wp, & !69-72
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !73-76
   &  0.0000000000_wp,  0.0000000000_wp,  2.0065590848_wp,  0.0000000000_wp, & !77-80
   &  2.4082419237_wp,  2.3308394398_wp,  1.7465130803_wp,  1.1185853327_wp, & !81-84
   &  1.6289727078_wp,  1.8723726163_wp,  0.0000000000_wp,  0.0000000000_wp, & !85-88
   &  0.0000000000_wp,  1.8217208128_wp,  1.8148675919_wp,  1.8224179744_wp, & !89-92
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !93-96
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !97-100
   &  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp] !101-103


   !> Specification of the GP3-xTB Hamiltonian
   type, public, extends(tb_h0spec) :: gp3_h0spec
      real(wp) :: kshell(0:3, 0:3)
      real(wp), allocatable :: kpair(:, :)
   contains
      !> Generator for the self energy / atomic levels of the Hamiltonian
      procedure :: get_selfenergy
      !> Generator for the coordination number dependent shift of the self energy
      procedure :: get_cnshift
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      procedure :: get_hscale
      !> Generator for the polynomial parameters for the distant dependent scaling
      procedure :: get_shpoly
      !> Generator for the reference occupation numbers of the atoms
      procedure :: get_reference_occ
   end type gp3_h0spec

   interface gp3_h0spec
      module procedure :: new_gp3_h0spec
   end interface gp3_h0spec

contains

subroutine new_gp3_calculator(calc, mol, error)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   ! Check if all atoms of mol%nat are supported (Z <= 103)
   if (any(mol%num > max_elem)) then
      call fatal_error(error, "No support for elements with Z >" // &
         & format_string(max_elem, '(i0)') // ".")
      return
   end if

   call add_basis(calc, mol)
   call add_ncoord(calc, mol)
   !call add_hamiltonian(calc, mol)
   call add_repulsion(calc, mol)
   !call add_dispersion(calc, mol)
   !call add_coulomb(calc, mol)

end subroutine new_gp3_calculator

subroutine add_basis(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   type(qvszp_basis_type), allocatable :: tmp
   type(qvszp_cgto_type), allocatable :: cgto(:, :)
   type(qvszp_cgto_type), allocatable :: cgto_scaled(:, :)
   integer, allocatable :: nsh_id(:)
   
   ! Setup qvszp basis
   call add_qvszp_cgtos(cgto, mol, nsh_id, .false.)
   ! Setup scaled qvszp basis for hamiltonian
   call add_qvszp_cgtos(cgto_scaled, mol, nsh_id, .false.) !, expscal())

   allocate(tmp)
   call new_basis(tmp, mol, nsh_id, cgto, 1.0_wp, cgto_scaled)
   call move_alloc(tmp, calc%bas)

end subroutine add_basis

subroutine add_ncoord(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call new_ncoord(calc%ncoord, mol, cn_type="erf")
end subroutine add_ncoord

subroutine add_hamiltonian(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call new_hamiltonian(calc%h0, mol, calc%bas, new_gp3_h0spec(mol))
end subroutine add_hamiltonian

subroutine add_dispersion(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   type(d4_dispersion), allocatable :: tmp

   allocate(tmp)
   call new_d4_dispersion(tmp, mol, s6=d4_s6, s8=d4_s8, a1=d4_a1, a2=d4_a2, s9=d4_s9)
   call move_alloc(tmp, calc%dispersion)
end subroutine add_dispersion

subroutine add_repulsion(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: alpha(:), zeff(:), rcov(:), kcn(:), kq(:), rcov_cn(:)

   zeff = rep_zeff(mol%num)
   alpha = rep_alpha(mol%num)
   kcn = rep_cn(mol%num)
   kq = rep_q(mol%num)
   rcov = rep_rcov(mol%num)
   rcov_cn = rep_rcov_cn(mol%num)
   call new_repulsion(calc%repulsion, mol, rep_type="gp3", zeff=zeff, alpha=alpha, &
      & rexp=rep_rexp, kcn=kcn, kq=kq, rcov_rep=rcov, rcov_cn=rcov_cn, exp_cn=rep_exp_cn)
end subroutine add_repulsion

subroutine add_coulomb(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: hardness(:, :), hubbard_derivs(:, :)
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)
   type(effective_coulomb), allocatable :: es2

   allocate(calc%coulomb)
   allocate(es2)
   call get_shell_hardness(mol, calc%bas, hardness)
   call new_effective_coulomb(es2, mol, gexp, hardness, arithmetic_average, &
      & calc%bas%nsh_id)
   call move_alloc(es2, calc%coulomb%es2)

   allocate(calc%coulomb%es3)
   call get_hubbard_derivs(mol, calc%bas, hubbard_derivs)
   call new_onsite_thirdorder(calc%coulomb%es3, mol, hubbard_derivs, calc%bas%nsh_id)

   !allocate(calc%coulomb%aes2)
   !dkernel = p_dkernel(mol%num)
   !qkernel = p_qkernel(mol%num)
   !rad = p_rad(mol%num)
   !vcn = p_vcn(mol%num)
   !!call new_damped_multipole(calc%coulomb%aes2, mol, kmp_ddamp, kmp_qdamp, dkernel, qkernel, &
   !   & mp_shift, mp_kexp, mp_rmax, rad, vcn)

end subroutine add_coulomb

subroutine get_shell_hardness(mol, bas, hardness)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Shell resolved hardness parameters
   real(wp), allocatable, intent(out) :: hardness(:, :)

   integer :: isp, izp, ish, il

   allocate(hardness(maxval(bas%nsh_id), mol%nid))
   hardness(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%ang
         hardness(ish, isp) = 1.0_wp !hubbard_parameter(izp) * shell_hubbard(il, izp)
      end do
   end do
end subroutine get_shell_hardness


pure function new_gp3_h0spec(mol) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Instance of the Hamiltonian specification
   type(gp3_h0spec) :: self

   integer :: il, jl

   allocate(self%kpair(mol%nid, mol%nid))
   self%kpair(:, :) = 1.0_wp

   do il = 0, 3
      do jl = 0, 3
         self%kshell(jl, il) = kshell(jl, il)
      end do
   end do

end function new_gp3_h0spec


!> Generator for the enhancement factor to for scaling Hamiltonian elements
subroutine get_hscale(self, mol, bas, hscale)
   !> Instance of the Hamiltonian specification
   class(gp3_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Scaling parameters for the Hamiltonian elements
   real(wp), intent(out) :: hscale(:, :, :, :)

   integer :: isp, jsp, izp, jzp, ish, jsh, il, jl
   real(wp) :: zi, zj, zij, den, enp, km

   hscale(:, :, :, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         den = (get_pauling_en(izp) - get_pauling_en(jzp))**2
         do ish = 1, bas%nsh_id(isp)
            il = bas%cgto(ish, isp)%ang
            do jsh = 1, bas%nsh_id(jsp)
               jl = bas%cgto(jsh, jsp)%ang
               !zi = slater_exponent(ish, izp)
               !zj = slater_exponent(jsh, jzp)
               !zij = (2*sqrt(zi*zj)/(zi+zj))!**wexp
               enp = 1.0_wp +  den !enscale *
               km = self%kpair(jsp, isp) * self%kshell(jl, il) * enp
               hscale(jsh, ish, jsp, isp) = 1.0_wp !zij * km
            end do
         end do
      end do
   end do
end subroutine get_hscale


!> Generator for the self energy / atomic levels of the Hamiltonian
subroutine get_selfenergy(self, mol, bas, selfenergy)
   !> Instance of the Hamiltonian specification
   class(gp3_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Self energy / atomic levels
   real(wp), intent(out) :: selfenergy(:, :)

   integer :: isp, izp, ish

   selfenergy(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         selfenergy(ish, isp) = 0.0_wp !p_selfenergy(ish, izp)
      end do
   end do
end subroutine get_selfenergy


!> Generator of the coordination number dependent shift of the self energy
subroutine get_cnshift(self, mol, bas, kcn)
   !> Instance of the Hamiltonian specification
   class(gp3_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Coordination number dependent shift
   real(wp), intent(out) :: kcn(:, :)

   integer :: isp, izp, ish

   kcn(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         kcn(ish, isp) = 0.0_wp !p_kcn(bas%cgto(ish, isp)%ang, izp)
      end do
   end do
end subroutine get_cnshift


!> Generator for the polynomial parameters for the distant dependent scaling
subroutine get_shpoly(self, mol, bas, shpoly)
   !> Instance of the Hamiltonian specification
   class(gp3_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Polynomial parameters for distant dependent scaleing
   real(wp), intent(out) :: shpoly(:, :)

   integer :: isp, izp, ish

   shpoly(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         shpoly(ish, isp) = 0.0_wp !p_shpoly(bas%cgto(ish, isp)%ang, izp)
      end do
   end do
end subroutine get_shpoly


subroutine get_reference_occ(self, mol, bas, refocc)
   !> Instance of the Hamiltonian specification
   class(gp3_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Reference occupation numbers
   real(wp), intent(out) :: refocc(:, :)

   integer :: isp, izp, ish

   refocc(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         refocc(ish, isp) = 0.0_wp !reference_occ(bas%cgto(ish, isp)%ang, izp)
      end do
   end do
end subroutine get_reference_occ


subroutine get_hubbard_derivs(mol, bas, hubbard_derivs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Shell resolved Hubbard derivatives
   real(wp), allocatable, intent(out) :: hubbard_derivs(:, :)

   integer :: isp, izp, ish, il

   allocate(hubbard_derivs(maxval(bas%nsh_id), mol%nid))
   hubbard_derivs(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%ang
         hubbard_derivs(ish, isp) = 0.0_wp !p_hubbard_derivs(izp) !* shell_hubbard_derivs(il)
      end do
   end do
end subroutine get_hubbard_derivs


subroutine export_gp3_param(param)
   type(param_record), intent(out) :: param

   integer :: izp, il, jl

   param%version = 1
   param%name = "GP3-xTB"
   param%reference = "tba"

   associate(par => param%hamiltonian)
      par%cn = "exp"
      !par%enscale = enscale
      !par%wexp = wexp
      par%lmax = 2
      do il = 0, 3
         do jl = 0, 3
            par%ksh(jl, il) = kshell(jl, il)
         end do
      end do
      par%kpol = 2.0_wp
      par%sym = to_symbol([(izp, izp=1, max_elem)])
      allocate(par%kpair(max_elem, max_elem), source=1.0_wp)
   end associate

   allocate(param%dispersion)
   associate(par => param%dispersion)
      par%s6 = d4_s6
      par%s8 = d4_s8
      par%a1 = d4_a1
      par%a2 = d4_a2
      par%s9 = d4_s9
      par%sc = .true.
      par%d3 = .false.
   end associate

   allocate(param%repulsion)
   associate(par => param%repulsion)
      par%rep_type = "gp3"
      par%exp_cn = rep_exp_cn
      par%rexp = rep_rexp
   end associate

   allocate(param%charge)
   associate(par => param%charge)
      par%kernel = coulomb_kernel%effective
      par%gexp = gexp
      par%average = "arithmetic"
   end associate

   allocate(param%thirdorder)
   associate(par => param%thirdorder)
      par%lmax = 2
      par%shell = .true.
      !par%ksh = shell_hubbard_derivs
   end associate

   allocate(param%multipole)
   associate(par => param%multipole)
      par%dmp3 = kmp_ddamp
      par%dmp5 = kmp_qdamp
      !par%kexp = mp_kexp
      !par%shift = mp_shift
      !par%rmax = mp_rmax
   end associate

   allocate(param%record(max_elem))
   do izp = 1, max_elem
      associate(par => param%record(izp))
         par%sym = to_symbol(izp)
         par%num = izp

         par%zeff = rep_zeff(izp)
         par%alpha = rep_alpha(izp)

         par%nsh = 1 !nshell(izp)
         allocate(par%lsh(par%nsh))
         allocate(par%pqn(par%nsh))
         allocate(par%ngauss(par%nsh))
         allocate(par%levels(par%nsh))
         allocate(par%slater(par%nsh))
         allocate(par%refocc(par%nsh))
         allocate(par%kcn(par%nsh))
         allocate(par%shpoly(par%nsh))
         allocate(par%lgam(par%nsh))

         !par%lsh = ang_shell(:par%nsh, izp)
         !par%pqn = principal_quantum_number(:par%nsh, izp)
         !par%ngauss = number_of_primitives(:par%nsh, izp)
         !par%levels = p_selfenergy(:par%nsh, izp)
         !par%slater = slater_exponent(:par%nsh, izp)
         !par%refocc = reference_occ(par%lsh, izp)
         !par%kcn = p_kcn(par%lsh, izp)
         !par%shpoly = p_shpoly(par%lsh, izp)

         !par%gam = hubbard_parameter(izp)
         !par%lgam = shell_hubbard(par%lsh, izp)

         !par%gam3 = p_hubbard_derivs(izp)

         !par%mpvcn = p_vcn(izp)
         !par%mprad = p_rad(izp)
         !par%dkernel = p_dkernel(izp)
         !par%qkernel = p_qkernel(izp)

         !par%en = get_pauling_en(izp)
      end associate
   end do
end subroutine export_gp3_param

elemental function kshell(k, l)
   integer, intent(in) :: k, l
   real(wp) :: kshell
   kshell = 2.0_wp !merge(2.0_wp, (kdiag(l)+kdiag(k))/2, k==2.and.any(l==[0,1]).or.l==2.and.any(k==[0,1]))
end function kshell

end module tblite_xtb_gp3
