
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/DynamicBuffer.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube {

void puncture_field_acc_0(
    gsl::not_null<Variables<tmpl::list<
        CurvedScalarWave::Tags::Psi, ::Tags::dt<CurvedScalarWave::Tags::Psi>,
        ::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                      Frame::Inertial>>> *>
        result,
    const tnsr::I<DataVector, 3, Frame::Inertial> &centered_coords,
    const tnsr::I<double, 3> &particle_position,
    const tnsr::I<double, 3> &particle_velocity,
    const tnsr::I<double, 3> &particle_acceleration, const double BH_mass) {
  const size_t grid_size = get<0>(centered_coords).size();
  result->initialize(grid_size);
  const double xp = particle_position[0];
  const double yp = particle_position[1];
  const double xpdot = particle_velocity[0];
  const double ypdot = particle_velocity[1];
  const double xpddot = particle_acceleration[0];
  const double ypddot = particle_acceleration[1];
  const double rp = get(magnitude(particle_position));
  const double rpdot = (xp * xpdot + yp * ypdot) / rp;

  const auto &Dx = get<0>(centered_coords);
  const auto &Dy = get<1>(centered_coords);
  const auto &z = get<2>(centered_coords);

  const double M = BH_mass;

  DynamicBuffer<DataVector> temps(57, grid_size);

  const double d_0 = a * a;
  const double d_1 = rp * rp;
  const double d_2 = d_0 + d_1;
  const double d_3 = 1.0 / d_2;
  const double d_4 = a * xp - rp * yp;
  const double d_5 = -d_4;
  const double d_6 = d_3 * d_5;
  const double d_7 = d_6 * ypdot;
  const double d_8 = a * yp;
  const double d_9 = rp * xp;
  const double d_10 = d_8 + d_9;
  const double d_11 = d_10 * d_3;
  const double d_12 = d_11 * xpdot;
  const double d_13 = 1.0 / rp;
  const double d_14 = d_13 * zp;
  const double d_15 = d_14 * zpdot;
  const double d_16 = d_15 + 1;
  const double d_17 = d_12 + d_16;
  const double d_18 = d_17 + d_7;
  const double d_19 = d_18 * d_18;
  const double d_20 = rp * rp * rp;
  const double d_21 = rp * rp * rp * rp;
  const double d_22 = zp * zp;
  const double d_23 = d_0 * d_22;
  const double d_24 = d_21 + d_23;
  const double d_25 = 1.0 / d_24;
  const double d_26 = M * d_25;
  const double d_27 = 2 * d_26;
  const double d_28 = d_20 * d_27;
  const double d_29 = -1 + xpdot * xpdot + ypdot * ypdot + zpdot * zpdot;
  const double d_30 = d_19 * d_28 + d_29;
  const double d_31 = 1.0 / d_30;
  const double d_32 = Dz * ypdot;
  const double d_33 = d_3 * d_4;
  const double d_34 = 2 * d_33;
  const double d_35 = d_20 * d_26;
  const double d_36 = d_34 * d_35;
  const double d_37 = Dz * xpdot;
  const double d_38 = 2 * d_11;
  const double d_39 = d_35 * d_38;
  const double d_40 = d_4 * d_4;
  const double d_41 = 1.0 / (d_2 * d_2);
  const double d_42 = d_28 * d_41;
  const double d_43 = d_40 * d_42;
  const double d_44 = Dz * zp;
  const double d_45 = d_1 * d_26;
  const double d_46 = 2 * d_45;
  const double d_47 = d_44 * d_46;
  const double d_48 = Dz * zpdot;
  const double d_49 = 2 * rp;
  const double d_50 = d_26 * d_49;
  const double d_51 = d_22 * d_50;
  const double d_52 = d_10 * d_10;
  const double d_53 = d_42 * d_52;
  const double d_54 = 4 * d_35;
  const double d_55 = d_44 * d_45;
  const double d_56 = 4 * d_55;
  const double d_57 = Dz * Dz;
  const double d_58 = d_22 * d_57;
  const double d_59 = Dz * d_14;
  const double d_60 = M_SQRT2;
  const double d_61 = xp * xp;
  const double d_62 = yp * yp;
  const double d_63 = d_22 + d_61 + d_62;
  const double d_64 = -d_0 + d_63;
  const double d_65 = sqrt(4 * d_23 + d_64 * d_64);
  const double d_66 = d_0 + d_63 + d_65;
  const double d_67 = d_60 * d_66;
  const double d_68 = sqrt(d_64 + d_65);
  const double d_69 =
      -3 * d_24 * d_60 * d_66 + 4 * rp * (d_0 * d_65 * d_68 + d_20 * d_67);
  const double d_70 = 1.0 / d_68;
  const double d_71 = d_70 * zp;
  const double d_72 = Dz * d_71;
  const double d_73 = 3 * d_23;
  const double d_74 = d_60 * d_68;
  const double d_75 = d_74 * (-d_21 + d_73);
  const double d_76 = d_75 * xp;
  const double d_77 = d_75 * yp;
  const double d_78 = 1.0 / d_65;
  const double d_79 = d_25 * d_78;
  const double d_80 = -d_38;
  const double d_81 = 1.0 / d_1;
  const double d_82 = d_68 * d_81;
  const double d_83 = d_82 * xp;
  const double d_84 = d_60 * d_78;
  const double d_85 = d_83 * d_84;
  const double d_86 = d_44 * d_85;
  const double d_87 = d_74 * d_78;
  const double d_88 = 2 * d_87;
  const double d_89 = d_88 * d_9;
  const double d_90 = d_10 * d_89 - d_2 * (d_49 + d_61 * d_87);
  const double d_91 = d_4 * d_89;
  const double d_92 = 2 * a;
  const double d_93 = -d_60 * d_68 * d_78 * xp * yp + d_92;
  const double d_94 = d_49 * d_8 - xp * (d_0 - d_1);
  const double d_95 = d_67 * d_78;
  const double d_96 = d_71 * d_95;
  const double d_97 = 2 * d_14;
  const double d_98 = -d_13 * d_22 * d_70 * d_95 + 2;
  const double d_99 = d_0 * yp - d_1 * yp + d_9 * d_92;
  const double d_100 = d_82 * yp;
  const double d_101 = d_100 * d_84;
  const double d_102 = d_101 * d_44;
  const double d_103 = d_88 * rp * yp;
  const double d_104 = d_87 * xp;
  const double d_105 = d_10 * d_103 - d_2 * (d_104 * yp + d_92);
  const double d_106 = d_103 * d_5 - d_2 * (d_49 + d_62 * d_87);
  const double d_107 = -d_33 * ypdot;
  const double d_108 = d_107 + d_16;
  const double d_109 = d_108 + d_12;
  const double d_110 = -d_69 * d_71;
  const double d_111 = -d_94;
  const double d_112 = -d_90;
  const double d_113 = -d_2 * d_93 + d_91;
  const double d_114 = d_113 * d_41;
  const double d_115 = 2 * d_6;
  const double d_116 = -d_106;
  const double d_117 = -d_105;
  const double d_118 = d_110 * zpdot + d_76 * xpdot + d_77 * ypdot;
  const double d_119 = d_101 * zp;
  const double d_120 = d_111 * d_41;
  const double d_121 = d_120 * d_96;
  const double d_122 = d_41 * d_99;
  const double d_123 = d_122 * d_96;
  const double d_124 = d_117 * d_41;
  const double d_125 = rp * zpdot - rpdot * zp;
  const double d_126 = d_125 * zp;
  const double d_127 = zp * zpdot;
  const double d_128 = M * 1.0 / (d_24 * d_24);
  const double d_129 =
      d_128 * (d_0 * d_127 * d_49 + d_21 * rpdot - d_73 * rpdot);
  const double d_130 = Dz * d_125 * d_50;
  const double d_131 = d_129 * d_44 * d_49;
  const double d_132 = d_1 * d_129;
  const double d_133 = d_132 * d_41;
  const double d_134 = d_49 * rpdot;
  const double d_135 =
      d_134 * d_5 - d_2 * (-a * xpdot + rp * ypdot + rpdot * yp);
  const double d_136 = 1.0 / (d_2 * d_2 * d_2);
  const double d_137 = d_136 * d_28;
  const double d_138 =
      d_10 * d_134 - d_2 * (a * ypdot + rp * xpdot + rpdot * xp);
  const double d_139 = d_10 * d_138;
  const double d_140 = d_135 * d_4;
  const double d_141 = d_45 * zp;
  const double d_142 = 2 * xpdot;
  const double d_143 = d_35 * d_41;
  const double d_144 = d_142 * d_143 * d_52 + xpdot;
  const double d_145 = d_34 * d_45;
  const double d_146 = 2 * ypdot;
  const double d_147 = d_143 * d_146;
  const double d_148 = Dz * zpddot;
  const double d_149 = 2 * d_129;
  const double d_150 = d_1 * d_149 * d_41;
  const double d_151 = d_136 * d_54;
  const double d_152 = d_135 * d_41;
  const double d_153 = d_96 * zpdot;
  const double d_154 = d_117 * ypdot;
  const double d_155 = d_41 * xpdot;
  const double d_156 = d_113 * xpdot;
  const double d_157 = d_18 * d_28;
  const double d_158 = -1.0 / d_30;
  const double d_159 = d_115 * d_35;
  const double d_160 = d_5 * d_5;
  const double d_161 = d_160 * d_42;
  const double d_162 = d_72 * d_95;
  const double d_163 = d_66 * d_70;
  const double d_164 = d_127 * d_84;
  const double d_165 = d_18 * rp;
  DataVector &dv_0 = temps.at(0);
  dv_0 = Dy * zpdot;
  DataVector &dv_1 = temps.at(1);
  dv_1 = Dy + d_14 * d_32 + d_14 * dv_0;
  DataVector &dv_2 = temps.at(2);
  dv_2 = Dx * ypdot;
  DataVector &dv_3 = temps.at(3);
  dv_3 = Dy * xpdot;
  DataVector &dv_4 = temps.at(4);
  dv_4 = dv_2 + dv_3;
  DataVector &dv_5 = temps.at(5);
  dv_5 = Dx - d_33 * dv_4;
  DataVector &dv_6 = temps.at(6);
  dv_6 = Dx * zpdot;
  DataVector &dv_7 = temps.at(7);
  dv_7 = d_14 * d_37 + d_14 * dv_6;
  DataVector &dv_8 = temps.at(8);
  dv_8 = Dy * ypdot;
  DataVector &dv_9 = temps.at(9);
  dv_9 = Dx * xpdot;
  DataVector &dv_10 = temps.at(10);
  dv_10 = d_48 + dv_8 + dv_9;
  DataVector &dv_11 = temps.at(11);
  dv_11 = d_47 + d_48 * d_51 + d_53 * dv_9 + dv_10;
  DataVector &dv_12 = temps.at(12);
  dv_12 = -d_36 * dv_1 + d_39 * (dv_5 + dv_7) + d_43 * dv_8 + dv_11;
  DataVector &dv_13 = temps.at(13);
  dv_13 = dv_12 * dv_12;
  DataVector &dv_14 = temps.at(14);
  dv_14 = Dy * d_33;
  DataVector &dv_15 = temps.at(15);
  dv_15 = -Dz * d_13 * zp + dv_14;
  DataVector &dv_16 = temps.at(16);
  dv_16 = -dv_15;
  DataVector &dv_17 = temps.at(17);
  dv_17 = Dx * d_11;
  DataVector &dv_18 = temps.at(18);
  dv_18 = d_54 * dv_17;
  DataVector &dv_19 = temps.at(19);
  dv_19 = Dy * Dy;
  DataVector &dv_20 = temps.at(20);
  dv_20 = d_40 * dv_19;
  DataVector &dv_21 = temps.at(21);
  dv_21 = Dx * Dx;
  DataVector &dv_22 = temps.at(22);
  dv_22 = d_52 * dv_21;
  DataVector &dv_23 = temps.at(23);
  dv_23 = d_42 * dv_22 + d_50 * d_58 + d_57 + dv_19 + dv_21;
  DataVector &dv_24 = temps.at(24);
  dv_24 = d_42 * dv_20 - d_56 * dv_14 + dv_23;
  DataVector &dv_25 = temps.at(25);
  dv_25 = -d_31 * dv_13 + dv_16 * dv_18 + dv_24;
  DataVector &dv_26 = temps.at(26);
  dv_26 = pow(dv_25, -1.0 / 2.0);
  DataVector &dv_27 = temps.at(27);
  dv_27 = d_59 - dv_14 + dv_17;
  DataVector &dv_28 = temps.at(28);
  dv_28 = Dx * d_76 + Dy * d_77;
  DataVector &dv_29 = temps.at(29);
  dv_29 = Dx * d_41;
  DataVector &dv_30 = temps.at(30);
  dv_30 = d_96 * dv_29;
  DataVector &dv_31 = temps.at(31);
  dv_31 = Dy * d_41;
  DataVector &dv_32 = temps.at(32);
  dv_32 = Dz * d_13 * d_98 + d_96 * d_99 * dv_31 + d_97;
  DataVector &dv_33 = temps.at(33);
  dv_33 = 2 * dv_14;
  DataVector &dv_34 = temps.at(34);
  dv_34 = 2 * dv_17;
  DataVector &dv_35 = temps.at(35);
  dv_35 = 2 * d_59 + dv_34;
  DataVector &dv_36 = temps.at(36);
  dv_36 = Dz * d_110 + dv_28;
  DataVector &dv_37 = temps.at(37);
  dv_37 = Dy * d_6;
  DataVector &dv_38 = temps.at(38);
  dv_38 = d_59 + dv_37;
  DataVector &dv_39 = temps.at(17);
  dv_39 = dv_17 + dv_38;
  DataVector &dv_40 = temps.at(39);
  dv_40 = d_79 * dv_39;
  DataVector &dv_41 = temps.at(40);
  dv_41 = dv_36 * dv_40;
  DataVector &dv_42 = temps.at(41);
  dv_42 = d_112 * dv_29;
  DataVector &dv_43 = temps.at(42);
  dv_43 = Dy * d_114 - d_86;
  DataVector &dv_44 = temps.at(43);
  dv_44 = d_116 * dv_31;
  DataVector &dv_45 = temps.at(44);
  dv_45 = -d_102 + d_117 * dv_29;
  DataVector &dv_46 = temps.at(45);
  dv_46 = Dx * (d_38 + dv_42 + dv_43) + Dy * (d_115 + dv_44 + dv_45) +
          Dz * (d_111 * dv_30 + dv_32);
  DataVector &dv_47 = temps.at(46);
  dv_47 = d_114 * (dv_2 - dv_3) + d_119 * d_32 - d_119 * dv_0 + 2 * d_12 +
          d_121 * d_37 - d_121 * dv_6 + d_123 * d_32 - d_123 * dv_0 -
          d_124 * dv_2 + d_124 * dv_3 + 2 * d_15 + 2 * d_7 -
          d_85 * zp * (-d_37 + dv_6) + 2;
  DataVector &dv_48 = temps.at(47);
  dv_48 = d_18 * d_49 * dv_46 + 2 * d_18 * dv_41 -
          dv_39 * (d_118 * d_25 * d_78 * dv_39 - d_49 * dv_47);
  DataVector &dv_49 = temps.at(15);
  dv_49 =
      -dv_13 * 1.0 / (d_28 * (d_109 * d_109) + d_29) - dv_15 * dv_18 + dv_24;
  DataVector &dv_50 = temps.at(24);
  dv_50 = d_135 * dv_31;
  DataVector &dv_51 = temps.at(48);
  dv_51 = d_28 * dv_29;
  DataVector &dv_52 = temps.at(12);
  dv_52 = d_31 * dv_12;
  DataVector &dv_53 = temps.at(49);
  dv_53 = Dx + d_52 * dv_51;
  DataVector &dv_54 = temps.at(50);
  dv_54 = d_39 * dv_16 - dv_52 * (d_108 * d_39 + d_144) + dv_53;
  DataVector &dv_55 = temps.at(51);
  dv_55 = d_10 * dv_51;
  DataVector &dv_56 = temps.at(52);
  dv_56 = 2 * Dy * M * d_20 * d_25 * d_40 * d_41 + Dy - d_145 * d_44 -
          d_4 * dv_55 - dv_52 * (d_147 * d_40 - d_17 * d_36 + ypdot);
  DataVector &dv_57 = temps.at(53);
  dv_57 = Dx * xpddot;
  DataVector &dv_58 = temps.at(54);
  dv_58 = Dy * ypddot;
  DataVector &dv_59 = temps.at(0);
  dv_59 = d_32 + dv_0;
  DataVector &dv_60 = temps.at(55);
  dv_60 = Dy + d_14 * dv_59;
  DataVector &dv_61 = temps.at(5);
  dv_61 = d_14 * (d_37 + dv_6) + dv_5;
  DataVector &dv_62 = temps.at(11);
  dv_62 = d_159 * dv_1 + d_161 * dv_8 + d_39 * (Dx + d_6 * dv_4 + dv_7) + dv_11;
  DataVector &dv_63 = temps.at(23);
  dv_63 = sqrt(d_158 * (dv_62 * dv_62) + d_161 * dv_19 + d_56 * dv_37 +
               dv_18 * dv_38 + dv_23);
  DataVector &dv_64 = temps.at(18);
  dv_64 = dv_39 * dv_39;
  DataVector &dv_65 = temps.at(36);
  dv_65 = d_79 * dv_36;
  DataVector &dv_66 = temps.at(37);
  dv_66 = dv_35 + 2 * dv_37 + dv_46;
  DataVector &dv_67 = temps.at(45);
  dv_67 = d_49 * dv_39;
  DataVector &dv_68 = temps.at(1);
  dv_68 = d_157 * dv_39 + dv_10;
  DataVector &dv_69 = temps.at(7);
  dv_69 = -d_158 * dv_48;
  DataVector &dv_70 = temps.at(4);
  dv_70 = (dv_64 * dv_65 + dv_66 * dv_67 - dv_68 * dv_69) * 1.0 / dv_63;
  DataVector &dv_71 = temps.at(6);
  dv_71 = 2 * d_45 * dv_70 - 4 * dv_63;
  DataVector &dv_72 = temps.at(11);
  dv_72 = d_158 * dv_62;
  DataVector &dv_73 = temps.at(18);
  dv_73 = d_79 * dv_64;
  DataVector &dv_74 = temps.at(42);
  dv_74 = d_117 * dv_31 + d_120 * d_162 + 2 * dv_42 + dv_43;
  DataVector &dv_75 = temps.at(37);
  dv_75 = d_49 * dv_66;
  DataVector &dv_76 = temps.at(41);
  dv_76 = d_18 * dv_40;
  DataVector &dv_77 = temps.at(39);
  dv_77 = d_118 * dv_40;
  DataVector &dv_78 = temps.at(36);
  dv_78 = d_18 * dv_65;
  DataVector &dv_79 = temps.at(17);
  dv_79 = dv_39 * rp;
  DataVector &dv_80 = temps.at(46);
  dv_80 = dv_47 * rp;
  DataVector &dv_81 = temps.at(1);
  dv_81 = 2 * d_158 * dv_68;
  DataVector &dv_82 = temps.at(56);
  dv_82 = (1.0 / 4.0) * dv_26 / pow(dv_49, 3.0 / 2.0);
  DataVector &dv_83 = temps.at(44);
  dv_83 = Dx * d_114 + d_122 * d_162 + 2 * dv_44 + dv_45;

  get(get<CurvedScalarWave::Tags::Psi>(*result)) =
      dv_26 * (-1.0 / 4.0 * d_45 * 1.0 / dv_25 *
                   (-d_31 * dv_48 * (d_109 * d_28 * dv_27 + dv_10) +
                    d_49 * dv_27 *
                        (Dx * (Dy * d_41 * (-d_2 * d_93 + d_91) - d_80 - d_86 -
                               d_90 * dv_29) -
                         Dy * (d_102 + d_105 * dv_29 + d_106 * dv_31 + d_34) +
                         Dz * (-d_94 * dv_30 + dv_32) - dv_33 + dv_35) +
                    d_79 * (-d_69 * d_72 + dv_28) * (dv_27 * dv_27)) +
               1);
  get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*result)) =
      dv_26 * 1.0 / dv_49 *
      (-d_126 * d_27 * d_57 + d_129 * d_58 + d_130 * dv_14 - d_131 * dv_14 +
       d_132 * dv_16 * dv_34 + d_133 * dv_20 + d_133 * dv_22 +
       d_137 * d_139 * dv_21 - d_137 * d_140 * dv_19 + d_138 * dv_16 * dv_51 -
       d_35 * dv_34 * (Dz * d_125 * d_81 - dv_50) + d_47 * dv_50 -
       1.0 / 2.0 * dv_13 *
           (d_1 * d_118 * d_128 * d_19 * d_78 + d_142 * xpddot +
            d_146 * ypddot +
            d_157 * (d_115 * ypddot +
                     d_13 * zpdot *
                         (-d_104 * d_14 * xpdot - d_14 * d_87 * yp * ypdot +
                          d_98 * zpdot) +
                     d_155 * (d_111 * d_153 + d_112 * xpdot + d_154) +
                     d_38 * xpddot +
                     d_41 * ypdot * (d_116 * ypdot + d_153 * d_99 + d_156) +
                     d_97 * zpddot) +
            2 * zpddot * zpdot) *
           1.0 / (d_30 * d_30) +
       dv_52 *
           (4 * d_126 * d_26 * d_48 + d_130 - d_131 + d_132 * d_34 * dv_60 -
            d_132 * d_38 * dv_61 - d_138 * d_42 * dv_61 - d_139 * d_151 * dv_9 +
            d_140 * d_151 * dv_8 -
            d_145 * (d_125 * d_13 * dv_59 + zp * (Dy * zpddot + Dz * ypddot)) +
            d_148 * d_51 + d_148 - d_149 * d_22 * d_48 - d_150 * d_40 * dv_8 -
            d_150 * d_52 * dv_9 - d_152 * d_28 * dv_60 +
            d_39 * (Dx * d_125 * d_81 * zpdot + Dz * d_125 * d_81 * xpdot +
                    d_13 * zp * (Dx * zpddot + Dz * xpddot) - d_152 * dv_2 -
                    d_152 * dv_3 - d_33 * (Dx * ypddot + Dy * xpddot)) +
            d_43 * dv_58 + d_53 * dv_57 + dv_57 + dv_58) +
       dv_54 * xpdot + dv_56 * ypdot +
       zpdot *
           (Dz * d_51 + Dz - d_141 * dv_33 + d_141 * dv_34 -
            dv_52 * (d_46 * zp * (d_107 + d_12 + 1) + d_51 * zpdot + zpdot)));
  get<0>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) =
      dv_82 *
      (d_45 * (-dv_63 * (d_11 * dv_75 + d_38 * dv_41 + d_76 * dv_73 +
                         dv_67 * dv_74 - dv_69 * (d_18 * d_39 + xpdot) +
                         dv_81 * (-d_11 * dv_77 + d_11 * dv_78 + d_11 * dv_80 +
                                  d_165 * (d_80 + dv_74) + d_76 * dv_76 +
                                  dv_79 * (d_113 * d_41 * ypdot - d_154 * d_41 -
                                           d_164 * (d_120 * d_163 + d_83)))) +
               dv_70 * (d_39 * dv_38 + dv_53 +
                        dv_72 * (d_144 + d_39 * (d_16 + d_7)))) +
       dv_54 * dv_71);
  get<1>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) =
      dv_82 *
      (d_45 * (-dv_63 * (d_115 * dv_41 + d_6 * dv_75 + d_77 * dv_73 +
                         dv_67 * dv_83 - dv_69 * (d_159 * d_18 + ypdot) +
                         dv_81 * (d_165 * (-d_115 + dv_83) - d_6 * dv_77 +
                                  d_6 * dv_78 + d_6 * dv_80 + d_77 * dv_76 -
                                  dv_79 * (-d_117 * d_155 + d_156 * d_41 +
                                           d_164 * (d_100 + d_122 * d_163)))) +
               dv_70 * (Dy + d_115 * d_55 + d_160 * d_28 * dv_31 + d_5 * dv_55 +
                        dv_72 * (d_147 * d_160 + d_159 * d_17 + ypdot))) +
       dv_56 * dv_71);
  get<2>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) = dzPsiP0;
}
} // namespace CurvedScalarWave::Worldtube
