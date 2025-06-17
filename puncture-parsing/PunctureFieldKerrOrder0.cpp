
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

namespace CurvedScalarWave::Worldtube
{

    void puncture_field_acc_0(
        gsl::not_null<Variables<tmpl::list<
            CurvedScalarWave::Tags::Psi, ::Tags::dt<CurvedScalarWave::Tags::Psi>,
            ::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                          Frame::Inertial>>> *>
            result,
        const tnsr::I<DataVector, 3, Frame::Inertial> &centered_coords,
        const tnsr::I<double, 3> &particle_position,
        const tnsr::I<double, 3> &particle_velocity,
        const tnsr::I<double, 3> &particle_acceleration, const double BH_mass,
        const std::array<double, 3> &BH_spin)
    {
        const size_t grid_size = get<0>(centered_coords).size();
        result->initialize(grid_size);
        const double xp = particle_position[0];
        const double yp = particle_position[1];
        const double zp = particle_position[2];
        const double xpdot = particle_velocity[0];
        const double ypdot = particle_velocity[1];
        const double zpdot = particle_velocity[2];

        const double xpddot = particle_acceleration[0];
        const double ypddot = particle_acceleration[1];
        const double zpddot = particle_acceleration[2];

        const double rp = get(magnitude(particle_position));
        const double rpdot = (xp * xpdot + yp * ypdot + zp * zpdot) / rp;

        const auto &Dx = get<0>(centered_coords);
        const auto &Dy = get<1>(centered_coords);
        const auto &Dz = get<2>(centered_coords);

        const double M = BH_mass;
        const double a = BH_spin[2];

        DynamicBuffer<DataVector> temps(38, grid_size);

        const double d_0 = a * a;
        const double d_1 = rp * rp;
        const double d_2 = d_0 + d_1;
        const double d_3 = 1.0 / d_2;
        const double d_4 = a * yp;
        const double d_5 = rp * xp;
        const double d_6 = d_4 + d_5;
        const double d_7 = d_3 * d_6;
        const double d_8 = d_7 * xpdot;
        const double d_9 = rp * yp;
        const double d_10 = a * xp - d_9;
        const double d_11 = d_10 * d_3;
        const double d_12 = -d_11 * ypdot;
        const double d_13 = 1.0 / rp;
        const double d_14 = d_13 * zp;
        const double d_15 = d_14 * zpdot + 1;
        const double d_16 = d_12 + d_15;
        const double d_17 = rp * rp * rp;
        const double d_18 = rp * rp * rp * rp;
        const double d_19 = zp * zp;
        const double d_20 = d_0 * d_19;
        const double d_21 = d_18 + d_20;
        const double d_22 = 1.0 / d_21;
        const double d_23 = M * d_22;
        const double d_24 = 2 * d_23;
        const double d_25 = d_17 * d_24;
        const double d_26 = -1 + xpdot * xpdot + ypdot * ypdot + zpdot * zpdot;
        const double d_27 = 1.0 / (d_25 * pow(d_16 + d_8, 2) + d_26);
        const double d_28 = d_11 * d_25;
        const double d_29 = d_25 * d_7;
        const double d_30 = d_24 * zp;
        const double d_31 = d_1 * d_30;
        const double d_32 = 2 * rp;
        const double d_33 = d_23 * d_32;
        const double d_34 = d_6 * d_6;
        const double d_35 = 1.0 / (d_2 * d_2);
        const double d_36 = d_25 * d_35;
        const double d_37 = d_34 * d_36;
        const double d_38 = d_10 * d_10;
        const double d_39 = d_36 * d_38;
        const double d_40 = 4 * d_23;
        const double d_41 = d_17 * d_40;
        const double d_42 = -d_10;
        const double d_43 = d_3 * d_42;
        const double d_44 = d_15 + d_8;
        const double d_45 = d_43 * ypdot + d_44;
        const double d_46 = d_45 * d_45;
        const double d_47 = d_25 * d_46 + d_26;
        const double d_48 = 1.0 / d_47;
        const double d_49 = rp * zpdot;
        const double d_50 = d_49 - rpdot * zp;
        const double d_51 = 3 * d_20;
        const double d_52 = M * 1.0 / (d_21 * d_21);
        const double d_53 = d_52 * (2 * d_0 * d_49 * zp + d_18 * rpdot - d_51 * rpdot);
        const double d_54 = d_1 * d_35;
        const double d_55 = d_53 * d_54;
        const double d_56 = d_32 * rpdot;
        const double d_57 = -d_2 * (-a * xpdot + rp * ypdot + rpdot * yp) + d_42 * d_56;
        const double d_58 = d_35 * d_57;
        const double d_59 = 1.0 / (d_2 * d_2 * d_2);
        const double d_60 = d_25 * d_59;
        const double d_61 = -d_2 * (a * ypdot + rp * xpdot + rpdot * xp) + d_56 * d_6;
        const double d_62 = d_6 * d_61;
        const double d_63 = d_10 * d_57;
        const double d_64 = 2 * d_53;
        const double d_65 = d_1 * d_64;
        const double d_66 = d_36 * d_61;
        const double d_67 = 1.0 / d_1;
        const double d_68 = d_19 * d_24 * d_49 + d_31 * (d_12 + d_8 + 1) + zpdot;
        const double d_69 = d_35 * xpdot;
        const double d_70 = d_16 * d_29 + d_25 * d_34 * d_69 + xpdot;
        const double d_71 = d_35 * ypdot;
        const double d_72 = d_40 * zp;
        const double d_73 = d_54 * d_64;
        const double d_74 = d_1 * d_11;
        const double d_75 = d_41 * d_59;
        const double d_76 = 2 * xpddot;
        const double d_77 = 2 * ypddot;
        const double d_78 = 2 * zpddot;
        const double d_79 = xp * xp;
        const double d_80 = yp * yp;
        const double d_81 = d_19 + d_79 + d_80;
        const double d_82 = -d_0 + d_81;
        const double d_83 = sqrt(4 * d_20 + d_82 * d_82);
        const double d_84 = 1.0 / d_83;
        const double d_85 = M_SQRT2;
        const double d_86 = sqrt(d_82 + d_83);
        const double d_87 = d_85 * d_86;
        const double d_88 = d_87 * (-d_18 + d_51);
        const double d_89 = yp * ypdot;
        const double d_90 = d_0 + d_81 + d_83;
        const double d_91 = d_85 * d_90;
        const double d_92 = 4 * rp;
        const double d_93 = 1.0 / d_86;
        const double d_94 = d_93 * zp * zpdot;
        const double d_95 = d_84 * d_87;
        const double d_96 = d_95 * xp;
        const double d_97 = d_84 * d_91;
        const double d_98 = d_94 * d_97;
        const double d_99 = 2 * d_95;
        const double d_100 = d_6 * d_99;
        const double d_101 = 2 * a;
        const double d_102 = d_1 * d_72;
        DataVector &dv_0 = temps.at(0);
        dv_0 = Dy * zpdot;
        DataVector &dv_1 = temps.at(1);
        dv_1 = Dz * ypdot;
        DataVector &dv_2 = temps.at(2);
        dv_2 = d_14 * dv_1;
        DataVector &dv_3 = temps.at(3);
        dv_3 = Dx * zpdot;
        DataVector &dv_4 = temps.at(4);
        dv_4 = Dz * xpdot;
        DataVector &dv_5 = temps.at(5);
        dv_5 = Dx * ypdot;
        DataVector &dv_6 = temps.at(6);
        dv_6 = Dy * xpdot;
        DataVector &dv_7 = temps.at(7);
        dv_7 = Dx - d_11 * (dv_5 + dv_6);
        DataVector &dv_8 = temps.at(8);
        dv_8 = Dx * xpdot;
        DataVector &dv_9 = temps.at(9);
        dv_9 = Dy * ypdot;
        DataVector &dv_10 = temps.at(10);
        dv_10 = Dz * zpdot;
        DataVector &dv_11 = temps.at(11);
        dv_11 = Dz * d_31;
        DataVector &dv_12 = temps.at(12);
        dv_12 = d_19 * dv_10;
        DataVector &dv_13 = temps.at(13);
        dv_13 = d_33 * dv_12 + d_37 * dv_8 + d_39 * dv_9 + dv_10 + dv_11 + dv_8 + dv_9;
        DataVector &dv_14 = temps.at(14);
        dv_14 = -d_28 * (Dy + d_14 * dv_0 + dv_2) + d_29 * (d_14 * dv_3 + d_14 * dv_4 + dv_7) + dv_13;
        DataVector &dv_15 = temps.at(15);
        dv_15 = dv_14 * dv_14;
        DataVector &dv_16 = temps.at(16);
        dv_16 = Dy * d_11;
        DataVector &dv_17 = temps.at(17);
        dv_17 = -Dz * d_14 + dv_16;
        DataVector &dv_18 = temps.at(18);
        dv_18 = Dx * d_7;
        DataVector &dv_19 = temps.at(19);
        dv_19 = d_41 * dv_18;
        DataVector &dv_20 = temps.at(20);
        dv_20 = Dx * Dx;
        DataVector &dv_21 = temps.at(21);
        dv_21 = Dy * Dy;
        DataVector &dv_22 = temps.at(22);
        dv_22 = Dz * Dz;
        DataVector &dv_23 = temps.at(23);
        dv_23 = d_19 * dv_22;
        DataVector &dv_24 = temps.at(24);
        dv_24 = Dz * d_23;
        DataVector &dv_25 = temps.at(25);
        dv_25 = -4 * d_1 * dv_16 * dv_24 * zp + d_33 * dv_23 + d_37 * dv_20 + d_39 * dv_21 + dv_20 + dv_21 + dv_22;
        DataVector &dv_26 = temps.at(26);
        dv_26 = -d_27 * dv_15 - dv_17 * dv_19 + dv_25;
        DataVector &dv_27 = temps.at(27);
        dv_27 = -dv_17;
        DataVector &dv_28 = temps.at(28);
        dv_28 = Dz * d_33;
        DataVector &dv_29 = temps.at(29);
        dv_29 = d_50 * dv_28;
        DataVector &dv_30 = temps.at(30);
        dv_30 = Dz * d_32 * d_53 * zp;
        DataVector &dv_31 = temps.at(31);
        dv_31 = Dy * d_58;
        DataVector &dv_32 = temps.at(32);
        dv_32 = d_48 * dv_14;
        DataVector &dv_33 = temps.at(33);
        dv_33 = Dx * d_37 + Dx;
        DataVector &dv_34 = temps.at(34);
        dv_34 = Dx * d_10 * d_36 * d_6 - 2 * Dy * M * d_17 * d_22 * d_35 * d_38 - Dy + d_11 * dv_11 + dv_32 * (d_25 * d_38 * d_71 - d_28 * d_44 + ypdot);
        DataVector &dv_35 = temps.at(35);
        dv_35 = Dx * xpddot;
        DataVector &dv_36 = temps.at(36);
        dv_36 = Dy * ypddot;
        DataVector &dv_37 = temps.at(37);
        dv_37 = Dz * zpddot;
        DataVector &dv_38 = temps.at(1);
        dv_38 = dv_0 + dv_1;
        DataVector &dv_39 = temps.at(0);
        dv_39 = Dy + d_14 * dv_38;
        DataVector &dv_40 = temps.at(4);
        dv_40 = d_14 * (dv_3 + dv_4) + dv_7;
        DataVector &dv_41 = temps.at(3);
        dv_41 = pow(dv_26, -3.0 / 2.0);

        get(get<CurvedScalarWave::Tags::Psi>(*result)) =
            pow(dv_26, -1.0 / 2.0);
        get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*result)) =
            1.0 / dv_26 * (-Dx * d_29 * (Dz * d_50 * d_67 - dv_31) + Dx * d_66 * dv_27 - d_30 * d_50 * dv_22 + d_34 * d_55 * dv_20 + d_38 * d_55 * dv_21 + d_53 * dv_23 + d_60 * d_62 * dv_20 - d_60 * d_63 * dv_21 + d_65 * dv_18 * dv_27 + dv_11 * dv_31 - 1.0 / 2.0 * dv_15 * (d_1 * d_46 * d_52 * d_84 * (d_88 * d_89 + d_88 * xp * xpdot + d_94 * (3 * d_21 * d_85 * d_90 - d_92 * (d_0 * d_83 * d_86 + d_17 * d_91))) + d_25 * d_45 * (d_13 * zpdot * (-d_14 * d_89 * d_95 - d_14 * d_96 * xpdot + zpdot * (-d_13 * d_19 * d_93 * d_97 + 2)) + d_14 * d_78 + d_43 * d_77 + d_69 * (d_98 * (-d_32 * d_4 + xp * (d_0 - d_1)) + xpdot * (-d_100 * d_5 + d_2 * (d_32 + d_79 * d_95)) + ypdot * (-d_100 * d_9 + d_2 * (d_101 + d_96 * yp))) + d_7 * d_76 + d_71 * (d_98 * (d_0 * yp - d_1 * yp + d_101 * d_5) + xpdot * (d_10 * d_5 * d_99 + d_2 * (-d_101 + d_84 * d_85 * d_86 * xp * yp)) + ypdot * (d_2 * (d_32 + d_80 * d_95) - d_42 * d_9 * d_99))) + d_76 * xpdot + d_77 * ypdot + d_78 * zpdot) * 1.0 / (d_47 * d_47) + dv_16 * dv_29 - dv_16 * dv_30 + dv_32 * (d_19 * d_33 * dv_37 - d_24 * d_74 * (d_13 * d_50 * dv_38 + zp * (Dy * zpddot + Dz * ypddot)) - d_25 * d_58 * dv_39 + d_29 * (Dx * d_50 * d_67 * zpdot + Dz * d_50 * d_67 * xpdot - d_11 * (Dx * ypddot + Dy * xpddot) + d_13 * zp * (Dx * zpddot + Dz * xpddot) - d_58 * dv_5 - d_58 * dv_6) - d_34 * d_73 * dv_8 + d_37 * dv_35 - d_38 * d_73 * dv_9 + d_39 * dv_36 + d_50 * d_72 * dv_10 - d_62 * d_75 * dv_8 + d_63 * d_75 * dv_9 + d_64 * d_74 * dv_39 - d_64 * dv_12 - d_65 * d_7 * dv_40 - d_66 * dv_40 + dv_29 - dv_30 + dv_35 + dv_36 + dv_37) - dv_34 * ypdot + xpdot * (d_29 * dv_27 - d_70 * dv_32 + dv_33) + zpdot * (Dz + d_19 * dv_28 - d_31 * dv_16 + d_31 * dv_18 - d_68 * dv_32)) /
            sqrt(-d_48 * dv_15 + dv_19 * dv_27 + dv_25);
        get<0>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>, Frame::Inertial>>(
            *result)) =
            -dv_41 * (-d_27 * d_70 * dv_14 - d_29 * dv_17 + dv_33);
        get<1>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>, Frame::Inertial>>(
            *result)) =
            dv_34 * dv_41;
        get<2>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>, Frame::Inertial>>(
            *result)) =
            -dv_41 * (Dz - d_102 * dv_16 + d_102 * dv_18 + d_19 * d_92 * dv_24 - d_27 * d_68 * (-d_28 * (Dy * d_15 + dv_2) + d_29 * (Dx * d_16 - dv_17 * xpdot) + dv_13));
    }
}
