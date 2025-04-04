#!/usr/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np
import matplotlib.pyplot as plt
import worldtube_postprocess as wt
import h5py
from spectre.SphericalHarmonics import SpherepackIterator


def ylm_no_normalization(data, l, m, l_max):
    iterator = SpherepackIterator(l_max, l_max)
    iterator.set(l, m)
    harmonic = data[:, iterator()]
    return harmonic


if __name__ == "__main__":
    # Use Latex
    plt.rcParams.update({"font.size": 20})
    plt.rcParams.update(
        {
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath} \usepackage{amsfonts}",
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica"],
            "lines.linewidth": 1.0,
        }
    )
    # Colorblind palettes
    # plt.style.use('tableau-colorblind10')

    sim_name = "r100_R4_n1_eps20_lev0_it3_final2"
    print(sim_name)

    plots_directory = (
        "/urania/ptmp/guilara/spectre/Worldtube/Test/" + sim_name + "/plots/"
    )
    sim_dict = wt.extract_sim_data(sim_name)

    # Print simulation data keys
    print(sim_dict.keys())

    # Plot trajectory

    fig = plt.figure(constrained_layout=True, figsize=(7.4, 6.8))
    widths = [1]
    heights = [1]
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=widths, height_ratios=heights)

    fig_ax = fig.add_subplot(gs[0, 0], aspect="equal")

    # plt.plot(sim_dict["times"], sim_dict["radii"], label = "radii")
    # plt.plot(sim_dict["times"], sim_dict["radial_vel"], label = "radial_vel")
    plt.plot(sim_dict["position"][:, 0], sim_dict["position"][:, 1], label="position")

    plt.xlim(-110, 110)
    plt.ylim(-110, 110)

    plt.savefig(plots_directory + "test.pdf")
    plt.close()

    #######################################
    #######################################

    # Modes
    # These are not the modes of the scalar field on the extraction sphere
    # These are the modes describing the shape of the Strahlkorper describing the
    # extraction sphere itself
    with h5py.File(f"{sim_name}/Reductions.h5", "r") as h5file:
        dat_file = h5file.get("/Spheres_Ylm.dat")
        legend = list(dat_file.attrs["Legend"])

        # print([legend.index(x) for x in legend])
        # print(legend)

        # Get time and coeffs
        time = dat_file[:, 0]
        coeffs = dat_file[:, :]

    fig = plt.figure(constrained_layout=True, figsize=(7.4, 6.8))
    widths = [1]
    heights = [1]
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=widths, height_ratios=heights)

    fig_ax = fig.add_subplot(gs[0, 0])

    lmodes = [1, 2]
    for l_temp in lmodes:
        lm_modes = [(l_temp, m) for m in range(-l_temp, l_temp + 1)]
        for l, m in lm_modes:
            label = "(" + str(l) + "," + str(m) + ")"
            index = legend.index("coef" + label)
            plt.plot(
                coeffs[:, 0],
                coeffs[:, index],
                label=(f"$ a^{{{label}}} $" if m >= 0 else f"$ b^{{{label}}} $"),
            )
            # print(label)
            # print("Min: ", np.min(coeffs[:, index]))
            # print("Max: ", np.max(coeffs[:, index]))

    plt.legend()

    plt.savefig(plots_directory + "test_strahlkorper.pdf")
    plt.close()

    # Inertial expansion center
    fig = plt.figure(constrained_layout=True, figsize=(7.4, 6.8))
    widths = [1]
    heights = [1]
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=widths, height_ratios=heights)

    fig_ax = fig.add_subplot(gs[0, 0])

    x = coeffs[:, legend.index("InertialExpansionCenter_x")]
    y = coeffs[:, legend.index("InertialExpansionCenter_y")]

    # plt.plot(x, y, label="Inertial Expansion Center")
    plt.plot(time, x, label="Inertial Expansion Center")
    plt.plot(time, y, label="Inertial Expansion Center")

    plt.legend()

    plt.savefig(plots_directory + "test_expansion_center.pdf")
    plt.close()

    #######################################
    #######################################

    # Modes from sphere data
    # These are not the actual modes. They are spherepack coefficients
    # with respect to a real spherical harmonic (cos and sin) decomposition
    # see Spherepack documentation

    # Get lmax from datafile
    with h5py.File(f"{sim_name}/Reductions.h5", "r") as h5file:
        dat_file = h5file.get("/Spheres_Ylm.dat")
        legend = list(dat_file.attrs["Legend"])
        lmax = int(dat_file[0, legend.index("Lmax")])
        print("lmax : ", lmax)

    times, all_spectral_data, radii = wt.extract_sphere_data(sim_name)

    fig = plt.figure(constrained_layout=True, figsize=(15.0, 5.0))
    widths = [1]
    heights = [1]
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=widths, height_ratios=heights)

    fig_ax = fig.add_subplot(gs[0, 0])

    lmodes = range(1, 2)
    for l_temp in lmodes:
        lm_modes = [(l_temp, m) for m in range(-l_temp, l_temp + 1)]
        for l, m in lm_modes:
            label = "(" + str(l) + "," + str(np.abs(m)) + ")"
            index = legend.index("coef" + label)
            radius_index = 0
            # multiply by radius
            plt.plot(
                times,
                radii[radius_index]
                * wt.ylm(data=all_spectral_data[radius_index], l=l, m=m, l_max=lmax),
                label=(f"$ a^{{{label}}} $" if m >= 0 else f"$ b^{{{label}}} $"),
            )

    plt.xlabel(r"$t$")
    plt.ylabel(r"$R \Psi^{(\ell, m)}$")
    plt.legend()

    plt.savefig(plots_directory + "test_modes_again.pdf")
    plt.close()

    #######################################
    #######################################

    # Plot the standard spherical harmonic coefficients

    fig = plt.figure(constrained_layout=True, figsize=(15.0, 5.0))
    widths = [1]
    heights = [1]
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=widths, height_ratios=heights)

    fig_ax = fig.add_subplot(gs[0, 0])

    lmodes = range(1, 2)
    for l_temp in lmodes:
        lm_modes = [(l_temp, m) for m in range(0, l_temp + 1)]
        for l, m in lm_modes:
            label = "(" + str(l) + "," + str(m) + ")"
            negative_m_mode_label = "(" + str(l) + "," + str(-m) + ")"
            index = legend.index("coef" + label)
            radius_index = 0

            # Assemble normalization factors explicitely
            sign = 1 if abs(m) % 2 == 0 else -1
            normalization_factor = np.sqrt(np.pi / 2.0)
            # Notice the different ylm function defined without normalizations
            alm = ylm_no_normalization(
                data=all_spectral_data[radius_index], l=l, m=m, l_max=lmax
            )
            blm = ylm_no_normalization(
                data=all_spectral_data[radius_index], l=l, m=-m, l_max=lmax
            )
            positive_m_mode = sign * normalization_factor * (alm + 1j * blm)
            # No sign here
            if m != 0:
                negative_m_mode = normalization_factor * (alm - 1j * blm)

            # Plot real part
            plt.plot(
                times,
                # multiply by radius
                radii[radius_index] * np.real(positive_m_mode),
                label=f"$ \mathrm{{Re}} \Psi^{{{label}}} $",
            )
            if m != 0:
                plt.plot(
                    times,
                    radii[radius_index] * np.real(negative_m_mode),
                    label=f"$ \mathrm{{Re}} \Psi^{{{negative_m_mode_label}}} $",
                )

            # Plot imaginary part
            # plt.plot(
            #     times,
            #     radii[radius_index] * np.imag(positive_m_mode),
            #     "--",
            #     label=f"$ \mathrm{{Im}} \Psi^{{{label}}} $",
            # )

            # if m != 0 :
            #     plt.plot(
            #         times,
            #         radii[radius_index] * np.imag(negative_m_mode),
            #         "--",
            #         label=f"$ \mathrm{{Im}} \Psi^{{{negative_m_mode_label}}} $",
            #     )

    plt.xlabel(r"$t$")
    plt.ylabel(r"$R \Psi^{(\ell, m)}$")
    plt.legend(loc="lower left")

    plt.savefig(plots_directory + "test_modes_ylm.pdf")
    plt.close()
