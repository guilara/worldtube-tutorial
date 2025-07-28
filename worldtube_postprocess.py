import numpy as np
import os
import yaml
import spectre.IO.H5 as spectre_h5
from spectre.SphericalHarmonics import Spherepack, SpherepackIterator


def spacetime_metric(position):
    radii = np.linalg.norm(position, axis=1)
    g00 = -1 + 2.0 / radii
    gi0 = np.einsum("...,...j", 2 / radii**2, position)
    gij = np.einsum("...,...j,...k->...jk", 2.0 / radii**3, position, position)
    gij += np.asarray([np.eye(3, 3) for _ in g00])
    g_mu_nu = np.empty((np.shape(g00)[0], 4, 4))
    g_mu_nu[:, 0, 0] = g00
    g_mu_nu[:, 1:4, 0] = gi0
    g_mu_nu[:, 0, 1:4] = gi0
    g_mu_nu[:, 1:4, 1:4] = gij
    return g_mu_nu


def lorentz_factor(position, velocity):
    metric = spacetime_metric(position)
    u0sq = -1 / (
        metric[:, 0, 0]
        + 2.0 * np.einsum("...i,...i->...", metric[:, 0, 1:], velocity)
        + np.einsum("...jk,...j,...k->...", metric[:, 1:, 1:], velocity, velocity)
    )
    return np.sqrt(u0sq)


def ylm(data, l, m, l_max):
    iterator = SpherepackIterator(l_max, l_max)
    iterator.set(l, m)
    harmonic = data[:, iterator()]
    sign = 1 if abs(m) % 2 == 0 else -1
    return np.sqrt(np.pi / 2.0) * sign * harmonic


def extract_sphere_data(file_name):
    with spectre_h5.H5File(f"{file_name}/Surface.h5") as f:
        input_data = list(yaml.full_load_all(f.input_source()))[1]

        l_max = input_data["InterpolationTargets"]["Spheres"]["LMax"]
        radii = input_data["InterpolationTargets"]["Spheres"]["Radius"]
        spherepack = Spherepack(l_max, l_max)
        physical_size = spherepack.physical_size
        spectral_size = spherepack.spectral_size
        sphere_file = f.get_vol("/Spheres")
        obs_ids = sphere_file.list_observation_ids()
        all_spectral_data = np.empty((len(radii), len(obs_ids), spectral_size))
        for k, id in enumerate(obs_ids):
            data_all_radii = np.asarray(
                sphere_file.get_tensor_component(id, "Psi").data
            )
            for i, radius in enumerate(radii):
                data_radius = data_all_radii[
                    i * physical_size : (i + 1) * physical_size
                ]
                spectral_data = spherepack.phys_to_spec(data_radius)
                all_spectral_data[i, k, :] = spectral_data
        times = np.asarray([sphere_file.get_observation_value(id) for id in obs_ids])

    restart_index = 1
    while os.path.exists(f"{file_name}/Restart{restart_index}/Reductions.h5"):
        with spectre_h5.H5File(f"{file_name}/Restart{restart_index}/Surface.h5") as f:
            sphere_file = f.get_vol("/Spheres")
            obs_ids = sphere_file.list_observation_ids()
            all_spectral_data_new = np.empty((len(radii), len(obs_ids), spectral_size))
            for k, id in enumerate(obs_ids):
                data_all_radii = np.asarray(
                    sphere_file.get_tensor_component(id, "Psi").data
                )
                for i, radius in enumerate(radii):
                    data_radius = data_all_radii[
                        i * physical_size : (i + 1) * physical_size
                    ]
                    spectral_data = spherepack.phys_to_spec(data_radius)
                    all_spectral_data_new[i, k, :] = spectral_data
            times_new = np.asarray(
                [sphere_file.get_observation_value(id) for id in obs_ids]
            )
        all_spectral_data = np.concatenate(
            (all_spectral_data, all_spectral_data_new), axis=1
        )
        times = np.concatenate((times, times_new), axis=0)
        restart_index += 1
    if os.path.exists(f"{file_name}/ringdown/Surfaces.h5"):
        with spectre_h5.H5File(f"{file_name}/ringdown/Surfaces.h5") as f:
            sphere_file = f.get_vol("/Spheres")
            obs_ids = sphere_file.list_observation_ids()
            all_spectral_data_new = np.empty((len(radii), len(obs_ids), spectral_size))
            for k, id in enumerate(obs_ids):
                data_all_radii = np.asarray(
                    sphere_file.get_tensor_component(id, "Psi").data
                )
                for i, radius in enumerate(radii):
                    data_radius = data_all_radii[
                        i * physical_size : (i + 1) * physical_size
                    ]
                    spectral_data = spherepack.phys_to_spec(data_radius)
                    all_spectral_data_new[i, k, :] = spectral_data
            times_new = np.asarray(
                [sphere_file.get_observation_value(id) for id in obs_ids]
            )
            times_new += times[-1]
        all_spectral_data = np.concatenate(
            (all_spectral_data, all_spectral_data_new), axis=1
        )
        times = np.concatenate((times, times_new), axis=0)
    return times, all_spectral_data, radii


def extract_sim_data(file_name):
    with spectre_h5.H5File(f"{file_name}/Reductions.h5") as f:
        input_data = list(yaml.full_load_all(f.input_source()))[1]
        dat_file = f.get_dat("/PsiTaylorCoefs")
        coef_data = np.asarray(dat_file.get_data())
        f.close_current_object()

    restart_index = 1
    while os.path.exists(f"{file_name}/Restart{restart_index}/Reductions.h5"):
        with spectre_h5.H5File(
            f"{file_name}/Restart{restart_index}/Reductions.h5"
        ) as f:
            dat_file = f.get_dat("/PsiTaylorCoefs")
            coef_data = np.concatenate(
                (coef_data, np.asarray(dat_file.get_data())), axis=0
            )
            f.close_current_object()
        restart_index += 1

    obs_interval = input_data["Worldtube"]["ObserveCoefficientsTrigger"]["Slabs"][
        "EvenlySpaced"
    ]["Interval"]
    times = coef_data[:, 0]
    time_steps = (times[1:] - times[:-1]) / obs_interval
    charge = input_data["Worldtube"]["Charge"]
    position = coef_data[:, 1:4]
    radii = np.linalg.norm(position, axis=1)
    velocity = coef_data[:, 4:7]
    acceleration = coef_data[:, 7:10]

    x = position[:, 0]
    y = position[:, 1]
    xdot = velocity[:, 0]
    ydot = velocity[:, 1]
    xddot = acceleration[:, 0]
    yddot = acceleration[:, 1]

    radial_vel = (x * xdot + y * ydot) / radii
    angular_vel = (x * ydot - y * xdot) / radii**2

    radial_acc = radii**-3 * (
        x**2 * (x * xddot + ydot**2)
        + x * y * (x * yddot - 2.0 * xdot * ydot)
        + y**2 * (x * xddot + xdot**2)
        + y**3 * yddot
    )

    angular_acc = radii**-4 * (
        x * y * (2.0 * xdot**2 + y * yddot - 2.0 * ydot**2)
        - x**2 * (y * xddot + 2.0 * xdot * ydot)
        + y**2 * (2.0 * xdot * ydot - y * xddot)
        + x**3 * yddot
    )

    phases = np.arctan2(y, x)
    orbits = 0
    for i, phase in enumerate(phases):
        phases[i] += orbits * 2.0 * np.pi
        if i + 1 < len(phases) and phases[i + 1] < phase:
            orbits += 1
    dx = coef_data[:, 11]
    dy = coef_data[:, 12]
    dz = coef_data[:, 13]
    dt = coef_data[:, 14]

    dr = np.cos(phases) * dx + np.sin(phases) * dy
    dphi = radii * (-np.sin(phases) * dx + np.cos(phases) * dy)

    time_killing = np.asarray([-1.0, 0.0, 0.0, 0.0])
    metric = spacetime_metric(position)
    u0 = lorentz_factor(position, velocity)
    ui = np.einsum("...i,...->...i", velocity, u0)
    u = np.concatenate((np.asarray([u0]).T, ui), axis=1)
    azimuthal_killing = np.zeros_like(u)
    azimuthal_killing[:, 1] = -position[:, 1]
    azimuthal_killing[:, 2] = position[:, 0]

    energy = np.einsum("...ij,...i,j->...", metric, u, time_killing)
    ang_mom = np.einsum("...ij,...i,...j->...", metric, u, azimuthal_killing)

    self_force_options = input_data["Worldtube"]["SelfForceOptions"]
    sim_info = {
        "times": coef_data[:, 0],
        "position": position,
        "radii": radii,
        "velocity": velocity,
        "angular_vel": angular_vel,
        "radial_vel": radial_vel,
        "acceleration": acceleration,
        "psi0": coef_data[:, 10],
        "time_steps": time_steps,
        "charge": charge,
        "radial_acc": radial_acc,
        "angular_acc": angular_acc,
        "phases": phases,
        "dr": dr,
        "dphi": dphi,
        "dx": dx,
        "dy": dy,
        "dz": dz,
        "dtpsi0": dt,
        "all_coef_data": coef_data,
        "energy": energy,
        "angular_momentum": ang_mom,
        "u0": u0,
        "turn_on_time": (
            self_force_options["TurnOnTime"] if self_force_options is not None else None
        ),
        "iterations": (
            self_force_options["Iterations"] if self_force_options is not None else None
        ),
        "turn_on_interval": (
            self_force_options["TurnOnInterval"]
            if self_force_options is not None
            else None
        ),
        "mass": self_force_options["Mass"] if self_force_options is not None else None,
        "wt_radius": input_data["DomainCreator"]["BinaryCompactObject"]["ObjectA"][
            "InnerRadius"
        ],
    }
    return sim_info


def get_osculating_elements(sim_data):
    energy = sim_data["energy"]
    ang_mom = sim_data["angular_momentum"]
    a = energy**2-1
    b = 2.
    c = -ang_mom**2
    d = 2*ang_mom**2

    bigp = (3 * a * c - b**2)/(3 * a**2)
    bigq = (2 * b**3 - 9 * a * b * c + 27 * a**2 * d)/(27 * a**3)

    bool_list=4 * bigp**3 + 27 * bigq**2 <= 0
    print('Fraction of data points with < 3 real roots:',np.sum(bool_list)/len(bool_list))
    bigp=bigp[bool_list]
    bigq=bigq[bool_list]
    a=a[bool_list]

    theta = np.arccos(3 * bigq / (2 * bigp) * np.sqrt(-3 / bigp))
    amp = 2 * np.sqrt(- bigp / 3)

    roots = np.stack([amp * np.cos(theta / 3) - b / (3 * a),
     amp * np.cos((theta - 2 * np.pi) / 3) - b / (3 * a),
     amp * np.cos((theta - 4 * np.pi) / 3) - b / (3 * a)], axis=1) # shape (N, 3)
    roots_sorted = np.sort(roots, axis=1) # shape (N, 3)
    _, rp, ra = roots_sorted.T # Unpack: smallest, middle, largest

    semi_osculating = 2.0 * ra * rp / (ra + rp)
    ecc_osculating = (ra - rp) / (ra + rp)

    # osculating_cutoff = (
    #     -1
    #     if not np.any(np.isnan(semi_osculating))
    #     else np.where(np.isnan(semi_osculating))[0][0]
    # )
    # cutoff_index = np.searchsorted(
    #     sim_data["times"],
    #     sim_data["turn_on_time"] + 2.0 * sim_data["turn_on_interval"],
    # )
    
    time=sim_data["times"][bool_list]
    return (ra,rp,semi_osculating,ecc_osculating,time)
    # return (
    #     ra[cutoff_index:osculating_cutoff],
    #     rp[cutoff_index:osculating_cutoff],
    #     semi_osculating[cutoff_index:osculating_cutoff],
    #     ecc_osculating[cutoff_index:osculating_cutoff],
    #     sim_data["times"][cutoff_index:osculating_cutoff],
    # )