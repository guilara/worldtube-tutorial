import os
from scipy.integrate import quad
import numpy as np
from scipy.optimize import root_scalar, root
from jinja2 import Template


def spacetime_metric(position):
    radii = np.linalg.norm(position)
    g00 = -1 + 2.0 / radii
    gi0 = np.einsum("...,...j", 2 / radii**2, position)
    gij = np.einsum("...,...j,...k->...jk", 2.0 / radii**3, position, position)
    gij += np.eye(3, 3)
    g_mu_nu = np.empty((4, 4))
    g_mu_nu[0, 0] = g00
    g_mu_nu[1:4, 0] = gi0
    g_mu_nu[0, 1:4] = gi0
    g_mu_nu[1:4, 1:4] = gij
    return g_mu_nu


def lorentz_factor(position, velocity):
    metric = spacetime_metric(position)
    u0sq = -1 / (
        metric[0, 0]
        + 2.0 * np.einsum("...i,...i->...", metric[0, 1:], velocity)
        + np.einsum("...jk,...j,...k->...", metric[1:, 1:], velocity, velocity)
    )
    return np.sqrt(u0sq)


# takes velocity and apoapsis, returns apoapsis, periapsis, semi-latus rectum, eccentricity, energy and angular momentum
def orbital_parameters(vel, rmax):
    time_killing = np.asarray([-1.0, 0.0, 0.0, 0.0])
    position = np.asarray([rmax, 0.0, 0.0])
    velocity = np.asarray([0.0, vel, 0.0])
    metric = spacetime_metric(position)
    u0 = lorentz_factor(position, velocity)
    ui = np.einsum("...i,...->...i", velocity, u0)
    u = np.concatenate((np.asarray([u0]), ui))
    azimuthal_killing = np.zeros_like(u)
    azimuthal_killing[1] = -position[1]
    azimuthal_killing[2] = position[0]
    energy = np.einsum("...ij,...i,j->...", metric, u, time_killing)
    ang_mom = np.einsum("...ij,...i,...j->...", metric, u, azimuthal_killing)
    a0 = 2.0 * ang_mom**2 / (energy**2 - 1)
    a1 = -(ang_mom**2) / (energy**2 - 1)
    a2 = 2.0 / (energy**2 - 1)
    bigq = (a2**2 - 3.0 * a1) / 9.0
    bigr = (2.0 * a2**3 - 9.0 * a2 * a1 + 27.0 * a0) / 54.0
    theta = np.arccos(bigr / np.sqrt(bigq**3))
    ra = -2.0 * np.sqrt(bigq) * np.cos((theta + 2.0 * np.pi) / 3.0) - a2 / 3.0
    rp = -2.0 * np.sqrt(bigq) * np.cos((theta - 2.0 * np.pi) / 3.0) - a2 / 3.0
    semi = 2.0 * ra * rp / (ra + rp)
    ecc = (ra - rp) / (ra + rp)
    return ra, rp, semi, ecc, energy, ang_mom


def p_from_vel(vel, rmax, target_p):
    _, _, semi, _, _, _ = orbital_parameters(vel, rmax)
    return semi - target_p


def radial_eq(r, energy, ang_mom):
    fact = (1.0 - 2.0 / r) * np.sqrt(
        1.0 - (1 - 2 / r) * (1 + ang_mom**2 / r**2) / energy**2
    )
    return -1 / fact


def orbital_period(ra, rp, energy, ang_mom):
    return 2.0 * quad(radial_eq, ra, rp, (energy, ang_mom))[0]


def determine_initial_geodesic(rmax, target_p):
    target_vel = root_scalar(
        p_from_vel, (rmax, target_p), x0=0.8 * rmax**-0.5, x1=0.9 * rmax**-0.5
    )
    ra, rp, semi, ecc, energy, ang_mom = orbital_parameters(target_vel.root, rmax)
    period = orbital_period(ra, rp, energy, ang_mom)
    print(
        f"Initial geodesic with apoapsis {ra:.4f} M, periapsis  {rp:.4f} M, semi-latus rectum {semi:.4f} M, eccentricity {ecc:.9f} and radial period {period:.4f} M"
    )
    return target_vel.root, period


def power_law(r, rb, exp, delta, amp):
    return amp * (r / rb) ** exp * (1 + (r / rb) ** (1 / delta)) ** (-exp * delta)


def power_law_root_find(rb, radius_at_isco, radius_at_inf, exponent, delta):
    return [
        power_law(6.0, rb[0], exponent, delta, radius_at_inf) - radius_at_isco,
    ]


def submit_job(
    input_template,
    launch_template,
    orbit_radius,
    initial_velocity,
    radial_period,
    worldtube_radius,
    r0,
    amp,
    r0_bh,
    amp_bh,
    particle_mass,
    particle_charge,
    expansion_order,
    lev,
    iterations,
    node_num,
    wt_radius_at_isco,
    num_geodesic_orbits,
):
    base_dir = f"r{int(orbit_radius)}_R{int(wt_radius_at_isco * 10)}_n{expansion_order}_eps{int(particle_charge * 1000)}_lev{lev}_it{iterations}_final2"
    os.mkdir(base_dir)

    launch_dict = {"sim_name": base_dir, "num_nodes": node_num}
    launch_template = Template(launch_template)
    launch_file = os.path.join(base_dir, "urania.sh")
    with open(launch_file, "w") as f:
        f.write(launch_template.render(**launch_dict))

    config_dict = {
        # the initial radius of the orbit
        "orbit_radius": orbit_radius,
        # the worldtube radius at apastron
        "worldtube_radius": worldtube_radius,
        # the angular velocity, needed by the rotation map
        "angular_vel": initial_velocity / orbit_radius,
        "particle_position": [orbit_radius, 0.0, 0.0],
        "particle_velocity": [0.0, initial_velocity, 0.0],
        # the slab interval where volume data is observed
        "observe_volume_interval": 20000,
        # the slab interval indicating how often data about the charge is observed
        "ylm_obs_interval": 100,
        # the slab interval indicating how often the waveform is observed
        "observe_spheres_interval": 100,
        "particle_mass": particle_mass,
        "particle_charge": particle_charge,
        # sets the resolution
        "P": lev,
        # expansion order of the worldtube scheme
        "expansion_order": expansion_order,
        # when the self force is turned on
        "turn_on_time": num_geodesic_orbits * radial_period - 100.0,
        # the interval over which the self force is turned on
        "turn_on_interval": 100.0,
        # how many iterations are done
        "iterations": iterations,
        # sets the parameters of the function that controls the excision sphere radii
        "rb": r0,
        "amp": amp,
        "exp": 1.5,
        "delta": 0.05,
        "rb_bh": r0_bh,
        "amp_bh": amp_bh,
        "exp_bh": 1.0,
        "delta_bh": 0.05,
    }
    input_template = Template(input_template)
    input_file = os.path.join(base_dir, "input_file.yaml")
    with open(input_file, "w") as f:
        f.write(input_template.render(**config_dict))
    os.chdir(base_dir)
    os.system("sbatch urania.sh")
    os.chdir("../")


if __name__ == "__main__":
    config_file = "input_template.yaml"
    with open(config_file, "r") as f:
        input_template = f.read()

    launch_file = "launch_template.sh"
    with open(launch_file, "r") as f:
        launch_template = f.read()

    # find r0 based on radius at isco and infinity
    delta = 0.05
    wt_radius_at_isco = 0.4
    wt_radius_at_inf = 3.0
    bh_radius_at_isco = 1.5
    bh_radius_at_inf = 1.9
    exp = 1.5
    exp_bh = 1.0

    power_law_solution = root(
        power_law_root_find,
        x0=(23.0),
        args=(wt_radius_at_isco, wt_radius_at_inf, exp, delta),
    )
    power_law_solution_bh = root(
        power_law_root_find,
        x0=(12.0),
        args=(bh_radius_at_isco, bh_radius_at_inf, exp_bh, delta),
    )
    r0 = power_law_solution.x[0]
    r0_bh = power_law_solution_bh.x[0]

    orbit_radius = 100.0
    epsilon = 0.02
    lev = 0
    expansion_order = 1
    semi_latus_rectum = 10.0
    initial_velocity, radial_period = determine_initial_geodesic(
        orbit_radius, semi_latus_rectum
    )

    node_num = 2
    iterations = 3
    num_geodesic_orbits = 4

    wt_radius_at_apostron = power_law(orbit_radius, r0, 1.5, delta, wt_radius_at_inf)

    submit_job(
        input_template,
        launch_template,
        orbit_radius,
        initial_velocity,
        radial_period,
        wt_radius_at_apostron,
        r0,
        wt_radius_at_inf,
        r0_bh,
        bh_radius_at_inf,
        epsilon,
        epsilon,
        expansion_order,
        lev,
        iterations,
        node_num,
        wt_radius_at_isco,
        num_geodesic_orbits,
    )
