#Import relevant packages
import os
from sys import argv
from scipy.integrate import quad
import numpy as np
np.set_printoptions(legacy='1.13')
from scipy.optimize import root_scalar,root
from jinja2 import Template

#Define useful functions for GR
def spacetime_metric(position):
    pass

def lorentz_factor(position, velocity):
    pass

#Define functions related to the the orbital parameters
def orbital_parameters(vel, rmax):
    pass

#Define functions to set the worldtube radii for BH and scalar charge
##It uses a broken power law series and root finds for the r0 values which is fixed base don values at r_p=6M
def power_law(r,r6,amp,exp,delta):
    return amp*(r/r6)**exp*(1+(r/r6)**(1/delta-1))**(-exp*delta)

def root_find_func(r6,r_isco,r_inf,exp,delta):
    return [power_law(6.0,r6[0],r_inf,exp,delta)-r_isco]

#Define the function to setup the simulation depending on the argument
##Thinking whether to use a class or not and whether to call for argv at this stage or in the main codeblock

#Run the actual code if this python file is run directly.
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
    num_geodesic_orbits,):

    base_dir = f"Quasicircular_rad{int(orbit_radius)}_eps0p{str(particle_charge).split('.')[1]}"
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
        "turn_on_time": num_geodesic_orbits * radial_period, #-100
        # the interval over which the self force is turned on
        "turn_on_interval": 1000.0,
        # how many iterations are done
        "iterations": iterations,
        # sets the parameters of the function that controls the excision sphere radii
        "r0": r0,
        "amp": amp,
        "exp": 1.5,
        "delta": 0.05,
        "r0_bh": r0_bh,
        "amp_bh": amp_bh,
        "exp_bh": 1.0,
        "delta_bh": 0.05,
    }
    # print(config_dict)
    input_template = Template(input_template)
    input_file = os.path.join(base_dir, "input_file.yaml")
    with open(input_file, "w") as f:
        f.write(input_template.render(**config_dict))
    os.chdir(base_dir)
    os.system("sbatch urania.sh")
    os.chdir("../")
    return

if __name__ == "__main__":
    #Read in the input and slurm launch files
    with open("input_template.yaml", "r") as f:
        input_template = f.read()
    with open("launch_template.sh", "r") as f:
        launch_template = f.read()

    #Set the orbital parameters
    orbit_radius=10.5
    initial_velocity=np.sqrt(1/orbit_radius)
    radial_period=orbit_radius**(-3/2)*np.sqrt((orbit_radius-6)/(orbit_radius-3)) #Temporarily setting for circular orbits only

    #Set up the wordtube and BH excision parameters
    delta=0.05

    exp=1.5 #Exponent for the power law
    amp=3.0 #WT radius @ infinity
    wt_radius_at_isco=0.4
    wt_radius_at_inf=3.0
    r0=root(root_find_func,x0=(23.0),
    args=(wt_radius_at_isco,wt_radius_at_inf,exp,delta)).x[0] #WT excision radius (Fixed parameters with value @ r_p=6M)
    worldtube_radius=power_law(orbit_radius,r0,wt_radius_at_inf,exp,delta) #WT radius at apastron

    exp_bh=1.0 #Exponent for the power law
    amp_bh=1.0 #BH excision @ infinity
    bh_radius_at_isco = 1.5
    bh_radius_at_inf = 1.9
    r0_bh=root(root_find_func,x0=(12.0),
    args=(bh_radius_at_isco,bh_radius_at_inf,exp_bh,delta)).x[0] #BH excision radius (Fixed parameters with value @ r_p=6M)

    #Strength of scalar field parameters
    epsilon=particle_mass=particle_charge=0.02 #Epsilon
    iterations=3 #The number of iterations for the self force to be updated
    num_geodesic_orbits=4 #The number of geodesic orbits to be run before the self force is turned on

    #Setting simulation parameters
    expansion_order=1 #Order of the puncture field?
    lev=0 #The level of refinement of the grid points
    node_num=1
    
    submit_job(
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
        num_geodesic_orbits
    )
       