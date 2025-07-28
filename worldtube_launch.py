#Import relevant packages
import os
import sys
from scipy.integrate import quad
import numpy as np
np.set_printoptions(legacy='1.13')
from scipy.optimize import root
from jinja2 import Template

#Define the radial roots of the orbit
def turning_points(p,e):
    ra=p/(1-e)
    rp=p/(1+e)
    return ra,rp

#Define the velocity at apoastron
def va(p,e):
    return (1-e)/p**0.5*np.sqrt((2-2*e-p)/(2+2*e-p))

#Define functions to set the worldtube radii for BH and scalar charge
##It uses a broken power law series and root finds for the r0 values which is fixed based on values at r_p=6M
def power_law(r,r6,amp,exp,delta):
    return amp*(r/r6)**exp*(1+(r/r6)**(1/delta))**(-exp*delta)

def root_find_func(r6,R_isco,R_inf,exp,delta): #Remember to set r6=6.0 for normal runs
    return [power_law(6.0,r6[0],R_inf,exp,delta)-R_isco]

#Define the function to setup the simulation depending on the argument
##Thinking whether to use a class or not and whether to call for argv at this stage or in the main codeblock

#Run the actual code if this python file is run directly.
def submit_job(
    input_template,
    launch_template,
    orbit_radius,
    initial_velocity,
    worldtube_radius,
    bh_radius,
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
    num_geodesic_orbits,
    turn_on_time,
    turn_on_interval,):


    circular_test='Quasicircular' if initial_velocity==1/orbit_radius**0.5 else 'Eccentric'
    base_dir = f"{circular_test}_r{int(orbit_radius)}_eps{str(particle_charge).split('.')[1]}"
    print(f"Running simulation in base directory: {base_dir}")
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
        #The BH radius at apastron
        "bh_radius": bh_radius,
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
        "turn_on_time": turn_on_time, #num_geodesic_orbits * 2*np.pi/(initial_velocity/orbit_radius) # for eccentric/hyperbolic radial_period, #-100
        # the interval over which the self force is turned on
        "turn_on_interval": turn_on_interval,
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

    if (len(sys.argv))<2:
        print("Please provide p,e as argument for running the .py file!")
        sys.exit()

    elif (len(sys.argv))>2:
        print("Extra arguments not allowed. Please provide p,e as a single string!")
        sys.exit()

    else:
        #Parse the string
        semi_lactus,ecc=[float(x) for x in sys.argv[1].split(',')]

        #Set the orbital parameters
        orbit_radius=turning_points(semi_lactus,ecc)[0]
        initial_velocity=va(semi_lactus,ecc)
        radial_period=None if ecc==0 else 2*np.pi/(orbit_radius**(-3/2)*np.sqrt((orbit_radius-6)/(orbit_radius-3)))
        print(f"Orbital radius:{orbit_radius}, Initial velocity:{initial_velocity}, Radial Period (if circular):{radial_period}")

    #Set up the wordtube and BH excision parameters
    delta=0.05

    exp=1.5 #Exponent for the power law
    amp=3.0 #WT radius @ infinity
    wt_radius_at_isco=0.8
    r0=root(root_find_func,x0=(23.0),args=(wt_radius_at_isco,amp,exp,delta)).x[0] #WT excision radius (Fixed parameters with value @ r_p=6M)
    worldtube_radius=power_law(orbit_radius,r0,amp,exp,delta) #WT radius at t=0
    fixed_worldtube_radius=power_law(6.0,r0,amp,exp,delta) #WT radius at 6M
    print(f"Initial WT radius:{worldtube_radius}, $r_0$:{r0}, WT radius at 6M is 0.8?:{fixed_worldtube_radius-0.8<0.0001}") #Verification of WT parameters

    exp_bh=1.0 #Exponent for the power law
    amp_bh=1.9 #BH excision @ infinity
    bh_radius_at_isco = 1.5
    r0_bh=root(root_find_func,x0=(23.0),
    args=(bh_radius_at_isco,amp_bh,exp_bh,delta)).x[0] #BH excision radius (Fixed parameters with value @ r_p=6M)
    bh_radius=power_law(orbit_radius,r0_bh,amp_bh,exp_bh,delta) #BH radius at t=0
    fixed_bh_radius=power_law(6.0,r0_bh,amp_bh,exp_bh,delta) #BH radius at 6M
    print(f'Initial BH radius:{bh_radius}, $r0_bh$:{r0_bh}, BH radius at 6M is 1.9?:{fixed_bh_radius-1.9<0.0001}') #Verification of BH parameters


    #Strength of scalar field parameters
    epsilon=particle_mass=particle_charge=0.08 #Epsilon
    print(f"Epsilon:{epsilon}")
    iterations=3 #The number of iterations for the self force to be updated
    num_geodesic_orbits=4 #The number of geodesic orbits to be run before the self force is turned on
    turn_on_interval= 1000. if ecc==0 else 100.
    turn_on_time=1500. if ecc==0 else num_geodesic_orbits * radial_period - turn_on_interval

    #Setting simulation parameters
    expansion_order=1 #Order of the puncture field?
    lev=0 #The level of refinement of the grid points
    node_num=4

    response=input("Do you want to run the simulation? (y/n): ")
    if response=='y':
        #Submit the job    
        submit_job(
            input_template,
            launch_template,
            orbit_radius,
            initial_velocity,
            worldtube_radius,
            bh_radius,
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
            num_geodesic_orbits,
            turn_on_time,
            turn_on_interval,
        )
    else:
        print('Batch Job was not sent. Rerun script with y')
    