# Compiled Notes for Worldtube Masters Project

## Git

We are working with two main repositories. One is the spectre repo which includes our executables and builds. Files are then run via batch jobs in the ptmp directory. We work on a fork of the develop branch of spectre. The goal is to work on multiple forks, with each fork representing versions of codes for each PR paper. Upon approval of each fork, we can get pull requests approved to merge these forks to the main develop repo.

The second repo is the worldtube-tutorial repo which contains the python files, instructions, questions, and notes (such as this file) for the project. The is no need for additional forking, but each branch exists for individuals to contribute to the main branch a their own pace.

Listed below are some commonly used git commands which is useful for working with git on the command line:

### 1. Status
````bash
git status
````
This commands specifies the current branch and any items in the branch that is currently staged.


### 2. Change Branches
```bash
git switch {Branch Name}
```
To change between branches. You can also use the older command,
```bash
git checkout {Branch Name}
```
Just be careful as that command can also create new branches with the "-o" flag.

### 3. Staging files
```bash
git add {File/Directory/} {File/Directory/2} ...
```
To stage items in the branch that you would like to record changes off and eventually (perhaps) commit as permanent changes to the remote branch.

To remove items from the staging area, one can use,
```bash
git restore --staged {File/Directory/1} {File/Directory/2} ...
```
These commands also work recursively for entire directories. By specifying only the directory, changes for files within will be recursively staged for their respective changes.

### 4. Committing staged changes
```bash
git commit -m "{Message for your commit}"
```
Commits all staged changes to the remote branch. Other people are not able to see your changes yet on the local repository for that push is required.

It is also highly advised that you tag your commits with messages as they help highlight what each commit was done for.

### 5. Pushing changes to the remote
```bash
git push
```
Commited changes only exist on the local repository. They are the "checkpoints" with which you can review changes as they are done, with documentation on what changed, and when.

## SpECTRE

### Notes for following: [A Hitchhiker's Guide to Running SpECTRE](https://spectre-code.org/beginners_guide.html)

#### 1. ExportCoordinates3D & InputTimeDependent3D

#### 2. EvolveScalarAdvection2D & Kuzmin2D

### Notes on running SpECTRE on Urania

#### The Input File (.yaml)

#### The Bash file (.sh)
The bash file is where you send commands via SLURM to put up jobs for the HPC to perform. They specify how many nodes, compute time, and commands to run on the HPC. There are a lot of prelude lines to ensure the actual files run work (e.g. setting up python environment, exporting important paths etc.). Nevertheless the important lines/options are as below:

1. BATCH commands. An example of the code block requires for a batch job to be realised is as follows.
```bash
#!/bin/bash -
#SBATCH -J Quasicircular_rad10_eps0p02_debug_restart1
#SBATCH --nodes 4
#SBATCH -t 24:00:00
#SBATCH -p p.urania
#SBATCH -o spectre.out
#SBATCH -e spectre.out
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=36
#SBATCH --no-requeue
#SBATCH --mail-type=all
#SBATCH --mail-user=clemens.dittmer@aei.mpg.de
# Distributed under the MIT License.
# See LICENSE.txt for details.
```

This commands block sets up the specifications of the job for SLURM to execute. Important command lines to observe are:
```bash
#SBATCH -J {Job Name}
```
which specifies the job name.
```bash
#SBATCH --nodes {Number of Nodes}
```
which specifies the number of cluster nodes to run the program on. The exists a maximum for 2 nodes for the debug queue and 82 nodes for the main queue.
```bash
#SBATCH -t {HH:MM:SS}
```
which specifies the amount of time the job is allowed to run on the nodes. For urania, we can only set times up to 24:00:00 and 04:00:00 for the main queue and debug queue respectively.
```bash
#SBATCH -p {Sub-cluster}
```
which specifies which cluster to run the job on. For urania, we have two sub-clusters. The main urania queue (labelled as p.urania), and the debug queue (labelled as p.debug). For tinkering, it is recommended to work on the debug queue until you figured out an appropriate configuration that works for your simulation before pushing for higher run times with greater number of nodes on the main queue.
```bash
#SBATCH --mail-type=all
#SBATCH --mail-user={E-Mail Address}
```

There are also optional lines of code that allows the cluster to notify you via email about the starts, terminations, and completion of your jobs. If one wishes not to be notified, they may simply omit those lines of code from the .sh file.

2. 

#### Python Bindings

## Theoretical Background

### The Worldtube

### The Self-Force

### Coordinate System

### Trajectories

### Scalar Fields

## Scalar Gauss-Bonnet Field

$S[\Psi, g_{a,b}]=\int -(\Nabla\Psi)^2 + f(\Psi)(Riemann^2)$
Shift Sym: $f(\Psi)=\lambda \Psi$
Spon. Scal: $f(\Psi)=\nu \Psi^2+\gamma \Psi^4$
Dilatonic: f(\Psi)=e^{\lambda\Psi}\approx \lambda +\lambda \Psi$

Such that,
$\box\phi=-f'(\curlypsi)\mathcal{G}$
