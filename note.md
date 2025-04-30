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
These commands also work recursively for entire directories. Bt specifying only the directory, changes for files within will be recursively staged for their respective changes.

### 4. Committing staged changes
```bash
git commit -m "{Message for your commit}"
```
Commits all staged changes to the remote branch. Other poeple are not abble to see your changes yet on the local repository for that push is required.

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

#### The bash file (.sh)

#### Python Bindings

## Theoretical Background

### The Worldtube

### The Self-Force

### Coordinate System

### Trajectories

### Scalar Fields

