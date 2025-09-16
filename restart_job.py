import os
import sys

rel_path = str(sys.argv[1])
cwd = os.getcwd()
abs_base_path = os.path.join(cwd, rel_path)
assert os.path.exists(abs_base_path), f"Path {abs_base_path} does not exist"

i = 1
while os.path.exists(os.path.join(abs_base_path, f"Restart{i}")):
    i += 1

checkpoint_folder = (
    os.path.join(abs_base_path, f"Restart{i-1}/Checkpoints/Checkpoint_0000")
    if i > 1
    else os.path.join(abs_base_path, "Checkpoints/Checkpoint_0000")
)
assert os.path.exists(
    checkpoint_folder
), f"Requested checkpoint folder {checkpoint_folder} does not exist"

restart_folder = os.path.join(abs_base_path, f"Restart{i}")
os.mkdir(restart_folder)

with open(os.path.join(abs_base_path, "urania.sh"), "r") as f:
    launch_script = f.readlines()
launch_script.append(f" +restart {checkpoint_folder}")

with open(os.path.join(restart_folder, "urania.sh"), "w") as f:
    f.writelines(launch_script)

os.chdir(restart_folder)
os.system("sbatch urania.sh")
