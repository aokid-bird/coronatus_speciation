# Building ngsDist Singularity image on Apple Silicon (via Lima x86_64)

This project runs `ngsDist` via Singularity/Apptainer. On Apple Silicon, build the SIF inside an x86_64 Linux VM using Lima, so it runs on the RedHat cluster.

Prerequisites on macOS:
- Homebrew installed
- Lima installed: `brew install lima`
- Apptainer inside the Lima VM (installed below)

Steps:

1) Start an x86_64 Lima VM
- Create and start a VM targeting amd64:
```
limactl start --name apptainer-x86 --arch x86_64 default
```
- Enter the VM:
```
limactl shell apptainer-x86
```

2) Install Apptainer in the VM
```
sudo apt-get update
sudo apt-get install -y wget git build-essential uuid-dev libseccomp-dev pkg-config squashfs-tools cryptsetup runc uidmap
export VER=1.3.1
wget https://github.com/apptainer/apptainer/releases/download/v${VER}/apptainer_${VER}_amd64.deb
sudo dpkg -i apptainer_${VER}_amd64.deb
apptainer --version
```

3) Build the SIF from the provided definition
```
cd /tmp
# Mount the project directory into the VM (if not mounted already)
# You can copy the def file into /tmp as an alternative.
mkdir -p singularity_out
apptainer build singularity_out/ngsdist.sif \
  /path/to/your/repo/workflow/singularity/ngsdist.def
```

4) Copy the SIF back to macOS host and place it under the directory configured by `singularity_dir` in your config YAML (default: `${HOME}/envs/singularity`). For example:
```
TARGET_DIR=${HOME}/envs/singularity
mkdir -p "$TARGET_DIR"
limactl copy apptainer-x86:/tmp/singularity_out/ngsdist.sif \
  "$TARGET_DIR/ngsdist.sif"
```

5) Configure the pipeline
- Ensure `singularity_dir` in `config/config_*.yaml` points to a directory synced to the cluster (e.g., `${HOME}/envs/singularity`).
- Place `ngsdist.sif` in that directory and sync to the cluster if needed.

Now the Snakemake rules will use `ngsdist.sif` on the cluster.
