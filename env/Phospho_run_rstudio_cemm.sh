#!/bin/bash
#SBATCH --job-name=rstudio_apptainer_interactive
#SBATCH --partition=interactiveq
#SBATCH --qos=interactiveq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --output=rstudio_apptainer_interactive_%j.log

set -ueo pipefail

# How to use the script - run it with sbatch command
# Other common SLURM variables https://docs.hpc.shef.ac.uk/en/latest/referenceinfo/scheduler/SLURM/SLURM-environment-variables.html#gsc.tab=0
echo "======================"
echo "Working directory:     $SLURM_SUBMIT_DIR"
echo "Job name:              $SLURM_JOB_NAME"
echo "Job id:                $SLURM_JOB_ID"
echo "Job queue (partition): $SLURM_JOB_PARTITION"
echo "Job nodes, tasks, CPUs: $SLURM_JOB_NUM_NODES, $SLURM_NTASKS, $SLURM_NPROCS"
echo "Job node name:         $SLURM_NODELIST"
echo "Job node address:      $(nslookup $(hostname) | grep Name: | cut -f2)"
echo "Job node IP address:   $(nslookup $(hostname) | grep Address: | tail -1 | cut -d' ' -f2)"
echo "Node allocated CPUs:   $SLURM_CPUS_ON_NODE"
echo "======================"

module load apptainer/1.1.9
workdir=$(pwd -P)
rstudio_server_config_dir="$(pwd)/.rstudio_server"
r_version="4.4"
r_apptainer_img="envs/ccribioinf_dockrstudio_4.4-java.sif"

mkdir -p -m 700 "${rstudio_server_config_dir}/run" "${rstudio_server_config_dir}/tmp" "${rstudio_server_config_dir}/var/lib/rstudio-server" \
	"${rstudio_server_config_dir}/R/${r_version}"

# R Session Configuration File https://docs.posit.co/ide/server-pro/reference/rsession_conf.html
cat >"${rstudio_server_config_dir}/rsession.conf" <<END
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with personal libraries from any R installation in the host environment
r-libs-user=R/$r_version
# Prevent R session from timeout
session-timeout-minutes=0

END

export APPTAINER_BIND="${rstudio_server_config_dir}/run:/run,\
  ${rstudio_server_config_dir}/tmp:/tmp,\
  ${rstudio_server_config_dir}/rsession.conf:/etc/rstudio/rsession.conf,\
  ${rstudio_server_config_dir}/var/lib/rstudio-server:/var/lib/rstudio-server,\
  ${rstudio_server_config_dir}/run:/var/run,\
  /nobackup:/nobackup,\
  /research:/research,\
  ${rstudio_server_config_dir}:${rstudio_server_config_dir},\
  ${workdir}/results/phosphoproteomics_analysis_R/renv:/renv/cache,\
  ${workdir}:/home/${USER}"

# Do not suspend idle sessions
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0
export APPTAINERENV_USER="${USER}"
export APPTAINERENV_PASSWORD="${RSTUDIO_PASSWORD:?Set RSTUDIO_PASSWORD before submitting the job}"
# Set renv cache path for the container
export APPTAINERENV_RENV_PATHS_CACHE="/renv/cache"

# Get unused socket between 8000 and 9000 (these are accessible within the CCRI network):
readonly port=$(python -c '
import socket
import random

def find_port_in_range(start=8000, end=9000):
    while True:
        port = random.randint(start, end)
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("", port))
                return port
            except OSError:
                continue

print(find_port_in_range())
')
readonly hostname=$(hostname)

cat 1>&2 <<END

To access the server, copy and paste this URL into your web browser (cmd + click for Mac users):
http://${hostname}.int.cemm.at:${port}

END

echo "======================"
echo "Job started at: $(date)"

apptainer exec \
	--cleanenv ${r_apptainer_img} \
	rserver --www-port ${port} \
	--server-user=${USER} \
	--auth-none=0 \
	--auth-pam-helper-path=pam-helper

echo "Job finished at: $(date)"

echo "Job stats:"
seff $SLURM_JOB_ID
echo "======================"
