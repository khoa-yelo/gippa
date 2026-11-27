#!/bin/bash
#SBATCH --job-name=gippa-campaign-build-trfp
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000M
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Builds the ThermoRawFileParser micromamba env (validation/campaign/envs/thermorawfileparser.yml,
# ADR-0077 decision 2) and writes its explicit lock next to the yml. Usage: sbatch build_trfp_env.sbatch
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
HERE="${CAMPAIGN_REPO:-$GIPPA_REPO/local/worktrees/campaign-ecoli}/validation/campaign"
YML="$HERE/envs/thermorawfileparser.yml"
PREFIX="$MAMBA_ROOT_PREFIX/envs/gippa-thermorawfileparser"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
rm -rf "$PREFIX"
micromamba create -y -p "$PREFIX" -f "$YML" --override-channels -c conda-forge -c bioconda
micromamba list -p "$PREFIX" --explicit --md5 >"$HERE/envs/locks/thermorawfileparser.explicit.txt"
micromamba run -p "$PREFIX" ThermoRawFileParser.sh --version || micromamba run -p "$PREFIX" thermorawfileparser --version
