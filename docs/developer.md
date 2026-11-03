# Developer guide

The local executor, the test layers, FarmShare notes and the [release process](release.md).

!!! note "Stub"
    Filled by M11-08 (the user guide and tutorials).

## FarmShare (the developers' site)

The `farmshare` and `farmshare_local` Nextflow profiles are the maintainers' site profiles for
Stanford FarmShare; `gippa run --profile auto` picks `farmshare,apptainer` when the Slurm cluster
name is `farmshare`. `farmshare` submits tasks to Slurm and uses the prebuilt tool environments
under `$GIPPA_ROOT`; `farmshare_local` uses the same environments but runs all tasks inside the
current Slurm allocation. The user pages do not describe them.

On FarmShare, GIPPA runs from a repository checkout with prebuilt environments
(`AGENTS.md › Setup`). Inside a Slurm allocation, never on a login node:

```bash
source scripts/farmshare_env.sh            # sets GIPPA_ROOT, NXF_APPTAINER_CACHEDIR=$GIPPA_ROOT/containers
bash scripts/bootstrap_farmshare.sh        # once: uv, Nextflow, nf-test, gh in user space
uv sync --frozen --all-extras
sbatch scripts/slurm/build_envs.sbatch <tool> [...]   # micromamba envs and Apptainer pulls
sbatch scripts/slurm/gippa_setup.sbatch     # the databases
sbatch scripts/slurm/gippa_run.sbatch <args>           # -profile farmshare,apptainer
```

`build_envs.sbatch` builds `$GIPPA_ROOT/micromamba/envs/gippa-<tool>/` with its `.gippa-built`
marker, including the conda-only tools and UniOP, and pulls every pinned image into
`$GIPPA_ROOT/containers`. Under `farmshare,apptainer`, the conda-only tools use their prebuilt
environments and the image-pinned tools use their images: `workflow/conf/container_common.config`,
shared by the `apptainer` and `docker` profiles, sets `conda = ''` for every image-pinned process,
so an environment never shadows an image. Under a `farmshare*` profile, `gippa setup` prints no
UniOP terms notice, since the prebuilt env is used.

**`$GIPPA_ROOT`.** When set, it is a fallback in three places: the Apptainer image cache is
`NXF_APPTAINER_CACHEDIR`, else `$GIPPA_ROOT/containers`, else `DIR/containers`; setup's
`ToolRunner` looks for `gippa-bakta` / `gippa-antismash` under `$MAMBA_ROOT_PREFIX/envs/`, default
`$GIPPA_ROOT/micromamba`; and preflight looks for the prebuilt envs under it. To reproduce a plain
user install on FarmShare, `unset NXF_APPTAINER_CACHEDIR GIPPA_ROOT` first.

The pilot (`docs/get-started/from-spectra.md`, run times) ran its steps with
`scripts/slurm/diann_pilot.sbatch` (DIA-NN), `scripts/slurm/mzml_pilot.sbatch` and
`scripts/slurm/denovo_pilot.sbatch` (de novo), `scripts/slurm/diann_convert.sbatch` and
`scripts/slurm/cascadia_convert.sbatch` (conversion), and `scripts/slurm/gippa_run.sbatch`; their
Slurm job IDs are in `validation/pilot/`. The runtime benchmark (`docs/user-guide/runtimes.md`)
comes from `scripts/slurm/benchmark_standard.sbatch` on FarmShare `wheat` nodes. The platform
matrix's `nightly` and `weekly` suites are `scripts/slurm/nightly.sbatch` and
`scripts/slurm/weekly.sbatch` (`-profile farmshare_local,apptainer`).
