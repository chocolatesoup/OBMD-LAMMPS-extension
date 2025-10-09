# Workflows

## compile LAMMPS

- compile_lammps.yml
Currently builds LAMMPS which is in the `code` subdir on top of `foss/2023b` EESSI version 2023.06.
If the build failes the files in the `/tmp/eb-*` will be coppied to the `RUNNER_TMP` so that they can be consulted after a failed run. This will make it easier to debug the builds.
