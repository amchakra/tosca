# Tosca container

| Image | Built from | Contents |
|---|---|---|
| `amclab/tosca:2.0.0dev` | `Dockerfile.base` + `environment.yml` (2025) | conda environment `tosca-2.0.0dev`; toscatools from the `dev` branch at build time |
| `amclab/tosca:2.0.1dev` | `Dockerfile` | 2.0.0dev (pinned by digest) + toscatools `v1.1.0.9001` |

## Known issue: the conda environment does not rebuild reproducibly

**What goes wrong.** `environment.yml` is a specification that conda re-solves at every build, not a lock file. Rebuilding it gives a different environment once the channels have moved on, for three reasons:

- **A single `=` is a fuzzy match.** `samtools=1.22` accepts 1.22.x; a rebuild in 2026 gave 1.22.1. An exact version needs `==`.
- **Versions are pinned, builds are not.** `star=2.7.11b` accepts any build of 2.7.11b; a rebuild picked `h5ca1c30_8` instead of `_6`. Under x86 emulation, all of `_8`'s SIMD variants crashed with "Illegal instruction".
- **Only the top level is pinned.** The file pins about 30 top-level packages, while the environment has about 375. The rest are solved fresh, e.g. htslib went from 1.22 to 1.24, and around 40 other packages changed.

**For now:** 2.0.1dev is layered on the published 2.0.0dev image instead of being rebuilt.

**To do at the next pipeline/container update:** make the environment reproducible.
- **Lock file:** commit one, from `conda list --explicit --md5` or `conda-lock`, exported from a known-good image, and build the environment from it.
- **Or exact pins:** pin builds as well as versions (`pkg==version=build`).
- **R packages:** pin the ones installed outside conda (toscatools and any of its CRAN dependencies).
- **Image tags:** keep them immutable.
