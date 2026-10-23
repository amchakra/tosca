# Tosca image 2.0.1dev: amclab/tosca:2.0.0dev plus toscatools v1.1.0.9001.
#
# Layered on the published 2.0.0dev image (pinned by digest) so that every tool and R package
# is byte-identical to 2.0.0dev and only toscatools changes (minus-strand convert_coordinates()
# fix; optional crosslink-distance filter in get_valid_hybrids()). Rebuilding the conda
# environment from Dockerfile.base would not reproduce 2.0.0dev: see docker/README.md.
#
# Build (the cluster is x86_64):
#   docker buildx build --platform linux/amd64 -t amclab/tosca:2.0.1dev --load docker/

FROM amclab/tosca@sha256:1a0d64eba4a2a8b30655b2c5251bf154de5241e6fb617dcb913ff8a02f385d52

LABEL authors="Anob Chakrabarti" \
      description="Docker image containing all software requirements for Tosca (2.0.0dev + toscatools v1.1.0.9001)"

# toscatools pinned to a release tag (not a moving branch); upgrade = "never" leaves the R
# packages already in the image alone; the build fails if the installed version does not match
ARG TOSCATOOLS_VERSION=1.1.0.9001
# (installed from the tag's source tarball, which needs no GitHub API calls and so is not rate-limited)
RUN R -e "remotes::install_url('https://github.com/amchakra/toscatools/archive/refs/tags/v${TOSCATOOLS_VERSION}.tar.gz', upgrade = 'never')" && \
    R -e "stopifnot(as.character(packageVersion('toscatools')) == '${TOSCATOOLS_VERSION}')"
