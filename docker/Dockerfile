# Tosca image 2.0.1dev: amclab/tosca:2.0.0dev plus toscatools v1.0.0.9002.
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
      description="Docker image containing all software requirements for Tosca (2.0.0dev + toscatools v1.0.0.9002)"

# toscatools pinned to a release tag (not a moving branch); upgrade = "never" leaves the R
# packages already in the image alone; the build fails if the installed version does not match
ARG TOSCATOOLS_VERSION=1.0.0.9002
RUN R -e "remotes::install_github('amchakra/toscatools', ref = 'v${TOSCATOOLS_VERSION}', upgrade = 'never')" && \
    R -e "stopifnot(as.character(packageVersion('toscatools')) == '${TOSCATOOLS_VERSION}')"
