FROM --platform=linux/amd64 continuumio/miniconda3

LABEL author="Filip Ramazan"
LABEL version="v1.0"
LABEL description="BWTandem — BWT-based tandem repeat finder"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# gcc: the Cython accelerator compiles at install time, and the four ctypes
# C libraries compile at first import
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential procps \
    && rm -rf /var/lib/apt/lists/*

# The exact pinned core the published numbers were measured on
COPY environment.core.lock.yml /tmp/environment.core.lock.yml
RUN conda env update -n base -f /tmp/environment.core.lock.yml && conda clean -a -y

# Install BWTandem itself (builds the accelerator into the package)
COPY . /opt/bwtandem/
RUN pip install /opt/bwtandem

# Pre-compile the ctypes C libraries into the installed package so read-only
# runtimes (e.g. Singularity converting this image) do not need to write at
# first import
RUN python -m bwtandem.c_extensions.build

ENTRYPOINT ["bwtandem"]
