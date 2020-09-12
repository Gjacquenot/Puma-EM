FROM debian:buster-slim

RUN apt-get update \
 && apt-get install --yes --no-install-recommends\
    autoconf \
    automake \
    cvs \
    g++ \
    gfortran \
    gmsh \
    libtool \
    make \
    openmpi-bin \
    python-dev \
    python-matplotlib \
    python-mpi4py \
    python-scipy \
    python-tk \
    libopenmpi-dev \
    dvipng \
    wget

RUN wget --no-check-certificate https://sourceforge.net/projects/blitz/files/blitz/Blitz%2B%2B%200.10/blitz-0.10.tar.gz/download -O blitz.tar.gz && \
    mkdir blitz_src && \
    cd blitz_src && \
    tar -xzf ../blitz.tar.gz --strip 1 && \
    ./configure && \
    make lib && \
    make install && \
    rm -rf  blitz_src && \
    rm -rf blitz.tar.gz

# #### ADD DEFAULT USER ####
# ARG USER=mpi
# ENV USER ${USER}
# RUN adduser ${USER} \
#       && echo "${USER}   ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
#
# ENV USER_HOME /home/${USER}
# RUN chown -R ${USER}:${USER} ${USER_HOME}
#
# WORKDIR /puma-em_src
# ADD . /puma-em_src
# RUN ls && make libs CFLAGS="-c -O3 -fPIC -pthread -march=native -mfpmath=both"
#
# USER ${USER}
# RUN make
# CMD ["/bin/bash"]