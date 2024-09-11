FROM debian:stable

RUN dpkg --add-architecture i386 && \
   apt-get update && \
   apt-get install -y \
      autoconf \
      git \
      libtool \
      libmpfr-dev \
      make \
      texinfo \
      gcc-multilib \
      libgmp-dev:i386 \
      libgmp10:i386 \
      libmpfr-dev:i386 \
      libmpfr6:i386

