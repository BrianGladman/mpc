FROM debian:stable

RUN apt-get update && apt-get install -y \
   autoconf \
   git \
   libtool \
   libmpfr-dev \
   make \
   texinfo

