#AlignQC Development Environment
FROM ubuntu:16.04
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
               python-pip \
               r-base \
               nano \
               wget \
               git \
    && apt-get autoremove \
    && apt-get clean

RUN pip install --upgrade pip
RUN pip install seq-tools==1.0.10

VOLUME /temp
VOLUME /root

RUN pip install AlignQC==2.0.5

ENV HOME /root
WORKDIR /root
