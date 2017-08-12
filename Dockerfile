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
RUN pip install seq-tools==1.0.4

VOLUME /temp
VOLUME /root

RUN mkdir /source

COPY . /source/AlignQC
RUN cd /source/AlignQC && pip install .

ENV HOME /root
WORKDIR /root

#CMD ["/usr/local/bin/alignqc"]
