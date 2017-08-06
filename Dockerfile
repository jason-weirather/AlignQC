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

VOLUME /temp
VOLUME /root

RUN mkdir /source
RUN cd /source/ \
    && git clone https://github.com/jason-weirather/py-seq-tools.git
RUN cd /source/py-seq-tools && pip install -e .

RUN cd /source/ && git clone https://github.com/jason-weirather/AlignQC.git
RUN cd /source/AlignQC && git checkout dev && pip install -e .

ENV HOME /root
WORKDIR /root

#CMD ["/usr/local/bin/alignqc"]
