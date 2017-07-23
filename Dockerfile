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

VOLUME /data
VOLUME /temp
Volume /output

ENV HOME /root

WORKDIR /root

RUN git clone https://github.com/jason-weirather/py-seq-tools.git
RUN cd /root/py-seq-tools && pip install -e .

RUN git clone https://github.com/jason-weirather/AlignQC.git
RUN cd /root/AlignQC && git checkout dev && pip install -e .

#CMD ["bash"]
