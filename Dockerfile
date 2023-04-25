
FROM ubuntu:20.04

# update yum
RUN apt-get -y update
RUN apt-get -y install libhts-dev python3-dev python3-pip git gcc


# install majiq

WORKDIR /tmp
RUN git clone https://bitbucket.org/biociphers/majiq_academic.git
WORKDIR majiq_academic
RUN pip3 install --no-cache-dir -U pip
RUN pip3 install --no-cache-dir ./voila
RUN pip3 install --no-cache-dir ./majiq
WORKDIR /tmp
RUN rm -rf majiq_academic


# remove these lines to go back to default majiq/voila usage (non CWL)
#COPY gen_majiq_cwl.py /opt/gen_majiq_cwl.py
#RUN chmod +x /opt/gen_majiq_cwl.py
#ENTRYPOINT ["/opt/gen_majiq_cwl.py"]
# end lines to remove

CMD []

# To build image: go to some temporary directory and $ docker build -t majiq_voila -f /path/to/this/Dockerfile .
# To run loaded image: docker run majiq_voila majiq --help
# To run loaded image: docker run majiq_voila -v /path/to/voila/data/files:/mnt -p 5010:5010 voila view /mnt --host 0.0.0.0 -p 5010 -j4
# note that you may not see the line "Serving on 0.0.0.0:5010" for voila view some reason, but it will work anyway
