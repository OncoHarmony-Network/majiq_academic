# This is the latest version of centos. Check the following link for other options:
# https://hub.docker.com/_/centos/
FROM centos:7.6.1810

# update yum
RUN yum -y update

# install Python3.6
# https://janikarhunen.fi/how-to-install-python-3-6-1-on-centos-7.html
RUN yum -y install https://centos7.iuscommunity.org/ius-release.rpm
RUN yum -y install python36u python36u-pip python36u-devel

# install needed python packages
RUN pip3.6 install --no-cache-dir -U pip
RUN pip3.6 install --no-cache-dir -U wheel setuptools
RUN pip3.6 install --no-cache-dir -U numpy cython

# install htslib
RUN yum -y install wget zlib-devel gcc make bzip2 bzip2-devel xz-devel libcurl-devel openssl-devel
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN tar -xvjf htslib-1.9.tar.bz2
WORKDIR htslib-1.9
RUN ./configure
RUN make && make install
WORKDIR /tmp
RUN rm -rf htslib-1.9
RUN rm htslib-1.9.tar.bz2

# install majiq
RUN yum -y install git gcc-c++
WORKDIR /tmp
RUN git clone https://bitbucket.org/biociphers/majiq.git
WORKDIR majiq
RUN python3.6 setup.py install
WORKDIR /tmp
RUN rm -rf majiq

# clean yum cache
RUN yum clean all

ENTRYPOINT ["majiq"]
CMD []

# To build image: docker build --squash -t majiq .
# To save image: docker save majiq -o majiq.tar
# To load image: docker load -i majiq.tar
# To run loaded image: docker run majiq -v
