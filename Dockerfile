# This is the latest version of centos. Check the following link for other options:
# https://hub.docker.com/_/centos/
FROM centos:7.6.1810

# update yum
RUN yum -y update

# install Python3.6
# https://janikarhunen.fi/how-to-install-python-3-6-1-on-centos-7.html
RUN yum -y install yum-utils
RUN yum -y groupinstall development
RUN yum -y install https://centos7.iuscommunity.org/ius-release.rpm
RUN yum -y install python36u python36u-pip python36u-devel

# install htslib
RUN yum -y install wget zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN tar -xvjf htslib-1.9.tar.bz2
WORKDIR htslib-1.9
RUN ./configure
RUN make && make install

# install majiq
RUN pip3.6 install -U pip
RUN pip3.6 install -U wheel setuptools
RUN pip3.6 install -U numpy cython

# TODO: This line is assuming majiq 2.0 has been moved to the majiq_stable repo.
RUN pip3.6 install git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq
# This would be the command if majiq is still in the majiq repo.
#RUN pip3.6 install git+https://<username>:<password>bitbucket.org/biociphers/majiq.git#egg=majiq

ENTRYPOINT ["majiq"]
CMD []

# To build image: docker build -t majiq .
# To save image: docker save majiq -o majiq.tar
# To load image: docker load -i majiq
# To run loaded image: docker run majiq -v
