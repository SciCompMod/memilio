FROM ubuntu:18.04
MAINTAINER Martin Siggel (martin.siggel@dlr.de)
RUN apt-get -qq update
RUN apt -qq -y install make wget gnupg cmake g++ git doxygen ninja-build python3.7-dev
RUN update-alternatives --install /usr/local/bin/python python /usr/bin/python3.7 3
RUN update-alternatives --install /usr/local/bin/python3 python3 /usr/bin/python3.7 3
