FROM ubuntu:18.04
MAINTAINER Martin Siggel (martin.siggel@dlr.de)
RUN apt-get -qq update
RUN apt -qq -y install wget gnupg
RUN sh -c "echo 'deb http://download.opensuse.org/repositories/science:/dlr/xUbuntu_18.04/ /' > /etc/apt/sources.list.d/science:dlr.list"
RUN wget -nv https://download.opensuse.org/repositories/science:dlr/xUbuntu_18.04/Release.key -O Release.key
RUN apt-key add - < Release.key
RUN apt-get -qq update
RUN apt -qq -y install make cmake g++ git doxygen ninja-build python3.7-dev
RUN update-alternatives --install /usr/local/bin/python python /usr/bin/python3.7 3
RUN update-alternatives --install /usr/local/bin/python3 python3 /usr/bin/python3.7 3
