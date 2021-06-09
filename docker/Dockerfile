# --------------------------disclaimer---------------------------------------------------------------------------
#This material was prepared as an account of work sponsored by an agency of the
#United States Government.  Neither the United States Government nor the United
#States Department of Energy, nor Battelle, nor any of their employees, nor any
#jurisdiction or organization that has cooperated in the development of these
#materials, makes any warranty, express or implied, or assumes any legal
#liability or responsibility for the accuracy, completeness, or usefulness or
#any information, apparatus, product, software, or process disclosed, or
#represents that its use would not infringe privately owned rights.
#
#Reference herein to any specific commercial product, process, or service by
#trade name, trademark, manufacturer, or otherwise does not necessarily
#constitute or imply its endorsement, recommendation, or favoring by the United
#States Government or any agency thereof, or Battelle Memorial Institute. The
#views and opinions of authors expressed herein do not necessarily state or
#reflect those of the United States Government or any agency thereof.
#
#                 PACIFIC NORTHWEST NATIONAL LABORATORY
#                              operated by
#                                BATTELLE
#                                for the
#                   UNITED STATES DEPARTMENT OF ENERGY
#                    under Contract DE-AC05-76RL01830
# ----------------------------------------------------------------------------------------------------------------
# -------------------------------------THANKS & CREDITS to:-------------------------------------------------------

#1. Proteowizard [Chambers, Matthew <matt.chambers42@gmail.com>]
#    ProteoWizard_description   :: "Convert MS RAW vendor files to open formats or analyze them with Skyline."
#    ProteoWizard_website       :: https://github.com/ProteoWizard/container
#    ProteoWizard_license       :: https://github.com/ProteoWizard/container
#    - NOTE: Use base Wine at least until the .NET version needs to be upgraded!
#               https://hub.docker.com/r/chambm/wine-dotnet
#               [ built of:https://github.com/ProteoWizard/container/blob/master/dotnet/Dockerfile]
#
#2. MoTrPAC  [David Jimenez-Morales <david.jimenez-morales@stanford.edu>]
#   https://github.com/MoTrPAC/motrpac-proteomics-pnnl-prototype
#3. Mono 5.16.0.220 [pgermishuys <https://hub.docker.com/u/pgermishuys>]
#   https://hub.docker.com/r/eventstore/eventstore-ci-ubuntu-16.04/dockerfile
#4. PNNL-Comp-Mass-Spec [Monroe, Matthew E <Matthew.Monroe@pnnl.gov>]
#   MSGFPlus                  :: https://github.com/MSGFPlus/msgfplus/releases/
#   MzidToTsvConverter        :: https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases/
#   PeptideHitResultsProcessor:: https://github.com/PNNL-Comp-Mass-Spec/PHRP/releases/
#   MASIC                     :: https://github.com/PNNL-Comp-Mass-Spec/MASIC/releases/
# ---------------------------------------------------------------------------------------------------------------
FROM chambm/wine-dotnet:4.7-x64

ARG BUILD_DATE
#ARG VCS_REF
ARG BUILD_VERSION

# Labels.
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.name="National Microbiome Data Collaborative: Meta-proteomics workflow"
LABEL org.label-schema.description="NMDC Meta-proteomics workflow \
                                    - This image is built on Ubuntu-16.04 and more information about the image and \
                                      how to run it's container can be found on ["org.label-schema.vcs-url"] \
                                    - Additionally, documentation[Examples, datasets etc] is at \
                                      https://nmdc-proteomics-workflow.readthedocs.io/en/latest/"
LABEL org.label-schema.url="https://microbiomedata.org/"
LABEL org.label-schema.vcs-url="https://github.com/PNNL-Comp-Mass-Spec/NMDC-Proteomics-Workflow"
#LABEL org.label-schema.vcs-ref=$VCS_REF
LABEL org.label-schema.vendor="Pacific Northwest National Laboratory-Integrative omics group"
LABEL org.label-schema.version=$BUILD_VERSION
LABEL maintainer="Anubhav <anubhav@pnnl.gov>"
LABEL tags="Proteomics,MassSpectrometry"

#                          ----- ProteoWizard :: msconvert :: Start -----
# .NET 4+ on Wine 64-bit with base image Ubuntu-16.04
# pwiz MSconvert :: Download Windows 64-bit tar.bz2 [ http://proteowizard.sourceforge.net/download.html ]
ADD pwiz-bin-windows-*.tar.bz2 /wineprefix64/drive_c/pwiz
ENV CONTAINER_GITHUB=https://github.com/ProteoWizard/container

ENV WINEDEBUG -all
ENV WINEPATH "C:\pwiz"
# sudo needed to run wine when container is run as a non-default user (e.g. -u 1234)
# wine*_anyuser scripts are convenience scripts that work like wine/wine64 no matter what user calls them
RUN apt-get update && \
    apt-get -y install sudo && \
    apt-get -y clean && \
    echo "ALL     ALL=NOPASSWD:  ALL" >> /etc/sudoers && \
    echo '#!/bin/sh\nsudo -E -u root wine64 "$@"' > /usr/bin/wine64_anyuser && \
    echo '#!/bin/sh\nsudo -E -u root wine "$@"' > /usr/bin/wine_anyuser && \
    chmod ugo+rx /usr/bin/wine*anyuser && \
    rm -rf \
      /var/lib/apt/lists/* \
      /usr/share/doc \
      /usr/share/doc-base \
      /usr/share/man \
      /usr/share/locale \
      /usr/share/zoneinfo

# create UIDs that Galaxy uses in default configs to launch docker containers; the UID must exist for sudo to work
RUN groupadd -r galaxy -g 1450 && \
    useradd -u 1450 -r -g galaxy -d /home/galaxy -c "Galaxy user" galaxy && \
    useradd -u 1000 -r -g galaxy -d /home/galaxy -c "Galaxy docker user" galaxy_docker && \
    useradd -u 2000 -r -g galaxy -d /home/galaxy -c "Galaxy Travis user" galaxy_travis && \
    useradd -u 999 -r -g galaxy -d /home/galaxy -c "usegalaxy.eu user" galaxy_eu

## Set up working directory and permissions to let user xclient save data
RUN mkdir /data
WORKDIR /data
CMD ["wine64_anyuser", "msconvert" ]

## If you need a proxy during build, don't put it into the Dockerfile itself:
## docker build --build-arg http_proxy=http://proxy.example.com:3128/  -t repo/image:version .

ADD mywine /usr/bin/
RUN chmod ugo+rx /usr/bin/mywine
#                          ----- ProteoWizard :: msconvert :: End -----

#                          ----- Java8 [ MSGFPlus ]-----
RUN apt-get update \
    && apt-get -y install gcc \
    && apt-get -y install -y git \
    && apt-get -y install openjdk-8-jre \
    && apt-get -y install tzdata \
    && apt-get -y install unixodbc unixodbc-dev freetds-dev freetds-bin tdsodbc \
    && apt-get clean

#                          ----- mono [ MASIC :: MzidToTsvConverter :: PeptideHitResultsProcessor ] -----
# install mono
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
    && apt install software-properties-common -y
RUN add-apt-repository ppa:git-core/ppa -y
RUN apt-get update \
    && apt-get install -y apt-transport-https curl unzip rake python sudo libcurl3
RUN curl -s -o packages-microsoft-prod.deb https://packages.microsoft.com/config/ubuntu/16.04/packages-microsoft-prod.deb \
    && sudo dpkg -i packages-microsoft-prod.deb
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
    && echo "deb http://download.mono-project.com/repo/ubuntu stable-xenial/snapshots/5.16.0.220 main" main | sudo tee /etc/apt/sources.list.d/mono-official.list \
    && apt-get update
RUN apt-get install -y mono-devel=5.16.0.220-0xamarin4+ubuntu1604b1 msbuild=1:16.0+xamarinxplat.2018.09.26.17.53-0xamarin3+ubuntu1604b1 dotnet-sdk-2.1

# MASIC layer
RUN wget -P /app/masic/ https://github.com/PNNL-Comp-Mass-Spec/MASIC/releases/download/v3.2.7465/MASIC_Console_Program.zip \
    && cd /app/masic/ \
    && unzip MASIC_Console_Program.zip \
    && rm MASIC_Console_Program.zip

# MSGFPlus layer
RUN wget -P /app/msgf https://github.com/MSGFPlus/msgfplus/releases/download/v2020.08.05/MSGFPlus_v20200805.zip \
    && cd /app/msgf/ \
    && unzip MSGFPlus_v20200805.zip \
    && rm MSGFPlus_v20200805.zip

# MzidToTSVConverter layer
RUN wget -P /app/mzid2tsv/ https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases/download/v1.4.3/MzidToTsvConverter.zip \
    && cd /app/mzid2tsv/ \
    && unzip MzidToTsvConverter.zip \
    && rm MzidToTsvConverter.zip

# install System.Data.SQLite
# donwload from https://system.data.sqlite.org/downloads/1.0.111.0/sqlite-netFx-full-source-1.0.111.0.zip
RUN wget -P /app/sqlite/ https://github.com/AshleyLab/motrpac-proteomics-pnnl-prototype/raw/master/step06/sqlite-netFx-full-source-1.0.111.0.zip \
    && cd /app/sqlite/ \
    && unzip sqlite-netFx-full-source-1.0.111.0.zip \
    && rm sqlite-netFx-full-source-1.0.111.0.zip \
    && chmod +x /app/sqlite/Setup/compile-interop-assembly-release.sh \
    && /app/sqlite/Setup/compile-interop-assembly-release.sh

# PeptideHitResultsProcessor layer
RUN wget -P /app/phrp/ https://github.com/PNNL-Comp-Mass-Spec/PHRP/releases/download/v1.5.7458/PeptideHitResultsProcessor.zip \
    && cd /app/phrp/ \
    && unzip PeptideHitResultsProcessor*.zip \
    && rm PeptideHitResultsProcessor*.zip \
    && cp /app/sqlite/bin/2013/Release/bin/libSQLite.Interop.so /app/phrp/ \
    && chmod -x /app/phrp/libSQLite.Interop.so

ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java
ENV MONO_HOME=/usr/bin/mono
ENV WINE_HOME=/usr/bin/wine

#                          ----- Download Meta-proteomics workflow -----
# to access DMS (PNNL-internal)
COPY freetds.conf /etc/freetds/freetds.conf
COPY odbc.ini /etc/odbc.ini
COPY odbcinst.ini /etc/odbcinst.ini

RUN cd /app/ \
    && git clone https://stash.pnnl.gov/scm/omcs/nmdc-proteomics-workflow.git

#                          ----- Miniconda3 :: Meta-proteomics workflow virtual environment :: Start -----
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENTRYPOINT [ "/usr/bin/tini", "--" ]
#                          ----- Miniconda3 :: Meta-proteomics workflow virtual environment :: END -----

RUN conda create --name proteomics_env python=3.7
SHELL ["conda", "run", "-n", "proteomics_env", "/bin/bash", "-c"]
RUN pip install pymssql pandas requests bs4 seaborn openpyxl
RUN echo "conda activate proteomics_env" >> ~/.bashrc
RUN cd /app/nmdc-proteomics-workflow/ \
    && pip install -e .

WORKDIR /app