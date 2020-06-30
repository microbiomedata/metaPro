#------------------MePro Processing------------------#

# .NET 4+ on Wine 64-bit with base image Ubuntu-16.04
FROM chambm/wine-dotnet:4.7-x64

# pwiz MSconvert
# download Windows 64-bit tar.bz2 from http://proteowizard.sourceforge.net/download.html
ADD pwiz-bin-windows-*.tar.bz2 /wineprefix64/drive_c/pwiz

ENV CONTAINER_GITHUB=https://github.com/ProteoWizard/container

LABEL description="Convert MS RAW vendor files to open formats or analyze them with Skyline."
LABEL website=https://github.com/ProteoWizard/container
LABEL documentation=https://github.com/ProteoWizard/container
LABEL license=https://github.com/ProteoWizard/container
LABEL tags="Metabolomics,Proteomics,MassSpectrometry"

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


RUN apt-get update \
    && apt-get -y install gcc \
    && apt-get -y install -y git \
    && apt-get -y install openjdk-8-jre \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
    && apt-get -y install apt-transport-https ca-certificates \
    && echo "deb https://download.mono-project.com/repo/ubuntu stable-xenial main" | sudo tee /etc/apt/sources.list.d/mono-official-stable.listn\
    && sudo apt-get update \
    && apt-get -y install mono-vbnc \
    && apt-get -y install mono-complete \
    && apt-get -y install unixodbc unixodbc-dev freetds-dev freetds-bin tdsodbc \
    && apt-get clean

# MASIC layer
RUN wget -P /app/masic/ https://github.com/PNNL-Comp-Mass-Spec/MASIC/releases/download/v3.0.7235/MASIC_Console_Program.zip \
    && cd /app/masic/ \
    && unzip MASIC_Console_Program.zip \
    && rm MASIC_Console_Program.zip

# MSGFPlus layer
RUN wget -P /app/msgf https://github.com/MSGFPlus/msgfplus/releases/download/v2019.06.28/MSGFPlus_v20190628.zip \
    && cd /app/msgf/ \
    && unzip MSGFPlus_v20190628.zip \
    && rm MSGFPlus_v20190628.zip

# MzidToTSVConverter layer
RUN wget -P /app/mzid2tsv/ https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases/download/v1.3.3/MzidToTsvConverter.zip \
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
RUN wget -P /app/phrp/ https://github.com/PNNL-Comp-Mass-Spec/PHRP/releases/download/v1.5.7130/PeptideHitResultsProcessor_Debian.zip \
    && cd /app/phrp/ \
    && unzip PeptideHitResultsProcessor*.zip \
    && rm PeptideHitResultsProcessor*.zip \
    && cp /app/sqlite/bin/2013/Release/bin/libSQLite.Interop.so /app/phrp/ \
    && chmod -x /app/phrp/libSQLite.Interop.so

ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java
ENV MONO_HOME=/usr/bin/mono
ENV WINE_HOME=/usr/bin/wine

#------------------mePro Analysis------------------#

# to access DMS
COPY freetds.conf /etc/freetds/freetds.conf
COPY odbc.ini /etc/odbc.ini
COPY odbcinst.ini /etc/odbcinst.ini

# copy entire analysis into the image.
#RUN wget -P /app/ https://stash.pnnl.gov/scm/omcs/nmdc-proteomics-workflow.git \
COPY ./nmdc-proteomics-workflow /app/nmdc-proteomics-workflow

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


RUN conda create --name proteomics_env python=3.7
SHELL ["conda", "run", "-n", "proteomics_env", "/bin/bash", "-c"]
RUN pip install pymssql pandas requests bs4 seaborn openpyxl
RUN echo "conda activate proteomics_env" >> ~/.bashrc
RUN cd /app/nmdc-proteomics-workflow/ \
    && pip install -e .

WORKDIR /app
