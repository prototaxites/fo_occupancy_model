FROM rocker/rstudio:4.3.3
ARG TARGETARCH

RUN apt-get update && apt-get install -y \ 
    libicu-dev \
    zlib1g-dev \
    make \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    pandoc \
    libxml2-dev \
    libgdal-dev \
    gdal-bin \
    libgeos-dev \
    libproj-dev \
    libsqlite3-dev \
    libssl-dev \
    libudunits2-dev \
    pandoc-citeproc \
    libglpk-dev \
    libgmp3-dev \
    curl

RUN mkdir /home/user
ARG QUARTO_VERSION="1.5.36"
RUN curl -o quarto.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb
RUN dpkg -i quarto.deb
COPY container/docker_build.R /home/user/docker_build.R
RUN Rscript '/home/user/docker_build.R'
RUN mkdir /opt/cmdstan
ENV CMDSTAN_VERSION="2.34.1"
RUN Rscript -e "cmdstanr::install_cmdstan(dir = '/opt/cmdstan', version = '$CMDSTAN_VERSION')"
RUN chmod -R 777 /opt/cmdstan
ENV CMDSTAN="/opt/cmdstan/cmdstan-$CMDSTAN_VERSION"
