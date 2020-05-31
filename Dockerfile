FROM rocker/verse

ARG R_VERSION
ARG BUILD_DATE
ARG CRAN
ENV BUILD_DATE ${BUILD_DATE:-2020-05-20}
ENV R_VERSION=${R_VERSION:-3.6.3} \
    CRAN=${CRAN:-https://cran.rstudio.com} \ 
    TERM=xterm

# Set the locale
#RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
#    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8   


RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN  apt-get update \
	&& DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends \
  apt-utils \
  pandoc-citeproc lmodern

RUN R -e 'BiocManager::install("Icens")'
  
RUN install2.r --error \
    magrittr \ 
    interval \ 
    cgam \ 
    icenReg \ 
    km.ci \ 
    flexsurv \ 
    caret \ 
    tableone \ 
    PropCIs
    
## add git clone and initial exec here - maybe a dummy compilation that pulls all the tex packages
WORKDIR /n95_refit

RUN git clone git@github.com:cryanking/n95_refit.git /n95_refit


ENTRYPOINT ["/usr/local/bin/R"]

CMD ["-e 'rmarkdown::render(\"n95_report.Rmd\") '"]
