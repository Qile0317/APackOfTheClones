FROM rocker/verse:4.5.1
COPY DESCRIPTION ./DESCRIPTION
RUN Rscript -e "remotes::install_deps(dependencies = TRUE)"
RUN apt-get update && apt-get install -y texlive-fonts-extra