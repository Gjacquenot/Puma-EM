# docker build . -t pumaemdoc
FROM debian:9-slim
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -q && \
    apt-get install -q --yes --no-install-recommends \
        doxygen \
        make \
        lmodern \
        texlive-base \
        texlive-fonts-recommended \
        texlive-latex-base \
        texlive-latex-extra \
        texlive-latex-recommended \
        && \
    rm -rf /var/lib/apt/lists/*
