FROM gitpod/workspace-full
RUN sudo apt-get install --yes libgsl0-dev

USER gitpod
ENV PATH=$PATH:$HOME/anaconda3
ENV PATH=$PATH:$HOME/anaconda3/bin
RUN wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
RUN bash Anaconda3-5.0.1-Linux-x86_64.sh -b
RUN rm Anaconda3-5.0.1-Linux-x86_64.sh
RUN conda create -c bioconda -c conda-forge -n rbt starcode bcftools