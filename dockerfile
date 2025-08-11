# Base image (Debian-based)
FROM python:3.8-slim

#  Forcer l'installation sans interaction
ENV DEBIAN_FRONTEND=noninteractive



# Installer les certificats et les outils système
RUN apt-get update && apt-get install -y \
    ca-certificates \
    curl \
    openssl \
    wget \
    sudo \
    libssl-dev \
    libbz2-dev \
    libexpat1-dev \
    libncurses5-dev \
    zlib1g-dev \
    bzip2 \
    build-essential \
    git \
    unzip \
    && rm -rf /var/lib/apt/lists/*





# Installer Miniconda
RUN curl  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda.sh \
    && bash Miniconda.sh -b -p /opt/conda \
    && rm -f Miniconda.sh

# Définir Conda dans le PATH
ENV PATH="/opt/conda/bin:$PATH"

# Mettre à jour Conda et les paquets de base
RUN conda update -n base -c defaults conda -y

# Créer un environnement Conda avec Python 3.10
RUN conda create -n bioinfo_env python=3.10 -y




# Installer Snakemake et Ensembl-VEP via Conda
RUN conda install --no-channel-priority --override-channels --insecure -n base -c conda-forge  \
    pandas \
    numpy \
    scipy \
    argparse\
	pysam \
    && conda clean --all -y || true


# Installer Snakemake et Ensembl-VEP via Conda
RUN conda install --no-channel-priority --override-channels --insecure -n base -c conda-forge  \
    scikit-learn \
	math \
	collections \
	pyro-ppl \
	pytorch \
	samtools \
    bcftools \
    jupyter \
    && conda clean --all -y || true





# Installer les plugins de VEP avec le bon `curl`
#RUN mkdir -p /workspace/vep_data \
#    && conda run -n base vep_install --AUTO a --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR /workspace/vep_data
#RUN mkdir -p /workspace/vep_data \
#    && bash -c "source /opt/conda/bin/activate base && vep_install --AUTO a --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR /workspace/vep_data --curl /usr/bin/curl" || true

COPY requirements.txt ./
RUN pip install --no-cache-dir \
    --trusted-host pypi.org \
    --trusted-host pypi.python.org \
    --trusted-host files.pythonhosted.org \
    -r requirements.txt
RUN pip install --no-cache-dir vcfpy
    


