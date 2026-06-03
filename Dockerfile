# CASPA — reproducible Linux image with the full Python + R/Bioconductor stack.
#
#   docker build -t caspa .
#   docker run --rm -v /host/MyExperiment:/work caspa --workdir /work --cores 8
#   docker run --rm -p 8501:8501 caspa gui        # the setup GUI at http://localhost:8501
#   docker run --rm caspa doctor                  # verify the install
#
# On Linux bioconda provides everything, so the image is just the conda
# environment from environment.yml. Windows users can run CASPA this way via
# Docker Desktop / WSL2 instead of installing the native stack.
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.source="https://github.com/vonkriegsheim/CASPA" \
      org.opencontainers.image.description="Context-Aware Single Cell Proteomic Analysis"

WORKDIR /opt/caspa

# Build the environment first so this expensive layer is cached across code changes.
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Put the env's python / snakemake / Rscript first on PATH (the Snakefile shells
# out to the bare commands), so no `conda activate` is needed at runtime.
ENV PATH=/opt/conda/envs/caspa/bin:$PATH

COPY . .

# Fail the build if any dependency is missing.
RUN python caspa/doctor.py

RUN chmod +x /opt/caspa/docker-entrypoint.sh
ENTRYPOINT ["/opt/caspa/docker-entrypoint.sh"]
CMD ["--help"]
