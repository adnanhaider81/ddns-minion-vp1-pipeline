FROM mambaorg/micromamba:1.5.10

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

WORKDIR /workspace
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace

ENV PATH="/opt/conda/bin:${PATH}"
CMD ["/bin/bash"]
