FROM mambaorg/micromamba:1.5

COPY env/environment.yml /tmp/environment.yml
RUN micromamba create -n bioport -f /tmp/environment.yml -y \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
