FROM condaforge/miniforge3:25.3.1-0 AS builder

# https://docs.astral.sh/uv/guides/integration/docker/#optimizations
ENV UV_COMPILE_BYTECODE=1
ENV UV_LINK_MODE=copy
ENV PYTHONUNBUFFERED=True
ENV PACKAGE_NAME=pgscatalog.utils
WORKDIR /workspace

COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/

COPY pyproject.toml uv.lock /workspace/

# dependencies for building bgen python package
RUN apt-get update && apt-get install -y cmake zlib1g-dev build-essential

# first use uv to grab dependencies
# --no-editable is important
RUN uv sync --frozen --no-install-workspace --package $PACKAGE_NAME --no-editable

COPY . /workspace

# now install the package itself
# --no-editable is important
RUN uv sync --frozen --package $PACKAGE_NAME --no-editable

FROM condaforge/miniforge3:25.3.1-0

# use a new runtime container to get rid of build time bloat
# just grabbing bgen apps + 7zip + the venv
RUN conda install -y -c bioconda bgen-cpp p7zip

COPY --from=builder /workspace/.venv /workspace/.venv

# nextflow dependency
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

ENV PATH="/workspace/.venv/bin:$PATH"