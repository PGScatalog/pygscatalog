FROM python:3.12-slim-bullseye

ENV UV_FROZEN=1

COPY ../../pgscatalog.utils /opt/pygscatalog/pgscatalog.utils

WORKDIR /opt/pygscatalog/pgscatalog.utils/

RUN pip install uv && uv sync --all-packages

RUN apt-get update && apt-get install -y procps

ENV PATH="/opt/pygscatalog/pgscatalog.utils/.venv/bin:$PATH"
