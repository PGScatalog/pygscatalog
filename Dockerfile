FROM python:3.11 as builder

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache 

WORKDIR /app

RUN pip install poetry

COPY pgscatalog.core /app/pgscatalog.core

COPY pgscatalog.calc /app/pgscatalog.calc

COPY pgscatalog.match /app/pgscatalog.match

COPY pgscatalog.utils /app/pgscatalog.utils

WORKDIR /app/pgscatalog.utils

RUN poetry install --no-root && rm -rf $POETRY_CACHE_DIR

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

ENV PATH="/app/pgscatalog.utils/.venv/bin:$PATH"


