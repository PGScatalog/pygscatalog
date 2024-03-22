FROM python:3.11 as builder

ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

WORKDIR /app

RUN pip install poetry

COPY poetry.lock pyproject.toml ./

RUN poetry install --no-root && rm -rf $POETRY_CACHE_DIR

FROM python:3.11-slim-bullseye

ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH"
    
COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}
    
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

ENV PATH="/venv/bin:${PATH}"
