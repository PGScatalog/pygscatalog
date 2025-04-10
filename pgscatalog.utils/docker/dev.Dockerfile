FROM python:3.12-slim-bullseye AS build

RUN apt-get update && apt install -y procps curl

COPY ../../pgscatalog.core /opt/pygscatalog/pgscatalog.core

COPY ../../pgscatalog.match /opt/pygscatalog/pgscatalog.match

COPY ../../pgscatalog.calc /opt/pygscatalog/pgscatalog.calc

COPY ../../pgscatalog.utils /opt/pygscatalog/pgscatalog.utils

COPY ../../pgscatalog.utils/docker/pyproject.dev /opt/pygscatalog/pgscatalog.utils/pyproject.toml

WORKDIR /opt/pygscatalog/pgscatalog.utils/

RUN pip install uv && uv build

FROM python:3.12-slim-bullseye

COPY --from=build /opt/pygscatalog/pgscatalog.utils/dist/ /opt/pgscatalog.utils

WORKDIR /opt/pgscatalog.utils 

RUN pip install *.whl
