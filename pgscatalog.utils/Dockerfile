FROM python:3.11-slim-bullseye

ARG VERSION=1.0.2

RUN apt-get update && apt install -y procps

RUN pip install pgscatalog-utils==${VERSION}
