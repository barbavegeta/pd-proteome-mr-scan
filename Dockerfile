FROM python:3.11-slim
LABEL org.opencontainers.image.title="pd-proteome-mr-scan"
LABEL org.opencontainers.image.description="Proteome-wide cis-MR and colocalisation scan against Parkinson's disease"
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt
COPY pyproject.toml /opt/pdmr/
COPY src /opt/pdmr/src
RUN pip install --no-cache-dir /opt/pdmr
