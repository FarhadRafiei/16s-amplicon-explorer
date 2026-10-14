FROM python:3.11-slim

WORKDIR /project
COPY . .
RUN useradd --create-home analyst && chown -R analyst:analyst /project
USER analyst

ENTRYPOINT ["python3", "src/analyze_16s.py"]
