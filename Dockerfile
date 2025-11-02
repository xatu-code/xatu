FROM ubuntu:22.04
RUN apt-get update && apt-get install -y \
    build-essential \
    libopenblas-dev \
    liblapack-dev \
    libarmadillo-dev \
    libarpack2-dev \ 
    gfortran && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/xatu
COPY . /opt/xatu
RUN make clean
RUN make build TEST=1

WORKDIR /opt/xatu/test
RUN make tests

CMD ["./bin/tests.x"]