FROM ubuntu:22.04
RUN apt-get update && apt-get install -y \
    build-essential \
    libopenblas-dev \
    liblapack-dev \
    libarmadillo-dev \
    libarpack2-dev \ 
    gfortran \
    libgtest-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/xatu
COPY . /opt/xatu
RUN make clean
RUN make build TEST=1

WORKDIR /opt/xatu/test
RUN make tests

ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV OPENBLAS_CORETYPE=generic

CMD ["./bin/tests.x", "--gtest_color=yes"]