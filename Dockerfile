FROM python:3.8

RUN apt-get update && apt-get install -y \
  cmake

RUN mkdir /app/
WORKDIR /app

# ----------------
# Install Gurobi
# ----------------
RUN wget https://packages.gurobi.com/9.0/gurobi9.0.2_linux64.tar.gz -O gurobi9.0.2_linux64.tar.gz
RUN tar xvzf gurobi9.0.2_linux64.tar.gz
RUN (cd gurobi902/linux64/src/build && make)
RUN (cd gurobi902/linux64/lib && ln -f -s ../src/build/libgurobi_c++.a libgurobi_c++.a)
ENV GUROBI_HOME /app/gurobi902

# ----------------
# Install SAMtools
# ----------------
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2/download -O samtools-1.7.tar.bz2
RUN tar xvjf samtools-1.7.tar.bz2
RUN (cd samtools-1.7 && ./configure && make)
ENV HATCHET_PATHS_SAMTOOLS /app/samtools-1.7

# ----------------
# Install BCFtools
# ----------------
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/bcftools-1.7.tar.bz2/download -O bcftools-1.7.tar.bz2
RUN tar xvjf bcftools-1.7.tar.bz2
RUN (cd bcftools-1.7 && ./configure && make)
ENV HATCHET_PATHS_BCFTOOLS /app/bcftools-1.7

# ----------------
# Copy source
# ----------------
COPY setup.py /app
COPY CMakeLists.txt /app
COPY FindGUROBI.cmake /app
COPY MANIFEST.in /app
COPY src/ /app/src/
COPY tests /app/tests/

# ----------------
# Install package
# ----------------
RUN CXXFLAGS=-pthread pip install .
