FROM julia:1.1.1

WORKDIR /polycomb

COPY packages.jl /polycomb/packages.jl

# Install Julia Dependencies
RUN julia packages.jl

# Install Python Dependencies (for upload script)
COPY requirements.txt /polycomb
RUN pip install -r requirements.txt

COPY upload.py /polycomb/upload.py

COPY PolycombCode/input/* /polycomb/input
COPY PolycombCode/src/* /polycomb/src

ENTRYPOINT ["julia"]