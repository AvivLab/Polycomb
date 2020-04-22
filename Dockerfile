FROM julia:1.1.1

WORKDIR /polycomb

COPY packages.jl /polycomb/packages.jl

# Install Julia Dependencies
RUN julia packages.jl

# Install Python Dependencies (for upload script)
RUN apt-get update && \
    apt-get install -y python3-pip

COPY requirements.txt /polycomb
RUN pip3 install -r requirements.txt

COPY upload.py /polycomb/upload.py

COPY PolycombCode/ /polycomb/

COPY gcp_service_account.json /polycomb/gcp_service_account.json

ENV GOOGLE_APPLICATION_CREDENTIALS /polycomb/gcp_service_account.json

ENTRYPOINT ["/bin/bash"]