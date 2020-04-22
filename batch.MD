Batch simulations were performed using Kubernetes.

# 1) Containerization

The Dockerfile specifies an image with all of the scripts in this repository.

# 2) Scripts

There is a python script called `make_scripts.py` in the `src/` folder. It
converts each line in `src/polycomb_var.txt` into a bash script which performs
the command in that line and then uploads the results to Google Cloud Storage.
See the `scripts/` folder to see what these look like.

# 3) Create Batch Jobs

There is a python script called `make_jobs.py` in the `src/` folder. It converts
each script in the `scripts/` folder into a kubernetes job-spec. These are
contained in the `jobs/` folder. Each of these jobs can then be submitted to a
kubernetes cluster to run in parallel.

    kubectl apply -f jobs/
