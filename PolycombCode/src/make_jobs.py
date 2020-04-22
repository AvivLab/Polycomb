# This script creates a set of kubernetes jobs for the scripts.
import os

with open("../../job_template.yaml") as f:
    template = f.read()

dirs = os.listdir("../scripts")

for dirname in dirs:
    job_name = dirname.rstrip(".sh")
    job = template.replace("YYY", job_name)
    clean_job_name = job_name.replace("_","-")
    job = job.replace("XXX", clean_job_name)
    with open("../../jobs/{0}.yaml".format(clean_job_name), "w") as j:
        j.write(job)
