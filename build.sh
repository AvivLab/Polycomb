python3 PolycombCode/src/make_scripts.py # Generates scripts
python3 PolycombCode/src/make_jobs.py # Generates job configurations
docker build -t gcr.io/kellylab/polycomb . # Builds Docker image
docker push gcr.io/kellylab/polycomb # Pushes Docker image to remote