# Imports the Google Cloud client library
from google.cloud import storage
from typing import List, Tuple, Union
import io
import os
import argparse
from tqdm import tqdm

# Instantiates a client
try:
    storage_client = storage.Client()
except Exception as e:
    print("{0}\n\n Try setting the environment variable GOOGLE_APPLICATION_CREDENTIALS to point to the file containing your service account key. \
        ie. try running 'export GOOGLE_APPLICATION_CREDENTIALS=(pwd)/gcp_service_account.json'".format(e))
    raise e

def download_folder_to_path(bucket_name, folder, path, suffix=None):
    """ Downloads a folder hosted in a bucket to the chosen path. """
    bucket = storage_client.get_bucket(bucket_name)
    blobs = list(bucket.list_blobs(prefix=folder))
    if suffix:
        blobs = [b for b in blobs if b.name.endswith(suffix)]
    for blob in tqdm(blobs):
        filename = blob.name.split('/')[-1]
        if not os.path.isdir(path):
            raise ValueError("You must first create a folder at {0} before running this command.".format(path))
        with open(os.path.join(path, filename), 'wb') as f:
            print("Downloading {0}".format(blob.name))
            storage_client.download_blob_to_file(blob, f)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Downloads files from a target folder in a bucket in GCS to a local path.")
    parser.add_argument('bucket', type=str, help="Name of bucket.")
    parser.add_argument('folder', type=str, help="Name of folder in bucket to download.")
    parser.add_argument('path', type=str, help="Local path to download the files into.")
    args = parser.parse_args()
    
    download_folder_to_path(args.bucket, args.folder, args.path)
    print("Download all files in {0}:{1} to {2}".format(args.bucket, args.folder, args.path))
