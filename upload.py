# Imports the Google Cloud client library
from google.cloud import storage
from typing import List, Tuple, Union
import io
import argparse
from tqdm import tqdm

# Instantiates a client
try:
    storage_client = storage.Client()
except Exception as e:
    print("{0}\n\n Try setting the environment variable GOOGLE_APPLICATION_CREDENTIALS to point to the file containing your service account key. \
        ie. try running 'export GOOGLE_APPLICATION_CREDENTIALS=(pwd)/gcp_service_account.json'".format(e))
    raise e

import os

def upload_folder_to_path(bucket_name, folder_name, path):
    """ Upload the content at the chosen path to the folder in the bucket. """
    bucket = storage_client.get_bucket(bucket_name)
    depth = len(path.split('/'))
    stripped_path = path.split('/')[-1]
    if os.path.isfile(path):
        # Upload just this file
        blob = bucket.blob(os.path.join(folder_name, stripped_path))
        blob.upload_from_filename(path)
    elif os.path.isdir(path):
        # Traverse folder and upload files
        for r, d, f in os.walk(path):
            for filename in f:
                full_filename = os.path.join(r, filename) # Path to file on disk
                base = os.path.join(*r.split('/')[depth-1:]) # Strip away preceding foldernames
                relative_filename = os.path.join(folder_name, base, filename) # Path to file in bucket
                print("Uploading {0}".format(full_filename))
                blob = bucket.blob(relative_filename)
                blob.upload_from_filename(full_filename)
    else:
        raise ValueError("The provided path does not point to a file or directory: {0}".format(path))

    print("Uploaded all files in {0} for bucket {1} under folder {2}".format(path, bucket_name, folder_name))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Uploads files from a local path to a bucket in GCS at target folder.")
    parser.add_argument('bucket', type=str, help="Name of bucket.")
    parser.add_argument('folder', type=str, help="Name of folder in bucket to download.")
    parser.add_argument('path', type=str, help="Local path to download the files into.")
    args = parser.parse_args()
    
    upload_folder_to_path(args.bucket, args.folder, args.path)