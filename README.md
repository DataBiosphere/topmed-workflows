# Topmed Workflows

## About
The original pipelines were assembled and written by Hyun Min Kang (hmkang@umich.edu) and Adrian Tan (atks@umich.edu) 
at the [Abecasis Lab at the University of Michigan](https://genome.sph.umich.edu/wiki/Abecasis_Lab)

See the [variant calling pipeline](https://github.com/statgen/topmed_freeze3_calling) and [alignment pipeline](https://github.com/statgen/docker-alignment) repositories

## Installing dependencies on your local system

### 1. Cloud SDK (`gcloud`, `gsutil`)
If you are on Debian / Ubuntu, follow the instructions on [Cloud SDK](https://cloud.google.com/sdk/downloads#apt-get). 
After you execute `gcloud init` the installer asks you to log in and you should respond with `Y`, head to the provided URL, copy the code and past it to the prompt. After that it will ask you for the cloud project you want to use. Pick the (platform-dev-178517). I picked `us-west1-b` as the region.

#### Configuration and credentials file
```bash
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)"
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
gcloud auth login
```
After that run `gcloud auth application-default --help` and follow the instructions. Briefly, run
```bash
gcloud iam service-accounts create <pick-a-username>
gcloud iam service-accounts keys create key.json --iam-account=<the-username-you-just-picked>@<your-service-account-name>.iam.gserviceaccount.com
```

That should print something like 
```bash
created key [<some long integer>] of type [json] as [key.json] for [<username-you-picked>@<your-service-account-name>.iam.gserviceaccount.com]
```
You can check in the Google Cloud Platform console under IAM Service Accounts. That account you just created should be in the list.

Next create an environment variable that points to the file `key.json`:
```bash
export GOOGLE_APPLICATION_CREDENTIALS=key.json
```

#### Providing credentials to your application
To run workflows of data stored on `gcloud` you need to set an environment variable `GOOGLE_APPLICATION_CREDENTIALS`, which holds the path to the credentials file.

### 2. Broad's execution engine `cromwell`
`cromwell` is a Java executable and requires a Java Runtime Engine. Follow the instruction [here](http://cromwell.readthedocs.io/en/develop/tutorials/FiveMinuteIntro/) for a complete installation.

### 3. Dockstore
For Dockstore to run you need to install the [Java Runtime Engine](https://www.digitalocean.com/community/tutorials/how-to-install-java-with-apt-get-on-ubuntu-16-04). Find installation instructions for Dockstore [here](https://dockstore.org/onboarding) (you need to be logged in to Dockstore).

## Running workflows

### Provisioning reference files
To copy contents of a SDK bucket to your local system (or a VM) use
```bash
gsutil -u [PROJECT_ID] cp gs://[BUCKET_NAME]/[OBJECT_NAME] [OBJECT_DESTINATION]
```

### Checker workflows


A WDL and a JSON file to test checker workflows are in the `test_data` directory. You need to adjust all paths in the JSON file to the paths on your system before running the checker. It has been tested with `cromwell-31.jar`. To run the checker workflow for the WDL aligner navigate to respective directory (usually it has _checker_ in its name) and run
```bash
java -Dconfig.file=<location_to_file> -jar ~/bin/<cromwell_version>.jar run <checker-workflow>.wdl -i  <checker-workflow>.json
```



