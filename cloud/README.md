# Running HATCHet in the cloud

HATCHet is a Docerizable application and comes with a Dockerfile for easy deployment. We have also made HATCHet
available as a publicly accessible Docker image at the [Google Cloud Container Registry](https://cloud.google.com/container-registry).
This facilitates running HATCHet in the cloud without worrying about downloading large BAM files, and without having to
build and install HATCHet locally.

This README provides details on how to run HATCHet entirely on the [Google Cloud Platform](https://cloud.google.com) (GCP)
on large datasets made available at [ISB-CGC](https://isb-cgc.appspot.com/).

## Running HATCHet on ISB-CGC Datasets

### Setting up access at ISB-CGC

The ISB-CGC [documentation](https://isb-cancer-genomics-cloud.readthedocs.io/) is the definitive guide to the preliminary
steps required to access *Controlled Data*. In particular, the section on [How to link your NIH/eRA & Google identities](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/controlled-access/Controlled-data-Interactive.html)
is a good place to get started on setting up your access. After associating your Google identity to your NIH or eRA
identity, proceed to the [Accessing Controlled Data](https://isb-cancer-genomics-cloud.readthedocs.io/en/latest/sections/Gaining-Access-To-Controlled-Access-Data.html)
section and follow the steps to register your Google project with ISB-CGC.

Note that your PI will most likely have to grant you access to one or more of these controlled datasets using
[dbGap](https://dbgap.ncbi.nlm.nih.gov/). The steps in the walk-throughs and tutorials on the ISB-CGC website will
verify that you do have the appropriate access you will need to programmatically read these datasets in HATCHet.

Also note that access to controlled datasets is typically granted only for 24 hours, so you will have to extend your
access period on the ISB-CGC website if it has expired.

### Setting up your environment to run HATCHet on GCP

You do not need to build or install HATCHet locally, either as a python package or a Docker image. The only pre-requisite
is that you have installed the [Google Cloud SDK](https://cloud.google.com/sdk/docs/quickstart).

This is most cleanly done by installing all required dependencies inside a new Python 3 Conda environment.
Specifically, the commands you will need to execute are:

```
conda create --name hatchet_cloud python=3.8
conda activate hatchet_cloud
conda install -c conda-forge google-cloud-sdk google-api-python-client
pip install oauth2client dsub
```

### Logging in to your GCP Account

After installing the required dependencies, make sure that you login to your Google account and set up your default
project. These are **one time steps** to make sure that HATCHet is able to correctly talk to your project.

```
gcloud auth application-default login
```

This will open up a browser window asking you to authorize the Google SDK to communicate with your GCP account.

```
gcloud auth application-default set-quota-project PROJECT_ID
```

Replace `PROJECT_ID` with the project ID (this is the project ID without spaces, as against the project name that can have spaces)
that you linked with your ISB-CGC account.

### Preparing a bucket for output files

In the Google project that you used in the steps above, use the following command to create a new bucket where the results
of your HATCHet analysis will be saved:

```
gsutil mb gs://BUCKET_NAME
```

Replace `BUCKET_NAME` with a globally-unique bucket name. This step can also be performed by logging in to the
[Google Cloud Console](https://console.cloud.google.com) and navigating to Home -> Storage -> Browser -> Create Bucket.

### Fine-tuning the HATCHet script

The `_run.sh` script provided with HATCHet is an end-end worflow of HATCHet. This will be familiar to you if you have
run HATCHet locally. You can comment out sections of this script to only run certain parts of HATCHet depending on your
needs, and specify the values of certain flags of the pipeline.

The part of the script that you will want to pay attention to is the `Reference Genome` section. Depending on the
reference genome used in your BAM files, you will want to uncomment one of the lines in this section. If none of the
provided links to reference genomes are applicable in your case, you can add a line on your own, pointing to a `.gz`
or `.fa` file available through `wget`.

### Running the scripts

The `cloud_run.sh` script provided with HATCHet is a single [dsub](https://github.com/DataBiosphere/dsub) command that
will run HATCHet in the cloud. This command leverages the [Google Life Sciences API](https://cloud.google.com/life-sciences/docs/reference/rest)
and internally performs the following series of steps:

<a name="cloud_steps"></a>
1. Create a Compute Engine virtual machine
2. Download the Docker image
3. Download input files from Google Storage (buckets) to the virtual machine
4. Run a new Docker container with the specified image and command
5. Upload the output files to Google Storage (a bucket of your choice)
6. Destroy the Compute Engine virtual machine

As an example, the command looks as follows:

```
dsub \
  --provider google-v2 \
  --project PROJECT_ID \
  --regions us-east1 \
  --logging "gs://BUCKET_NAME/RUN_NAME/logging" \
  --disk-size 1000 \
  --name "RUN_NAME" \
  --image gcr.io/durable-tracer-294016/hatchet \
  --input NORMALBAM="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/b9774dd35c320f70de8f2b81c15d5a98.bam" \
  --input NORMALBAI="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/b9774dd35c320f70de8f2b81c15d5a98.bam.bai" \
  --input TUMORBAM1="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/2258e57e8e0af9db6969a1da86177ca7.bam" \
  --input TUMORBAI1="gs://gdc-tcga-phs000178-controlled/BRCA/DNA/WGS/WUGSC/ILLUMINA/2258e57e8e0af9db6969a1da86177ca7.bam.bai" \
  --output-recursive OUTPUT_FOLDER="gs://BUCKET_NAME/RUN_NAME/output" \
  --script "_run.sh" \
  --wait
```

In the above command, you will want to replace `PROJECT_ID` with your project id, `BUCKET_NAME` with the bucket name that
you created above, `RUN_NAME` with any unique name (no spaces!) that identifies your HATCHet run. In addition:

- The `NORMALBAM` parameter should be replaced with the `gs://..` path to the matched-normal sample of the patient.

- The `NORMALBAI` parameter should be replaced with the `gs://..` path to the index file (with `.bai` extension) for the `NORMALBAM` file.
This is typically just a file with `.bai` appended to the filename.

- The `TUMORBAM1` parameter should be replaced with the `gs://..` path to the first tumor sample of the patient.

- The `TUMORBAI1` parameter should be replaced with the `gs://..` path to the index file (with `.bai` extension) for the `TUMORBAM1` file.
This is typically just a file with `.bai` appended to the filename.

  For additional DNA sequencing reads for tumor samples, use input names like `TUMORBAM2`/`TUMORBAI2`, `TUMORBAM3`/`TUMORBAI3` etc.

- Replace the `--disk-size` with an appropriate disk size (in GB) of your virtual machine (this will typically be dictated by the size of the BAMs),
and replace `us-east1` with the appropriate [region](https://cloud.google.com/compute/docs/regions-zones) for the virtual machine (this will depend on your geographic location).

> :exclamation: The most likely problem you will encounter in this entire scenario is to make sure that
you have proper access to the above BAM/BAI files. You can save yourself a lot of time by executing a `gsutil stat gs://..`
command on all your input files at this point, which will query these cloud files for details, and unearth any access issues in
the process.

Running the script `cloud_run.sh` will trigger a cloud operation that starts the process outlined in the section [above](#cloud_steps).

**This command will wait (block) the command line till the process is complete.
Depending on the size of your BAMs the complexity of your analysis, this may take several hours.**

Note that pressing *Ctrl-C* on the terminal will stop the command from waiting for the completion of your operation, but your
operation on the cloud will continue to run. To query the status of the operation in the cloud, follow the instructions you see on the screen
when you ran the command. This will look something like:

```
dstat --provider google-v2 --project ...
```

The folder you specify in the `--logging` flag in the above command will contain live logging output on a minute-by-minute basis, which is helpful in diagnosing potential issues.
