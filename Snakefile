# type: ignore
### ISSUES
# for some reason the Snakemake file is run several times so variable NOW changes
# each time when formatted up to seconds
# ._ files from MacOS in input directory will interfere with the run

import pandas as pd
import glob
import datetime
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.7.0")

##### load config and sample sheets #####
include: "snakemake/rules/common.smk"
include: "snakemake/rules/functions.smk"
configfile: "config.yaml"
# validate(config, schema="snakemake/schema/config.schema.yaml")

Inputs = list(glob_wildcards(INPUT_DIR+"{contrast}.txt")[0])
NOW = str(datetime.datetime.now().strftime("%Y%m%d"))

rule all:
    input:
        expand(REPORT_OUTDIR + "{contrast}/" + "report.html", contrast=Inputs),
        expand(REPORT_OUTDIR + "{contrast}/" + "canonical_pathways.txt", contrast=Inputs),
        expand(REPORT_OUTDIR + "{contrast}/" + "upstream_regulators.txt", contrast=Inputs),
        expand(REPORT_OUTDIR + "{contrast}/" + "diseases_functions.txt", contrast=Inputs),
        expand(ARCHIVE_OUTDIR + NOW + "_" + "{contrast}" + "_result_archive.tar.gz", contrast=Inputs)

include: "snakemake/rules/ipa.smk"
include: "snakemake/rules/result_archive.smk"
        
