import pandas as pd
import glob
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.7.0")

##### load config and sample sheets #####
include: "snakemake/rules/common.smk"
include: "snakemake/rules/functions.smk"
configfile: "config.yaml"
# validate(config, schema="snakemake/schema/config.schema.yaml")

# Metadata = pd.read_table(config["metadata"])
# validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")
Inputs = list(glob_wildcards(INPUT_DIR+"{contrast}.txt"))[0]


rule all:
    input:
        expand(REPORT_OUTDIR + "{contrast}/" + "ipa.html", contrast=Inputs),
        expand(REPORT_OUTDIR + "{contrast}/" + "canonical_pathways.txt", contrast=Inputs),
        expand(REPORT_OUTDIR + "{contrast}/" + "upstream_regulators.txt", contrast=Inputs),
        expand(REPORT_OUTDIR + "{contrast}/" + "diseases_functions.txt", contrast=Inputs),
        expand(ARCHIVE_OUTDIR + NOW + "_" + "{contrast}" + "_result_archive.tar.gz")

rule ipa:
    input:
        INPUT_DIR+"{contrast}.txt"
    output:
        ipa=REPORT_OUTDIR + "{contrast}/" + "ipa.html",
        canonical=REPORT_OUTDIR + "{contrast}/" + "canonical_pathways.txt",
        regulators=REPORT_OUTDIR + "{contrast}/" + "upstream_regulators.txt",
        functions=REPORT_OUTDIR + "{contrast}/" + "diseases_functions.txt"
    params:
        pandoc_path=config["pandoc_path"],
        canonical_string=config["parse_strings"]["canonical"],
        upstream_string=config["parse_strings"]["upstream"],
        disease_string=config["parse_strings"]["disease"],
        tox_string=config["parse_strings"]["tox"],
        regulator_string=config["parse_strings"]["regulator"],
        networks_string=config["parse_strings"]["networks"],
        molecules_string=config["parse_strings"]["molecules"],
        graph_nitems=config["canonical"]["graph_nitems"],
        upstream_pval_threshold=config["upstream"]["graph_pval_threshold"],
        upstream_zscore_threshold=config["upstream"]["graph_zscore_threshold"],
        disease_pval_threshold=config["disease"]["graph_pval_threshold"],
        disease_zscore_threshold=config["disease"]["graph_zscore_threshold"],
        self_contained=config['self_contained']
    script:
        "snakemake/scripts/ipa.R"

rule result_archive:
    input:
        rules.ipa.output
    output:
        ARCHIVE_OUTDIR + NOW + "_" + "{contrast}" + "_result_archive.tar.gz"
    params:
        ARCHIVE_OUTDIR
    shell:
        input_string = " ".join(input)
        "tar -czf {output} {input_string} -C {params}"
