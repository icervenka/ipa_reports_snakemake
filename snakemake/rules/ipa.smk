rule ipa:
    input:
        INPUT_DIR+"{contrast}.txt"
    output:
        ipa=REPORT_OUTDIR + "{contrast}/" + "report.html",
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
        "../../snakemake/scripts/ipa.R"