rule ipa:
    input:
        get_ipa_file
    output:
        ipa = REPORT_OUTDIR + "{contrast}/" + "ipa.html",
        canonical = REPORT_OUTDIR + "{contrast}/" + "canonical_pathways.txt",
        regulators = REPORT_OUTDIR + "{contrast}/" + "upstream_regulators.txt",
        functions = REPORT_OUTDIR + "{contrast}/" + "diseases_functions.txt"
    params:
        canonical_string=config["ipa"]["canonical_string"],
        upstream_string=config["ipa"]["upstream_string"],
        disease_string=config["ipa"]["disease_string"],
        tox_string=config["ipa"]["tox_string"],
        regulator_string=config["ipa"]["regulator_string"],
        networks_string=config["ipa"]["networks_string"],
        molecules_string=config["ipa"]["molecules_string"],
        graph_nitems=config["ipa"]["canonical"]["graph_nitems"],
        upstream_pval_threshold=config["ipa"]["upstream"]["graph_pval_threshold"],
        upstream_zscore_threshold=config["ipa"]["upstream"]["graph_zscore_threshold"],
        disease_pval_threshold=config["ipa"]["disease"]["graph_pval_threshold"],
        disease_zscore_threshold=config["ipa"]["disease"]["graph_zscore_threshold"]
    shell:
        "echo {input}"

      # "snakemake/scripts/ipa.R"

# ipa:
#   canonical_string: "Canonical Pathways for"
#   upstream_string: "Upstream Regulators"
#   disease_string: "Diseases and Bio Functions"
#   tox_string: "Tox Functions"
#   regulator_string: "Regulator Effects"
#   networks_string: "Networks"
#   molecules_string: "Analysis Ready Molecules"
#   canonical:
#     graph_nitems: 20
#     select_rows: NULL
#     exclude_rows: NULL
#   upstream:
#     graph_pval_threshold: 0.05
#     graph_zscore_threshold: 1.5
#   disease:
#     graph_pval_threshold: 0.05
#     graph_zscore_threshold: 1.5
