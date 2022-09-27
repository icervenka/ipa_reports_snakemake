# rule ipa:
#     input:
#         get_ipa_file
#     output:
#         ipa = REPORT_OUTDIR + "{contrast}/" + "ipa.html",
#         canonical = REPORT_OUTDIR + "{contrast}/" + "canonical_pathways.txt",
#         regulators = REPORT_OUTDIR + "{contrast}/" + "upstream_regulators.txt",
#         functions = REPORT_OUTDIR + "{contrast}/" + "diseases_functions.txt"
#     params:
#         canonical_string=config["ipa"]["canonical_string"],
#         upstream_string=config["ipa"]["upstream_string"],
#         disease_string=config["ipa"]["disease_string"],
#         tox_string=config["ipa"]["tox_string"],
#         regulator_string=config["ipa"]["regulator_string"],
#         networks_string=config["ipa"]["networks_string"],
#         molecules_string=config["ipa"]["molecules_string"],
#         graph_nitems=config["ipa"]["canonical"]["graph_nitems"],
#         upstream_pval_threshold=config["ipa"]["upstream"]["graph_pval_threshold"],
#         upstream_zscore_threshold=config["ipa"]["upstream"]["graph_zscore_threshold"],
#         disease_pval_threshold=config["ipa"]["disease"]["graph_pval_threshold"],
#         disease_zscore_threshold=config["ipa"]["disease"]["graph_zscore_threshold"],
#         self_contained=config['self_contained']
#     shell:
#         "snakemake/scripts/ipa.R"
