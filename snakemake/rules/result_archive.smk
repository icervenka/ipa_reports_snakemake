rule result_archive:
    input:
        rules.ipa.output
    output:
        ARCHIVE_OUTDIR + NOW + "_" + "{contrast}" + "_result_archive.tar.gz"
    params:
        ARCHIVE_OUTDIR
    run:
        shell(
            "tar -C reports -czf {output} {wildcards.contrast}"
        )