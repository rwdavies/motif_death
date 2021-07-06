rule reference_all:
    input:
        REF_DIR + REFNAME + ".fa.sa",
        REF_DIR + REFNAME + ".fa.fai",
        REF_DIR + REFNAME + ".dict",
        REF_DIR + REFNAME + ".stidx"

rule download_ref:
    input:
    output:
        ref = REF_DIR + REFNAME + ".fa.gz"
    params: N='make_ref', threads=1, queue = "short.qc"
    shell:
        'mkdir -p {REF_DIR} && cd {REF_DIR} && '
        'wget {REF_URL}'

rule unzip_ref:
    input:
        ref = REF_DIR + REFNAME + ".fa.gz"
    output:
        ref = REF_DIR + REFNAME + ".fa"
    params: N='unzip_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
	'gunzip -c {REFNAME}.fa.gz > {REFNAME}.fa'

rule bwa_mem_ref:
    input:
        ref = REF_DIR + REFNAME + ".fa"
    output:
        REF_DIR + REFNAME + ".fa.sa"
    params: N='bwa_mem_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        'bwa index {REFNAME}.fa '

rule faidx_ref:
    input:
        ref = REF_DIR + REFNAME + ".fa"
    output:
        ref = REF_DIR + REFNAME + ".fa.fai"
    params: N='faidx_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        'samtools faidx {REFNAME}.fa'

rule picard_ref:
    input:
        ref = REF_DIR + REFNAME + ".fa"
    output:
        ref = REF_DIR + REFNAME + ".dict"
    params: N='picard_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        '${{JAVA}} -Xmx12G -jar ${{PICARD}} CreateSequenceDictionary R={REFNAME}.fa O={REFNAME}.dict'

rule stampy_ref:
    input:
        ref = REF_DIR + REFNAME + ".fa"
    output:
        REF_DIR + REFNAME + ".sthash",
        REF_DIR + REFNAME + ".stidx"
    params: N='stampy_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        '${{PYTHON_278}} ${{STAMPY}} -G {REFNAME} {REFNAME}.fa && '
        '${{PYTHON_278}} ${{STAMPY}} -g {REFNAME} -H {REFNAME}'