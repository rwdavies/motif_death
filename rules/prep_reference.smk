rule reference_all:
    input:
        f"{REF_DIR}/{REF_NAME}.fa.sa",
        f"{REF_DIR}/{REF_NAME}.fa.fai",
        f"{REF_DIR}/{REF_NAME}.dict",
        f"{REF_DIR}/{REF_NAME}.stidx"

rule unzip_ref:
    input:
        ref = f"{REF_DIR}/{REF_NAME}.fa.gz"
        # ref = ancient(f"{REF_DIR}/{REF_NAME}.fa.gz")
    output:
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    params: N='unzip_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
	'gunzip -c {REF_NAME}.fa.gz > {REF_NAME}.fa'

rule bwa_mem_ref:
    input:
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    output:
        f"{REF_DIR}/{REF_NAME}.fa.sa",
        f"{REF_DIR}/{REF_NAME}.fa.amb",
        f"{REF_DIR}/{REF_NAME}.fa.ann"
    params: N='bwa_mem_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        'bwa index {REF_NAME}.fa '

rule faidx_ref:
    input:
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    output:
        ref = f"{REF_DIR}/{REF_NAME}.fa.fai"
    params: N='faidx_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        'samtools faidx {REF_NAME}.fa'

rule picard_ref:
    input:
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    output:
        ref = f"{REF_DIR}/{REF_NAME}.dict"
    params: N='picard_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        '${{JAVA}} -Xmx12G -jar ${{PICARD}} CreateSequenceDictionary R={REF_NAME}.fa O={REF_NAME}.dict'

rule stampy_ref:
    input:
        ref = f"{REF_DIR}/{REF_NAME}.fa"
    output:
        f"{REF_DIR}/{REF_NAME}.sthash",
        f"{REF_DIR}/{REF_NAME}.stidx"
    params: N='stampy_ref', threads=1, queue = "short.qc"
    shell:
        'cd {REF_DIR} && '
        '${{PYTHON_278}} ${{STAMPY}} -G {REF_NAME} {REF_NAME}.fa && '
        '${{PYTHON_278}} ${{STAMPY}} -g {REF_NAME} -H {REF_NAME}'