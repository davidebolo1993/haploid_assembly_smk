rule medaka_wtdbg2_raw:
    input:
        reads=rules.fastcat_ont_reads.output.reads,
        assembly=rules.wtdbg2_consensus.output
    output:
        'results/{sample}/assembly/medaka/wtdbg2/raw/consensus.fasta'
    threads:
        config['medaka']['threads']
    conda:
        '../envs/medaka.yml'
    resources:
        mem_mb=config['medaka']['mem'],
        time=config['medaka']['time']
    benchmark:
        'benchmarks/{sample}.medaka_wtdbg2_raw.tsv'
    params:
        prefix='results/{sample}/assembly/medaka/wtdbg2/raw',
        model=config['medaka']['model']
    shell:
        '''
        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.prefix} \
            -t {threads} \
            -m {params.model}
        '''

rule medaka_wtdbg2_cnsp:
    input:
        reads=rules.fastcat_ont_reads.output.reads,
        assembly=rules.wtdbg2_polish_consensus.output
    output:
        'results/{sample}/assembly/medaka/wtdbg2/cnsp/consensus.fasta'
    threads:
        config['medaka']['threads']
    conda:
        '../envs/medaka.yml'
    resources:
        mem_mb=config['medaka']['mem'],
        time=config['medaka']['time']
    benchmark:
        'benchmarks/{sample}.medaka_wtdbg2_cnsp.tsv'
    params:
        prefix='results/{sample}/assembly/medaka/wtdbg2/cnsp',
        model=config['medaka']['model']
    shell:
        '''
        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.prefix} \
            -t {threads} \
            -m {params.model}
        '''

rule medaka_canu:
    input:
        reads=rules.fastcat_ont_reads.output.reads,
        assembly=rules.canu_assemble.output
    output:
        'results/{sample}/assembly/medaka/canu/consensus.fasta'
    threads:
        config['medaka']['threads']
    conda:
        '../envs/medaka.yml'
    resources:
        mem_mb=config['medaka']['mem'],
        time=config['medaka']['time']
    benchmark:
        'benchmarks/{sample}.medaka_canu.tsv'
    params:
        prefix='results/{sample}/assembly/medaka/canu',
        model=config['medaka']['model']
    shell:
        '''
        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.prefix} \
            -t {threads} \
            -m {params.model}
        '''