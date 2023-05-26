from glob import glob

rule wtdbg2_assemble:
    input:
        lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/wtdbg2/assembly.ctg.lay.gz'
    threads:
        config['wtdbg2_assemble']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_assemble']['mem'],
        time=config['wtdbg2_assemble']['time']
    params:
        prefix='results/{sample}/wtdbg2/assembly'
    benchmark:
        'benchmarks/{sample}.wtdbg2_assemble.tsv'
    shell:
        '''
        zcat {input} | wtdbg2 \
            -p 0 -k 15 -AS 2 -s 0.05 -L 5000 \
            -g 120m \
            -t {threads} \
            -fo {params.prefix} \
            -i -
        '''

rule wtdbg2_consensus:
    input:
        rules.wtdbg2_assemble.output
    output:
        'results/{sample}/wtdbg2/assembly.raw.fa'
    threads:
        config['wtdbg2_consensus']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_consensus']['mem'],
        time=config['wtdbg2_consensus']['time']
    benchmark:
        'benchmarks/{sample}.wtdbg2_consensus.tsv'
    shell:
        '''
        wtpoa-cns \
            -t {threads} \
            -i {input} \
            -fo {output}
        '''

rule minimap2_align:
    input:
        reference=rules.wtdbg2_consensus.output,
        reads=lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/minimap2/alignment.sam'
    threads:
        config['minimap2']['threads']
    conda:
        '../envs/minimap2.yml'
    resources:
        mem_mb=config['minimap2']['mem'],
        time=config['minimap2']['time']
    benchmark:
        'benchmarks/{sample}.minimap2.tsv'
    shell:
        '''
        minimap2 -t {threads} -ax map-ont {input.reference} {input.reads} > {output}
        '''

rule samtools_filter_sort:
    input:
        rules.minimap2_align.output
    output:
        'results/{sample}/samtools/fltrd.srt.bam'
    threads:
        config['samtools']['threads']
    conda:
        '../envs/samtools.yml'
    resources:
        mem_mb=config['samtools']['mem'],
        time=config['samtools']['time']
    benchmark:
        'benchmarks/{sample}.samtools.tsv'
    shell:
        '''
        samtools view -b -F 0x900 -@ {threads} {input} | samtools sort -@ {threads} --write-index -o {output} -
        '''

rule wtdbg2_polish_consensus:
    input:
        reads=rules.samtools_filter_sort.output,
        assembly=rules.wtdbg2_consensus.output
    output:
        'results/{sample}/wtdbg2/assembly.cnsp.fa'
    threads:
        config['wtdbg2_consensus']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_consensus']['mem'],
        time=config['wtdbg2_consensus']['time']
    benchmark:
        'benchmarks/{sample}.wtdbg2_consensus_polish.tsv'
    shell:
        '''
        samtools view {input.reads} | wtpoa-cns \
            -t {threads} \
            -d {input.assembly} \
            -i - \
            -fo {output}
        '''
     