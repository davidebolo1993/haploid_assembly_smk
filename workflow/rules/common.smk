rule get_ont_reads:
    input:
        lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))        
    output:
        fastq='results/{sample}/ont.fq.gz',
        rsummary='results/{sample}/ont.per-read.tsv',
        fsummary='results/{sample}/ont.per-file.tsv'
    threads:
        1
    conda:
        '../envs/fastcat.yml'
    benchmark:
        'benchmarks/{sample}.fastcat.tsv'
    params:
        prefix='resources/reads/{sample}/ont'
    shell:
        '''
        fastcat -x \
            -f {output.fsummary} \
            -r {output.rsummary} \
            {params.prefix} > {output}
        '''
