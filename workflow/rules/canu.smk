rule canu_assemble_ont:
    input:
        lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/canu/canu.contigs.fasta'
    threads:
        config['canu']['threads']
    conda:
        '../envs/canu.yml'
    resources:
        mem_mb=config['canu']['mem'],
        time=config['canu']['time']
    params:
        working_directory='results/{sample}/canu',
    benchmark:
        'benchmarks/{sample}.canu_assemble.tsv',
    shell:
        '''
        canu \
            -p canu \
            -d {params.working_directory} \
            genomeSize=120m \
            maxInputCoverage=100 \
            -nanopore \
            useGrid=false \
            maxThreads={threads} \
            {input}
        '''