rule canu_assemble:
    input:
        rules.fastcat_ont_reads.output.reads
    output:
        'results/{sample}/assembly/canu/assembly.contigs.fasta'
    threads:
        config['canu']['threads']
    conda:
        '../envs/canu.yml'
    resources:
        mem_mb=config['canu']['mem'],
        time=config['canu']['time']
    params:
        working_directory='results/{sample}/assembly/canu',
    benchmark:
        'benchmarks/{sample}.canu_assemble.tsv',
    shell:
        '''
        canu \
            -p assembly \
            -d {params.working_directory} \
            genomeSize=120m \
            maxInputCoverage=100 \
            -nanopore \
            useGrid=false \
            maxThreads={threads} \
            {input}
        '''