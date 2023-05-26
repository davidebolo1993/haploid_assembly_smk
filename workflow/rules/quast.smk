rule cat_reads:
    input:
        lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/quast/cns/tmp.fq.gz'
    threads:
        1
    shell:
        '''
        cat {input} > {output}
        '''

rule quast_on_cns:
    input:
        assembly=rules.wtdbg2_polish_consensus.output,
        reference=config['reference'],
        genes=config['genes'],
        reads=rules.cat_reads.output
    output:
        'results/{sample}/quast/cns/latest/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_cns.tsv'
    params:
        prefix='results/{sample}/quast/cns'
    shell:
        '''
        quast.py \
            {input.assembly} \
            -r {input.reference} \
            -g {input.genes} \
            -o {params.prefix} \
            --threads {threads} \
            --large \
            --k-mer-stats \
            --circos \
            --conserved-genes-finding \
            --glimmer \
            --nanopore {input.reads}
        '''    