rule quast:
    input:
        assembly=''
        reference=config['reference'],
        genes=config['genes'],
        reads=''
    output:
        ''
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        ''
    params:
        prefix=''
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
            --conserved-genes-finding \
            --glimmer
        '''

