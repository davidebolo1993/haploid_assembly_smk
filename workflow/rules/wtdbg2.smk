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
        zcat {input} | wtdbg2 -x ont \
            -g 120m \
            -t {threads} \
            -fo {params.prefix} \
            -i -
        '''
