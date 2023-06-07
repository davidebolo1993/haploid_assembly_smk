import os

rule repeatmasker_reference:
    input:
        config['reference']
    output:
        fa='results/common/repeatmasker/' + os.path.basename(config['reference']) + '.masked',
        out='results/common/repeatmasker/' + os.path.basename(config['reference']) + '.out'
    threads:
        config['repeatmasker']['threads']
    conda:
        '../envs/repeatmasker.yml'
    resources:
        mem_mb=config['repeatmasker']['mem'],
        time=config['repeatmasker']['time']
    benchmark:
        'benchmarks/repeatmasker_reference.tsv'
    params:
        prefix='results/common/repeatmasker',
        parallel='4'
    shell:
        '''
        RepeatMasker \
            -e ncbi \
            -pa {params.parallel} \
            -q \
            -nocut \
            -gff \
            -dir {params.prefix} \
            -species "Chlamydomonas reinhardtii" \
            {input}
        '''

rule repeatmasker_assemblies:
    input:
        rules.ragtag_scaffold_best_assembly.output
    output:
        fa='results/{sample}/assembly/repeatmasker/ragtag.scaffold.fasta.masked',
        out='results/{sample}/assembly/repeatmasker/ragtag.scaffold.fasta.out'
    threads:
        config['repeatmasker']['threads']
    conda:
        '../envs/repeatmasker.yml'
    resources:
        mem_mb=config['repeatmasker']['mem'],
        time=config['repeatmasker']['time']
    benchmark:
        'benchmarks/{sample}.repeatmasker_assemblies.tsv'
    params:
        prefix='results/{sample}/assembly/repeatmasker',
        parallel='4'
    shell:
        '''
        RepeatMasker \
            -e ncbi \
            -pa {params.parallel} \
            -q \
            -nocut \
            -gff \
            -dir {params.prefix} \
            -species "Chlamydomonas reinhardtii" \
            {input}
        '''
