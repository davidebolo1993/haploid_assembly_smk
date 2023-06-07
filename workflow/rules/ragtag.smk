rule ragtag_correct_best_assembly:
    input:
        target=config['reference'],
        query=rules.medaka_canu.output,
        reads=rules.fastcat_ont_reads.output.reads
    output:
        'results/{sample}/assembly/ragtag/correct/ragtag.correct.fasta'
    threads:
        config['ragtag']['threads']
    conda:
        '../envs/ragtag.yml'
    resources:
        mem_mb=config['ragtag']['mem'],
        time=config['ragtag']['time']
    benchmark:
        'benchmarks/{sample}.ragtag_correct_best_assembly.tsv'
    params:
        prefix='results/{sample}/assembly/ragtag/correct'
    shell:
        '''
        ragtag.py correct \
            {input.target} \
            {input.query} \
            --remove-small \
            -o {params.prefix} \
            -w \
            -u \
            -t {threads} \
            -R {input.reads} \
            -T ont
        '''

rule ragtag_scaffold_best_assembly:
    input:
        target=config['reference'],
        query=rules.ragtag_correct_best_assembly.output
    output:
        'results/{sample}/assembly/ragtag/scaffold/ragtag.scaffold.fasta'
    threads:
        config['ragtag']['threads']
    conda:
        '../envs/ragtag.yml'
    resources:
        mem_mb=config['ragtag']['mem'],
        time=config['ragtag']['time']
    benchmark:
        'benchmarks/{sample}.ragtag_scaffold_best_assembly.tsv'
    params:
        prefix='results/{sample}/assembly/ragtag/scaffold'
    shell:
        '''
        ragtag.py scaffold \
            {input.target} \
            {input.query} \
            --remove-small \
            -o {params.prefix} \
            -w \
            -u \
            -t {threads} \
            -r
        '''
    
