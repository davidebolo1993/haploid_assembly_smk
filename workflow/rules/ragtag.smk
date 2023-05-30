rule ragtag_correct_wtdbg2:
    input:
        query=rules.wtdbg2_polish_illumina.output,
        reference=config['reference'],
        reads=rules.get_ont_reads.output.fastq
    output:
        'results/{sample}/ragtag/correct/wtdbg2/ragtag.correct.fasta'
    threads:
        config['ragtag']['threads']
    conda:
        '../envs/ragtag.yml'
    resources:
        mem_mb=config['ragtag']['mem'],
        time=config['ragtag']['time']
    benchmark:
        'benchmarks/{sample}.ragtag_correct_wtdbg2.tsv'
    params:
        prefix='results/{sample}/ragtag/correct/wtdbg2'
    shell:
        '''
        ragtag.py correct \
            -w \
            -u \
            -R {input.reads} \
            -t {threads} \
            -T ont \
            -o {params.prefix} \
            {input.reference} \
            {input.query}            
        '''

rule ragtag_scaffold_wtdbg2:
    input:
        query=rules.ragtag_correct_wtdbg2.output,
        reference=config['reference']
    output:
        'results/{sample}/ragtag/scaffold/wtdbg2/ragtag.scaffold.fasta'
    threads:
        config['ragtag']['threads']
    conda:
        '../envs/ragtag.yml'
    resources:
        mem_mb=config['ragtag']['mem'],
        time=config['ragtag']['time']
    benchmark:
        'benchmarks/{sample}.ragtag_scaffold_wtdbg2.tsv'
    params:
        prefix='results/{sample}/ragtag/scaffold/wtdbg2'
    shell:
        '''
        ragtag.py scaffold \
            -w \
            -u \
            -t {threads} \
            -i 0.3 \
            -a 0.3 \
            -s 0.3 \
            -o {params.prefix} \
            {input.reference} \
            {input.query}            
        '''

rule ragtag_patch_wtdbg2:
    input:
        target=rules.ragtag_scaffold_wtdbg2.output,
        query=config['reference']
    output:
        'results/{sample}/ragtag/patch/wtdbg2/ragtag.patch.fasta'
    threads:
        config['ragtag']['threads']
    conda:
        '../envs/ragtag.yml'
    resources:
        mem_mb=config['ragtag']['mem'],
        time=config['ragtag']['time']
    benchmark:
        'benchmarks/{sample}.ragtag_patch_wtdbg2.tsv'
    params:
        prefix='results/{sample}/ragtag/patch/wtdbg2'
    shell:
        '''
        ragtag.py patch \
            {input.target} \
            {input.query} \
            -o {params.prefix} \
            -w \
            -u \
            -t {threads}
        '''    
