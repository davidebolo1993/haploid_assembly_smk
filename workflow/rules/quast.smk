rule get_ont_reads:
    input:
        lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))        
    output:
        'results/{sample}/quast/ont.fq.gz'
    threads:
        1
    benchmark:
        'benchmarks/{sample}.cat_fastq.tsv'
    shell:
        '''
        cat {input} > {output}
        '''

rule quast_on_cns_nanopore:
    input:
        assembly=rules.wtdbg2_polish_consensus.output,
        reference=config['reference'],
        genes=config['genes'],
        reads=rules.get_ont_reads.output
    output:
        'results/{sample}/quast/cns/ont/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_cns_ont.tsv'
    params:
        prefix='results/{sample}/quast/cns/ont'
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

rule quast_on_cns_illumina:
    input:
        assembly=rules.wtdbg2_polish_consensus.output,
        reference=config['reference'],
        genes=config['genes'],
        reads_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
        reads_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/quast/cns/illumina/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_cns_illumina.tsv'
    params:
        prefix='results/{sample}/quast/cns/illumina'
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
            --pe1 {input.reads_1} \
            --pe2 {input.reads_2}
        '''

rule quast_on_srp_nanopore:
    input:
        assembly=rules.wtdbg2_polish_illumina.output,
        reference=config['reference'],
        genes=config['genes'],
        reads=rules.get_ont_reads.output
    output:
        'results/{sample}/quast/srp/ont/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_srp_ont.tsv'
    params:
        prefix='results/{sample}/quast/srp/ont'
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

rule quast_on_srp_illumina:
    input:
        assembly=rules.wtdbg2_polish_illumina.output,
        reference=config['reference'],
        genes=config['genes'],
        reads_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
        reads_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/quast/srp/illumina/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_srp_illumina.tsv'
    params:
        prefix='results/{sample}/quast/srp/illumina'
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
            --pe1 {input.reads_1} \
            --pe2 {input.reads_2}
        '''