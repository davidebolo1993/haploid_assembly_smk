rule bwa_index_wtdbg2_medaka_cnsp:
    input:
        rules.medaka_wtdbg2_cnsp.output
    output:
        multiext('results/{sample}/assembly/medaka/wtdbg2/cnsp/consensus.fasta', '.amb', '.ann', '.bwt', '.pac', '.sa')
    threads:
        1
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_index']['mem'],
        time=config['bwa_index']['time']
    benchmark:
        'benchmarks/{sample}.bwa_index_wtdbg2_medaka_cnsp.tsv'
    shell:
        '''
        bwa index {input}
        '''

rule bwa_illumina_read1_to_wtdbg2_medaka_cnsp_assembly:
    input:
        reference=rules.medaka_wtdbg2_cnsp.output,
        index=rules.bwa_index_wtdbg2_medaka_cnsp.output,
        read_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
    output:
        'results/{sample}/assembly/polypolish/wtdbg2_medaka/cnsp/bwa_read1/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_read1_to_wtdbg2_medaka_cnsp_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} -a {input.reference} {input.read_1} > {output}
        '''

rule bwa_illumina_read2_to_wtdbg2_medaka_cnsp_assembly:
    input:
        reference=rules.medaka_wtdbg2_cnsp.output,
        index=rules.bwa_index_wtdbg2_medaka_cnsp.output,
        read_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample)),
    output:
        'results/{sample}/assembly/polypolish/wtdbg2_medaka/cnsp/bwa_read2/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_read2_to_wtdbg2_medaka_cnsp_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} -a {input.reference} {input.read_2} > {output}
        '''

rule polypolish_insert_filter_wtdbg2_medaka_cnsp:
    input:
        alignment1=rules.bwa_illumina_read1_to_wtdbg2_medaka_cnsp_assembly.output,
        alignment2=rules.bwa_illumina_read2_to_wtdbg2_medaka_cnsp_assembly.output
    output:
        alignment1_filtered='results/{sample}/assembly/polypolish/wtdbg2_medaka/cnsp/bwa_read1/alignment.fltrd.sam',
        alignment2_filtered='results/{sample}/assembly/polypolish/wtdbg2_medaka/cnsp/bwa_read2/alignment.fltrd.sam'
    threads:
        1
    conda:
        '../envs/polypolish.yml'
    resources:
        mem_mb=config['polypolish']['mem'],
        time=config['polypolish']['time']
    benchmark:
        'benchmarks/{sample}.polypolish_insert_filter_wtdbg2_medaka_cnsp.tsv'
    shell:
        '''
        polypolish_insert_filter.py \
            --in1 {input.alignment1} \
            --in2 {input.alignment2} \
            --out1 {output.alignment1_filtered}\
            --out2 {output.alignment2_filtered}
        '''
    

rule polypolish_polish_wtdbg2_medaka_cnsp:
    input:
        alignment1=rules.polypolish_insert_filter_wtdbg2_medaka_cnsp.output.alignment1_filtered,
        alignment2=rules.polypolish_insert_filter_wtdbg2_medaka_cnsp.output.alignment2_filtered,
        assembly=rules.medaka_wtdbg2_cnsp.output
    output:
        'results/{sample}/assembly/polypolish/wtdbg2_medaka/cnsp/polished.fasta'
    threads:
        1
    conda:
        '../envs/polypolish.yml'
    resources:
        mem_mb=config['polypolish']['mem'],
        time=config['polypolish']['time']
    benchmark:
        'benchmarks/{sample}.polypolish_polish_wtdbg2_medaka_cnsp.tsv'
    shell:
        '''
        polypolish \
            {input.assembly} \
            {input.alignment1} \
            {input.alignment2} > {output}
        '''

rule bwa_index_wtdbg2_medaka_raw:
    input:
        rules.medaka_wtdbg2_raw.output
    output:
        multiext('results/{sample}/assembly/medaka/wtdbg2/raw/consensus.fasta', '.amb', '.ann', '.bwt', '.pac', '.sa')
    threads:
        1
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_index']['mem'],
        time=config['bwa_index']['time']
    benchmark:
        'benchmarks/{sample}.bwa_index_wtdbg2_medaka_raw.tsv'
    shell:
        '''
        bwa index {input}
        '''

rule bwa_illumina_read1_to_wtdbg2_medaka_raw_assembly:
    input:
        reference=rules.medaka_wtdbg2_raw.output,
        index=rules.bwa_index_wtdbg2_medaka_raw.output,
        read_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
    output:
        'results/{sample}/assembly/polypolish/wtdbg2_medaka/raw/bwa_read1/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_read1_to_wtdbg2_medaka_raw_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} -a {input.reference} {input.read_1} > {output}
        '''

rule bwa_illumina_read2_to_wtdbg2_medaka_raw_assembly:
    input:
        reference=rules.medaka_wtdbg2_raw.output,
        index=rules.bwa_index_wtdbg2_medaka_raw.output,
        read_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample)),
    output:
        'results/{sample}/assembly/polypolish/wtdbg2_medaka/raw/bwa_read2/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_read2_to_wtdbg2_medaka_raw_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} -a {input.reference} {input.read_2} > {output}
        '''

rule polypolish_insert_filter_wtdbg2_medaka_raw:
    input:
        alignment1=rules.bwa_illumina_read1_to_wtdbg2_medaka_raw_assembly.output,
        alignment2=rules.bwa_illumina_read2_to_wtdbg2_medaka_raw_assembly.output
    output:
        alignment1_filtered='results/{sample}/assembly/polypolish/wtdbg2_medaka/raw/bwa_read1/alignment.fltrd.sam',
        alignment2_filtered='results/{sample}/assembly/polypolish/wtdbg2_medaka/raw/bwa_read2/alignment.fltrd.sam'
    threads:
        1
    conda:
        '../envs/polypolish.yml'
    resources:
        mem_mb=config['polypolish']['mem'],
        time=config['polypolish']['time']
    benchmark:
        'benchmarks/{sample}.polypolish_insert_filter_wtdbg2_medaka_raw.tsv'
    shell:
        '''
        polypolish_insert_filter.py \
            --in1 {input.alignment1} \
            --in2 {input.alignment2} \
            --out1 {output.alignment1_filtered}\
            --out2 {output.alignment2_filtered}
        '''
    

rule polypolish_polish_wtdbg2_medaka_raw:
    input:
        alignment1=rules.polypolish_insert_filter_wtdbg2_medaka_raw.output.alignment1_filtered,
        alignment2=rules.polypolish_insert_filter_wtdbg2_medaka_raw.output.alignment2_filtered,
        assembly=rules.medaka_wtdbg2_raw.output
    output:
        'results/{sample}/assembly/polypolish/wtdbg2_medaka/raw/polished.fasta'
    threads:
        1
    conda:
        '../envs/polypolish.yml'
    resources:
        mem_mb=config['polypolish']['mem'],
        time=config['polypolish']['time']
    benchmark:
        'benchmarks/{sample}.polypolish_polish_wtdbg2_medaka_raw.tsv'
    shell:
        '''
        polypolish \
            {input.assembly} \
            {input.alignment1} \
            {input.alignment2} > {output}
        '''


rule bwa_index_canu_medaka:
    input:
        rules.medaka_canu.output
    output:
        multiext('results/{sample}/assembly/medaka/canu/consensus.fasta', '.amb', '.ann', '.bwt', '.pac', '.sa')
    threads:
        1
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_index']['mem'],
        time=config['bwa_index']['time']
    benchmark:
        'benchmarks/{sample}.bwa_index_canu_medaka.tsv'
    shell:
        '''
        bwa index {input}
        '''

rule bwa_illumina_read1_to_canu_medaka_assembly:
    input:
        reference=rules.medaka_canu.output,
        index=rules.bwa_index_canu_medaka.output,
        read_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
    output:
        'results/{sample}/assembly/polypolish/canu_medaka/bwa_read1/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_read1_to_canu_medaka_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} -a {input.reference} {input.read_1} > {output}
        '''

rule bwa_illumina_read2_to_canu_medaka_assembly:
    input:
        reference=rules.medaka_canu.output,
        index=rules.bwa_index_canu_medaka.output,
        read_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample)),
    output:
        'results/{sample}/assembly/polypolish/canu_medaka/bwa_read2/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_read2_to_canu_medaka_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} -a {input.reference} {input.read_2} > {output}
        '''

rule polypolish_insert_filter_canu_medaka:
    input:
        alignment1=rules.bwa_illumina_read1_to_canu_medaka_assembly.output,
        alignment2=rules.bwa_illumina_read2_to_canu_medaka_assembly.output
    output:
        alignment1_filtered='results/{sample}/assembly/polypolish/canu_medaka/bwa_read1/alignment.fltrd.sam',
        alignment2_filtered='results/{sample}/assembly/polypolish/canu_medaka/bwa_read2/alignment.fltrd.sam'
    threads:
        1
    conda:
        '../envs/polypolish.yml'
    resources:
        mem_mb=config['polypolish']['mem'],
        time=config['polypolish']['time']
    benchmark:
        'benchmarks/{sample}.polypolish_insert_filter_canu_medaka.tsv'
    shell:
        '''
        polypolish_insert_filter.py \
            --in1 {input.alignment1} \
            --in2 {input.alignment2} \
            --out1 {output.alignment1_filtered}\
            --out2 {output.alignment2_filtered}
        '''
    
rule polypolish_polish_canu_medaka:
    input:
        alignment1=rules.polypolish_insert_filter_canu_medaka.output.alignment1_filtered,
        alignment2=rules.polypolish_insert_filter_canu_medaka.output.alignment2_filtered,
        assembly=rules.medaka_canu.output
    output:
        'results/{sample}/assembly/polypolish/canu_medaka/polished.fasta'
    threads:
        1
    conda:
        '../envs/polypolish.yml'
    resources:
        mem_mb=config['polypolish']['mem'],
        time=config['polypolish']['time']
    benchmark:
        'benchmarks/{sample}.polypolish_polish_canu_medaka.tsv'
    shell:
        '''
        polypolish \
            {input.assembly} \
            {input.alignment1} \
            {input.alignment2} > {output}
        '''
