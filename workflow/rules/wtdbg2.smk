from glob import glob

rule wtdbg2_assemble:
    input:
       rules.fastcat_ont_reads.output.reads
    output:
        'results/{sample}/assembly/wtdbg2/raw/assembly.ctg.lay.gz'
    threads:
        config['wtdbg2_assemble']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_assemble']['mem'],
        time=config['wtdbg2_assemble']['time']
    params:
        prefix='results/{sample}/assembly/wtdbg2/raw/assembly'
    benchmark:
        'benchmarks/{sample}.wtdbg2_assemble.tsv'
    shell:
        '''
        wtdbg2 \
            -p 0 -k 15 -AS 2 -s 0.05 -L 5000 \
            -g 120m \
            -t {threads} \
            -fo {params.prefix} \
            -i {input}
        '''

rule wtdbg2_consensus:
    input:
        rules.wtdbg2_assemble.output
    output:
        'results/{sample}/assembly/wtdbg2/raw/assembly.raw.fa'
    threads:
        config['wtdbg2_consensus']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_consensus']['mem'],
        time=config['wtdbg2_consensus']['time']
    benchmark:
        'benchmarks/{sample}.wtdbg2_consensus.tsv'
    shell:
        '''
        wtpoa-cns \
            -t {threads} \
            -i {input} \
            -fo {output}
        '''

rule minimap2_ont_to_assembly:
    input:
        reference=rules.wtdbg2_consensus.output,
        reads=rules.fastcat_ont_reads.output.reads
    output:
        'results/{sample}/assembly/minimap2/ont.sam'
    threads:
        config['minimap2']['threads']
    conda:
        '../envs/minimap2.yml'
    resources:
        mem_mb=config['minimap2']['mem'],
        time=config['minimap2']['time']
    benchmark:
        'benchmarks/{sample}.minimap2_ont_to_assembly.tsv'
    shell:
        '''
        minimap2 -t {threads} -ax map-ont {input.reference} {input.reads} > {output}
        '''

rule samtools_ont_assembly_filter_sort:
    input:
        rules.minimap2_ont_to_assembly.output
    output:
        'results/{sample}/assembly/samtools/ont.srt.bam'
    threads:
        config['samtools']['threads']
    conda:
        '../envs/samtools.yml'
    resources:
        mem_mb=config['samtools']['mem'],
        time=config['samtools']['time']
    benchmark:
        'benchmarks/{sample}.samtools_ont_assembly_filter_sort.tsv'
    shell:
        '''
        samtools view -b -F 0x900 -@ {threads} {input} | samtools sort -@ {threads} --write-index -o {output} -
        '''

rule wtdbg2_polish_consensus:
    input:
        reads=rules.samtools_ont_assembly_filter_sort.output,
        assembly=rules.wtdbg2_consensus.output
    output:
        'results/{sample}/assembly/wtdbg2/cnsp/assembly.cnsp.fa'
    threads:
        config['wtdbg2_consensus']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_consensus']['mem'],
        time=config['wtdbg2_consensus']['time']
    benchmark:
        'benchmarks/{sample}.wtdbg2_polish_consensus.tsv'
    shell:
        '''
        samtools view {input.reads} | wtpoa-cns \
            -t {threads} \
            -d {input.assembly} \
            -i - \
            -fo {output}
        '''

rule bwa_index_assembly:
    input:
        rules.wtdbg2_polish_consensus.output
    output:
        multiext('results/{sample}/assembly/wtdbg2/cnsp/assembly.cnsp.fa', '.amb', '.ann', '.bwt', '.pac', '.sa')
    threads:
        1
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_index']['mem'],
        time=config['bwa_index']['time']
    benchmark:
        'benchmarks/{sample}.bwa_index_assembly.tsv'
    shell:
        '''
        bwa index {input}
        '''

rule bwa_illumina_to_assembly:
    input:
        reference=rules.wtdbg2_polish_consensus.output,
        index=rules.bwa_index_assembly.output,
        read_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
        read_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/assembly/bwa/alignment.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_illumina_to_assembly.tsv'
    shell:
        '''
        bwa mem -t {threads} {input.reference} {input.read_1} {input.read_2} > {output}
        '''

rule samtools_illumina_assembly_sort_sam:
    input:
        rules.bwa_illumina_to_assembly.output
    output:
        'results/{sample}/assembly/samtools/alignment.srt.sam'
    threads:
        config['samtools']['threads']
    conda:
        '../envs/samtools.yml'
    resources:
        mem_mb=config['samtools']['mem'],
        time=config['samtools']['time']
    benchmark:
        'benchmarks/{sample}.samtools_illumina_assembly_sort_sam.tsv'
    shell:
        '''
        samtools sort -@ {threads} -o {output} -O SAM {input}
        '''

rule wtdbg2_polish_illumina:
    input:
        reads=rules.samtools_illumina_assembly_sort_sam.output,
        assembly=rules.wtdbg2_polish_consensus.output
    output:
        'results/{sample}/assembly/wtdbg2/srp/assembly.srp.fa'
    threads:
        config['wtdbg2_consensus']['threads']
    conda:
        '../envs/wtdbg2.yml'
    resources:
        mem_mb=config['wtdbg2_consensus']['mem'],
        time=config['wtdbg2_consensus']['time']
    benchmark:
        'benchmarks/{sample}.wtdbg2_polish_illumina.tsv'
    shell:
        '''
        samtools view {input.reads} | wtpoa-cns \
            -t {threads} \
            -x sam-sr \
            -d {input.assembly} \
            -i - \
            -fo {output}
        '''