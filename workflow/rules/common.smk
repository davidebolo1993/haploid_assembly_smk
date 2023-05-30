rule fastcat_ont_reads:
    input:
        lambda wildcards: glob('resources/reads/{sample}/ont/*.gz'.format(sample=wildcards.sample))        
    output:
        reads='results/{sample}/common/fastcat/ont.fq.gz',
        rsummary='results/{sample}/common/fastcat/ont.per-read.tsv',
        fsummary='results/{sample}/common/fastcat/ont.per-file.tsv'
    threads:
        1
    conda:
        '../envs/fastcat.yml'
    benchmark:
        'benchmarks/{sample}.get_ont_reads.tsv'
    params:
        prefix='resources/reads/{sample}/ont'
    shell:
        '''
        fastcat -x \
            -f {output.fsummary} \
            -r {output.rsummary} \
            {params.prefix} > {output}
        '''

rule minimap2_ont_to_reference:
    input:
        reads=rules.fastcat_ont_reads.output.reads,
        reference=config['reference']
    output:
        'results/{sample}/common/minimap2/ont.sam'
    threads:
        config['minimap2']['threads']
    conda:
        '../envs/minimap2.yml'
    resources:
        mem_mb=config['minimap2']['mem'],
        time=config['minimap2']['time']
    benchmark:
        'benchmarks/{sample}.minimap2_ont_to_reference.tsv'
    shell:
        '''
        minimap2 \
            -t {threads} \
            -ax map-ont \
            --MD \
            {input.reference} {input.reads} > {output}
        '''

rule samtools_ont_sort_index:
    input:
       rules.minimap2_ont_to_reference.output
    output:
        'results/{sample}/common/samtools/ont.srt.bam'
    threads:
        config['samtools']['threads']
    conda:
        '../envs/samtools.yml'
    resources:
        mem_mb=config['samtools']['mem'],
        time=config['samtools']['time']
    benchmark:
        'benchmarks/{sample}.ont_sort_index.tsv'
    shell:
        '''
        samtools sort \
            -@ {threads} \
            --write-index \
            -o {output} \
            {input}
        '''

rule mosdepth_ont:
    input:
        rules.samtools_ont_sort_index.output
    output:
        'results/{sample}/common/mosdepth/ont.mosdepth.global.dist.txt'
    threads:
        config['mosdepth']['threads'],
    conda:
        '../envs/mosdepth.yml'
    resources:
        mem_mb=config['mosdepth']['mem'],
        time=config['mosdepth']['time']
    benchmark:
        'benchmarks/{sample}.mosdepth_ont.tsv'
    params:
        prefix='results/{sample}/common/mosdepth/ont'
    shell:
        '''
        mosdepth \
            -n \
            -x \
            -b 500 \
            -t {threads} \
            {params.prefix} \
            {input}
        '''

rule alfred_ont_qc:
    input:
        reference=config['reference'],
        reads=rules.samtools_ont_sort_index.output
    output:
        tsv='results/{sample}/common/alfred/ont.qc.tsv.gz',
        json='results/{sample}/common/alfred/ont.qc.json.gz'
    threads:
        1
    conda:
        '../envs/alfred.yml'
    resources:
        mem_mb=config['alfred']['mem'],
        time=config['alfred']['time']
    benchmark:
        'benchmarks/{sample}.alfred_ont_qc.tsv'
    shell:
        '''
        alfred qc \
            -r {input.reference} \
            -o {output.tsv} \
            -j {output.json} \
            {input.reads}
        '''

rule bwa_index_reference:
    input:
        config['reference']
    output:
        multiext(config['reference'], '.amb', '.ann', '.bwt', '.pac', '.sa')
    threads:
        1
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_index']['mem'],
        time=config['bwa_index']['time']
    benchmark:
        'benchmarks/bwa_index_reference.tsv'
    shell:
        '''
        bwa index {input}
        '''

rule bwa_mem_illumina_to_reference:
    input:
        reference=config['reference'],
        index=rules.bwa_index_reference.output,
        read_1=lambda wildcards: glob('resources/reads/{sample}/illumina/*R1*.gz'.format(sample=wildcards.sample)),
        read_2=lambda wildcards: glob('resources/reads/{sample}/illumina/*R2*.gz'.format(sample=wildcards.sample))
    output:
        'results/{sample}/common/bwa/illumina.sam'
    threads:
        config['bwa_mem']['threads']
    conda:
        '../envs/bwa.yml'
    resources:
        mem_mb=config['bwa_mem']['mem'],
        time=config['bwa_mem']['time']
    benchmark:
        'benchmarks/{sample}.bwa_mem_illumina_to_reference.tsv'
    shell:
        '''
        bwa mem \
            -t {threads} \
            {input.reference} {input.read_1} {input.read_2} > {output}
        '''

rule samtools_illumina_sort_index:
    input:
       rules.bwa_mem_illumina_to_reference.output
    output:
        'results/{sample}/common/samtools/illumina.srt.bam'
    threads:
        config['samtools']['threads']
    conda:
        '../envs/samtools.yml'
    resources:
        mem_mb=config['samtools']['mem'],
        time=config['samtools']['time']
    benchmark:
        'benchmarks/{sample}.illumina_sort_index.tsv'
    shell:
        '''
        samtools sort \
            -@ {threads} \
            --write-index \
            -o {output} \
            {input}
        '''

rule mosdepth_illumina:
    input:
        rules.samtools_illumina_sort_index.output
    output:
        'results/{sample}/common/mosdepth/illumina.mosdepth.global.dist.txt'
    threads:
        config['mosdepth']['threads'],
    conda:
        '../envs/mosdepth.yml'
    resources:
        mem_mb=config['mosdepth']['mem'],
        time=config['mosdepth']['time']
    benchmark:
        'benchmarks/{sample}.mosdepth_illumina.tsv'
    params:
        prefix='results/{sample}/common/mosdepth/illumina'
    shell:
        '''
        mosdepth \
            -n \
            -x \
            -b 500 \
            -t {threads} \
            {params.prefix} \
            {input}
        '''