from glob import glob
import pandas as pd

df=(pd.read_table(config['samples'], dtype={'sample': str, 'path': str})
	.set_index('sample', drop=False)
	.sort_index()
)

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
        'benchmarks/{sample}.fastcat_ont_reads.tsv'
    params:
        prefix='resources/reads/{sample}/ont'
    shell:
        '''
        fastcat -x \
            -f {output.fsummary} \
            -r {output.rsummary} \
            {params.prefix} | pigz > {output}
        '''


rule seqtk_ont_reads:
    input:
        rules.fastcat_ont_reads.output.reads
    output:
        'results/{sample}/common/seqtk/ont.fasta'
    threads:
        1
    conda:
        '../envs/seqtk.yml'
    benchmark:
        'benchmarks/{sample}.seqtk_ont_reads.tsv'
    shell:
        '''
        seqtk seq \
            -a {input} > {output}
        '''


rule fastcat_ont_reads_to_plot:
    input:
        rules.fastcat_ont_reads.output.rsummary
    output:
        'results/{sample}/common/fastcat/ont.per-read.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var}}' {input} | tail -n +2 > {output}
        '''
    
rule fastcat_ont_reads_combine_to_plot:
    input:
        expand('results/{sample}/common/fastcat/ont.per-read.mod.tsv', sample=df['sample'].tolist())
    output:
        'results/all.ont.per-read.mod.tsv'
    threads:
        1
    shell:
        '''
        cat {input} > {output}
        '''

rule fastcat_ont_reads_plot:
    input:
        rules.fastcat_ont_reads_combine_to_plot.output
    output:
        'results/all.ont.per-read.mod.pdf'
    threads:
        1
    conda:
        '../envs/r.yml'
    shell:
        '''
        Rscript workflow/scripts/plotstats.r {input} {output}
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

rule mosdepth_ont_to_plot:
    input:
        rules.mosdepth_ont.output
    output:
        'results/{sample}/common/mosdepth/ont.mosdepth.global.dist.mod.txt'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var , "ont"}}' {input} > {output}
        '''

rule mosdepth_illumina_to_plot:
    input:
        rules.mosdepth_illumina.output
    output:
        'results/{sample}/common/mosdepth/illumina.mosdepth.global.dist.mod.txt'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var , "illumina"}}' {input} > {output}
        '''

rule mosdepth_combine_to_plot:
    input:
        expand('results/{sample}/common/mosdepth/{platform}.mosdepth.global.dist.mod.txt', sample=df['sample'].tolist(),platform=['ont', 'illumina'])
    output:
        'results/all.mosdepth.global.dist.mod.txt'
    threads:
        1
    shell:
        '''
        cat {input} > {output}
        '''
rule mosdepth_plot:
    input:
        rules.mosdepth_combine_to_plot.output
    output:
        'results/all.mosdepth.global.dist.mod.pdf'
    threads:
        1
    conda:
        '../envs/r.yml'
    shell:
        '''
        Rscript workflow/scripts/plotcov.r {input} {output}
        '''