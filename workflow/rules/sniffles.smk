
from glob import glob
import pandas as pd

df=(pd.read_table(config['samples'], dtype={'sample': str, 'path': str})
	.set_index('sample', drop=False)
	.sort_index()
)

rule sniffles_mapping_variants:
    input:
        bam=rules.samtools_ont_sort_index.output,
        reference=config['reference'],
        bed=rules.bedops_makebed_repeatmasker_reference.output
    output:
        snf='results/{sample}/variants/mapping/sniffles/snf.snf',
        vcf='results/{sample}/variants/mapping/sniffles/snf.vcf.gz'       
    threads:
        config['sniffles']['threads']
    conda:
        '../envs/sniffles.yml'
    resources:
        mem_mb=config['sniffles']['mem'],
        time=config['sniffles']['time']
    benchmark:
        'benchmarks/{sample}.sniffles_mapping_variants.tsv'
    params:
        samplename='{sample}'
    shell:
        '''
        sniffles --input {input.bam} \
            --reference {input.reference} \
            -t {threads} \
            -v {output.vcf} \
            --snf {output.snf} \
            --minsvlen 50 \
            --sample-id {params.samplename} \
            --tandem-repeats {input.bed}
        '''

rule sniffles_mapping_combined:
    input:
        expand('results/{sample}/variants/mapping/sniffles/snf.snf', sample=df['sample'].tolist())
    output:
        'results/mapping.vcf.gz'
    threads:
        config['sniffles']['threads']
    conda:
        '../envs/sniffles.yml'
    resources:
        mem_mb=config['sniffles']['mem'],
        time=config['sniffles']['time']
    benchmark:
        'benchmarks/sniffles_mapping_combined.tsv'
    shell:
        '''
        sniffles \
            --input {input} \
            -v {output} \
            -t {threads} \
            --minsvlen 50
        '''


rule sniffles_assemblies:
    input:
        reference='results/wild_type/assembly/repeatmasker/ragtag.scaffold.fasta.masked',
        bam=rules.samtools_assemblies_bam.output,
        bed='results/wild_type/assembly/bedops/ragtag.scaffold.fasta.bed'
    output:
        'results/common/variants/assembly/sniffles/variants.vcf.gz'
    threads:
        config['sniffles']['threads']
    conda:
        '../envs/sniffles.yml'
    resources:
        mem_mb=config['sniffles']['mem'],
        time=config['sniffles']['time']
    benchmark:
        'benchmarks/sniffles_assemblies.tsv'
    shell:
        '''
        sniffles --input {input.bam} \
            --reference {input.reference} \
            -t {threads} \
            -v {output} \
            --minsvlen 50 \
            --sample-id mutant \
            --tandem-repeats {input.bed}
        '''