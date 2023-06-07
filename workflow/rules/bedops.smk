rule bedops_makebed_repeatmasker_reference:
    input:
        rules.repeatmasker_reference.output.out
    output:
        'results/common/bedops/' + os.path.basename(config['reference']) + '.bed'
    threads:
        1
    conda:
        '../envs/bedops.yml'
    resources:
        mem_mb=config['bedops']['mem'],
        time=config['bedops']['time']
    benchmark:
        'benchmarks/bedops_makebed_repeatmasker_reference.tsv'
    shell:
        '''
        rmsk2bed < {input} > {output}
        '''

rule bedops_makebed_repeatmasker_assemblies:
    input:
        rules.repeatmasker_assemblies.output.out
    output:
        'results/{sample}/assembly/bedops/ragtag.scaffold.fasta.bed'
    threads:
        1
    conda:
        '../envs/bedops.yml'
    resources:
        mem_mb=config['bedops']['mem'],
        time=config['bedops']['time']
    benchmark:
        'benchmarks/{sample}.bedops_makebed_repeatmasker_assemblies.tsv'
    shell:
        '''
        rmsk2bed < {input} > {output}
        '''