rule minimap2_assemblies_sam:
    input:
        ref='results/wild_type/assembly/repeatmasker/ragtag.scaffold.fasta.masked',
        alt='results/mutant/assembly/repeatmasker/ragtag.scaffold.fasta.masked'
    output:
        'results/common/variants/assembly/assemblies.sam'
    threads:
        config['minimap2']['threads']
    conda:
        '../envs/minimap2.yml'
    resources:
        mem_mb=config['minimap2']['mem'],
        time=config['minimap2']['time']
    benchmark:
        'benchmarks/minimap2_assemblies_sam.tsv'
    shell:
        '''
        minimap2 \
            -a \
            -x asm5 \
            -t {threads} \
            --cs \
            --MD \
            -r2k \
            {input.ref} {input.alt} > {output}
        '''
   
rule samtools_assemblies_bam:
    input:
        rules.minimap2_assemblies_sam.output
    output:
         'results/common/variants/assembly/assemblies.srt.bam'
    threads:
        config['samtools']['threads']
    conda:
        '../envs/samtools.yml'
    resources:
        mem_mb=config['samtools']['mem'],
        time=config['samtools']['time']
    benchmark:
        'benchmarks/samtools_assemblies_bam.tsv'
    shell:
        '''
        samtools sort \
            -@ {threads} \
            -o {output} \
            --write-index \
            {input}
        '''

rule svim_asm_assemblies:
    input:
        ref='results/wild_type/assembly/repeatmasker/ragtag.scaffold.fasta.masked',
        bam=rules.samtools_assemblies_bam.output
    output:
        'results/common/variants/assembly/svim-asm/variants.vcf'
    threads:
        1
    conda:
        '../envs/svim-asm.yml'
    resources:
        mem_mb=config['svim-asm']['mem'],
        time=config['svim-asm']['time']
    benchmark:
        'benchmarks/svim_asm_assemblies.tsv'
    params:
        prefix='results/common/variants/assembly/svim-asm'
    shell:
        '''
        svim-asm haploid \
            {params.prefix} {input.bam} {input.ref}
        '''
