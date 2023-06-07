rule seqtk_filter_for_paf:
    input:
        reference=config['reference'],
        assembly=rules.repeatmasker_assemblies.output.fa
    output:
        'results/{sample}/assembly/seqtk/ragtag.scaffold.fasta.masked.clean.fa'
    threads:
        1
    conda:
        '../envs/seqtk.yml'
    benchmark:
        'benchmarks/{sample}.filter_for_paf.tsv'
    params:
        tmp_file='results/{sample}/assembly/seqtk/names.txt',
        samplename='{sample}'
    shell:
        '''
        grep "^>" {input.reference} | \
            cut -d " " -f 1 | \
            sed 's/>//g' | \
            sed 's/$/_RagTag/g' > {params.tmp_file} \
        && seqtk subseq \
            {input.assembly} {params.tmp_file} | \
            awk -v var="{params.samplename}" '/^>/{{$1=$1"_"var}} 1' > {output}
        '''

rule minimap2_assemblies_paf:
    input:
        expand('results/{sample}/assembly/seqtk/ragtag.scaffold.fasta.masked.clean.fa', sample=df['sample'].tolist())
    output:
        'results/common/paf/assemblies.paf'
    threads:
        config['minimap2']['threads']
    conda:
        '../envs/minimap2.yml'
    resources:
        mem_mb=config['minimap2']['mem'],
        time=config['minimap2']['time']
    benchmark:
        'benchmarks/minimap2_assemblies_paf.tsv'
    shell:
        '''
        minimap2 \
            -x asm5 \
            -t {threads} \
            {input} > {output}
        '''

rule plot_assemblies_paf:
    input:
        rules.minimap2_assemblies_paf.output
    output:
        'results/common/paf/assemblies.dot.pdf'
    threads:
        1
    conda:
        '../envs/r.yml'
    params:
        prefix='results/common/paf'
    shell:
        '''
        Rscript workflow/scripts/plotpaf.r {input} {params.prefix}
        '''