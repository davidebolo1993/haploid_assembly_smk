import pandas as pd

df=(pd.read_table(config['samples'], dtype={'sample': str, 'path': str})
	.set_index('sample', drop=False)
	.sort_index()
)


rule quast_wtdbg2_cnsp:
    input:
        assembly=rules.wtdbg2_polish_consensus.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/wtdbg2/cnsp/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_wtdbg2_cns.tsv'
    params:
        prefix='results/{sample}/assembly/quast/wtdbg2/cnsp'
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
            --conserved-genes-finding \
            --glimmer
        '''

rule quast_wtdbg2_srp:
    input:
        assembly=rules.wtdbg2_polish_illumina.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/wtdbg2/srp/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_wtdbg2_srp.tsv'
    params:
        prefix='results/{sample}/assembly/quast/wtdbg2/srp'
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
            --conserved-genes-finding \
            --glimmer
        '''

rule quast_wtdbg2_medaka_raw:
    input:
        assembly=rules.medaka_wtdbg2_raw.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka/raw/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_wtdbg2_medaka_raw.tsv'
    params:
        prefix='results/{sample}/assembly/quast/wtdbg2_medaka/raw'
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
            --conserved-genes-finding \
            --glimmer
        '''

rule quast_wtdbg2_medaka_cnsp:
    input:
        assembly=rules.medaka_wtdbg2_cnsp.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka/cnsp/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_wtdbg2_medaka_cnsp.tsv'
    params:
        prefix='results/{sample}/assembly/quast/wtdbg2_medaka/cnsp'
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
            --conserved-genes-finding \
            --glimmer
        '''


rule quast_wtdbg2_medaka_polypolish_raw:
    input:
        assembly=rules.polypolish_polish_wtdbg2_medaka_raw.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/raw/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_wtdbg2_medaka_polypolish_raw.tsv'
    params:
        prefix='results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/raw'
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
            --conserved-genes-finding \
            --glimmer
        '''


rule quast_wtdbg2_medaka_polypolish_cnsp:
    input:
        assembly=rules.polypolish_polish_wtdbg2_medaka_cnsp.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/cnsp/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_wtdbg2_medaka_polypolish_cnsp.tsv'
    params:
        prefix='results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/cnsp'
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
            --conserved-genes-finding \
            --glimmer
        '''

rule quast_canu_medaka:
    input:
        assembly=rules.medaka_canu.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/canu_medaka/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_canu_medaka.tsv'
    params:
        prefix='results/{sample}/assembly/quast/canu_medaka'
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
            --conserved-genes-finding \
            --glimmer
        '''

rule quast_canu_medaka_polypolish:
    input:
        assembly=rules.polypolish_polish_canu_medaka.output,
        reference=config['reference'],
        genes=config['genes']
    output:
        'results/{sample}/assembly/quast/canu_medaka_polypolish/report.txt'
    threads:
        config['quast']['threads']
    conda:
        '../envs/quast.yml'
    resources:
        mem_mb=config['quast']['mem'],
        time=config['quast']['time']
    benchmark:
        'benchmarks/{sample}.quast_canu_medaka.tsv'
    params:
        prefix='results/{sample}/assembly/quast/canu_medaka_polypolish'
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
            --conserved-genes-finding \
            --glimmer
        '''

rule quast_wtdbg2_cnsp_to_plot:
    input:
        rules.quast_wtdbg2_cnsp.output
    output:
        'results/{sample}/assembly/quast/wtdbg2/cnsp/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "wtdbg2-cnsp"}}' > {output}
        '''

rule quast_wtdbg2_srp_to_plot:
    input:
        rules.quast_wtdbg2_srp.output
    output:
        'results/{sample}/assembly/quast/wtdbg2/srp/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "wtdbg2-cnsp-srp"}}' > {output}
        '''

rule quast_wtdbg2_medaka_raw_to_plot:
    input:
        rules.quast_wtdbg2_medaka_raw.output
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka/raw/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "wtdbg2-raw_medaka-cnsp"}}' > {output}
        '''

rule quast_wtdbg2_medaka_cnsp_to_plot:
    input:
        rules.quast_wtdbg2_medaka_cnsp.output
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka/cnsp/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "wtdbg2-cnsp_medaka-cnsp"}}' > {output}
        '''


rule quast_wtdbg2_medaka_polypolish_raw_to_plot:
    input:
        rules.quast_wtdbg2_medaka_polypolish_raw.output
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/raw/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "wtdbg2-raw_medaka-cnsp_polypolish-srp"}}' > {output}
        '''

rule quast_wtdbg2_medaka_polypolish_cnsp_to_plot:
    input:
        rules.quast_wtdbg2_medaka_polypolish_cnsp.output
    output:
        'results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/cnsp/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "wtdbg2-cnsp_medaka-cnsp_polypolish-srp"}}' > {output}
        '''


rule quast_canu_medaka_to_plot:
    input:
        rules.quast_canu_medaka.output
    output:
        'results/{sample}/assembly/quast/canu_medaka/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "canu-raw_medaka-cnsp"}}' > {output}
        '''

rule quast_canu_medaka_polypolish_to_plot:
    input:
        rules.quast_canu_medaka_polypolish.output
    output:
        'results/{sample}/assembly/quast/canu_medaka_polypolish/report.mod.tsv'
    threads:
        1
    params:
        samplename='{sample}'
    shell:
        '''
        tail -n+4 {input} | sed 's/ \+ /\t/g' | sed 's/\t$//' | awk -v var="{params.samplename}" '{{OFS=FS="\\t"}}{{print $0 , var, "canu-raw_medaka-cnsp_polypolish-srp"}}' > {output}
        '''

rule quast_combine_to_plot:
    input:
        expand('results/{sample}/assembly/quast/wtdbg2/{type}/report.mod.tsv', sample=df['sample'].tolist(),type=['cnsp', 'srp']),
        expand('results/{sample}/assembly/quast/wtdbg2_medaka/{type}/report.mod.tsv', sample=df['sample'].tolist(),type=['raw', 'cnsp']),
        expand('results/{sample}/assembly/quast/wtdbg2_medaka_polypolish/{type}/report.mod.tsv', sample=df['sample'].tolist(),type=['raw', 'cnsp']),
        expand('results/{sample}/assembly/quast/{type}/report.mod.tsv', sample=df['sample'].tolist(),type=['canu_medaka', 'canu_medaka_polypolish'])
    output:
        'results/all.quast.tsv'
    threads:
        1
    shell:
        '''
        cat {input} > {output}
        '''

rule quast_plot:
    input:
        rules.quast_combine_to_plot.output
    output:
        'results/all.quast.pdf'
    threads:
        1
    conda:
        '../envs/r.yml'
    shell:
        '''
        Rscript workflow/scripts/plotquast.r {input} {output}
        '''