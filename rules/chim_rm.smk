rule uchime:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.fasta'
    output:
        uchime_out = 'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.uchime.txt',
        chim = 'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.chimera.fasta',
        nonchim = 'results/finalData/{sample}_{unit}/{sample}_{unit}.clustered100.nonchimera.fasta'
    params:
        minh = config['chim']['minh'],
        mindiffs = config['chim']['mindiffs'],
        mindiv = config['chim']['mindiv'],
        beta = config['chim']['beta'],
        pseudo_c = config['chim']['pseudo_count'],
        abskew = config['chim']['abskew']
    threads: config['general']['cores']
    shell:
        './bin/usearch7 -uchime_denovo {input} -minh {params.minh}'
        ' -mindiffs {params.mindiffs} -mindiv {params.mindiv} -xn {params.beta}'
        ' -dn {params.pseudo_c} -abskew {params.abskew} -uchimeout'
        ' {output.uchime_out} -chimeras {output.chim} -nonchimeras'
        ' {output.nonchim} 2>&1'
