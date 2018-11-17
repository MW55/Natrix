rule vsearch:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.fasta'
    output:
        uchime_out = 'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.uchime.txt',
        chim = 'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.chimera.fasta',   
        nonchim = 'results/finalData/{sample}_{unit}.clustered100.nonchimera.fasta' 
    params:
        beta = config['chim']['beta'],
        pseudo_c = config['chim']['pseudo_count'],
        abskew = config['chim']['abskew']
    threads: config['general']['cores']                                                         
    conda:                                                                                      
         '../envs/vsearch.yml'                                                                  
    shell:
        'vsearch --uchime3_denovo {input} -uchimeout {output.uchime_out}'                       
        ' -chimeras {output.chim} -nonchimeras {output.nonchim} -xn {params.beta}'              
        '  -dn {params.pseudo_c} -abskew {params.abskew} 2>&1'