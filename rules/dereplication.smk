import dinopy
from glob import glob
from collections import Counter

rule dereplication:
    input:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.fasta'
    output:
        'results/assembly/{sample}_{unit}/{sample}_{unit}.clustered100.fasta'
    run:
        seqs = dinopy.FastaReader(str(input))
        sorted_counter = Counter(seqs.reads()).most_common()
        sample_name = wildcards.sample + '_' + wildcards.unit
        with dinopy.FastaWriter(str(output), line_width=1000) as clust:
            for i in range(len(sorted_counter)):
                name = '{}_{};size={};'.format(sample_name, i+1, sorted_counter[i][1])
                clust.write_entry((sorted_counter[i][0],bytes(name, 'utf-8')))
            clust.close()


