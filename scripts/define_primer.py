import pandas

# Script to define the primer used by pandaseq. The defined
# primer depends on the existance off the barcode, if the reads are
# single or paired end and if pandaseq should use an offset of
# the primer length instead of the sequence. Using an offset can
# be useful if the sequence has many uncalled bases in the primer
# region, preventing a nucleotide primer from matching.

primertable = pd.read_csv(snakemake.input, index_col = 'Probe')

if snakemake.params.paired_end:
    if snakemake.params.offset:
        if snakemake.params.bar_removed:
            # I took it like it was in the old pipeline, should't the poly N for the forward primer also get in there???
            primertable['f_primer'] = (
                primertable['specific_forward_primer'].str.len())
            primertable['r_primer'] = (
                primertable[['poly_N_rev',
                'specific_reverse_primer']].sum(axis=1).str.len())
        else:
            primertable['f_primer'] = (
                primertable[['poly_N', 'Barcode_forward',
                'specific_forward_primer']].sum(axis=1).str.len())
            primertable['r_primer'] = (
                primertable[['poly_N_rev','Barcode_reverse',
                'specific_reverse_primer']].sum(axis=1).str.len())
    else:
        if snakemake.params.bar_removed:
            primertable['f_primer'] = (
                primertable['specific_forward_primer'])
            primertable['r_primer'] = (
                primertable['specific_reverse_primer'])
        else:
            primertable['f_primer'] = (
                primertable[['Barcode_forward',
                'specific_forward_primer']].sum(axis=1))
            primertable['r_primer'] = (
                primertable[['Barcode_reverse',
                'specific_reverse_primer']].sum(axis=1))
else:
    if snakemake.params.offset:
        if snakemake.params.bar_removed:
            # I took it like it was in the old pipeline, should't the poly N for the forward primer also get in there???
            primertable['f_primer'] = (
                primertable['specific_forward_primer'].str.len())
        else:
            primertable['f_primer'] = (
                primertable[['poly_N', 'Barcode_forward',
                'specific_forward_primer']].sum(axis=1).str.len())
    else:
        if snakemake.params.bar_removed:
            primertable['f_primer'] = (
                primertable['specific_forward_primer'])
        else:
            primertable['f_primer'] = (
                primertable[['Barcode_forward',
                'specific_forward_primer']].sum(axis=1))

primertable.to_csv(snakemake.output)
