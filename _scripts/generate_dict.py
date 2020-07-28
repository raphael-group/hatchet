import gzip
from Bio import SeqIO

with gzip.open("/media/vineetb/t5-vineetb/raphael-group/hg19.fa.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        chromosome = str(record.id)
        if not chromosome.upper().startswith('CHR'):
            chromosome = 'chr' + chromosome
        print('@SQ SN:' + chromosome + ' LN:' + str(len(record)))
