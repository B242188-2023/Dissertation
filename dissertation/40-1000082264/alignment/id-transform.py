import re
gff3_file = 'AeUmbellulata_TA1851_v1.gff3'
output_file = 'rename.txt'

contig_dict = {}
ID = re.compile(r'ID=([^;]+)')
with open(gff3_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        seqid = fields[0]
        attributes = fields[8]
        match = ID.search(attributes)
        if match:
            contig_id = match.group(1)
            if contig_id not in contig_dict:
                contig_dict[contig_id] = seqid
with open(output_file, 'w') as f:
    for contig_id, seqid in contig_dict.items():
        f.write(f'{contig_id}\t{seqid}\n')