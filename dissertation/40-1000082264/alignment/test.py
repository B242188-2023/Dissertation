def read_chain_file(chain_file_path):
    with open(chain_file_path, 'r') as file:
        chain_data = file.read()
    return chain_data

def parse_chain(chain_data):
    mapping = {}
    for line in chain_data.strip().split("\n"):
        parts = line.split()
        genome_chr = parts[0]
        genome_start = int(parts[1])
        genome_end = int(parts[2])
        genome_strand = parts[3]
        cds_name = parts[4]
        cds_start = int(parts[5])
        cds_end = int(parts[6])
        cds_strand = parts[7]
        
        if cds_name not in mapping:
            mapping[cds_name] = []
        mapping[cds_name].append((genome_chr, genome_start, genome_end, cds_start, cds_end, genome_strand, cds_strand))
    
    return mapping

def convert_coordinates(cds_name, cds_position, mapping):
    if cds_name not in mapping:
        return None
    
    for genome_chr, genome_start, genome_end, cds_start, cds_end, genome_strand, cds_strand in mapping[cds_name]:
        if cds_start <= cds_position <= cds_end:
            if cds_strand == '+' and genome_strand == '+':
                genome_position = genome_start + (cds_position - cds_start)
            elif cds_strand == '-' and genome_strand == '-':
                genome_position = genome_end - (cds_position - cds_start)
            elif cds_strand == '+' and genome_strand == '-':
                genome_position = genome_end - (cds_position - cds_start)
            elif cds_strand == '-' and genome_strand == '+':
                genome_position = genome_start + (cds_position - cds_start)
            return (genome_chr, genome_position)
    
    return None

def read_vcf(vcf_file_path):
    vcf_header = []
    mutations = []
    with open(vcf_file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                vcf_header.append(line)
            else:
                parts = line.strip().split("\t")
                cds_name = parts[0]
                cds_position = int(parts[1])
                mutations.append((cds_name, cds_position, parts))
    return vcf_header, mutations

def write_vcf(vcf_file_path, vcf_header, converted_mutations, include_orig_id=False):
    with open(vcf_file_path, 'w') as file:
        for line in vcf_header:
            file.write(line)
        for mutation in converted_mutations:
            parts = mutation[2]
            if include_orig_id:
                parts.insert(0, mutation[3]) 
                parts[1] = mutation[0] 
                parts[2] = str(mutation[1])
            else:
                parts[0] = mutation[0]  
                parts[1] = str(mutation[1])
            file.write("\t".join(parts) + "\n")

# 读取链文件内容
chain_file_path = 'cdsToGenome.txt'
chain_data = read_chain_file(chain_file_path)

# 解析链文件
mapping = parse_chain(chain_data)

samples = ["LD1", "LD2", "LL1", "LL2", "SD1", "SD2", "SL1", "SL2"]
for sample in samples:
    # 读取 VCF 文件中的突变信息
    vcf_file_path = f'{sample}_snp.vcf'
    vcf_header, mutations = read_vcf(vcf_file_path)

    # 转换突变坐标
    converted_mutations = []
    for cds_name, cds_position, original_parts in mutations:
        result = convert_coordinates(cds_name, cds_position, mapping)
        if result:
            genome_chr, genome_position = result
            converted_mutations.append((genome_chr, genome_position, original_parts, cds_name))
        else:
            converted_mutations.append((cds_name, cds_position, original_parts, cds_name))

    # 输出转换后的结果到两个文件
    output_vcf_file_path = f'{sample}_output.vcf'
    output_vcf_with_id_file_path = f'{sample}_output_with_orig_id.vcf'
    write_vcf(output_vcf_file_path, vcf_header, converted_mutations, include_orig_id=False)
    write_vcf(output_vcf_with_id_file_path, vcf_header, converted_mutations, include_orig_id=True)

