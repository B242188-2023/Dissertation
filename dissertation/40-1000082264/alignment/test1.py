import subprocess
import os

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

def read_output_file(output_file_path):
    data = []
    with open(output_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            read_id = parts[0]
            cds_name = parts[1]
            cds_start = int(parts[2])
            cds_end = int(parts[3])
            strand = parts[4]
            data.append((read_id, cds_name, cds_start, cds_end, strand))
    return data

def write_converted_output(file_path, converted_data):
    with open(file_path, 'w') as file:
        for entry in converted_data:
            file.write("\t".join(map(str, entry)) + "\n")

# 读取链文件内容
chain_file_path = 'cdsToGenome.txt'
chain_data = read_chain_file(chain_file_path)

# 解析链文件
mapping = parse_chain(chain_data)

# 读取 output.txt 文件中的信息
samples = ["LD1", "LD2", "LL1", "LL2", "SD1", "SD2", "SL1", "SL2"]
for sample in samples:
    output_file_path = f'{sample}_output.txt'
    #print(output_file_path)
    awk_script = r"""
    {
        read_id = $1;
        chrom = $3;
        start = $4;
        cigar = $6;
        flag = $2;
        len = 0;
        while (match(cigar, /[0-9]+[MIDNSHPX=]/)) {
            if (substr(cigar, RSTART + RLENGTH - 1, 1) ~ /[MDN=X]/) {
                len += substr(cigar, RSTART, RLENGTH - 1);
            }
            cigar = substr(cigar, RSTART + RLENGTH);
        }
        end = start + len - 1;
        strand = (and(flag, 16) ? "-" : "+");
        print read_id, chrom, start, end, strand;
    }
    """
    subprocess.run(f"samtools view {sample}_cds_sorted.bam | awk '{awk_script}' > {output_file_path}", shell=True)


    # 读取 output.txt 文件中的信息
    data = read_output_file(output_file_path)

    converted_data = []
    new_references = set()
    for read_id, cds_name, cds_start, cds_end, strand in data:
        genome_start = convert_coordinates(cds_name, cds_start, mapping)
        genome_end = convert_coordinates(cds_name, cds_end, mapping)

        if genome_start and genome_end:
            converted_data.append((read_id, genome_start[0], genome_start[1], genome_end[1], strand))
            new_references.add(genome_start[0])
        else:
            converted_data.append((read_id, cds_name, cds_start, cds_end, strand))

    # 输出转换后的结果到文件
    converted_output_file_path = f'{sample}_converted_output.txt'
    write_converted_output(converted_output_file_path, converted_data)
    print(f"Coordinates conversion for {sample} completed successfully.")


    def convert_bam_coordinates(bam_file_path, output_bam_file_path, converted_file_path, new_references):
        temp_sam_file = "temp.sam"
        temp_modified_sam_file = "temp_modified.sam"

        # 提取 BAM 文件中的数据到临时 SAM 文件
        with open(temp_sam_file, 'w') as sam_output:
            subprocess.run(["samtools", "view", "-h", bam_file_path], stdout=sam_output)

        # 读取转换后的数据
        converted_data = {}
        with open(converted_file_path, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                read_id = parts[0]
                genome_chr = parts[1]
                genome_start = int(parts[2])
                genome_end = int(parts[3])
                strand = parts[4]
                converted_data[read_id] = (genome_chr, genome_start, genome_end, strand)

        # 修改临时 SAM 文件中的坐标，并更新头部信息
        header_lines = []
        alignment_lines = []
        with open(temp_sam_file, 'r') as sam_input:
            for line in sam_input:
                if line.startswith('@'):
                    header_lines.append(line)
                else:
                    parts = line.strip().split('\t')
                    read_id = parts[0]
                    if read_id in converted_data:
                        genome_chr, genome_start, genome_end, strand = converted_data[read_id]
                        parts[2] = genome_chr
                        parts[3] = str(genome_start + 1)  # SAM 格式中的位置是1-based
                        read_length = len(parts[9])  # 读取序列的长度
                        parts[5] = f"{read_length}M"  # 简单的CIGAR字符串
                        alignment_lines.append('\t'.join(parts) + '\n')
                    else:
                        alignment_lines.append(line)
    
        # 确保所有参考名称都在头部中定义
        existing_refs = {line.split()[1].split(':')[1] for line in header_lines if line.startswith('@SQ')}
        for ref in new_references:
            if ref not in existing_refs:
                header_lines.append(f"@SQ\tSN:{ref}\tLN:1000000\n")  # 假设参考长度为1000000

        # 将修改后的头部和对齐信息写入临时 SAM 文件
        with open(temp_modified_sam_file, 'w') as sam_output:
            sam_output.writelines(header_lines)
            sam_output.writelines(alignment_lines)
    
        # 使用 samtools view 将修改后的临时 SAM 文件转换回 BAM 文件
        with open(output_bam_file_path, 'wb') as bam_output:
            subprocess.run(["samtools", "view", "-bS", temp_modified_sam_file], stdout=bam_output)
    
        # 清理临时文件
        os.remove(temp_sam_file)
        os.remove(temp_modified_sam_file)

    bam_file_path = f'{sample}_cds_sorted.bam'  # 使用提取的 BAM 文件进行测试
    output_bam_file_path = f'{sample}_genome_sorted.bam'
    convert_bam_coordinates(bam_file_path, output_bam_file_path, converted_output_file_path, new_references)
    print(f"Conversion for {sample} BAM file completed successfully.")