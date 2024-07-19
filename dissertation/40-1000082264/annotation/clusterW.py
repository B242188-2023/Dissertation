import pysam
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

gene_vcf_mapping = {
    'LD_SD_genes.csv': ['LD1_snp.vcf', 'LD2_snp.vcf', 'SD1_snp.vcf', 'SD2_snp.vcf'],
    'LL_SL_genes.csv': ['LL1_snp.vcf', 'LL2_snp.vcf', 'SL1_snp.vcf', 'SL2_snp.vcf'],
    'SD_SL_genes.csv': ['SD1_snp.vcf', 'SD2_snp.vcf', 'SL1_snp.vcf', 'SL2_snp.vcf'],
    'LD_LL_genes.csv': ['LD1_snp.vcf', 'LD2_snp.vcf', 'LL1_snp.vcf', 'LL2_snp.vcf']
}
#print(gene_vcf_mapping.keys())

vcf_dir = '../alignment/'

def strip_version(gene_id):
    parts = gene_id.split('.')
    if len(parts) > 4 and parts[-1] == '1':
        return gene_id  # 保留包含.1版本号的ID
    elif len(parts) > 4 and parts[-1] == '2':
        return gene_id
    elif len(parts) == 4:
        return gene_id  # 保留不带版本号的ID
    return None  # 如果不是.1的转录本，则返回None

def get_main_id(gene_id):
    parts = gene_id.split('.')
    return '.'.join(parts[:4])  # 返回前四部分

gene_sequences = {}
missing_gene_ids = {}

for gene_file in gene_vcf_mapping.keys():
    with open(gene_file) as f:
        gene_ids = set(get_main_id(line.strip().replace('"', '')) for line in f)
        print(f"{gene_file} 中提取的基因ID: {len(gene_ids)} 个")  # 打印提取的基因ID数量
        
        sequences = {}
        main_id_to_record = {}
        for record in SeqIO.parse('../fastq/AeUmbellulata_TA1851_v1-cds.fasta', 'fasta'):
            main_id = get_main_id(record.id)
            if main_id in gene_ids:
                if main_id not in main_id_to_record:
                    main_id_to_record[main_id] = record
                elif record.id.endswith(".1"):
                    main_id_to_record[main_id] = record

        # 将记录添加到序列中
        for main_id, record in main_id_to_record.items():
            sequences[record.id] = record.seq
        
        # 找出没有匹配到的基因ID
        missing_ids = gene_ids - set(main_id_to_record.keys())
        # 尝试使用 .2 版本来补充缺失的基因ID
        for record in SeqIO.parse('../fastq/AeUmbellulata_TA1851_v1-cds.fasta', 'fasta'):
            main_id = get_main_id(record.id)
            if main_id in missing_ids and record.id.endswith(".2"):
                sequences[record.id] = record.seq
                missing_ids.remove(main_id)
        
        missing_gene_ids[gene_file] = missing_ids
        gene_sequences[gene_file] = sequences
        #print(f"{gene_file} 中提取的序列: {list(sequences.keys())}")

        # 打印提取的序列ID数量和缺失的基因ID数量
        print(f"{gene_file} 中提取的序列: {len(sequences)} 个")
        #print(f"{gene_file} 中缺失的基因ID: {len(missing_ids)} 个")
'''
# 输出缺失的基因ID
for gene_file, missing_ids in missing_gene_ids.items():
    if missing_ids:
        print(f"{gene_file} 中缺失的基因ID:")
        for missing_id in missing_ids:
            print(missing_id)
'''
    # 打印从CSV文件中提取的基因ID
    #print(f"{gene_file} 中提取的基因ID: ")
    #for gene_id in gene_ids:
        #print(gene_id)
# 输出缺失的基因ID



# 保存突变前的蛋白质序列
for gene_file in gene_vcf_mapping.keys():
    protein_sequences = {}
    sequences = gene_sequences[gene_file]
    for gene_id, dna_seq in sequences.items():
        protein_seq = dna_seq.translate(to_stop=True)
        protein_sequences[gene_id] = protein_seq

    # 保存蛋白质序列到FASTA文件
    protein_records = [SeqRecord(Seq(seq), id=gene_id + "_original", description="") for gene_id, seq in protein_sequences.items()]
    output_file = gene_file.replace('genes.csv', 'original_protein_sequences.fasta')
    SeqIO.write(protein_records, output_file, 'fasta')
print("突变前的蛋白质序列已保存到对应的FASTA文件中。")


for gene_file, vcf_files in gene_vcf_mapping.items():
    sequences = gene_sequences[gene_file]
    for vcf_file in vcf_files:
        vcf_path = os.path.join(vcf_dir, vcf_file)
        vcf_reader = pysam.VariantFile(vcf_path)
        #print(f"处理 {vcf_file} 中的变异：")
        for record in vcf_reader.fetch():
            #print(record)
            gene_id = strip_version(record.chrom)
            #print(gene_id)  # 打印从VCF文件中提取的基因ID
            if gene_id in sequences:
                seq = sequences[gene_id]
                for alt in record.alts:
                    if len(record.ref) == 1 and len(alt) == 1:
                        # SNP
                        seq = seq[:record.pos-1] + alt + seq[record.pos:]
                sequences[gene_id] = seq
    gene_sequences[gene_file] = sequences

for gene_file in gene_vcf_mapping.keys():
    protein_sequences = {}
    sequences = gene_sequences[gene_file]
    for gene_id, dna_seq in sequences.items():
        protein_seq = dna_seq.translate(to_stop=True)
        protein_sequences[gene_id] = protein_seq

    # 保存蛋白质序列到FASTA文件
    protein_records = [SeqRecord(Seq(seq), id=gene_id, description="") for gene_id, seq in protein_sequences.items()]
    output_file = gene_file.replace('genes.csv', 'mutated_protein_sequences.fasta')
    SeqIO.write(protein_records, output_file, 'fasta')
print("DNA序列已提取、修改并翻译成蛋白质序列, 结果已分别保存到对应的FASTA文件中。")





import subprocess
input_files = {
    'LD_SD_original': 'LD_SD_original_protein_sequences.fasta',
    'LD_SD_mutated': 'LD_SD_mutated_protein_sequences.fasta',
    'LL_SL_original': 'LL_SL_original_protein_sequences.fasta',
    'LL_SL_mutated': 'LL_SL_mutated_protein_sequences.fasta'
}
output_dir = 'clustalw_results'

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

# 读取序列文件
sequences = {}
for label, file in input_files.items():
    sequences[label] = {record.id: record for record in SeqIO.parse(file, 'fasta')}

# 去掉_original后缀以便于匹配
def strip_suffix(gene_id):
    if gene_id.endswith('_original'):
        #print(gene_id[:-9])
        return gene_id[:-9]
    #print(gene_id)
    return gene_id

# 获取LD_SD的基因ID
ld_sd_ids = set(strip_suffix(gene_id) for gene_id in sequences['LD_SD_original'].keys()).union(
             strip_suffix(gene_id) for gene_id in sequences['LD_SD_mutated'].keys())


# 获取LL_SL的基因ID
ll_sl_ids = set(strip_suffix(gene_id) for gene_id in sequences['LL_SL_original'].keys()).union(
             strip_suffix(gene_id) for gene_id in sequences['LL_SL_mutated'].keys())


print(f"LD_SD_original: {len(sequences['LD_SD_original'])} IDs")
print(f"LD_SD_mutated: {len(sequences['LD_SD_mutated'])} IDs")
print(f"LL_SL_original: {len(sequences['LL_SL_original'])} IDs")
print(f"LL_SL_mutated: {len(sequences['LL_SL_mutated'])} IDs")
print(f"LD_SD unique IDs: {len(ld_sd_ids)}")
print(f"LL_SL unique IDs: {len(ll_sl_ids)}")

log_file = 'clustalw.log'
# 将ClustalW的输出重定向到日志文件
with open(log_file, 'w') as log:

    # 对LD_SD进行比对
    for gene_id in ld_sd_ids:
        original_seq = sequences['LD_SD_original'].get(gene_id + '_original')
        mutated_seq = sequences['LD_SD_mutated'].get(gene_id)
        
        # 打印调试信息
        log.write(f"Processing gene_id: {gene_id}\n")
        if original_seq:
            log.write(f"Original: {original_seq.id}\n")
        else:
            log.write("Original: None\n")
        
        if mutated_seq:
            log.write(f"Mutated: {mutated_seq.id}\n")
        else:
            log.write("Mutated: None\n")
        
        if original_seq and mutated_seq:
            # 创建LD_SD的比对输入文件
            alignment_input_file = os.path.join(output_dir, f'{gene_id}_LD_SD.fasta')
            with open(alignment_input_file, 'w') as f:
                SeqIO.write([original_seq, mutated_seq], f, 'fasta')
            
            # 运行ClustalW
            alignment_output_file = os.path.join(output_dir, f'{gene_id}_LD_SD.aln')
            clustalw_command = f'clustalw -INFILE={alignment_input_file} -OUTFILE={alignment_output_file} -OUTPUT=FASTA'
            log.write(f'运行命令: {clustalw_command}\n')
            result = subprocess.run(clustalw_command, shell=True, capture_output=True, text=True)
            log.write(result.stdout)
            if result.returncode != 0:
                log.write(f'错误: ClustalW命令执行失败，基因ID: {gene_id}\n')
                log.write(result.stderr)
            else:
                log.write(f'成功: ClustalW命令执行成功，基因ID: {gene_id}\n')
            
            # 删除临时输入文件
            os.remove(alignment_input_file)

    # 对LL_SL进行比对
    for gene_id in ll_sl_ids:
        original_seq = sequences['LL_SL_original'].get(gene_id + '_original')
        mutated_seq = sequences['LL_SL_mutated'].get(gene_id)

        # 打印调试信息
        log.write(f"Processing gene_id: {gene_id}\n")
        if original_seq:
            log.write(f"Original: {original_seq.id}\n")
        else:
            log.write("Original: None\n")
        
        if mutated_seq:
            log.write(f"Mutated: {mutated_seq.id}\n")
        else:
            log.write("Mutated: None\n")

        if original_seq and mutated_seq:
            # 创建LL_SL的比对输入文件
            alignment_input_file = os.path.join(output_dir, f'{gene_id}_LL_SL.fasta')
            with open(alignment_input_file, 'w') as f:
                SeqIO.write([original_seq, mutated_seq], f, 'fasta')
            
            # 运行ClustalW
            alignment_output_file = os.path.join(output_dir, f'{gene_id}_LL_SL.aln')
            clustalw_command = f'clustalw -INFILE=' + alignment_input_file + ' -OUTFILE=' + alignment_output_file + ' -OUTPUT=FASTA'
            log.write(f'运行命令: {clustalw_command}\n')
            result = subprocess.run(clustalw_command, shell=True, capture_output=True, text=True)
            log.write(result.stdout)
            if result.returncode != 0:
                log.write(f'错误: ClustalW命令执行失败，基因ID: {gene_id}\n')
                log.write(result.stderr)
            else:
                log.write(f'成功: ClustalW命令执行成功，基因ID: {gene_id}\n')
            
            # 删除临时输入文件
            os.remove(alignment_input_file)

print("所有比对已完成。")

# 检查输出文件数量
ld_sd_output_files = [f for f in os.listdir(output_dir) if f.endswith('_LD_SD.aln')]
ll_sl_output_files = [f for f in os.listdir(output_dir) if f.endswith('_LL_SL.aln')]

print(f"LD_SD 比对文件数量: {len(ld_sd_output_files)}")
print(f"LL_SL 比对文件数量: {len(ll_sl_output_files)}")
