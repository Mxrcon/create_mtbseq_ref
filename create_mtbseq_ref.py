
from Bio import GenBank
from Bio import SeqIO
import sys
import os

#Determinando os inputs e outputs
gbk_file = sys.argv[1]


def create_fasta(output_name,sample_sequence):
    # Function to create a fasta output as required by mtbseq
    with open(output_name+".fasta", "w") as output_fasta:
        output_fasta.write(">{0}\n{1}\n".format(output_name, sample_sequence))
        output_fasta.close()

def create_genes_txt(output_name, feature_list):
    with open(output_name+"_genes.txt","w") as output_txt:
            genes_txt_header = "# ID      name    start   stop    frame   product description     function        cogcats status_region   status_function type    region_number   function_number"
            output_txt.write(genes_txt_header)
            output_txt.write("\n".join(sorted(feature_list)))
            output_txt.close()
def get_frame(start_pos, end_pos, strand):
    frame = min(start_pos,end_pos)%3
    if frame == 0:
        frame = 3
    return str(frame*strand)


#getting sample information
with open(gbk_file) as handle:
    record =  GenBank.read(handle)
    sample_date = record.date
    sample_name = "_".join(str(record.organism).split(" "))
    sample_sequence = record.sequence
    handle.close()
    
output_name = sample_name+sample_date

create_fasta(output_name, sample_sequence)


#Parsing objects
feature_list = []
for seq_record in SeqIO.parse(gbk_file , 'genbank'):
    for feature in seq_record.features:
        if feature.type == ["CDS","rRNA","tRNA"]:
            if "gene" in feature.qualifiers:
                feature_gene = feature.qualifiers['gene'][0]
            else:
                feature_gene = '-'
            frame = get_frame(feature.location.start+1,feature.location.end,feature.location.strand)
            if int(frame) <0:
                feature_start = feature.location.end
                feature_end = feature.location.start+1
            else: 
                feature_start,feature_end = feature.location.start+1, feature.location.end
            
            if feature.type == "CDS":
                feature_product = feature.qualifiers['product'][0]
            else: 
                feature_product = feature.type
            
            if feature.type =="rRNA":
                feature_function= "rRNA"
            else: 
                feature_function = "-"
            plain_feature = [feature.qualifiers['locus_tag'][0],
                             feature_gene,
                             str(feature_start),
                             str(feature_end),
                             frame,
                             feature_product,
                             feature_function,
                             '-',
                             '-',
                             'status 3',
                             'annotated',
                             feature.type,
                             '5',
                             '6'                                                
                    ]
            feature_list.append("   ".join(plain_feature))
        
create_genes_txt(output_name,feature_list)


