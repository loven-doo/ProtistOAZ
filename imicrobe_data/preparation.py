import os
import json
from collections import defaultdict
import argparse

import numpy as np
import pandas as pd

import eagledb
import eagle
from eagledb.scheme import GenomeInfo
from eagle.lib.seqs import SeqsDict


SOURCES_DIR = "source"
FNA_PATH = os.path.join(SOURCES_DIR, "CAM_P_0001000.nt.fa")
FNA_META_PATH = os.path.join(SOURCES_DIR, "fna_meta.csv")
FNA_META_SEP = ","
RNA_18S_PATH = os.path.join(SOURCES_DIR, "18s.fasta")
PREPARED_DIR = "prepared"
PREPARED_18S_PATH = os.path.join(PREPARED_DIR, "18s.fasta")
PREPARED_FNA_PATH = os.path.join(PREPARED_DIR, "transcriptomes.fna")
TRANSCRIPTOMES_PATH = os.path.join(PREPARED_DIR, "transcriptomes.json")


def prepare_data_cmd():
    parser = argparse.ArgumentParser()

    parser.add_argument("-fna",
                        "--fna-path",
                        help="Path to source Imicrobe project fna",
                        required=False,
                        default=FNA_PATH)
    parser.add_argument("-meta",
                        "--fna-meta-path",
                        help="Path to source Imicrobe project meta table",
                        required=False,
                        default=FNA_META_PATH)
    parser.add_argument("-rna",
                        "--rna-18s-path",
                        help="Path to source Imicrobe project 18S rRNA sequences",
                        required=False,
                        default=RNA_18S_PATH)
    parser.add_argument("-tres",
                        "--transcriptomes-path",
                        help="Path to transcriptomes info json (the result of preparation)",
                        required=False,
                        default=TRANSCRIPTOMES_PATH)
    parser.add_argument("-pfna",
                        "--prepared-fna-path",
                        help="Path to prepared fna",
                        required=False,
                        default=PREPARED_FNA_PATH)
    parser.add_argument("-prna",
                        "--prepared-18s-path",
                        help="Path to prepared 18S rRNA sequences",
                        required=False,
                        default=PREPARED_18S_PATH)
    
    cmd_args = parser.parse_args()

    prepare_data(**cmd_args.__dict__)
    

def prepare_data(fna_path=FNA_PATH, 
                 rna_18s_path=RNA_18S_PATH, 
                 fna_meta_path=FNA_META_PATH, 
                 fna_meta_sep=FNA_META_SEP,
                 transcriptomes_path=TRANSCRIPTOMES_PATH,
                 prepared_fna_path=PREPARED_FNA_PATH,
                 prepared_18s_path=PREPARED_18S_PATH):
    
    samples_dict = get_samples_dict(fna_meta_path=fna_meta_path, sep=fna_meta_sep)

    rna_seqs = SeqsDict.load_from_file(rna_18s_path, seqs_format="fasta", low_memory=False)
    print("18S rRNA sequences read")
    rna_names_conv = convert_rna_names(rna_seqs, samples_dict)
    print("\n")
    fna_names_conv = convert_fna_names(in_fna_path=fna_path, 
                                       sample_names=rna_names_conv, 
                                       out_fna_path=prepared_fna_path)
    print("Transcriptome sequences read")

    transcriptomes = list()
    fna_ids_conv = dict()
    for sample_name in fna_names_conv:
        transcriptomes.append(GenomeInfo(genome_id=sample_name, 
                                         org_name=samples_dict[sample_name],
                                         fna_path=prepared_fna_path,
                                         fna_id_list=list(fna_names_conv[sample_name])).get_json())
    with open(transcriptomes_path, "w") as transcriptomes_f:
        json.dump(transcriptomes, transcriptomes_f, indent=2)
    print("%s trascriptomes prepared" % len(transcriptomes))

    rna_seqs.rename_seqs({rna_names_conv[t_name]: t_name for t_name in fna_names_conv})
    rna_seqs.get_sample(list(fna_names_conv.keys()), low_memory=False).dump(prepared_18s_path, seqs_format="fasta")

def get_samples_dict(fna_meta_path, sep):
    fna_meta_df = pd.read_csv(fna_meta_path, sep=sep)[["BIOMATERIAL_NAME", "SAMPLE_DESCRIPTION"]]
    return dict(filter(lambda s: s[0], fna_meta_df.apply(prepare_sample_meta, axis=1)))


def prepare_sample_meta(sample_meta):
    s_name = sample_meta["BIOMATERIAL_NAME"]
    s_descr = sample_meta["SAMPLE_DESCRIPTION"]
    if s_name[:7] == "(MMETSP":
        return s_name[:12].strip("()"), s_descr.split("(")[0].strip().replace(" ", "_")
    else:
        return None, None
    

def convert_rna_names(in_rna_seqs, sample_names):   
    transformed_names = dict()
    for seq_name in in_rna_seqs:
        transformed_name = None
        transformed_name = seq_name.split("|")[0]
        if transformed_name in sample_names:
            transformed_names[transformed_name] = seq_name
        else:
            print("sample '%s' is absent in sample_names" % seq_name)
    return transformed_names


def convert_fna_names(in_fna_path, sample_names, out_fna_path):
    sample_fna_ids = defaultdict(set)
    with open(out_fna_path, "w") as out_fna_f:
        with open(in_fna_path) as in_fna_f:
            seq_list = list()
            transformed_name = None
            sample_name = None
            for line_ in in_fna_f:
                line = None
                line = line_.strip()
                if not line:
                    continue
                if line[0] == ">":
                    if transformed_name is not None:
                        if sample_name in sample_names:
                            if transformed_name in sample_fna_ids[sample_name]:
                                l = len(sample_fna_ids[sample_name])
                                sample_fna_ids[sample_name].add(transformed_name+"_"+str(l))
                                out_fna_f.write(">"+transformed_name+"_"+str(l)+"\n")
                            else:
                                sample_fna_ids[sample_name].add(transformed_name)
                                out_fna_f.write(">"+transformed_name+"\n")
                            out_fna_f.write("".join(seq_list)+"\n")
                        else:
                            print("sample '%s' is absent in sample_names" % transformed_name)
                        transformed_name = None
                        sample_name = None
                        seq_list = list()
                    transformed_name_list = line.split(" ")[1].split("=")[-1].split("-")
                    transformed_name = transformed_name_list[0][:10] + "-" + transformed_name_list[1]
                    sample_name = transformed_name.split("-")[0]
                else:
                    seq_list.append(line)
                    
            if transformed_name is not None:
                if sample_name in sample_names:
                    if transformed_name in sample_fna_ids[sample_name]:
                        l = len(sample_fna_ids[sample_name])
                        sample_fna_ids[sample_name].add(transformed_name+"_"+str(l))
                        out_fna_f.write(">"+transformed_name+"_"+str(l)+"\n")
                    else:
                        sample_fna_ids[sample_name].add(transformed_name)
                        out_fna_f.write(">"+transformed_name+"\n")
                    out_fna_f.write("".join(seq_list)+"\n")
                else:
                    print("sample '%s' is absent in sample_names" % transformed_name)
                transformed_name = None
                sample_name = None
                seq_list = list()
        
    return sample_fna_ids
