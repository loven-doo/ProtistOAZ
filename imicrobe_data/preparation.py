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
    
    prepare_data()
    

def prepare_data(fna_path=FNA_PATH, 
                 rna_18s_path=RNA_18S_PATH, 
                 fna_meta_path=FNA_META_PATH, 
                 fna_meta_sep=FNA_META_SEP,
                 transcriptomes_path=TRANSCRIPTOMES_PATH,
                 prepared_fna_path=PREPARED_FNA_PATH,
                 prepared_18s_path=PREPARED_18S_PATH):
    
    samples_dict = get_samples_dict(fna_meta_path=fna_meta_path, sep=fna_meta_sep)

    rna_seqs = SeqsDict.load_from_file(rna_18s_path, seqs_format="fasta", low_memory=False)
    rna_names_conv = convert_rna_names(rna_seqs, samples_dict)
    print("\n")
    fna_seqs = SeqsDict.load_from_file(fna_path, seqs_format="fasta", low_memory=True)
    fna_names_conv = convert_fna_names(fna_seqs, rna_names_conv)

    transcriptomes = list()
    fna_ids_conv = dict()
    for sample_name in fna_names_conv:
        transcriptomes.append(GenomeInfo(genome_id=sample_name, 
                                         org_name=samples_dict[sample_name],
                                         fna_path=prepared_fna_path,
                                         fna_id_list=list(fna_names_conv[sample_name].keys())).get_json())
        fna_ids_conv.update(fna_names_conv[sample_name])
    with open(transcriptomes_path, "w") as transcriptomes_f:
        json.dump(transcriptomes, transcriptomes_f)

    rna_seqs.rename_seqs({rna_names_conv[t_name]: t_name for t_name in fna_names_conv})
    rna_seqs.get_sample(list(fna_names_conv.keys()), low_memory=False).dump(prepared_18s_path, seqs_format="fasta")
    fna_seqs.rename_seqs({fna_ids_conv[t_name]: t_name for t_name in fna_ids_conv})
    fna_seqs.get_sample(list(fna_ids_conv.keys()), low_memory=True).dump(prepared_fna_path, seqs_format="fasta")


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


def convert_fna_names(in_fna_seqs, sample_names):
    transformed_names = defaultdict(dict)
    for seq_name in in_fna_seqs:
        transformed_name = None
        sample_name = None
        transformed_name = seq_name.split(" ")[1].split("=")[-1]
        sample_name = transformed_name.split("-")[0]
        if sample_name in sample_names:
            if transformed_name in transformed_names[sample_name]:
                l = len(transformed_names[sample_name])
                transformed_names[sample_name][transformed_name+"_"+str(l)] = seq_name
            else:
                transformed_names[sample_name][transformed_name] = seq_name
        else:
            print("sample '%s' is absent in sample_names" % transformed_name)
    return transformed_names
