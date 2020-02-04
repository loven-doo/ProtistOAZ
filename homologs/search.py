import os
import sys
import json

import eagledb
import eagle
from eagledb.scheme import SeqProfileInfo, GenomeInfo
from eagle.lib.alignment import construct_mult_aln, MultAln
from eagle.lib.orthan import hom_search_profile


NUM_THREADS = 4

REPO_DATA_DIR = os.path.join("..", "repo_data")
OAZ_REPR_FASTA = os.path.join(REPO_DATA_DIR, "oaz_repr.fasta")
OAZ_REPR_ALN = os.path.join(REPO_DATA_DIR, "oaz_repr_aln.fasta")
OAZ_REPR_HMM = os.path.join(REPO_DATA_DIR, "oaz_repr.hmm")
OAZ_REPR_INFO = os.path.join(REPO_DATA_DIR, "oaz_repr_info.json")
IMICROBE_DATA_DIR = os.path.join("..", "imicrobe_data", "prepared")
RNA_FASTA = os.path.join(IMICROBE_DATA_DIR, "18s.fasta")
REF_TREE_NWK = os.path.join(REPO_DATA_DIR, "18S_rRNA_tree.nwk")
TRANSCRIPTOMES_PATH = os.path.join(IMICROBE_DATA_DIR, "transcriptomes.json")

DATA_DIR = "data"
RESULT_DIR = "result"


def build_ref_tree(ref_aln_list, seqs2orgs=None, ref_tree_nwk=None):
    mult_aln = ref_aln_list[0]
    assert isinstance(mult_aln, MultAln), "ERROR: items of 'ref_aln_list' should be an objects of " \
                                          "class 'eagle.lib.alignment.MultAln'"
    mult_aln.improve_aln(inplace=True)
    if seqs2orgs is not None:
        mult_aln.remove_paralogs(seq_ids_to_orgs=seqs2orgs, inplace=True)
        mult_aln.improve_aln(inplace=True)
    phylo_tree = mult_aln.build_tree(tree_name=mult_aln.aln_name, method="FastME", options={"-b": 100})
    if seqs2orgs is not None:
        phylo_tree.full_seq_names = seqs2orgs
        phylo_tree.set_full_names(inplace=True) 
    if ref_tree_nwk is not None:
        return phylo_tree.dump_tree(tree_path=ref_tree_nwk, tree_format="newick")
    else:
        return phylo_tree.newick


def build_oaz_repr_profile(oaz_repr_fasta, oaz_repr_hmm, oaz_repr_info=None, oaz_repr_aln=None, 
                           num_threads=NUM_THREADS):
    mult_aln = construct_mult_aln(fasta_path=oaz_repr_fasta, 
                                  method="MSAProbs", 
                                  aln_name="OAZ_repr", 
                                  aln_type="prot",
                                  num_threads=NUM_THREADS)
    mult_aln.improve_aln(inplace=True)
    if oaz_repr_aln is not None:
        mult_aln.dump_alignment(aln_path=oaz_repr_aln, aln_format="fasta")
    mult_aln.get_hmm_profile(profile_path=oaz_repr_hmm, method="hmmer")
    profile_info = SeqProfileInfo(name="OAZ_repr", path=oaz_repr_hmm, seq_type="prot")
    if oaz_repr_info is not None:
        with open(oaz_repr_info, "w") as profile_info_f:
            json.dump(profile_info.get_json(), profile_info_f, indent=2)
        return oaz_repr_info
    else:
        return profile_info