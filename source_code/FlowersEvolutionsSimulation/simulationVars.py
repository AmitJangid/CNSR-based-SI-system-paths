from enum import Enum
from typing import Dict, Any, Tuple, Union
import numpy as np

class Allele_evolution_item:
    def __init__(self, allele_ind_List, allele_origin_List, birth_iteration_List):
        self.Allele_ind_List = allele_ind_List
        self.Allele_origin_List = allele_origin_List
        self.Birth_iteration_List = birth_iteration_List

class PartnerDetails:
    def __init__(self, iteration, partners, allele_freqs, partners_freqs, partners_haplotype, partners_haplotype_freq, allele_in_haplotypes):
        self.Iteration = iteration
        self.Partners = partners
        self.AlleleFreqs = allele_freqs
        self.PartnersFreqs = partners_freqs
        self.PartnersHaplotype = partners_haplotype
        self.PartnersHaplotypeFreq = partners_haplotype_freq
        self.AlleleInHaplotypes = allele_in_haplotypes

class AA_letter(Enum):
    # 1 - 4 AA letters:
    # 0 - (H) - Hydrophbic,
    # 1 - (P) - neural polar,
    # 2 - (+) - (Plus) - positively charged
    # 3 - (-) - (Minus) - negatively charged
    # 4 - (S) - (Stop Codon) - signals the termination of the translation process of the current protein

    H = 0
    P = 1
    Plus = 2
    Minus = 3
    Stop = 4

class AA:
    def __init__(self, letter):
        self.Letter = letter
        # self.Codon = codon
        # self.Item = item

class Allele:
    '''
    an array of L = 10 AAs
    '''
    def __init__(self, allele, indexInRNasesPool, indexInSLFsPool, slf_allele_count, rnase_allele_count,
                 birthIter, birthReason, deathIter, deathReason, previousAllele, interactWith):
        self.Allele = allele
        self.IndexInRNasesPool = indexInRNasesPool
        self.IndexInSLFsPool = indexInSLFsPool
        #self.SLF_allele_count = slf_allele_count
        #self.RNase_allele_count = rnase_allele_count
        #self.BirthIter = birthIter
        #self.BirthReason = birthReason #Initiation, duplication, mutation
        #self.DeathIter = deathIter
        #self.DeathReason = deathReason  # mutation, deletion
        self.PreviousAllele = previousAllele
        #self.InteractWith = interactWith # a list of current active alleles that specific interact with this allele

class Haplotype:
    def __init__(self, female_part, male_part, haplotype_index_in_pool, Ancestor_hap_index,
             RNasesIndices, SLFsIndices, PreviousHaplotype,
             birthIter, birthReason, is_SC, SLFs_len, RNases_len, SLFs_unique_len, SLFsIndices_unique,
             SLFsAncestralIndices, RNasesAncestralIndices):
        self.Female_part = female_part
        self.Male_part = male_part
        self.Haplotype_index_in_pool = haplotype_index_in_pool
        # self.Haplotype_count = haplotype_count
        self.RNasesIndices = RNasesIndices
        self.SLFsIndices = SLFsIndices
        self.PreviousHaplotype = PreviousHaplotype
        # self.FollowingHaplotypesLst = FollowingHaplotypesLst
        self.BirthIter = birthIter
        self.BirthReason = birthReason
        self.IS_SC = is_SC
        # self.SC_type = sc_type
        self.SLFs_len = SLFs_len
        self.RNases_len = RNases_len
        self.SLFs_unique_len = SLFs_unique_len
        self.SLFsIndices_unique = SLFsIndices_unique
        self.SLFsAncestralIndices = SLFsAncestralIndices
        self.RNasesAncestralIndices = RNasesAncestralIndices
        self.Ancestor_hap_index = Ancestor_hap_index

# 1/0 True/False buttons
ToyExample = False #True #
enable_mutations = 1

# initiate haplotypes set compatibility:
one_on_one = True #False #
init_constrained_RNases = True#False #

if init_constrained_RNases:
    initiation_RNases_costrained_to_log_file = "First, whole SI haplotypes set was initiated, only the RNases were kept"
else:
    initiation_RNases_costrained_to_log_file = "First, only the RNases set was initiated"
if one_on_one:
    initiation_text_to_log_file = "Finally, SLFs that were interacting with one and only one RNase were added to the SI haplotypes"
else:
    initiation_text_to_log_file = "Finally, SLFs set was initiated, non-compatible SLFs were added to SI haplotypes, untill all the haplotypes were fully cross-compatible"

# If I want to use the same population again
use_seed = 1
generate_new_population = 1
test_freqs = 0
debug_mode = 1
plot_stuff = 0
initiate_bioclasses_freq = "bioclasses_by_UNIPROT" #"bioclasses_by_BLOSUM" #"Uniform" #"input" #

#general population parameters
num_of_iters = 100000
restart_from_iter = 0
popu_size = 500  # population size, needs to be eve
seq_len = 18  # length of each allele
num_types_of_haplotypes = 10
slf_amount_in_haploid = num_types_of_haplotypes - 1 #
added_slfs = 0
alphabet_size = 4  # size of alphabet for positions in alleles
num_of_codons = 64 # (= 4^3) is the amount of all possible combinations of 3 letters length sequence out of 4 possible letters (A,C,G,U)

# SC parameters
alpha = 0.95 #np.arange(0, 1, 1.1) #0.8 # proportion of self pollination
delta = [1] #np.arange(0,1.1,0.1) #  a fraction of the self fertilized offspring that won't survive (a fraction of 1 - delta will survive)
SC_type = "SI" # "SCf" # "SCp" #
#is_SCp = 0
#is_SCf = 1
SC_freq_thresh = 0.8
#interaction parameters
interaction_thresh = -6
num_of_fertilizations_trials = 5
#mutation probabilities
p_mutation = 1e-4 #1e-10 #0.0001#
p_duplicate = 1e-6 #0.01 #1e-6 #
p_delete = 1e-6
disable_RNase_dup = True # dont let RNases to be duplicated
enable_dup = str("disabled")
disable_RNase_del = True # dont let RNases to be deleted
enable_del = str("disabled")

write_to_txt_file = 0

'''
% H -        Hydrophbic
% P -        neural polar
% + (Plus)-  positively charged
% - (Minus)- negatively charged
% eij is organised as follows:

This classification is taken from a study of protein folding [REF].

how connections between amino acids in a peptide chain allow for different protein folds.

% thus the potential matrix is condensed into a 4 X 4 symetric matrix

%    H P + - S
%  H _________
%  P _________
%  + _________
%  - _________
%  S _________
'''
eij = [[-1, 0, 0, 0],
       [0, -0.75, -0.25, -0.25],
       [0, -0.25, 1, -1.25],
       [0, -0.25, -1.25, 1]]

eij = np.asarray(eij)
