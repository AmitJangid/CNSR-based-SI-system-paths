from typing import Dict, Any, Tuple, Union
import numpy as np
from numpy import random
# import numpy.matlib
# import pandas as pd
from copy import deepcopy
# from itertools import combinations
import cProfile, pstats #, StringIO #io
import simulationVars as Model
import allele_duplication_module as Duplicate
import collections
import analyze_mutations_supportive_functions as supportFuncs


def mutate_population(iter_ind, bioclass_freqs, haplotypes_dict, diploids_count_dict, haplotypes_count_dict,
          AA_rnases_dict, AA_slfs_dict):

#ver4 - 18/1/21 - stop codon is a legitimate mutants.
#                 alleles with stop codons are not added to the list of RNases and SLFs
#                 rather they are marked in diploids_obj and haplotypes_obj as "non-functional alleles"
# (replacesed by the earlier ver (ver3) - while the mutated codon is a STOP CODON - keep mutating
# ver9 - 28/2/22 -


    # input:
    #iter_ind                      - for forcing the pseudo randon data to be the same as before
    #diploids_obj              - list of flower objects (contains only existing diploids)
    #haplotypes_obj            - list of haplotype objects (contains only existing haplotypes)
    #diploids_count_dict       - dictionary of diploids amount (contains only existing diploids)
    #haplotypes_count_dict     - dictionary of haplotypes amount (contains only existing haplotypes)

    #all_RNase_ar              - array of all the RNases in the population
    #all_SLF_ar                - array of all the SLFs in the population

    # output :
    # diploids_mutated_obj        - list of flower objects after mutations (excluding empty diploids)
    # haplotypes_obj              - list of haplotypes objects after mutations (including empty haplotypes)
    # diploids_count_mutated_dict - dictionary of diploids amounts (excluding empty diploids)
    # haplotypes_count_dict       - dictionary of haplotypes amounts (including empty haplotypes)
    # all_RNase_ar                - array of unique sequences of female AA's labeled by their bioclasses
    # all_SLF_ar                  - array of unique sequences of male AA's labeled by their bioclasses

    # pr = cProfile.Profile()
    # pr.enable()

    num_mut_haplotypes = 0
    AA_RNases_count_dict, AA_SLFs_count_dict = \
        count_classified_AAs(haplotypes_count_dict, haplotypes_dict, AA_rnases_dict, AA_slfs_dict)

    mutated_bioClass0_indices_bool = np.random.choice([False, True],
                             AA_SLFs_count_dict[0] + AA_RNases_count_dict[0], replace=True, p = [1-Model.p_mutation, Model.p_mutation])
    mutated_bioClass1_indices_bool = np.random.choice([False, True],
                                 AA_SLFs_count_dict[1] + AA_RNases_count_dict[1], replace=True, p = [1-Model.p_mutation, Model.p_mutation])
    mutated_bioClass2_indices_bool = np.random.choice([False, True],
                                 AA_SLFs_count_dict[2] + AA_RNases_count_dict[2], replace=True, p = [1-Model.p_mutation, Model.p_mutation])
    mutated_bioClass3_indices_bool = np.random.choice([False, True],
                                 AA_SLFs_count_dict[3] + AA_RNases_count_dict[3], replace=True, p = [1-Model.p_mutation, Model.p_mutation])

    mutated_bioClass0_indices_tmp = np.random.choice([0, 1, 2, 3],
                             AA_SLFs_count_dict[0] + AA_RNases_count_dict[0], replace=True, p = bioclass_freqs)
    mutated_bioClass1_indices_tmp = np.random.choice([0, 1, 2, 3],
                             AA_SLFs_count_dict[1] + AA_RNases_count_dict[1], replace=True, p = bioclass_freqs)
    mutated_bioClass2_indices_tmp = np.random.choice([0, 1, 2, 3],
                             AA_SLFs_count_dict[2] + AA_RNases_count_dict[2], replace=True, p = bioclass_freqs)
    mutated_bioClass3_indices_tmp = np.random.choice([0, 1, 2, 3],
                             AA_SLFs_count_dict[3] + AA_RNases_count_dict[3], replace=True, p = bioclass_freqs)

    mutated_bioClass0_indices = np.zeros((AA_SLFs_count_dict[0] + AA_RNases_count_dict[0]))
    mutated_bioClass1_indices = np.ones((AA_SLFs_count_dict[1] + AA_RNases_count_dict[1]))
    mutated_bioClass2_indices = 2*np.ones((AA_SLFs_count_dict[2] + AA_RNases_count_dict[2]))
    mutated_bioClass3_indices = 3*np.ones((AA_SLFs_count_dict[3] + AA_RNases_count_dict[3]))

    if any(mutated_bioClass0_indices_bool):
        for ii, bool_val in enumerate(mutated_bioClass0_indices_bool):
            if bool_val:
                #num_muts_from_to[0,mutated_bioClass0_indices_tmp[ii]] += 1
                mutated_bioClass0_indices[ii] = mutated_bioClass0_indices_tmp[ii]
    if any(mutated_bioClass1_indices_bool):
        for ii, bool_val in enumerate(mutated_bioClass1_indices_bool):
            if bool_val:
                #num_muts_from_to[1, mutated_bioClass1_indices_tmp[ii]] += 1
                mutated_bioClass1_indices[ii] = mutated_bioClass1_indices_tmp[ii]
    if any(mutated_bioClass2_indices_bool):
        for ii, bool_val in enumerate(mutated_bioClass2_indices_bool):
            if bool_val:
                #num_muts_from_to[2, mutated_bioClass2_indices_tmp[ii]] += 1
                mutated_bioClass2_indices[ii] = mutated_bioClass2_indices_tmp[ii]
    if any(mutated_bioClass3_indices_bool):
        for ii, bool_val in enumerate(mutated_bioClass3_indices_bool):
            if bool_val:
                #num_muts_from_to[3, mutated_bioClass3_indices_tmp[ii]] += 1
                mutated_bioClass3_indices[ii] = mutated_bioClass3_indices_tmp[ii]
    b=6
    if all(mutated_bioClass0_indices == 0) and all(mutated_bioClass1_indices == 1) and all(mutated_bioClass2_indices == 2) and \
        all(mutated_bioClass3_indices == 3):
        diploids_count_mutated_dict = deepcopy(diploids_count_dict)
        haplotypes_count_mutated_dict = deepcopy(haplotypes_count_dict)

    else:
        mutated_bioClass_indices_dict = {}
        mutated_bioClass_indices_dict[0] = mutated_bioClass0_indices
        mutated_bioClass_indices_dict[1] = mutated_bioClass1_indices
        mutated_bioClass_indices_dict[2] = mutated_bioClass2_indices
        mutated_bioClass_indices_dict[3] = mutated_bioClass3_indices
        haplotypes_count_mutated_dict = deepcopy(haplotypes_count_dict)
        diploids_count_mutated_dict = deepcopy(diploids_count_dict)
        ind_bioclass = np.zeros((4), dtype=int)

        for keys, values in diploids_count_dict.items():  # one diploid index
            AA_SLFs1_bioClass = make_AAs_male_list(haplotypes_dict[keys[0]], AA_slfs_dict)
            AA_SLFs2_bioClass = make_AAs_male_list(haplotypes_dict[keys[1]], AA_slfs_dict)
            for value in range(values):  # number of such diploid in the population
                ###
                AA_SLFs1_mutated_bioClass, ind_bioclass = make_mut_AAs_male_arr\
                    (haplotypes_dict[keys[0]].SLFs_len, 1, AA_SLFs1_bioClass, mutated_bioClass_indices_dict, ind_bioclass)
                # try:
                is_mut_hap1_male = any(AA_SLFs1_bioClass != AA_SLFs1_mutated_bioClass)
                # except:
                #     b=6
                AA_RNases_bioClass = make_AAs_female_list(haplotypes_dict[keys[0]], AA_rnases_dict)
                AA_RNases_mutated_bioClass, ind_bioclass = make_mut_AAs_female_arr\
                    (haplotypes_dict[keys[0]].RNases_len, 1, AA_RNases_bioClass, mutated_bioClass_indices_dict, ind_bioclass)
                is_mut_hap1_female = any(AA_RNases_bioClass != AA_RNases_mutated_bioClass)
                if any([is_mut_hap1_female, is_mut_hap1_male]):
                    num_mut_haplotypes += 1
                    Haplotype1, haplotypes_dict, haplotypes_count_mutated_dict, AA_rnases_dict, AA_slfs_dict = \
                        mutate_haplotype_ver7(iter_ind, AA_SLFs1_bioClass, AA_RNases_bioClass,
                        AA_SLFs1_mutated_bioClass, AA_RNases_mutated_bioClass, haplotypes_dict[keys[0]],
                        haplotypes_count_mutated_dict, haplotypes_dict,
                                              AA_rnases_dict, AA_slfs_dict)
                    b=6
                else:
                    Haplotype1 = haplotypes_dict[keys[0]]
                ###
                AA_SLFs2_mutated_bioClass, ind_bioclass = make_mut_AAs_male_arr \
                    (haplotypes_dict[keys[1]].SLFs_len, 1, AA_SLFs2_bioClass, mutated_bioClass_indices_dict,
                     ind_bioclass)
                #try:
                is_mut_hap2_male = any(AA_SLFs2_bioClass != AA_SLFs2_mutated_bioClass)
                # except:
                #     b=6
                AA_RNases_bioClass = make_AAs_female_list(haplotypes_dict[keys[1]], AA_rnases_dict)
                AA_RNases_mutated_bioClass, ind_bioclass = make_mut_AAs_female_arr \
                    (haplotypes_dict[keys[1]].RNases_len, 1, AA_RNases_bioClass, mutated_bioClass_indices_dict,
                     ind_bioclass)
                is_mut_hap2_female = any(AA_RNases_bioClass != AA_RNases_mutated_bioClass)
                if any([is_mut_hap2_female, is_mut_hap2_male]):
                    num_mut_haplotypes += 1
                    Haplotype2, haplotypes_dict, haplotypes_count_mutated_dict, AA_rnases_dict, AA_slfs_dict = \
                        mutate_haplotype_ver7(iter_ind, AA_SLFs2_bioClass, AA_RNases_bioClass,
                                              AA_SLFs2_mutated_bioClass, AA_RNases_mutated_bioClass,
                                              haplotypes_dict[keys[1]],
                                              haplotypes_count_mutated_dict, haplotypes_dict,
                                              AA_rnases_dict, AA_slfs_dict)
                    b = 6
                else:
                    Haplotype2 = haplotypes_dict[keys[1]]
                ###
                if any([is_mut_hap1_female, is_mut_hap1_male, is_mut_hap2_female,
                        is_mut_hap2_male]):  # the flower was deleted
                    diploids_count_mutated_dict[keys] = diploids_count_mutated_dict[keys] - 1
                    diploid_tuple1 = tuple([Haplotype1.Haplotype_index_in_pool, Haplotype2.Haplotype_index_in_pool])
                    diploid_tuple2 = tuple([Haplotype2.Haplotype_index_in_pool, Haplotype1.Haplotype_index_in_pool])
                    if tuple(diploid_tuple1) in diploids_count_mutated_dict:
                        # already axist in the diploids pool
                        diploids_count_mutated_dict[diploid_tuple1] = diploids_count_mutated_dict[
                                                                          diploid_tuple1] + 1
                    elif tuple(diploid_tuple2) in diploids_count_mutated_dict:
                        # already axist in the diploids pool
                        diploids_count_mutated_dict[diploid_tuple2] = diploids_count_mutated_dict[
                                                                          diploid_tuple2] + 1
                    else:
                        # it is a new diploid
                        if diploid_tuple1[0] <= diploid_tuple1[1]:
                            diploids_count_mutated_dict[diploid_tuple1] = 1
                        else:
                            diploids_count_mutated_dict[diploid_tuple2] = 1

        for key in list(diploids_count_mutated_dict.keys()):  # one diploid index
            if diploids_count_mutated_dict[key] == 0:
                del diploids_count_mutated_dict[key]

    return haplotypes_dict, diploids_count_mutated_dict, haplotypes_count_mutated_dict, \
           AA_rnases_dict, AA_slfs_dict

def mutate_haplotype_ver7(iter_ind, AA_SLFs_bioClass, AA_RNases_bioClass,
                    AA_SLFs_mutated_bioClass, AA_RNases_mutated_bioClass, haplotype,
                    haplotypes_count_mutated_dict, haplotypes_dict,
                    AA_rnases_dict, AA_slfs_dict):
    # pr = cProfile.Profile()
    # pr.enable()

    # haplotypes_count_mutated_dict = deepcopy(haplotypes_count_dict)
    haplotype_orig_ind = haplotype.Haplotype_index_in_pool

    #-----------------------------------------------------
    #                 mutate the RNases
    #-----------------------------------------------------
    haplotype_RNases_mutated_list = []
    if any(AA_RNases_mutated_bioClass != AA_RNases_bioClass):  # mutate the RNases
        ind_start_allele = 0
        for ii_rnase, rnase_ind in enumerate(haplotype.RNasesIndices):  # (Model.rnase_amount_in_haploid):
            AA_RNase_mutated_bioClass = AA_RNases_mutated_bioClass[
                np.arange(ind_start_allele, ind_start_allele + Model.seq_len)]
            AA_RNase_bioClass = AA_RNases_bioClass[
                np.arange(ind_start_allele, ind_start_allele + Model.seq_len)]
            ind_start_allele = ind_start_allele + Model.seq_len
            if any(AA_RNase_mutated_bioClass != AA_RNase_bioClass):  # mutate this allele
                all_RNase_ind = mutate_allele_ver7(AA_RNase_mutated_bioClass, AA_rnases_dict)
                haplotype_RNases_mutated_list.append(all_RNase_ind)
                AA_rnases_dict[all_RNase_ind] = AA_RNase_mutated_bioClass
            else:  # this RNase was not mutated
                haplotype_RNases_mutated_list.append(rnase_ind)

    else:  # non of these rnases were mutated
        haplotype_RNases_mutated_list = haplotype.RNasesIndices

    # -----------------------------------------------------
    #                 mutate the SLFs
    # -----------------------------------------------------

    haplotype_SLFs_mutated_list = []
    if any(AA_SLFs_mutated_bioClass != AA_SLFs_bioClass):  # mutate the SLFs
        ind_start_allele = 0
        for ii_slf, slf_ind in enumerate(haplotype.SLFsIndices):  # (Model.slf_amount_in_haploid):
            AA_SLF_mutated_bioClass = AA_SLFs_mutated_bioClass[
                    np.arange(ind_start_allele, ind_start_allele + Model.seq_len)]
            AA_SLF_bioClass = AA_SLFs_bioClass[
                np.arange(ind_start_allele, ind_start_allele + Model.seq_len)]
            ind_start_allele = ind_start_allele + Model.seq_len
            if any(AA_SLF_mutated_bioClass != AA_SLF_bioClass):  # mutate this allele
                all_SLF_ind = mutate_allele_ver7(AA_SLF_mutated_bioClass, AA_slfs_dict)
                haplotype_SLFs_mutated_list.append(all_SLF_ind)
                AA_slfs_dict[all_SLF_ind] = AA_SLF_mutated_bioClass
            else: # this SLF was not mutated
                haplotype_SLFs_mutated_list.append(slf_ind)

    else:  # non of these slfs were mutated
        haplotype_SLFs_mutated_list = haplotype.SLFsIndices

    haplotype_RNases_tuple = tuple(haplotype_RNases_mutated_list)
    haplotype_SLFs_tuple = tuple(haplotype_SLFs_mutated_list)

    # ----------------------------------------------------------------------
    #       update haplotypes_obj, haplotypes_count_mutated_dict
    # ----------------------------------------------------------------------
    # if testing is required: small population, p_duplicate = 1,
    # the output after looping over the whole population, is the haplotypes originated in haplotypes_count_dict are empty
    # and are written in haplotypes_count_mutated_dict, a new list of haplotypes with the same amount is added at the end of this dict
    haplotype = deepcopy(haplotype)
    haplotype.RNasesIndices = haplotype_RNases_mutated_list
    haplotype.SLFsIndices = haplotype_SLFs_mutated_list
    haplotype.SLFsIndices_unique = np.unique(haplotype_SLFs_mutated_list)
    haplotype.SLFs_unique_len = len(np.unique(haplotype.SLFsIndices))
    haplotype.PreviousHaplotype = haplotype_orig_ind

    hap_ii = 0
    is_same_RNases = False  # collections.Counter(haplotype_RNases_duplicated_tuple) == collections.Counter(haplotypes_obj[haplotype_ind].RNasesIndices)
    is_same_SLFs = False  # collections.Counter(haplotype_SLFs_duplicated_tuple) == collections.Counter(haplotypes_obj[haplotype_ind].SLFsIndices)
    haplotypes_count_list = list(haplotypes_count_mutated_dict.keys())
    len_haplotypes_count_list = len(haplotypes_count_list)
    while (not is_same_RNases or not is_same_SLFs) and hap_ii < len_haplotypes_count_list:
        haplotype_ind = haplotypes_count_list[hap_ii]
        is_same_RNases = collections.Counter(haplotype_RNases_tuple) == collections.Counter(
            haplotypes_dict[haplotype_ind].RNasesIndices)
        is_same_SLFs = collections.Counter(haplotype_SLFs_tuple) == collections.Counter(
            haplotypes_dict[haplotype_ind].SLFsIndices)
        hap_ii = hap_ii + 1
    if is_same_RNases and is_same_SLFs:
        # already axist in the haplotypes pool
        hap_ii = hap_ii - 1
        haplotype_ind = haplotypes_count_list[hap_ii]
        haplotype_mutated_ind = haplotype_ind
        if haplotype_mutated_ind in haplotypes_count_mutated_dict:
            haplotypes_count_mutated_dict[haplotype_mutated_ind] = haplotypes_count_mutated_dict[
                                                                       haplotype_mutated_ind] + 1
        else:
            haplotypes_count_mutated_dict[haplotype_mutated_ind] = 1

        haplotypes_count_mutated_dict[haplotype_orig_ind] = haplotypes_count_mutated_dict[haplotype_orig_ind] - 1
        if haplotypes_count_mutated_dict[haplotype_orig_ind] == 0:
            del [haplotypes_count_mutated_dict[haplotype_orig_ind]]
        haplotype.Haplotype_index_in_pool = haplotype_mutated_ind

    else:  # it is a new haplotype
        Haplotype_index_in_pool = np.max(list(haplotypes_dict.keys())) + 1
        haplotypes_count_mutated_dict[Haplotype_index_in_pool] = 1
        haplotypes_count_mutated_dict[haplotype_orig_ind] = haplotypes_count_mutated_dict[haplotype_orig_ind] - 1
        if haplotypes_count_mutated_dict[haplotype_orig_ind] == 0:
            del [haplotypes_count_mutated_dict[haplotype_orig_ind]]
        haplotype.Haplotype_index_in_pool = Haplotype_index_in_pool
        # add amount the new haplotype (without the deleted alleles)
        # add birth iteration and reason
        haplotype.BirthIter = iter_ind
        haplotype.BirthReason = 'mutation'
        # copy the Haplotype object to haplotypes object
        haplotypes_dict[Haplotype_index_in_pool] = haplotype

    # pr.disable()
    # sortby = 'cumulative'
    #ps = pstats.Stats(pr).sort_stats(sortby)  # , stream=s
    #ps.print_stats()
    #print("##########################################################################")

    return haplotype, haplotypes_dict, haplotypes_count_mutated_dict, AA_rnases_dict, AA_slfs_dict


def mutate_allele_ver7(mutated_allele_ar, AA_dict):

    #all_ind            - the index of the mutated allele list in the alleles pool (all_ar)
    #                     if this mutated allele is not part of the allele pool the index equals n+1 (n number of alleles in the pool)

    not_in_all_ar = True
    for key, val in AA_dict.items():
        not_in_all_ar = any(val != mutated_allele_ar)
        if not not_in_all_ar:
            all_ind = key
            break
    if not_in_all_ar:
        all_ind = key+1

    return all_ind


def count_classified_AAs(haplotypes_count_dict, haplotypes_dict, AA_rnases_dict, AA_slfs_dict):
    RNases_count_dict = {}
    for keys, values in haplotypes_count_dict.items():
        if values > 0:
            RNases_hap = haplotypes_dict[keys].RNasesIndices
            for ii_rnase in RNases_hap:
                if ii_rnase not in RNases_count_dict:
                    RNases_count_dict[ii_rnase] = 0
                RNases_count_dict[ii_rnase] += values
    AA_RNases_count_dict = {}
    for keys, values in RNases_count_dict.items():
        aa_list = AA_rnases_dict[keys]
        for ii_aa in aa_list:
            if ii_aa not in AA_RNases_count_dict:
                AA_RNases_count_dict[ii_aa] = 0
            AA_RNases_count_dict[ii_aa] += values

    for ii_aa in range(5):
        if ii_aa not in AA_RNases_count_dict:
            AA_RNases_count_dict[ii_aa] = 0

    SLFs_count_dict = {}
    for keys, values in haplotypes_count_dict.items():
        if values > 0:
            SLFs_hap = haplotypes_dict[keys].SLFsIndices

            for ii_slf in SLFs_hap:
                if ii_slf not in SLFs_count_dict:
                    SLFs_count_dict[ii_slf] = 0
                SLFs_count_dict[ii_slf] += values

    AA_SLFs_count_dict = {}
    for keys, values in SLFs_count_dict.items():
        aa_list = AA_slfs_dict[keys]
        for ii_aa in aa_list:
            if ii_aa not in AA_SLFs_count_dict:
                AA_SLFs_count_dict[ii_aa] = 0
            AA_SLFs_count_dict[ii_aa] += values

    for ii_aa in range(5):
        if ii_aa not in AA_SLFs_count_dict:
            AA_SLFs_count_dict[ii_aa] = 0

    return AA_RNases_count_dict, AA_SLFs_count_dict

def make_AAs_male_list(haplotype, AA_slfs_dict):
    aa_ind = 0
    SLFs_indices = haplotype.SLFsIndices
    SLFs_len = haplotype.SLFs_len
    AA_SLFs_bioClass = np.zeros(SLFs_len * Model.seq_len, dtype=int)

    aa_ind = 0
    for slf_ind in SLFs_indices:
        aa_end = aa_ind + Model.seq_len
        AA_SLFs_bioClass[aa_ind:aa_end] = AA_slfs_dict[slf_ind]
        aa_ind = aa_end

    return AA_SLFs_bioClass

def make_mut_AAs_male_arr(SLFs_len, values, AA_SLFs_bioClass,
                  mutated_bioClass_indices_dict, ind_bioclass):

    AA_SLFs_mutated_bioClass = np.zeros(SLFs_len * Model.seq_len, dtype=int)
    for aa_ind in range(SLFs_len * Model.seq_len):
        AA_bioClass = AA_SLFs_bioClass[aa_ind]
        AA_SLFs_mutated_bioClass[aa_ind] = mutated_bioClass_indices_dict[AA_bioClass][
            ind_bioclass[AA_bioClass]]
        ind_bioclass[AA_bioClass] += 1

    return AA_SLFs_mutated_bioClass, ind_bioclass

def make_AAs_female_list(haplotype,AA_rnases_dict):

    RNases_indices = haplotype.RNasesIndices
    RNases_len = haplotype.RNases_len
    AA_RNases_bioClass = np.zeros(RNases_len * Model.seq_len, dtype=int)

    aa_ind = 0
    for rnase_ind in RNases_indices:
        aa_end = aa_ind + Model.seq_len
        AA_RNases_bioClass[aa_ind:aa_end] = AA_rnases_dict[rnase_ind]
        aa_ind = aa_end

    return AA_RNases_bioClass

def make_mut_AAs_female_arr(RNases_len, values, AA_RNases_bioClass,
                            mutated_bioClass_indices_dict, ind_bioclass):

    AA_RNases_mutated_bioClass = np.zeros(RNases_len * Model.seq_len, dtype=int)

    for aa_ind in range(RNases_len * Model.seq_len):
        AA_bioClass = AA_RNases_bioClass[aa_ind]
        AA_RNases_mutated_bioClass[aa_ind] = mutated_bioClass_indices_dict[AA_bioClass][
                        ind_bioclass[AA_bioClass]]
        ind_bioclass[AA_bioClass] += 1

    return AA_RNases_mutated_bioClass, ind_bioclass
