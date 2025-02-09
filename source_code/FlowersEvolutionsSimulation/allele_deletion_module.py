from typing import Dict, Any, Tuple, Union
import numpy as np
from numpy import random
#import numpy.matlib
#import pandas as pd
from copy import deepcopy
#from itertools import combinations
#import cProfile, pstats #, StringIO #io
import simulationVars as Model
import collections
import allele_duplication_module as Duplicate

def delete_population(iter_ind, haplotypes_dict, diploids_count_dict, haplotypes_count_dict):
    #choose randomly alleles to be deleted

    # FIX total_num_of_alleles calculation for different amount of alleles !!!!
    #total_num_of_alleles = Model.popu_size * 2 * (1 + Model.slf_amount_in_haploid)
    total_num_of_alleles = 0
    for keys, values in haplotypes_count_dict.items():
        num_type_RNases = len(haplotypes_dict[keys].RNasesIndices)
        num_type_SLFs = len(haplotypes_dict[keys].SLFsIndices)
        num_Alleles_hap = num_type_SLFs
        total_num_of_alleles = total_num_of_alleles + values * num_Alleles_hap
    prob_binary_tobe_deleted = np.random.choice([0, 1], total_num_of_alleles, replace=True,
                                                p=[1 - Model.p_delete, Model.p_delete])
    if sum(prob_binary_tobe_deleted)>0:
        ind_end = 0
        diploids_count_deleted_dict = deepcopy(diploids_count_dict)
        for keys, values in diploids_count_dict.items():  # one diploid index
            for value in range(values):  # number of such diploid in the population
                num_type_RNases1 = haplotypes_dict[keys[0]].RNases_len
                num_type_SLFs1 = haplotypes_dict[keys[0]].SLFs_len
                num_type_RNases2 = haplotypes_dict[keys[1]].RNases_len
                num_type_SLFs2 = haplotypes_dict[keys[1]].SLFs_len
                num_Alleles_hap1 = num_type_SLFs1
                num_Alleles_hap2 = num_type_SLFs2
                ind_start = ind_end
                ind_end = ind_start + num_Alleles_hap1 + num_Alleles_hap2

                flower_prob_binary = prob_binary_tobe_deleted[np.arange(ind_start, ind_end)]
                if sum(flower_prob_binary)>0:
                    is_deleted_hap1_female = False
                    is_deleted_hap1_male = False
                    is_deleted_hap2_female = False
                    is_deleted_hap2_male = False
                    # Haplotype1 = flower.Haplotype1
                    # Haplotype2 = flower.Haplotype2
                    haplotype1_prob_binary = prob_binary_tobe_deleted[np.arange(ind_start, ind_start + num_Alleles_hap1)]
                    haplotype2_prob_binary = prob_binary_tobe_deleted[np.arange(ind_start + num_Alleles_hap1, ind_start + num_Alleles_hap1 + num_Alleles_hap2)]

                    if sum(haplotype1_prob_binary)>0:
                        is_deleted_hap1_female, is_deleted_hap1_male, Haplotype1, haplotypes_count_dict = \
                            delete_haplotype (iter_ind, prob_binary_tobe_deleted, ind_start,
                                              haplotypes_dict[keys[0]],
                                                 haplotypes_count_dict, haplotypes_dict)
                        # previousHap_ind = Haplotype1.PreviousHaplotype
                        # currentHap_ind = Haplotype1.Haplotype_index_in_pool
                        # haplotypes_dict[previousHap_ind].FollowingHaplotypesLst.append(currentHap_ind)
                        Haplotype_index_in_pool = np.max(list(haplotypes_dict.keys())) + 1

                        if is_deleted_hap1_female or is_deleted_hap1_male:
                            haplotypes_dict[Haplotype_index_in_pool] = Haplotype1
                        else:
                            b=6
                    else:
                        Haplotype1 = haplotypes_dict[keys[0]]
                    if sum(haplotype2_prob_binary)>0:
                        is_deleted_hap2_female, is_deleted_hap2_male, Haplotype2, haplotypes_count_dict = \
                            delete_haplotype(iter_ind, prob_binary_tobe_deleted, ind_start + num_Alleles_hap1,
                                             haplotypes_dict[keys[1]],
                                                haplotypes_count_dict, haplotypes_dict)
                        previousHap_ind = Haplotype2.PreviousHaplotype
                        # currentHap_ind = Haplotype2.Haplotype_index_in_pool
                        # haplotypes_dict[previousHap_ind].FollowingHaplotypesLst.append(currentHap_ind)
                        Haplotype_index_in_pool = np.max(list(haplotypes_dict.keys())) + 1
                        if is_deleted_hap2_female or is_deleted_hap2_male:
                            haplotypes_dict[Haplotype_index_in_pool] = Haplotype2
                        else:
                            b=6
                    else:
                        Haplotype2 = haplotypes_dict[keys[1]]

                    if any([is_deleted_hap1_female, is_deleted_hap1_male, is_deleted_hap2_female, is_deleted_hap2_male]):  # the flower was deleted
                        diploids_count_deleted_dict[keys] = diploids_count_deleted_dict[keys] - 1
                        diploid_tuple1 = tuple([Haplotype1.Haplotype_index_in_pool, Haplotype2.Haplotype_index_in_pool])
                        diploid_tuple2 = tuple([Haplotype2.Haplotype_index_in_pool, Haplotype1.Haplotype_index_in_pool])
                        if tuple(diploid_tuple1) in diploids_count_deleted_dict:
                            # already axist in the diploids pool
                            diploids_count_deleted_dict[diploid_tuple1] = diploids_count_deleted_dict[diploid_tuple1] + 1
                        elif tuple(diploid_tuple2) in diploids_count_deleted_dict:
                            # already axist in the diploids pool
                            diploids_count_deleted_dict[diploid_tuple2] = diploids_count_deleted_dict[diploid_tuple2] + 1
                        else:
                            # it is a new diploid
                            if diploid_tuple1[0] <= diploid_tuple1[1]:
                                diploids_count_deleted_dict[diploid_tuple1] = 1
                            else:
                                diploids_count_deleted_dict[diploid_tuple2] = 1

        for key in list(diploids_count_deleted_dict.keys()):  # one diploid index
            if diploids_count_deleted_dict[key] == 0:
                del diploids_count_deleted_dict[key]
    else:
        diploids_count_deleted_dict = diploids_count_dict
    return haplotypes_dict, diploids_count_deleted_dict, haplotypes_count_dict


def delete_haplotype (iter_ind, prob_binary_tobe_deleted, ind, haplotype, haplotypes_count_dict, haplotypes_dict):
    haplotypes_count_deleted_dict = deepcopy(haplotypes_count_dict)

    haplotype_orig_ind = haplotype.Haplotype_index_in_pool

    # -----------------------------------------------------
    #                delete the RNases
    # -----------------------------------------------------
    haplotype_RNases_list = []
    haplotype_Ancestral_RNases_list = []
    is_deleted_hap_female = False
    if not Model.disable_RNase_del:

        for ii_rnase, rnase_ind in enumerate(haplotype.RNasesIndices):
            allele_prob_binary = prob_binary_tobe_deleted[ind]
            ind = ind + 1
            if not allele_prob_binary:  # do not delete this SLF
                haplotype_RNases_list.append(rnase_ind)
                haplotype_Ancestral_RNases_list.append(haplotype.RNasesAncestralIndices[ii_rnase])
        if len(haplotype_RNases_list) != len(haplotype.RNasesIndices):
            is_deleted_hap_female = True
    else:
        haplotype_RNases_list = haplotype.RNasesIndices
        haplotype_Ancestral_RNases_list = haplotype.RNasesAncestralIndices
    # -----------------------------------------------------
    #                delete the SLFs
    # -----------------------------------------------------
    is_deleted_hap_male = False
    haplotype_SLFs_list = []
    haplotype_Ancestral_SLFs_list = []
    for ii_slf, slf_ind in enumerate (haplotype.SLFsIndices):
        allele_prob_binary = prob_binary_tobe_deleted[ind]
        ind = ind + 1
        if not allele_prob_binary:  # do not delete this SLF
            haplotype_SLFs_list.append(slf_ind)
            haplotype_Ancestral_SLFs_list.append(haplotype.SLFsAncestralIndices[ii_slf])

    if len(haplotype_SLFs_list) != haplotype.SLFs_len:
        is_deleted_hap_male = True

    haplotype_RNases_tuple = tuple(haplotype_RNases_list)
    haplotype_SLFs_tuple = tuple(haplotype_SLFs_list)

    #----------------------------------------------------------------------
    #       update haplotypes_dict, haplotypes_count_deleted_dict
    #----------------------------------------------------------------------

    # if testing is required: small population, p_duplicate = 1,
    # the output after looping over the whole population, is the haplotypes originated in haplotypes_count_dict are empty
    # and are written in haplotypes_count_duplicated_dict, a new list of haplotypes with the same amount is added at the end of this dict

    if is_deleted_hap_female or is_deleted_hap_male:
        haplotype = deepcopy(haplotype)
        haplotype.RNasesIndices = haplotype_RNases_list
        haplotype.SLFsIndices = haplotype_SLFs_list
        haplotype.SLFsIndices_unique = np.unique(haplotype_SLFs_list)
        haplotype.SLFs_unique_len = len(np.unique(haplotype_SLFs_list))
        haplotype.SLFs_len = len(haplotype_SLFs_list)
        haplotype.RNases_len = len(haplotype_RNases_list)
        haplotype.PreviousHaplotype = haplotype_orig_ind
        haplotype.SLFsAncestralIndices = haplotype_Ancestral_SLFs_list
        haplotype.RNasesAncestralIndices = haplotype_Ancestral_RNases_list


        hap_ii = 0
        is_same_RNases = False
        is_same_SLFs = False
        haplotypes_count_list = list(haplotypes_count_dict.keys())
        while (not is_same_RNases or not is_same_SLFs) and hap_ii < len(haplotypes_count_list):
            haplotype_ind = haplotypes_count_list[hap_ii]
            is_same_RNases = collections.Counter(haplotype_RNases_tuple) == collections.Counter(haplotypes_dict[haplotype_ind].RNasesIndices)
            is_same_SLFs = collections.Counter(haplotype_SLFs_tuple) == collections.Counter(haplotypes_dict[haplotype_ind].SLFsIndices)
            hap_ii = hap_ii + 1
        if is_same_RNases and is_same_SLFs:# haplotype_ind < len(haplotypes_count_duplicated_dict):
            # already axist in the haplotypes pool
            hap_ii = hap_ii - 1
            haplotype_ind = haplotypes_count_list[hap_ii]
            haplotype_deleted_ind = haplotype_ind
            if haplotype_deleted_ind in haplotypes_count_deleted_dict:
                haplotypes_count_deleted_dict[haplotype_deleted_ind] = haplotypes_count_deleted_dict[
                                                                           haplotype_deleted_ind] + 1
            else:
                haplotypes_count_deleted_dict[haplotype_deleted_ind] = 1

            haplotypes_count_deleted_dict[haplotype_orig_ind] = haplotypes_count_deleted_dict[haplotype_orig_ind] - 1
            if haplotypes_count_deleted_dict[haplotype_orig_ind] == 0:
                del[haplotypes_count_deleted_dict[haplotype_orig_ind]]
            haplotype.Haplotype_index_in_pool = haplotype_deleted_ind
            # haplotype.Haplotype_count = haplotypes_count_deleted_dict[haplotype_deleted_ind]
        else:  # it is a new haplotype
            Haplotype_index_in_pool = np.max(list(haplotypes_dict.keys())) + 1
            haplotypes_count_deleted_dict[Haplotype_index_in_pool] = 1
            haplotypes_count_deleted_dict[haplotype_orig_ind] = haplotypes_count_deleted_dict[haplotype_orig_ind] - 1
            if haplotypes_count_deleted_dict[haplotype_orig_ind] == 0:
                del[haplotypes_count_deleted_dict[haplotype_orig_ind]]
            # update the index of the deleted haplotype
            haplotype.Haplotype_index_in_pool = Haplotype_index_in_pool
            # add amount the new haplotype (without the deleted alleles)
            # haplotype.Haplotype_count = 1
            # add birth iteration and reason
            haplotype.BirthIter = iter_ind
            haplotype.BirthReason = 'Deletion'
            # copy the Haplotype object to haplotypes object

    return is_deleted_hap_female, is_deleted_hap_male, haplotype, haplotypes_count_deleted_dict