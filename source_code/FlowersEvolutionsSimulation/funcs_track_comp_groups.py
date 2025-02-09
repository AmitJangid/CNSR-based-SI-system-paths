import numpy as np
import simulationVars as Model
def parent_hap_func(haplotypes_dict, haplotypes_dict_previous, current_hap, parents_haps, iter_ind, d_iters, n):

    hap = current_hap
    try:
        birth_iter = haplotypes_dict[hap].BirthIter
    except:
        b=6
        parents_hap = np.nan
        return parents_hap

    if birth_iter <= iter_ind - n*d_iters:
        parents_hap = hap
    else:
        parent_in_previous_iter = True
        while parent_in_previous_iter:
            try:
                if parents_haps[hap] in haplotypes_dict:
                    birth_iter = haplotypes_dict[parents_haps[hap]].BirthIter
                    # if birth_iter <= iter_ind - n * d_iters:
                    parents_hap = parents_haps[hap]
                    parent_in_previous_iter = False
                    # else:
                    #     hap = parents_haps[hap]
                elif parents_haps[hap] in haplotypes_dict_previous:
                    birth_iter = haplotypes_dict_previous[parents_haps[hap]].BirthIter
                    if birth_iter <= iter_ind - n * d_iters:
                        parents_hap = parents_haps[hap]
                        parent_in_previous_iter = False
                    else:
                        hap = parents_haps[hap]
                else:
                    hap = parents_haps[hap]
            except:
                b=6
                parents_hap = np.nan
                parent_in_previous_iter = False

    return parents_hap

def parent_hap_n_iters_before(haplotypes_dict, current_hap, parents_haps):

    hap = current_hap
    if hap in haplotypes_dict:
        #current hap was born before more than 10 generations
        parents_hap = hap
        return parents_hap
    else:
        # current hap was born before less than 10 generations
        if parents_haps[hap] in haplotypes_dict:
            parents_hap = parents_haps[hap]
            return parents_hap
        else:
            parent_not_in_previous_iter = True
            parents_hap = parents_haps[hap]
            while parent_not_in_previous_iter:
                if parents_hap in haplotypes_dict:
                    parent_not_in_previous_iter = False
                    return parents_hap
                else:
                    try:
                        parents_haps[parents_hap]
                    except:
                        return np.nan
                    parents_hap = parents_haps[parents_hap]
                    parent_not_in_previous_iter = True

def parent_hap_current_iter(haplotypes_dict, current_hap, parents_haps):

    hap = current_hap
    if parents_haps[hap] in haplotypes_dict:
        parents_hap = parents_haps[hap]
        return parents_hap
    else:
        parent_not_in_previous_iter = True
        parents_hap = parents_haps[hap]
        while parent_not_in_previous_iter:
            if parents_hap in haplotypes_dict:
                parent_not_in_previous_iter = False
                return parents_hap
            else:
                parents_hap = parents_haps[parents_hap]
                parent_not_in_previous_iter = True

def parent_haps(haplotypes_dict, haplotypes_dict_previous, current_haps_lst, parents_haps, iter_ind, d_iters, n):
    parents_haps_arr = np.zeros(len(current_haps_lst), dtype = int)
    for curr_hap_ii, current_hap in enumerate(current_haps_lst):
        hap = current_hap
        try:
            birth_iter = haplotypes_dict[hap].BirthIter
        except:
            b=6
            parents_haps_arr[curr_hap_ii] = np.nan
            break
        if birth_iter <= iter_ind - n*d_iters:
            parents_haps_arr[curr_hap_ii] = hap
        else:
            parent_in_previous_iter = True
            while parent_in_previous_iter:
                if parents_haps[hap] in haplotypes_dict:
                    birth_iter = haplotypes_dict[parents_haps[hap]].BirthIter
                    if birth_iter <= iter_ind - d_iters:
                        parents_haps_arr[curr_hap_ii] = parents_haps[hap]
                        parent_in_previous_iter = False
                    else:
                        hap = parents_haps[hap]
                elif parents_haps[hap] in haplotypes_dict_previous:
                    birth_iter = haplotypes_dict_previous[parents_haps[hap]].BirthIter
                    if birth_iter <= iter_ind - d_iters:
                        parents_haps_arr[curr_hap_ii] = parents_haps[hap]
                        parent_in_previous_iter = False
                    else:
                        hap = parents_haps[hap]
                else:
                    hap = parents_haps[hap]
    return parents_haps_arr

def label_parent_hap(parent_hap, copy_mis_groups_hap_all_iters_dict, iter_ind, d_iters, n):
    label = np.nan
    for former_corrected_group_ind, former_haps_lst in \
            copy_mis_groups_hap_all_iters_dict[iter_ind - n*d_iters].items():
        if parent_hap in former_haps_lst:
            label = former_corrected_group_ind
            break
    return label
def parent_non_grouped_haps(haplotypes_dict, haplotypes_dict_previous, mis_non_grouped_haps_dict, parents_haps, iter_ind, d_iters, n):
    mis_non_grouped_parent_haps_dict = {}
    for current_hap_ind, group_lst in mis_non_grouped_haps_dict.items():
        hap = current_hap_ind
        birth_iter = haplotypes_dict[hap].BirthIter
        if birth_iter <= iter_ind - d_iters:
            parent_hap = hap
        else:
            parent_in_previous_iter = True
            while parent_in_previous_iter:
                if parents_haps[hap] in haplotypes_dict:
                    birth_iter = haplotypes_dict[parents_haps[hap]].BirthIter
                    if birth_iter <= iter_ind - d_iters:
                        parent_hap = parents_haps[hap]
                        parent_in_previous_iter = False
                    else:
                        hap = parents_haps[hap]
                elif parents_haps[hap] in haplotypes_dict_previous:
                    birth_iter = haplotypes_dict_previous[parents_haps[hap]].BirthIter
                    if birth_iter <= iter_ind - d_iters:
                        parent_hap = parents_haps[hap]
                        parent_in_previous_iter = False
                    else:
                        hap = parents_haps[hap]
                else:
                    hap = parents_haps[hap]
        mis_non_grouped_parent_haps_dict[parent_hap] = group_lst
    return mis_non_grouped_parent_haps_dict

def set_current_to_new_group_labeling(parents_haps_arr, copy_mis_groups_hap_all_iters_dict, iter_ind, d_iters, n):
    parents_haps_arr_unique = np.unique(parents_haps_arr)
    current_to_new_group_labeling = \
        np.zeros(len(parents_haps_arr_unique))
    current_to_new_group_labeling[:] = np.nan

    for par_hap_ii, parent_hap in enumerate(parents_haps_arr_unique):
        for former_corrected_group_ind, former_haps_lst in \
                copy_mis_groups_hap_all_iters_dict[iter_ind - n*d_iters].items():
            if parent_hap in former_haps_lst:
                current_to_new_group_labeling[par_hap_ii] \
                    = former_corrected_group_ind
                break
    return current_to_new_group_labeling


def check_comp_with_groups \
                (hap_lst, haplotypes_dict, AA_rnases_dict, AA_slfs_dict, E_Ri_Fj_dict, hap_partner):
    comp_hap_with_group_dict = {}

    partner_rnase_ind = haplotypes_dict[hap_partner].RNasesIndices[0]
    for hap_ii, hap in enumerate(hap_lst):
        comp_cond1 = False
        # the mutant haplotype is a male
        # the group haplotypes are the females
        rnase_ind = haplotypes_dict[hap].RNasesIndices[0]
        rnase_AA = AA_rnases_dict[rnase_ind]
        for slf_ind in haplotypes_dict[hap_partner].SLFsIndices:
            slf_AA = AA_slfs_dict[slf_ind]
            if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
            else:
                E = E_Ri_Fj_dict[rnase_ind, slf_ind]
            if E < Model.interaction_thresh:
                # interacting
                comp_cond1 = True
                break
        comp_cond2 = False
        # the mutant haplotype is a female
        # the group haplotypes are the males
        rnase_AA = AA_rnases_dict[partner_rnase_ind]
        for slf_ind in haplotypes_dict[hap].SLFsIndices:
            slf_AA = AA_slfs_dict[slf_ind]
            if tuple([partner_rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
            else:
                E = E_Ri_Fj_dict[partner_rnase_ind, slf_ind]
            if E < Model.interaction_thresh:
                # interacting
                comp_cond2 = True
                break
        if comp_cond1 and comp_cond2:
            # symetric compatibility
            comp_hap_with_group_dict[hap_ii] = [2]
        elif not comp_cond1 and not comp_cond2:
            comp_hap_with_group_dict[hap_ii] = [0]
        elif comp_cond1 and not comp_cond2: # cond1-True and cond2-False
            txt = 'interact: mutant male-group female, non interact: mutant female-group male'
            comp = [2, 0]
            comp_hap_with_group_dict[hap_ii] = [1, comp]
        elif comp_cond2 and not comp_cond1: # cond2-True and cond1-False
            txt = 'interact: mutant female-group male, non interact: mutant male-group female'
            comp = [0, 2]
            comp_hap_with_group_dict[hap_ii] = [1, comp]
        else:
            b=6


    # for group_ind in comp_haps_in_group_dict.keys():
    #     if all(comp_haps_in_group_dict[group_ind] == 2):
    #         comp_group[group_ind] = 2  # comp
    #
    #     elif all(comp_haps_in_group_dict[group_ind] == 0):
    #         comp_group[group_ind] = 0  # non_comp
    #     else:
    #         comp_group[group_ind] = 1

    return comp_hap_with_group_dict