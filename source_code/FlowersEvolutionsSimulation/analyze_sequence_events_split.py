import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import simulationVars as Model
import funcs_track_comp_groups as funcs
import copy
min_duration = 10
d_iters = 10
thresh_num_hap_in_group = 10
os.chdir("..")

interaction_thresh = -6
os.chdir("FlowersEvolutionsSimulation_ver13_data")
# directory = "output_data_10haplotypes_k5_18AAs_E_thresh" + str(interaction_thresh) + "_p10-4"
# # with open(directory + '/diversity_mut_dict_' + str(5 * 10) + '_iters.pkl', 'rb') as input:
# #     diversity_mut_dict = pickle.load(input)
# # subdirectories = list(diversity_mut_dict.keys())
# subdirectories = list(set(next(os.walk(directory))[1]) - set(['2023-04-19-15-41-27-33']) - set(['2023-04-19-15-41-27-36'])                      -set(['2023-03-19-11-56-28-30'])
#     - set(['2023-03-19-11-56-28-25']))
# # subdirectories = ['2022-11-12-15-57-45-11']

# os.chdir("FlowersEvolutionsSimulation_ver13_data3") # ... data2 -> n_haps = 250, data3 -> n_haps = 750
# directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_alpha_08_delta_09"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_alpha_095_delta_085"

subdirectories = list(set(next(os.walk(directory))[1]) - set(['2023-06-07-10-54-24-41']) \
      -set(['2023-06-07-10-54-24-13']) - set(['2023-05-28-06-52-45-50']))
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p2times10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_phalf_times10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_no_dup_no_del"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-4_p10-4_rerun"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-4_p10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_2000haps"
subdirectories = list(set(next(os.walk(directory))[1])
              - set(['2023-06-28-10-19-50-3']) - set(['2023-06-28-10-19-51-22']))

num_files = len(subdirectories)
analyze_func = False

diversity_mut_dict = {}
diversity_orig_dict = {}
if analyze_func:
    frac_func_orig_slf_from_orig_group = np.zeros(len(subdirectories))
    frac_func_orig_slf_from_orig_group[:] = np.nan
    frac_func_orig_slf_from_mut_group = np.zeros(len(subdirectories))
    frac_func_orig_slf_from_mut_group[:] = np.nan

slf_orig_rnase_slf_mut_split_SC_slf_orig = 0
slf_orig_rnase_slf_mut_split_SC_slf_mut = 0
slf_orig_slf_mut_rnase_split_SC_slf_orig = 0
slf_orig_slf_mut_rnase_split_SC_slf_mut = 0
slf_mut_slf_orig_rnase_split_SC_slf_mut = 0
slf_mut_slf_orig_rnase_split_SC_slf_orig = 0
slf_mut_rnase_slf_orig_split_SC_slf_orig = 0
slf_mut_rnase_slf_orig_split_SC_slf_mut = 0
rnase_slf_mut_slf_orig_split_SC_rnase = 0
rnase_slf_mut_slf_orig_split_SC_slf_orig = 0
rnase_slf_mut_slf_orig_split_SC_slf_mut = 0
rnase_slf_orig_slf_mut_split_SC_rnase = 0
rnase_slf_orig_slf_mut_split_SC_slf_orig = 0
rnase_slf_orig_slf_mut_split_SC_slf_mut = 0

rnase_slf_mut_slf_orig_split = 0
rnase_slf_orig_slf_mut_split = 0
slf_mut_rnase_slf_orig_split = 0
slf_orig_rnase_slf_mut_split = 0
slf_orig_slf_mut_rnase_split = 0
slf_mut_slf_orig_rnase_split = 0

slf_mut_rnase_slf_orig_split_NSC_more_than_one_rnase_in_orig_gro = 0
slf_mut_slf_orig_rnase_split_NSC_more_than_one_rnase_in_orig_gro = 0
slf_orig_slf_mut_rnase_split_NSC_more_than_one_rnase_in_orig_gro = 0
for ii_dir, subdirectory in enumerate(subdirectories):
    output_dir = directory + "/" + subdirectory
    if not os.path.isfile(output_dir + '/analyze_split_group_all_iters_dict_' + str(min_duration * 10) + '_iters.pkl'):
        if os.path.isfile(output_dir + '/mut_rnases_' + str(min_duration * 10) + '_iters_dict.pkl') and \
            os.path.isfile(output_dir + '/split_one_group_into_two_all_iters_' + str(min_duration*10) + '_iters_dict.pkl'):
            print(os.path.join(subdirectory))
            with open(output_dir + '/AA_rnases_dict.pkl', 'rb') as input:
                AA_rnases_dict = pickle.load(input)
            with open(output_dir + '/AA_slfs_dict.pkl', 'rb') as input:
                AA_slfs_dict = pickle.load(input)
            with open(output_dir + '/E_Ri_Fj_dict.pkl', 'rb') as input:
                E_Ri_Fj_dict = pickle.load(input)
            with open(output_dir + '/split_one_group_into_two_all_iters_' + str(min_duration*10) + '_iters_dict.pkl', 'rb') as input:
                split_one_group_into_two_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_strict_groups_hap_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_strict_groups_hap_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_strict_groups_rnase_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_strict_groups_rnase_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_hap_strict_groups_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_hap_strict_groups_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_rnase_strict_groups_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_rnase_strict_groups_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_hap_strict_group_extension_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_hap_strict_group_extension_all_iters_dict = pickle.load(input)
            with open(output_dir + '/reorder_mis_strict_group_extension_hap_all_iters_dict.pkl', 'rb') as input:
                reorder_mis_strict_group_extension_hap_all_iters_dict = pickle.load(input)
            with open(output_dir + '/mut_rnases_' + str(min_duration*10) + '_iters_dict.pkl', 'rb') as input:
                mut_rnases_dict = pickle.load(input)
            with open(output_dir + '/parents_haps.pkl', 'rb') as input:
                parents_haps = pickle.load(input)

            mut_iters_lst = list(mut_rnases_dict.keys())
            mut_iters_lst.sort(reverse=True)
            mut_iters_arr = np.array(mut_iters_lst)
            count_mut_in_orig_in_cause_split = 0
            diversity_mut = np.zeros(len(split_one_group_into_two_all_iters_dict))
            diversity_mut[:] = np.nan
            diversity_orig = np.zeros(len(split_one_group_into_two_all_iters_dict))
            diversity_orig[:] = np.nan
            analyze_split_group_all_iters_dict = {}
            analyze_split_rnase_cros_slfs_info = {}
            if analyze_func:
                num_func_orig_slf_from_mut_group = 0
                num_func_orig_slf_from_orig_group = 0
                num_non_func_orig_slf_from_mut_group = 0
                num_non_func_orig_slf_from_orig_group = 0

            for ii_iter, [iter_ind, tmp_dict] in enumerate(split_one_group_into_two_all_iters_dict.items()):
                mut_iters_arr_before_split = mut_iters_arr[mut_iters_arr<= iter_ind]
                if len(mut_iters_arr_before_split) == 0:
                    continue
                gro = list(tmp_dict[0].keys())[0]
                gro1 = list(split_one_group_into_two_all_iters_dict[iter_ind][0][gro][0].keys())[0]
                gro2 = list(split_one_group_into_two_all_iters_dict[iter_ind][0][gro][1].keys())[0]
                rnases1_lst = reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind][gro1]
                rnases2_lst = reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind][gro2]
                if len(list(set(rnases1_lst).intersection(set(rnases2_lst)))) > 0:
                    b=6
                rnases_in_gro_after_split = np.unique([item for sublist in [rnases1_lst,rnases2_lst] for item in sublist])
                ii = 0
                iter_first = 0
                sec_iter_ind = mut_iters_arr_before_split[ii]
                x = range(len(mut_rnases_dict[mut_iters_arr_before_split[0]])) # num of rnase mutations occoured at iteration 0
                for ii_mut_rnase in x:
                    rnases_in_gro_after_mut = [mut_rnases_dict[mut_iters_arr_before_split[0]][ii_mut_rnase]['orig RNase'], \
                        mut_rnases_dict[mut_iters_arr_before_split[0]][ii_mut_rnase]['mutant RNase']]
                    # same_group = len(list(set(rnases_in_gro_after_split) - set(rnases_in_gro_after_mut))) \
                    #              == len(rnases_in_gro_after_split) - 2
                    same_group = len(list(set(np.unique(rnases1_lst)).intersection(set(rnases_in_gro_after_mut)))) > 0 \
                        and len(list(set(np.unique(rnases2_lst)).intersection(set(rnases_in_gro_after_mut)))) > 0


                    cond1 = mut_iters_arr_before_split[0] <= iter_ind #true if the rnase mutation took place befor the split
                    cond2 = same_group #true if the splited group and the mutated group with the neutral rnase mutation are the same groups
                    cond = cond1 and cond2 #the two conditions fulfil
                    if cond: #the same mut rnase and its orig are those in the two splitted groups
                        count_mut_in_orig_in_cause_split += 1
                        continue
                if len(mut_iters_arr_before_split) == 1:
                    continue
                cond = True
                while cond:
                    ii += 1
                    sec_iter_ind = mut_iters_arr_before_split[ii]
                    x = range(len(mut_rnases_dict[sec_iter_ind]))
                    for ii_mut_rnase in x:
                        # for cases of more than one RNase mutant in a certain iteration
                        rnases_in_gro_after_mut = [mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['orig RNase'], \
                                                   mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['mutant RNase']]
                        same_group = len(list(set(rnases_in_gro_after_split) - set(rnases_in_gro_after_mut))) \
                                     == len(rnases_in_gro_after_split) - 2

                        cond1 = not sec_iter_ind <= iter_ind #true if not the rnase mutation took place befor the split
                        cond2 = not same_group #True only if the original RNase and its mutant RNase
                        # in the two splitted groups
                        cond3 = ii < len(mut_iters_arr_before_split) - 1 # true if before going through all the neutral rnase mutants
                        cond = (cond1 or cond2) and cond3
                        b = 6
                        if not cond:
                            break
                if (not cond1 and not cond2):
                    count_mut_in_orig_in_cause_split += 1
                    with open(output_dir + '/haplotypes_dict_' + str(sec_iter_ind) + '.pkl', 'rb') as input:
                        haplotypes_dict = pickle.load(input)
                    # mut_hap = list(appendix_mut_in_orig_in_all_iters_dict[sec_iter_ind].keys())[0]
                    mut_hap = mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['mutant hap']
                    mut_rnase = mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['mutant RNase']
                    orig_rnase = mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['orig RNase']
                    mut_rnase_iter = mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['mutant hap birth']
                    if orig_rnase in rnases1_lst:
                        diversity_orig[ii_iter] = len(np.unique(rnases1_lst))
                        diversity_mut[ii_iter] = len(np.unique(rnases2_lst))
                        orig_gro = gro1
                        mut_gro = gro2
                        orig_gro_ii = 0
                        mut_gro_ii = 1
                        b=6
                    else:
                        diversity_orig[ii_iter] = len(np.unique(rnases2_lst))
                        diversity_mut[ii_iter] = len(np.unique(rnases1_lst))
                        orig_gro = gro2
                        mut_gro = gro1
                        orig_gro_ii = 1
                        mut_gro_ii = 0
                        b=6

                    if orig_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind]:
                        with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_dict_ = pickle.load(input)
                        haps_orig_gro = [reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][orig_gro], \
                                         reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][orig_gro]]
                        haps_orig_gro = [item for sublist in haps_orig_gro for item in sublist]
                        # rnases_orig_gro = reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind][orig_gro]
                        # for hap in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][orig_gro]:
                        #     rnases_orig_gro.append(haplotypes_dict_[hap].RNasesIndices)
                        # rnases_orig_gro = [item for sublist in rnases_orig_gro for item in sublist]
                        # rnases_orig_gro = np.unique(rnases_orig_gro)
                    else:
                        haps_orig_gro = reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][orig_gro]
                    rnases_orig_gro = np.unique(reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind][orig_gro])

                    if mut_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind]:
                        with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_dict_ = pickle.load(input)
                        haps_mut_gro = [reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][mut_gro], \
                                         reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][mut_gro]]
                        haps_mut_gro = [item for sublist in haps_mut_gro for item in sublist]
                        # rnases_mut_gro = reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind][mut_gro]
                        # for hap in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][mut_gro]:
                        #     rnases_mut_gro.append(haplotypes_dict_[hap].RNasesIndices)
                        # rnases_mut_gro = [item for sublist in rnases_mut_gro for item in sublist]
                        # rnases_mut_gro = np.unique(rnases_mut_gro)
                    else:
                        haps_mut_gro = reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][mut_gro]
                    rnases_mut_gro = np.unique(reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind][mut_gro])

                    hap_in_orig_pollinate_rnase_in_mut = {}
                    for hap_ii, hap_ind in enumerate(haps_orig_gro):
                        with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_dict = pickle.load(input)
                        with open(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_count_dict = pickle.load(input)
                        interact_hap_count = haplotypes_count_dict[hap_ind]
                        interact_hap_birth = haplotypes_dict[hap_ind].BirthIter
                        slfs_ind = haplotypes_dict[hap_ind].SLFsIndices

                        for slf_ii, slf_ind in enumerate(slfs_ind):
                            slf_ii_correct = slf_ii
                            with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                                haplotypes_dict = pickle.load(input)
                            interact_hap = hap_ind

                            slf_AA = AA_slfs_dict[slf_ind]
                            E = np.zeros(len(rnases_mut_gro))
                            for rnase_ii, rnase_ind in enumerate(rnases_mut_gro):
                                rnase_AA = AA_rnases_dict[rnase_ind]
                                if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                    E[rnase_ii] = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                else:
                                    E[rnase_ii] = E_Ri_Fj_dict[rnase_ind, slf_ind]
                            if all(E < interaction_thresh):
                                # interacting SLF
                                n = 1
                                interact_slf = slf_ind
                                slf_already_exist = True
                                while slf_already_exist:
                                    try:
                                        with open(output_dir + '/haplotypes_dict_' + str(iter_ind - n * d_iters) + '.pkl',
                                            'rb') as input:
                                            haplotypes_dict_previous = pickle.load(input)
                                    except:
                                        slf_already_exist = False
                                    parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous,
                                       interact_hap, parents_haps)

                                    if parent_hap == interact_hap:
                                        orig_slfs_lst = haplotypes_dict[parent_hap].SLFsIndices
                                        orig_slf_ind = orig_slfs_lst[slf_ii_correct]
                                        interact_slf = orig_slf_ind
                                        haplotypes_dict = copy.deepcopy(haplotypes_dict_previous)
                                        n += 1
                                    else:
                                        orig_slfs_lst = haplotypes_dict_previous[parent_hap].SLFsIndices
                                        if len(haplotypes_dict[interact_hap].SLFsIndices) != \
                                            len(haplotypes_dict_previous[parent_hap].SLFsIndices):
                                            tmp_slf_ii = 0
                                            same_slfs = True
                                            while same_slfs and tmp_slf_ii <= slf_ii_correct:
                                                if haplotypes_dict[interact_hap].SLFsIndices[tmp_slf_ii] != \
                                                    haplotypes_dict_previous[parent_hap].SLFsIndices[tmp_slf_ii]:
                                                    same_slfs = False
                                                else:
                                                    tmp_slf_ii +=1
                                            if not same_slfs:
                                                b=6
                                                if len(haplotypes_dict[interact_hap].SLFsIndices)> \
                                                    len(haplotypes_dict_previous[parent_hap].SLFsIndices):
                                                    if len(list(set(np.unique(haplotypes_dict[interact_hap].SLFsIndices)) \
                                                        - set(np.unique(haplotypes_dict_previous[parent_hap].SLFsIndices)))) == 0:
                                                        slf_ii_correct -=1
                                                    else:
                                                        b=6
                                                if len(haplotypes_dict[interact_hap].SLFsIndices)< \
                                                    len(haplotypes_dict_previous[parent_hap].SLFsIndices):
                                                    if len(list(set(np.unique(haplotypes_dict[interact_hap].SLFsIndices)) \
                                                        - set(np.unique(haplotypes_dict_previous[parent_hap].SLFsIndices)))) == 0:
                                                        slf_ii_correct +=1
                                                    else:
                                                        b=6
                                        # what if SLFs where deleted or duplicated in addition to the mutation
                                        # check if the number of total SLFs remain the same
                                        # although it can still be equal number of deletions and duplications
                                        # but this is a rare case -so it is not taken into consideration
                                        orig_slf_ind = orig_slfs_lst[slf_ii_correct]
                                        slf_AA = AA_slfs_dict[orig_slf_ind]
                                        E = np.zeros(len(rnases_mut_gro))
                                        for rnase_ii, rnase_ind in enumerate(rnases_mut_gro):
                                            rnase_AA = AA_rnases_dict[rnase_ind]
                                            if tuple([rnase_ind, orig_slf_ind]) not in E_Ri_Fj_dict:
                                                E[rnase_ii] = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                            else:
                                                E[rnase_ii] = E_Ri_Fj_dict[rnase_ind, orig_slf_ind]
                                        if all(E < interaction_thresh):
                                            # interacting- the mutation in SLF was nonfunctional
                                            interact_hap = parent_hap
                                            interact_hap_birth = haplotypes_dict_previous[parent_hap].BirthIter
                                            interact_slf = orig_slf_ind
                                            haplotypes_dict = copy.deepcopy(haplotypes_dict_previous)
                                            n += 1
                                        else:
                                            if parents_haps[hap_ind] in haplotypes_dict and not parent_hap == parents_haps[hap_ind]:
                                                slf_exists_in_parent_hap_alive = True
                                                tmp_hap = hap_ind
                                                while slf_exists_in_parent_hap_alive:
                                                    tmp_hap = parents_haps[tmp_hap]
                                                    orig_slfs_lst = haplotypes_dict[tmp_hap].SLFsIndices
                                                    # if len(haplotypes_dict[tmp_hap].SLFsIndices) != \
                                                    #         len(haplotypes_dict[hap_ind].SLFsIndices):
                                                    #     b = 6
                                                    # what if SLFs where deleted or duplicated in addition to the mutation
                                                    # check if the number of total SLFs remain the same
                                                    # although it can still be equal number of deletions and duplications
                                                    # but this is a rare case -so it is not taken into consideration
                                                    orig_slf_ind = orig_slfs_lst[slf_ii_correct]
                                                    slf_AA = AA_slfs_dict[orig_slf_ind]
                                                    E = np.zeros(len(rnases_mut_gro))
                                                    for rnase_ii, rnase_ind in enumerate(rnases_mut_gro):
                                                        rnase_AA = AA_rnases_dict[rnase_ind]
                                                        if tuple([rnase_ind, orig_slf_ind]) not in E_Ri_Fj_dict:
                                                            E[rnase_ii] = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                                        else:
                                                            E[rnase_ii] = E_Ri_Fj_dict[rnase_ind, orig_slf_ind]
                                                    if all(E < interaction_thresh):
                                                        # interacting- the mutation in SLF didnt chane the functionality
                                                        interact_hap = parent_hap
                                                        interact_hap_birth = haplotypes_dict[tmp_hap].BirthIter
                                                        interact_slf = orig_slf_ind
                                                        if parents_haps[tmp_hap] not in haplotypes_dict:
                                                            slf_exists_in_parent_hap_alive = False
                                                            slf_already_exist = False
                                                    else:
                                                        slf_exists_in_parent_hap_alive = False
                                                        slf_already_exist = False
                                                b=6
                                            else:
                                                slf_already_exist = False
                                # end of while
                                if hap_ind not in hap_in_orig_pollinate_rnase_in_mut:
                                    hap_in_orig_pollinate_rnase_in_mut[hap_ind] = []
                                if analyze_func:
                                    # check if the SLF before mutation was functional
                                    # (decompose any of the rest of the RNases in the population))
                                    # at least one RNase

                                    slf_before_mut = haplotypes_dict_previous[parent_hap].SLFsIndices[slf_ii_correct]
                                    RNases_before_mut = []
                                    for h_keys, h_vals in haplotypes_dict_previous.items():
                                        RNases_before_mut.append(h_vals.RNasesIndices[0])
                                    RNases_before_mut = np.unique(RNases_before_mut)
                                    slf_AA_ = AA_slfs_dict[slf_before_mut]
                                    func_slf_cond = False
                                    for r_ind in RNases_before_mut:
                                        rnase_AA_ = AA_rnases_dict[r_ind]
                                        if tuple([r_ind, slf_before_mut]) not in E_Ri_Fj_dict:
                                            E_ = np.sum(Model.eij[rnase_AA_, slf_AA_], axis=0)
                                        else:
                                            E_ = E_Ri_Fj_dict[r_ind, slf_before_mut]
                                        if E_ < interaction_thresh:
                                            # interacting- the mutation in SLF was nonfunctional
                                            b = 6
                                            num_func_orig_slf_from_orig_group += 1
                                            func_slf_cond = True
                                            break
                                    if not func_slf_cond:
                                        num_non_func_orig_slf_from_orig_group += 1
                                    b = 6

                                tmp_dict = {}
                                if analyze_func:
                                    tmp_dict['is orig SLF func'] = E_ < interaction_thresh # true if func, false if not func
                                tmp_dict['interacting SLF mutant'] = interact_slf
                                tmp_dict['interacting hap birth'] = interact_hap_birth
                                tmp_dict['interacting hap count'] = interact_hap_count
                                hap_in_orig_pollinate_rnase_in_mut[hap_ind].append(tmp_dict)

                            else:  # non-interacting SLF
                                b = 6
                        b = 6  # for slfs in hap
                        b = 6  # for rnases in the orig group

                    if len(hap_in_orig_pollinate_rnase_in_mut) == 0:
                        continue
                    hap_in_mut_pollinate_rnase_in_orig = {}
                    for hap_ii, hap_ind in enumerate(haps_mut_gro):
                        interact_hap = hap_ind
                        with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_dict = pickle.load(input)
                        with open(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_count_dict = pickle.load(input)
                        interact_hap_count = haplotypes_count_dict[hap_ind]
                        interact_hap_birth = haplotypes_dict[hap_ind].BirthIter
                        slfs_ind = haplotypes_dict[hap_ind].SLFsIndices
                        for slf_ii, slf_ind in enumerate(slfs_ind):
                            with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                                haplotypes_dict = pickle.load(input)
                            interact_hap = hap_ind
                            slf_AA = AA_slfs_dict[slf_ind]
                            E = np.zeros(len(rnases_orig_gro))
                            for rnase_ii, rnase_ind in enumerate(rnases_orig_gro):
                                rnase_AA = AA_rnases_dict[rnase_ind]
                                if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                    E[rnase_ii] = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                else:
                                    E[rnase_ii] = E_Ri_Fj_dict[rnase_ind, slf_ind]
                            if all(E < interaction_thresh):
                                # interacting SLF
                                n=1
                                interact_slf = slf_ind
                                slf_already_exist = True
                                while slf_already_exist:
                                    try:
                                        with open(output_dir + '/haplotypes_dict_' + str(iter_ind - n * d_iters) + '.pkl',
                                            'rb') as input:
                                            haplotypes_dict_previous = pickle.load(input)
                                    except:
                                        slf_already_exist = False

                                    parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous,
                                                                   interact_hap, parents_haps)

                                    if parent_hap == interact_hap:
                                        orig_slfs_lst = haplotypes_dict[parent_hap].SLFsIndices
                                        orig_slf_ind = orig_slfs_lst[slf_ii]
                                        interact_slf = orig_slf_ind
                                        haplotypes_dict = copy.deepcopy(haplotypes_dict_previous)
                                        n += 1

                                    else:
                                        orig_slfs_lst = haplotypes_dict_previous[parent_hap].SLFsIndices
                                        if len(haplotypes_dict[interact_hap].SLFsIndices) != \
                                            len(haplotypes_dict_previous[parent_hap].SLFsIndices):
                                            tmp_slf_ii = 0
                                            same_slfs = True
                                            while same_slfs and tmp_slf_ii <= slf_ii:
                                                if haplotypes_dict[interact_hap].SLFsIndices[tmp_slf_ii] != \
                                                    haplotypes_dict_previous[parent_hap].SLFsIndices[tmp_slf_ii]:
                                                    same_slfs = False
                                                else:
                                                    tmp_slf_ii +=1
                                            if not same_slfs:
                                                b=6

                                        # what if SLFs where deleted or duplicated in addition to the mutation
                                        # check if the number of total SLFs remain the same
                                        # although it can still be equal number of deletions and duplications
                                        # but this is a rare case -so it is not taken into consideration

                                        orig_slfs_lst = haplotypes_dict_previous[parent_hap].SLFsIndices
                                        orig_slf_ind = orig_slfs_lst[slf_ii]
                                        slf_AA = AA_slfs_dict[orig_slf_ind]
                                        E = np.zeros(len(rnases_orig_gro))
                                        for rnase_ii, rnase_ind in enumerate(rnases_orig_gro):
                                            rnase_AA = AA_rnases_dict[rnase_ind]
                                            if tuple([rnase_ind, orig_slf_ind]) not in E_Ri_Fj_dict:
                                                E[rnase_ii] = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                            else:
                                                E[rnase_ii] = E_Ri_Fj_dict[rnase_ind, orig_slf_ind]
                                        if all(E < interaction_thresh):
                                            #interacting- the mutation in SLF was nonfunctional
                                            interact_hap = parent_hap
                                            interact_hap_birth = haplotypes_dict_previous[parent_hap].BirthIter
                                            interact_slf = orig_slf_ind
                                            haplotypes_dict = copy.deepcopy(haplotypes_dict_previous)
                                            n += 1
                                        else:
                                            b=6
                                            if parents_haps[hap_ind] in haplotypes_dict and not parent_hap == parents_haps[hap_ind]:
                                                slf_exists_in_parent_hap_alive = True
                                                tmp_hap = hap_ind
                                                while slf_exists_in_parent_hap_alive:
                                                    tmp_hap = parents_haps[tmp_hap]
                                                    orig_slfs_lst = haplotypes_dict[tmp_hap].SLFsIndices
                                                    # if len(haplotypes_dict[tmp_hap].SLFsIndices) != \
                                                    #         len(haplotypes_dict[hap_ind].SLFsIndices):
                                                    #     b = 6
                                                    # what if SLFs where deleted or duplicated in addition to the mutation
                                                    # check if the number of total SLFs remain the same
                                                    # although it can still be equal number of deletions and duplications
                                                    # but this is a rare case -so it is not taken into consideration
                                                    orig_slf_ind = orig_slfs_lst[slf_ii]
                                                    slf_AA = AA_slfs_dict[orig_slf_ind]
                                                    E = np.zeros(len(rnases_orig_gro))
                                                    for rnase_ii, rnase_ind in enumerate(rnases_orig_gro):
                                                        rnase_AA = AA_rnases_dict[rnase_ind]
                                                        if tuple([rnase_ind, orig_slf_ind]) not in E_Ri_Fj_dict:
                                                            E[rnase_ii] = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                                                        else:
                                                            E[rnase_ii] = E_Ri_Fj_dict[rnase_ind, orig_slf_ind]
                                                    if all(E < interaction_thresh):
                                                        # interacting- the mutation in SLF was nonfunctional
                                                        interact_hap = parent_hap
                                                        interact_hap_birth = haplotypes_dict[tmp_hap].BirthIter
                                                        interact_slf = orig_slf_ind
                                                        if parents_haps[tmp_hap] not in haplotypes_dict:
                                                            slf_exists_in_parent_hap_alive = False
                                                            slf_already_exist = False
                                                    else:
                                                        slf_exists_in_parent_hap_alive = False
                                                        slf_already_exist = False
                                                b = 6
                                            else:
                                                slf_already_exist = False

                                #end of while
                                if hap_ind not in hap_in_mut_pollinate_rnase_in_orig:
                                    hap_in_mut_pollinate_rnase_in_orig[hap_ind] = []
                                if analyze_func:
                                    # check if the SLF before mutation was functional
                                    # (decompose any of the rest of the RNases in the population))

                                    slf_before_mut = haplotypes_dict_previous[parent_hap].SLFsIndices[
                                        slf_ii]
                                    RNases_before_mut = []
                                    for h_keys, h_vals in haplotypes_dict_previous.items():
                                        RNases_before_mut.append(h_vals.RNasesIndices[0])
                                    RNases_before_mut = np.unique(RNases_before_mut)
                                    slf_AA_ = AA_slfs_dict[slf_before_mut]
                                    func_slf_cond = False
                                    for r_ind in RNases_before_mut:
                                        rnase_AA_ = AA_rnases_dict[r_ind]
                                        if tuple([r_ind, slf_before_mut]) not in E_Ri_Fj_dict:
                                            E_ = np.sum(Model.eij[rnase_AA_, slf_AA_], axis=0)
                                        else:
                                            E_ = E_Ri_Fj_dict[r_ind, slf_before_mut]
                                        if E_ < interaction_thresh:
                                            # interacting- the mutation in SLF didnt add functionality
                                            b = 6
                                            num_func_orig_slf_from_mut_group += 1
                                            func_slf_cond = True
                                            break
                                    if not func_slf_cond:
                                        num_non_func_orig_slf_from_mut_group += 1
                                    b = 6

                                tmp_dict = {}
                                if analyze_func:
                                    tmp_dict['is orig SLF func'] = E_ < interaction_thresh  # true if was func, false if was not func
                                tmp_dict['interacting SLF mutant'] = interact_slf
                                tmp_dict['interacting hap birth'] = interact_hap_birth
                                tmp_dict['interacting hap count'] = interact_hap_count
                                hap_in_mut_pollinate_rnase_in_orig[hap_ind].append(tmp_dict)

                            else: #non-interacting SLF
                                b = 6
                            b=6 # for slfs in hap
                        b=6 # for rnases in the orig group
                    if len(hap_in_mut_pollinate_rnase_in_orig) == 0:
                        continue

                    analyze_split_group_all_iters_dict[iter_ind] = {}
                    analyze_split_group_all_iters_dict[iter_ind]['hap_in_mut_pollinate_rnase_in_orig'] = \
                        hap_in_mut_pollinate_rnase_in_orig
                    analyze_split_group_all_iters_dict[iter_ind]['hap_in_orig_pollinate_rnase_in_mut'] = \
                        hap_in_orig_pollinate_rnase_in_mut
                    slf_mut_birth_arr = np.zeros(len(hap_in_mut_pollinate_rnase_in_orig), dtype=int)
                    slf_mut_ind_arr = np.zeros(len(slf_mut_birth_arr), dtype=int)
                    SC_slf_in_mut = 0
                    iter_first_slf_mut = 1000000
                    for hap_ii, [hap_ind, slf_mut_info] in enumerate(
                            hap_in_mut_pollinate_rnase_in_orig.items()):
                        for ii_slf, slf_in_hap in enumerate(slf_mut_info):
                             iter_first_slf_mut = np.min([iter_first_slf_mut, slf_in_hap['interacting hap birth']])
                    for hap_ii, [hap_ind, slf_mut_info] in enumerate(
                            hap_in_mut_pollinate_rnase_in_orig.items()):
                        slfs_hap_birth_arr = np.zeros(len(slf_mut_info), dtype=int)
                        slfs_hap_orig_ind_arr = np.zeros(len(slf_mut_info), dtype=int)
                        for ii_slf, slf_in_hap in enumerate(slf_mut_info):
                            slfs_hap_birth_arr[ii_slf] = slf_in_hap['interacting hap birth']
                            slfs_hap_orig_ind_arr[ii_slf] = slf_in_hap['interacting SLF mutant']
                            if hap_ind in haplotypes_dict:
                                rnase = haplotypes_dict[hap_ind].RNasesIndices[0]
                            else:
                                with open(output_dir + '/haplotypes_dict_' + str(int(np.ceil(slfs_hap_birth_arr[ii_slf]/d_iters)*d_iters)) + '.pkl', 'rb') as input:
                                    haplotypes_dict1 = pickle.load(input)
                                if hap_ind in haplotypes_dict1:
                                    rnase = haplotypes_dict1[hap_ind].RNasesIndices[0]
                                else:
                                    continue
                            rnase_AA = AA_rnases_dict[rnase]
                            slf_AA = AA_slfs_dict[slfs_hap_orig_ind_arr[ii_slf]]
                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                            if E < interaction_thresh:
                                # interacting- the mutation in SLF caused a SC
                                if rnase in rnases_in_gro_after_split and slfs_hap_birth_arr[ii_slf]==iter_first_slf_mut:
                                    SC_slf_in_mut += 1
                        slf_in_orig = slfs_hap_orig_ind_arr[np.argmin(slfs_hap_birth_arr)]
                        slf_mut_birth_arr[hap_ii] = np.min(slfs_hap_birth_arr)
                        slf_mut_ind_arr[hap_ii] = slf_in_orig
                        b=6

                    slf_orig_birth_arr = np.zeros(len(hap_in_orig_pollinate_rnase_in_mut), dtype=int)
                    slf_orig_ind_arr = np.zeros(len(slf_orig_birth_arr), dtype=int)
                    SC_slf_in_orig = 0
                    iter_first_slf_mut = 1000000
                    for hap_ii, [hap_ind, slf_mut_info] in enumerate(
                            hap_in_orig_pollinate_rnase_in_mut.items()):
                        for ii_slf, slf_in_hap in enumerate(slf_mut_info):
                            iter_first_slf_mut = np.min([iter_first_slf_mut, slf_in_hap['interacting hap birth']])
                    for hap_ii, [hap_ind, slf_mut_info] in enumerate(
                            hap_in_orig_pollinate_rnase_in_mut.items()):
                        slfs_hap_birth_arr = np.zeros(len(slf_mut_info), dtype=int)
                        slfs_hap_orig_ind_arr= np.zeros(len(slf_mut_info), dtype=int)
                        for ii_slf, slf_in_hap in enumerate (slf_mut_info):
                            slfs_hap_birth_arr[ii_slf] = slf_in_hap['interacting hap birth']
                            slfs_hap_orig_ind_arr[ii_slf] = slf_in_hap['interacting SLF mutant']
                            if hap_ind in haplotypes_dict:
                                rnase = haplotypes_dict[hap_ind].RNasesIndices[0]
                            else:
                                with open(output_dir + '/haplotypes_dict_' + str(int(np.ceil(slfs_hap_birth_arr[ii_slf]/d_iters)*d_iters)) + '.pkl', 'rb') as input:
                                    haplotypes_dict1 = pickle.load(input)
                                if hap_ind in haplotypes_dict1:
                                    rnase = haplotypes_dict1[hap_ind].RNasesIndices[0]
                                else:
                                    continue
                            rnase_AA = AA_rnases_dict[rnase]
                            slf_AA = AA_slfs_dict[slfs_hap_orig_ind_arr[ii_slf]]
                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                            if E < interaction_thresh:
                                # interacting- the mutation in SLF caused a SC
                                if rnase in rnases_in_gro_after_split and slfs_hap_birth_arr[ii_slf]==iter_first_slf_mut:
                                    SC_slf_in_orig += 1 #the first SLF mut cause a SC (its RNase in one of the splitted groups)
                        slf_in_orig = slfs_hap_orig_ind_arr[np.argmin(slfs_hap_birth_arr)]
                        slf_orig_birth_arr[hap_ii] = np.min(slfs_hap_birth_arr)
                        slf_orig_ind_arr[hap_ii] = slf_in_orig
                        b=6

                    b = 6
                    tmp_dict = {}
                    split_proc_ind_lst = [0, 1, 2]
                    tmp_dict['rnase_mut'] = mut_rnase_iter
                    tmp_dict['split'] = iter_ind
                    tmp_dict['first_iter_slf_mut_in_mut_gro_full_comp'] = np.min(slf_mut_birth_arr)
                    tmp_dict['first_slf_mut_in_mut_gro'] = slf_mut_ind_arr[np.argmin(slf_mut_birth_arr)]
                    tmp_dict['first_iter_slf_mut_in_orig_gro_full_comp'] = np.min(slf_orig_birth_arr)
                    tmp_dict['first_slf_mut_in_orig_gro'] = slf_orig_ind_arr[np.argmin(slf_orig_birth_arr)]
                    first_cros_mut_hap_in_mut_gro  = list(hap_in_mut_pollinate_rnase_in_orig.keys())[
                        np.argmin(slf_mut_birth_arr)]
                    first_cros_mut_hap_in_orig_gro = list(hap_in_orig_pollinate_rnase_in_mut.keys())[
                        np.argmin(slf_orig_birth_arr)]
                    if iter_ind not in analyze_split_rnase_cros_slfs_info:
                        analyze_split_rnase_cros_slfs_info[iter_ind] = []

                    split_cascade_dict = {}
                    split_cascade_dict['orig_gro'] = gro
                    split_cascade_dict['orig_rnase'] = orig_rnase
                    split_cascade_dict['mut_rnase'] = mut_rnase
                    split_cascade_dict['mut_rnase_hap'] = mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['mutant hap']
                    split_cascade_dict['orig_rnase_hap'] = mut_rnases_dict[sec_iter_ind][ii_mut_rnase]['orig hap']
                    split_cascade_dict['rnase_mut_iter'] = sec_iter_ind
                    split_cascade_dict['iter_cros_mut_slf_mut_gro'] = np.min(slf_mut_birth_arr)
                    split_cascade_dict['cros_mut_slf_mut_gro'] = slf_mut_ind_arr[np.argmin(slf_mut_birth_arr)]
                    split_cascade_dict['hap_cros_mut_slf_mut_gro'] = first_cros_mut_hap_in_mut_gro
                    split_cascade_dict['iter_cros_mut_slf_orig_gro'] = np.min(slf_orig_birth_arr)
                    split_cascade_dict['cros_mut_slf_orig_gro'] = slf_orig_ind_arr[np.argmin(slf_orig_birth_arr)]
                    split_cascade_dict['hap_cros_mut_slf_orig_gro'] = first_cros_mut_hap_in_orig_gro

                    split_proc_iters = [tmp_dict['rnase_mut'], tmp_dict['first_iter_slf_mut_in_mut_gro_full_comp'],
                                        tmp_dict['first_iter_slf_mut_in_orig_gro_full_comp']]
                    if all(np.array(split_proc_ind_lst)[np.argsort(split_proc_iters)] ==
                           np.array([2, 0, 1])): #SLF_orig-> RNase-> SLF_mut-> SPLIT
                        split_cascade_dict['mut cascade'] = 'SLF_orig-> RNase-> SLF_mut-> SPLIT'
                        split_cascade_dict['mut iters'] = split_proc_iters
                        slf_orig_rnase_slf_mut_split += 1
                        if SC_slf_in_orig > 0:
                            slf_orig_rnase_slf_mut_split_SC_slf_orig += 1
                        if SC_slf_in_mut > 0:
                            slf_orig_rnase_slf_mut_split_SC_slf_mut += 1

                    if all(np.array(split_proc_ind_lst)[np.argsort(split_proc_iters)] ==
                            np.array([2, 1, 0])): #SLForig->SLFmut->RNase
                        split_cascade_dict['mut cascade'] = 'SLF_orig-> SLF_mut-> RNase-> SPLIT'
                        split_cascade_dict['mut iters'] = split_proc_iters
                        slf_orig_slf_mut_rnase_split += 1
                        if SC_slf_in_orig > 0:
                            slf_orig_slf_mut_rnase_split_SC_slf_orig += 1
                        if SC_slf_in_mut > 0:
                            slf_orig_slf_mut_rnase_split_SC_slf_mut += 1
                        else:
                            # check if there are more than one RNase in the orig group
                            if len(rnases_orig_gro) > 1:
                                slf_orig_slf_mut_rnase_split_NSC_more_than_one_rnase_in_orig_gro += 1

                    if all(np.array(split_proc_ind_lst)[np.argsort(split_proc_iters)] ==
                                 np.array([1, 0, 2])):#SLFmut->RNase->SLForig
                        split_cascade_dict['mut cascade'] = 'SLF_mut-> RNase-> SLF_orig-> SPLIT'
                        split_cascade_dict['mut iters'] = split_proc_iters
                        slf_mut_rnase_slf_orig_split += 1
                        if SC_slf_in_orig > 0:
                            slf_mut_rnase_slf_orig_split_SC_slf_orig += 1
                        if SC_slf_in_mut > 0:
                            slf_mut_rnase_slf_orig_split_SC_slf_mut += 1
                        else:
                            # check if there are more than one RNase in the orig group
                            if len(rnases_orig_gro) > 1:
                                slf_mut_rnase_slf_orig_split_NSC_more_than_one_rnase_in_orig_gro += 1

                    if all(np.array(split_proc_ind_lst)[np.argsort(split_proc_iters)] ==
                             np.array([1, 2, 0])):#SLFmut-> SLForig-> RNase
                        split_cascade_dict['mut cascade'] = 'SLF_mut-> SLF_orig-> RNase-> SPLIT'
                        split_cascade_dict['mut iters'] = split_proc_iters
                        slf_mut_slf_orig_rnase_split += 1
                        if SC_slf_in_orig > 0:
                            slf_mut_slf_orig_rnase_split_SC_slf_orig += 1
                        if SC_slf_in_mut > 0:
                            slf_mut_slf_orig_rnase_split_SC_slf_mut += 1
                        else:
                            #check if there are more than one RNase in the orig group
                            if len(rnases_orig_gro) > 1:
                                slf_mut_slf_orig_rnase_split_NSC_more_than_one_rnase_in_orig_gro += 1

                    if all(np.array(split_proc_ind_lst)[np.argsort(split_proc_iters)] ==
                           np.array([0, 1, 2])):# RNase-> SLFmut-> SLForig
                        split_cascade_dict['mut cascade'] = 'RNase-> SLF_mut-> SLF_orig-> SPLIT'
                        split_cascade_dict['mut iters'] = split_proc_iters
                        if SC_slf_in_orig > 0:
                            rnase_slf_mut_slf_orig_split_SC_slf_orig += 1
                        if SC_slf_in_mut > 0:
                            rnase_slf_mut_slf_orig_split_SC_slf_mut += 1
                        rnase_slf_mut_slf_orig_split += 1
                        #check SC between the RNase mutant and its SLFs
                        RNase_mut_iter = split_proc_iters[0]
                        with open(output_dir + '/haplotypes_dict_' + str(sec_iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_dict = pickle.load(input)
                        rnase_AA = AA_rnases_dict[mut_rnase]
                        SLFs_in_mut_hap = haplotypes_dict[mut_hap].SLFsIndices
                        is_SC = False
                        for slf in SLFs_in_mut_hap:
                            slf_AA = AA_slfs_dict[slf]
                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                            if E < interaction_thresh:
                                # interacting- the mutation in RNase caused a SC
                                is_SC = True
                                break
                        if is_SC:
                            rnase_slf_mut_slf_orig_split_SC_rnase += 1

                    if all(np.array(split_proc_ind_lst)[np.argsort(split_proc_iters)] ==
                           np.array([0, 2, 1])): ## RNase-> SLForig-> SLFmut
                        split_cascade_dict['mut cascade'] = 'RNase-> SLF_orig-> SLF_mut-> SPLIT'
                        split_cascade_dict['mut iters'] = split_proc_iters
                        if SC_slf_in_orig > 0:
                            rnase_slf_orig_slf_mut_split_SC_slf_orig += 1
                        if SC_slf_in_mut > 0:
                            rnase_slf_orig_slf_mut_split_SC_slf_mut += 1
                        rnase_slf_orig_slf_mut_split += 1
                        # check SC between the RNase mutant and its SLFs
                        RNase_mut_iter = split_proc_iters[0]
                        with open(output_dir + '/haplotypes_dict_' + str(sec_iter_ind) + '.pkl', 'rb') as input:
                            haplotypes_dict = pickle.load(input)
                        rnase_AA = AA_rnases_dict[mut_rnase]
                        SLFs_in_mut_hap = haplotypes_dict[mut_hap].SLFsIndices
                        is_SC = False
                        for slf in SLFs_in_mut_hap:
                            slf_AA = AA_slfs_dict[slf]
                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                            if E < interaction_thresh:
                                # interacting- the mutation in RNase caused a SC
                                is_SC = True
                                break
                        if is_SC:
                            rnase_slf_orig_slf_mut_split_SC_rnase += 1

                    analyze_split_rnase_cros_slfs_info[iter_ind].append(split_cascade_dict)
                else:
                    continue

            diversity_mut_dict[subdirectory] = diversity_mut
            diversity_orig_dict[subdirectory] = diversity_orig
            if analyze_func:
                if num_func_orig_slf_from_mut_group != 0:
                    frac_func_orig_slf_from_mut_group[ii_dir] = num_func_orig_slf_from_mut_group/\
                            (num_func_orig_slf_from_mut_group + num_non_func_orig_slf_from_mut_group)
                if num_func_orig_slf_from_orig_group != 0:
                    frac_func_orig_slf_from_orig_group[ii_dir] = num_func_orig_slf_from_orig_group / \
                            (num_func_orig_slf_from_orig_group + num_non_func_orig_slf_from_orig_group)

            if analyze_func:
                with open(directory + '/frac_func_orig_slf_from_mut_group_ver2.pkl', 'wb') as output:
                    pickle.dump(frac_func_orig_slf_from_mut_group, output, pickle.HIGHEST_PROTOCOL)
                with open(directory + '/frac_func_orig_slf_from_orig_group_ver2.pkl', 'wb') as output:
                    pickle.dump(frac_func_orig_slf_from_orig_group, output, pickle.HIGHEST_PROTOCOL)

            with open(output_dir + '/analyze_split_group_all_iters_dict_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
                pickle.dump(analyze_split_group_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(output_dir + '/analyze_split_rnase_cros_slfs_info_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
                pickle.dump(analyze_split_rnase_cros_slfs_info, output, pickle.HIGHEST_PROTOCOL)

b=6
with open(directory + '/diversity_mut_dict_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(diversity_mut_dict, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/diversity_orig_dict_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(diversity_orig_dict, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/slf_orig_rnase_slf_mut_split_SC_slf_orig_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_rnase_slf_mut_split_SC_slf_orig, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_orig_rnase_slf_mut_split_SC_slf_mut_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_rnase_slf_mut_split_SC_slf_mut, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_orig_slf_mut_rnase_split_SC_slf_orig_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_slf_mut_rnase_split_SC_slf_orig, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_orig_slf_mut_rnase_split_SC_slf_mut_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_slf_mut_rnase_split_SC_slf_mut, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_mut_rnase_slf_orig_split_SC_slf_orig_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_mut_rnase_slf_orig_split_SC_slf_orig, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_mut_rnase_slf_orig_split_SC_slf_mut_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_mut_rnase_slf_orig_split_SC_slf_mut, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_mut_slf_orig_rnase_split_SC_slf_mut_' + str(min_duration * 10) + '_iters.pkl',
          'wb') as output:
    pickle.dump(slf_mut_slf_orig_rnase_split_SC_slf_mut, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/slf_mut_slf_orig_rnase_split_SC_slf_orig_' + str(min_duration * 10) + '_iters.pkl',
          'wb') as output:
    pickle.dump(slf_mut_slf_orig_rnase_split_SC_slf_orig, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/rnase_slf_mut_slf_orig_split_SC_rnase_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(rnase_slf_mut_slf_orig_split_SC_rnase, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/rnase_slf_orig_slf_mut_split_SC_rnase_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(rnase_slf_orig_slf_mut_split_SC_rnase, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/slf_mut_slf_orig_rnase_split_NSC_more_than_one_rnase_in_orig_gro_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_mut_slf_orig_rnase_split_NSC_more_than_one_rnase_in_orig_gro, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/slf_mut_rnase_slf_orig_split_NSC_more_than_one_rnase_in_orig_gro_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_mut_rnase_slf_orig_split_NSC_more_than_one_rnase_in_orig_gro, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/slf_orig_slf_mut_rnase_split_NSC_more_than_one_rnase_in_orig_gro_' + str(min_duration*10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_slf_mut_rnase_split_NSC_more_than_one_rnase_in_orig_gro, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/rnase_slf_mut_slf_orig_split_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(rnase_slf_mut_slf_orig_split, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/rnase_slf_orig_slf_mut_split_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(rnase_slf_orig_slf_mut_split, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_mut_rnase_slf_orig_split_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_mut_rnase_slf_orig_split, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_orig_rnase_slf_mut_split_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_rnase_slf_mut_split, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_orig_slf_mut_rnase_split_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_orig_slf_mut_rnase_split, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/slf_mut_slf_orig_rnase_split_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
    pickle.dump(slf_mut_slf_orig_rnase_split, output, pickle.HIGHEST_PROTOCOL)

b=6
