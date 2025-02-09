import numpy as np
import pickle
import os
import copy
import funcs_track_comp_groups as funcs
import calc_weighted_num_partners_of_allele as weighted_num_partners
min_duration = 10 #10 #5 #1 #
d_iters = 10
thresh_num_hap_in_group = 10
os.chdir("..")
os.chdir("FlowersEvolutionsSimulation_ver13_data")
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4"  # output_data_10haplotypes_k5_" + str(L) + "AAs_E_thresh-6_p2_2_10-5_mut_2" #test2" #test" #
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_alpha_08_delta_09"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_alpha_095_delta_085"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p2times10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_phalf_times10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_no_dup_no_del"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_alpha_06_delta_09"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-4_p10-4_rerun"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-4_p10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_2000haps"
subdirectories = list(set(next(os.walk(directory))[1]) - set(['2023-07-03-10-00-56-2']))#
# subdirectories = ['2023-06-25-15-01-01-11']
# subdirectories = ['2023-06-25-15-01-01-7']
# subdirectories = ['2023-06-25-22-40-35-5']

num_files = len(subdirectories)
for ii_dir, subdirectory in enumerate(subdirectories):
    output_dir = directory + "/" + subdirectory
    if os.path.isfile(output_dir + '/mis_groups_hap_all_iters_dict.pkl') \
            and os.path.isfile(output_dir + '/parents_haps.pkl'): #\
            #and not os.path.isfile(output_dir + '/reorder_mis_hap_strict_group_extension_all_iters_dict.pkl'):
            print(os.path.join(subdirectory))

            with open(output_dir + '/RNasesInHaplotypes_all_iters_Dict.pkl', 'rb') as input:
                RNasesInHaplotypes_all_iters_Dict = pickle.load(input)

            reorder_mis_strict_groups_hap_all_iters_dict = {}
            reorder_mis_strict_groups_rnase_all_iters_dict = {}
            reorder_mis_rnase_strict_groups_all_iters_dict = {}
            reorder_mis_hap_strict_groups_all_iters_dict = {}
            reorder_mis_strict_group_extension_hap_all_iters_dict = {}
            reorder_mis_hap_strict_group_extension_all_iters_dict = {}

            group_naming_current_to_new_all_iters_dict = {}

            mut_rnases_dict = {}

            # rnase_reordered_groups_dict = {}
            split_one_group_into_two_all_iters_dict = {}
            extinct_group_all_iters_dict = {}
            new_group_all_iters_dict = {}
            new_to_current_group_labeling_all_iters_dict = {}
            current_to_new_group_labeling_all_iters_dict = {}
            mis_strict_groups_former_parent_hap_all_iters_dict = {}
            last_iter = 0

            #-------------------------------------------
            '''
            last_iter = 0
            not_able_to_load = True
            while not_able_to_load:
                try:
                    open(output_dir + '/haplotypes_dict_' + str(last_iter) + '.pkl', 'rb')
                    last_iter += 500
                except:
                    not_able_to_load = False
                    last_iter -= 500
            not_able_to_load = True
            while not_able_to_load:
                try:
                    open(output_dir + '/haplotypes_dict_' + str(last_iter) + '.pkl', 'rb')
                    last_iter += 10
                except:
                    not_able_to_load = False
                    last_iter -= 10
            # detect the first iter used (in smaller saving files rate )
            first_iter = last_iter
            not_able_to_load = True
            while not_able_to_load:
                try:
                    open(output_dir + '/haplotypes_dict_' + str(first_iter) + '.pkl', 'rb')
                    first_iter -= 10
                except:
                    not_able_to_load = False
                    first_iter += 10
            
            if first_iter != last_iter:
                # last_iter = 10000
                if first_iter < 0:
                    # -------------------------------------------
                    files_lst = os.listdir(output_dir)
                    iters_lst = []
                    for file_str in files_lst:
                        if 'diploid_count_dict_' in file_str:
                            iters_lst.append(int(file_str[19:-4]))
                    if len(iters_lst) > 0:
                        iters_arr = np.array(iters_lst)
                        iters_arr = iters_arr[iters_arr % d_iters == 0]
                        last_iter = np.max(iters_arr)
                        first_iter = np.min(iters_arr)
                        # -------------------------------------------
            # '''
            #-------------------------------------------
            with open(output_dir + '/mis_groups_hap_all_iters_dict.pkl', 'rb') as input:
                mis_groups_hap_all_iters_dict = pickle.load(input)
            with open(output_dir + '/mis_groups_rnase_all_iters_dict.pkl', 'rb') as input:
                mis_groups_rnase_all_iters_dict = pickle.load(input)
            with open(output_dir + '/mis_rnase_groups_all_iters_dict.pkl', 'rb') as input:
                mis_rnase_groups_all_iters_dict = pickle.load(input)
            with open(output_dir + '/mis_hap_groups_all_iters_dict.pkl', 'rb') as input:
                mis_hap_groups_all_iters_dict = pickle.load(input)

            with open(output_dir + '/parents_haps.pkl', 'rb') as input:
                parents_haps = pickle.load(input)

            first_iter = list(mis_groups_hap_all_iters_dict.keys())[0]
            last_iter = list(mis_groups_hap_all_iters_dict.keys())[-1]
            try:
                with open(output_dir + '/haplotypes_dict_' + str(first_iter) + '.pkl', 'rb') as input:
                    haplotypes_dict = pickle.load(input)
            except:
                b=6
            reorder_mis_strict_groups_hap_all_iters_dict[first_iter] = copy.deepcopy(mis_groups_hap_all_iters_dict[first_iter])
            reorder_mis_strict_groups_rnase_all_iters_dict[first_iter] = copy.deepcopy(mis_groups_rnase_all_iters_dict[first_iter])
            reorder_mis_rnase_strict_groups_all_iters_dict[first_iter] = copy.deepcopy(mis_rnase_groups_all_iters_dict[first_iter])
            reorder_mis_hap_strict_groups_all_iters_dict[first_iter] = copy.deepcopy(mis_hap_groups_all_iters_dict[first_iter])
            reorder_mis_strict_group_extension_hap_all_iters_dict[first_iter] = {}
            reorder_mis_hap_strict_group_extension_all_iters_dict[first_iter] = {}

            group_naming_current_to_new_dict = {}
            # group_naming_new_to_current_dict = {}
            for ii in range(len(mis_groups_hap_all_iters_dict[first_iter])):
                group_naming_current_to_new_dict[ii] = ii
            group_naming_current_to_new_all_iters_dict[first_iter] = group_naming_current_to_new_dict
            new_any_rnases_lst = []
            t_new_any_RNase = []
            max_group_ind = np.max(list(reorder_mis_strict_groups_hap_all_iters_dict[first_iter].keys()))

            for iter_ind in range(first_iter + d_iters, last_iter - d_iters*min_duration, d_iters):
                # print(iter_ind)
                with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                    haplotypes_dict = pickle.load(input)
                with open(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                    haplotypes_count_dict = pickle.load(input)

                group_naming_current_to_new_dict = {}
                # group_naming_new_to_current_dict = {}
                regroup_current_dict = {}

                previous_iter = iter_ind - d_iters
                with open(output_dir + '/haplotypes_dict_' + str(previous_iter) + '.pkl', 'rb') as input:
                    haplotypes_dict_previous = pickle.load(input)

                #--------------------------------------------------------
                # determine groups re-naming based on the last iteration
                #--------------------------------------------------------
                current_to_new_group_labeling_dict = {}
                n = 1

                merge_groups_to_one_dict = {}
                for current_group_ind, current_haps_lst in mis_groups_hap_all_iters_dict[iter_ind].items():

                    parents_haps_arr = np.zeros(len(current_haps_lst), dtype=int)
                    current_to_new_group_labeling = np.zeros(len(current_haps_lst))
                    current_to_new_group_labeling[:] = np.nan
                    for curr_hap_ii, current_hap in enumerate(current_haps_lst):

                        parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous, current_hap,
                                                                     parents_haps)

                        if parent_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters]:
                            label_by_grouping = reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters][parent_hap]
                        else:
                            label_by_grouping = np.nan
                        if parent_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind - n * d_iters]:
                            label_by_extenction = reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind - n * d_iters][
                                parent_hap]
                        else:
                            label_by_extenction = np.nan
                        if not np.isnan(label_by_grouping) and not np.isnan(label_by_extenction):
                            if label_by_grouping == label_by_extenction:
                                label = label_by_grouping
                                current_to_new_group_labeling[curr_hap_ii] = label
                            else:
                                b = 6
                        elif not np.isnan(label_by_grouping) and np.isnan(label_by_extenction):
                            label = label_by_grouping
                            current_to_new_group_labeling[curr_hap_ii] = label
                        elif np.isnan(label_by_grouping) and not np.isnan(label_by_extenction):
                            label = label_by_extenction
                            current_to_new_group_labeling[curr_hap_ii] = label
                        else:
                            b=6
                            # print('partner is neither in a group nor in an extension')

                    current_to_new_group_labeling_dict[current_group_ind] = list(current_to_new_group_labeling)

                    group_label = np.unique(current_to_new_group_labeling
                                            [~ np.isnan(current_to_new_group_labeling)])

                    if len(group_label) == 1:
                        # all the haps parents are labeled by the same group
                        group_naming_current_to_new_dict[current_group_ind] = \
                            int(group_label[0])

                    elif len(group_label) > 1: # more than one labeling within one group
                        group_naming_current_to_new_dict[current_group_ind] = []  # 'unite'
                        merge_groups_to_one_dict[current_group_ind] = np.unique(
                            group_label)
                        b=6
                    else:
                        b=6
                        group_naming_current_to_new_dict[current_group_ind] = max_group_ind + 1
                        max_group_ind += 1
                #---------------------------------------------------
                #  merge groups
                #---------------------------------------------------
                if len(merge_groups_to_one_dict)>0:
                    for current_group_ind, new_group_lst in merge_groups_to_one_dict.items():
                        is_group_exist_arr = np.zeros(len(new_group_lst), dtype = int)
                        #all the united groups are also groups of their own
                        for gro_ii, gro_ind in enumerate(new_group_lst):
                            if gro_ind in list(group_naming_current_to_new_dict.values()):
                                is_group_exist_arr[gro_ii] = 1

                        if all(is_group_exist_arr == 1):
                            # both the original groups are still alive so the unite group of offspring haplotype,
                            # originated in the two groups, should be called by a new index
                            b=6
                            gro_ind = max_group_ind + 1

                            # print('except for merging, both legitimate groups survived')
                        elif any(is_group_exist_arr == 0):
                            if len(is_group_exist_arr[is_group_exist_arr==0]) == 1:
                                # one of the original merged groups still exist and one does not exist anymore.
                                # give the merged group the name of the one that does not exist anymore
                                gro_ind = int(new_group_lst[np.where(is_group_exist_arr == 0)][0])
                                # print('except for merging, one legitimate group survived')
                            else:
                                # both the groups does not exist anymore, give the merged group the name of the bigger one
                                max_tmp = 0
                                group_size = np.zeros(len(new_group_lst), dtype=int)
                                for ii_gro, gro_ind in enumerate(new_group_lst):
                                    gro_ind = int(gro_ind)
                                    hap_lst_in_current = mis_groups_hap_all_iters_dict[iter_ind][current_group_ind]
                                    hap_lst = np.array(hap_lst_in_current)[np.array(current_to_new_group_labeling_dict[current_group_ind]) == gro_ind]
                                    # hap_lst = np.array(current_haps_lst)[np.array(np.where(
                                    #     current_to_new_group_labeling == gro_ind)[0])]
                                    for hap in hap_lst:
                                        group_size[ii_gro] += haplotypes_count_dict[hap]
                                    if group_size[ii_gro] > max_tmp:
                                        tmp_gro = gro_ind
                                        max_tmp = group_size[ii_gro]
                                    b = 6
                                b = 6
                                gro_ind = tmp_gro
                                # print('except for merging, none of the legitimate groups survived')
                                #the current group index should be the bigger among the two groups in the iteration before
                        # not_chosen_gro = int(list(set(new_group_lst) - set([gro_ind]))[0])

                        group_naming_current_to_new_dict[current_group_ind] = gro_ind
                        max_group_ind += 1

                        for ii_gro, gro_ind in enumerate(new_group_lst):
                            gro_ind = int(gro_ind)
                            hap_lst_in_current = mis_groups_hap_all_iters_dict[iter_ind][current_group_ind]
                            hap_lst = np.array(hap_lst_in_current)[
                                np.array(current_to_new_group_labeling_dict[current_group_ind]) == gro_ind]
                            is_origin_group_arr = np.zeros(len(hap_lst), dtype=int)
                            is_same_as_orig_arr = np.zeros(len(hap_lst), dtype=int)
                            is_orig_alive_arr = np.zeros(len(hap_lst))
                            is_orig_alive_arr[:] = np.nan
                            for hap_ii, hap in enumerate(hap_lst):
                                # parent_hap_ = funcs.parent_hap_func(haplotypes_dict, haplotypes_dict_previous,
                                #                                    hap, parents_haps, iter_ind, d_iters, n)
                                parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous, hap, parents_haps)
                                # if parent_hap_ != parent_hap:
                                #     b=6
                                if parent_hap == hap:
                                    is_same_as_orig_arr[hap_ii] = 1
                                    is_orig_alive_arr[hap_ii] = 1
                                else:
                                    # check if the parent is still alive
                                    b = 6
                                    if parent_hap in list(haplotypes_dict.keys()):
                                        is_orig_alive_arr[hap_ii] = 1
                                    else:
                                        is_orig_alive_arr[hap_ii] = 0
                                if parent_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters]:
                                    label_by_grouping = \
                                        reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters][parent_hap]
                                    is_origin_group_arr[hap_ii] = 1
                                else:
                                    label_by_grouping = np.nan
                                if parent_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[
                                    iter_ind - n * d_iters]:
                                    label_by_extenction = \
                                        reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind - n * d_iters][
                                            parent_hap]
                                else:
                                    label_by_extenction = np.nan


                    b=6
                #---------------------------------------------------------
                #    group splitting
                #---------------------------------------------------------
                uniq, c = np.unique(list(group_naming_current_to_new_dict.values()), return_counts=True)
                if any(c > 1):
                    # one group in the previous iter was splitted into two groups
                    # save the iter and the groups indices
                    b = 6
                    for gro in uniq[np.where(c > 1)]:
                        current_gro = []
                        for key, val in group_naming_current_to_new_dict.items():
                            if val == gro:
                                current_gro.append(key)

                        min_tmp = 1000
                        group_size = np.zeros(len(current_gro), dtype=int)
                        for ii_gro, gro_ind in enumerate(current_gro):
                            hap_lst = mis_groups_hap_all_iters_dict[iter_ind][gro_ind]
                            for hap in hap_lst:
                                group_size[ii_gro] += haplotypes_count_dict[hap]
                            if group_size[ii_gro] < min_tmp:
                                curr_gro = gro_ind
                                min_tmp = group_size[ii_gro]
                            b=6
                        b=6

                        # replace the double labeling of the smaller group(among the two) by a new labeling
                        group_naming_current_to_new_dict[curr_gro] = \
                            max_group_ind + 1
                        tmp_dict = {}
                        tmp_dict[gro] = []
                        for ii_gro, cur_gro in enumerate(current_gro):
                            hap_lst = mis_groups_hap_all_iters_dict[iter_ind][cur_gro]
                            is_origin_group_arr = np.zeros(len(hap_lst), dtype=int)
                            is_same_as_orig_arr = np.zeros(len(hap_lst), dtype=int)
                            is_orig_alive_arr = np.zeros(len(hap_lst))
                            is_orig_alive_arr[:] = np.nan
                            for hap_ii, hap in enumerate(hap_lst):
                                # parent_hap_ = funcs.parent_hap_func(haplotypes_dict, haplotypes_dict_previous,
                                #                                    hap, parents_haps, iter_ind, d_iters, n)

                                parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous, hap,
                                                                             parents_haps)
                                # if parent_hap_ != parent_hap:
                                #     b=6
                                if parent_hap == hap:
                                    is_same_as_orig_arr[hap_ii] = 1
                                else:
                                    #check if the parent is still alive
                                    b=6
                                    if parent_hap in list(haplotypes_dict.keys()):
                                        is_orig_alive_arr[hap_ii] = 1
                                    else:
                                        is_orig_alive_arr[hap_ii] = 0
                                if parent_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters]:
                                    label_by_grouping = \
                                    reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters][parent_hap]
                                    is_origin_group_arr[hap_ii] = 1

                                else:
                                    label_by_grouping = np.nan
                                if parent_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[
                                    iter_ind - n * d_iters]:
                                    label_by_extenction = \
                                    reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind - n * d_iters][
                                        parent_hap]

                                else:
                                    label_by_extenction = np.nan
                            b=6
                            is_origin_group_dict = {}
                            is_origin_group_dict['is_origin_grouped'] = is_origin_group_arr
                            is_same_as_orig_dict = {}
                            is_same_as_orig_dict['is_same_as_orig'] = is_same_as_orig_arr
                            is_orig_alive_dict = {}
                            is_orig_alive_dict['is_orig_still_alive'] = is_orig_alive_arr
                            num_haps_dict = {}
                            num_haps_dict['num_haps'] = group_size[ii_gro]
                            tmp_dict2 = {}
                            tmp_dict2[group_naming_current_to_new_dict[cur_gro]] = [is_origin_group_dict,
                                                                                    is_same_as_orig_dict,
                                                                                    is_orig_alive_dict,
                                                                                    num_haps_dict]
                            tmp_dict[gro].append(tmp_dict2)
                            b=6
                        if iter_ind not in split_one_group_into_two_all_iters_dict:
                            split_one_group_into_two_all_iters_dict[iter_ind] = []
                        split_one_group_into_two_all_iters_dict[iter_ind].append(tmp_dict)
                            # ([gro, [gro, group_naming_current_to_new_dict[curr_gro]]])
                        max_group_ind += 1
                reorder_mis_groups_rnase_dict = {}
                reorder_mis_groups_hap_dict = {}
                for group_ind, new_group_ind in group_naming_current_to_new_dict.items():
                    reorder_mis_groups_rnase_dict[new_group_ind] = copy.deepcopy(
                        mis_groups_rnase_all_iters_dict[iter_ind][group_ind])
                    reorder_mis_groups_hap_dict[new_group_ind] = copy.deepcopy(
                        mis_groups_hap_all_iters_dict[iter_ind][group_ind])
                    # for rnase_ind in np.unique(mis_strict_groups_rnase_all_iters_dict[iter_ind][group_ind]):
                    #     if rnase_ind not in rnase_reordered_groups_dict:
                    #         rnase_reordered_groups_dict[rnase_ind] = []
                    #     if isinstance(rnase_reordered_groups_dict[rnase_ind], int):
                    #         if rnase_reordered_groups_dict[rnase_ind] != new_group_ind:
                    #             b = 6
                    #     rnase_reordered_groups_dict[rnase_ind] = new_group_ind
                if len(reorder_mis_groups_hap_dict) != len(mis_groups_hap_all_iters_dict[iter_ind]):
                    b = 6
                reorder_mis_hap_groups_dict = {}
                for hap, current_group_name in mis_hap_groups_all_iters_dict[iter_ind].items():
                    former_group_name = group_naming_current_to_new_dict[current_group_name]
                    reorder_mis_hap_groups_dict[hap] = former_group_name

                reorder_mis_rnase_groups_dict = {}
                for rnase, current_group_name in mis_rnase_groups_all_iters_dict[iter_ind].items():
                    former_group_name = group_naming_current_to_new_dict[current_group_name]
                    reorder_mis_rnase_groups_dict[rnase] = former_group_name

                #----------------------------------------------------------
                #   group extentions - non grouped haplotypes
                #----------------------------------------------------------
                reorder_mis_strict_group_extension_hap_dict = {}
                non_grouped_haps = list(set(haplotypes_count_dict.keys())
                                        - set(list(mis_hap_groups_all_iters_dict[iter_ind].keys())))
                for curr_hap_ii, current_hap in enumerate(non_grouped_haps):
                    n = 1
                    previous_iter = iter_ind - d_iters
                    if previous_iter >= first_iter:
                        with open(output_dir + '/haplotypes_dict_' + str(previous_iter) + '.pkl', 'rb') as input:
                            haplotypes_dict_previous = pickle.load(input)
                        try:
                            aa = parents_haps[current_hap]
                        except:
                            b = 6
                        # parent_hap_ = funcs.parent_hap_func(haplotypes_dict, haplotypes_dict_previous,
                        #                                    current_hap, parents_haps, iter_ind, d_iters, n)
                        parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous, current_hap, parents_haps)
                        # if parent_hap_ != parent_hap:
                        #     b=6
                        if parent_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters]:
                            label_by_grouping = reorder_mis_hap_strict_groups_all_iters_dict[iter_ind - n * d_iters][parent_hap]
                        else:
                            label_by_grouping = np.nan
                        if parent_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind - n * d_iters]:
                            label_by_extenction = reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind - n * d_iters][
                                parent_hap]
                        else:
                            label_by_extenction = np.nan
                        if not np.isnan(label_by_grouping) and not np.isnan(label_by_extenction):
                            if label_by_grouping == label_by_extenction:
                                label = label_by_grouping
                                if label not in reorder_mis_strict_group_extension_hap_dict:
                                    reorder_mis_strict_group_extension_hap_dict[label] = []
                                reorder_mis_strict_group_extension_hap_dict[label].append(current_hap)
                            else:
                                b = 6
                        elif not np.isnan(label_by_grouping) and np.isnan(label_by_extenction):
                            label = label_by_grouping
                            if label not in reorder_mis_strict_group_extension_hap_dict:
                                reorder_mis_strict_group_extension_hap_dict[label] = []
                            reorder_mis_strict_group_extension_hap_dict[label].append(current_hap)
                        elif np.isnan(label_by_grouping) and not np.isnan(label_by_extenction):
                            label = label_by_extenction
                            if label not in reorder_mis_strict_group_extension_hap_dict:
                                reorder_mis_strict_group_extension_hap_dict[label] = []
                            reorder_mis_strict_group_extension_hap_dict[label].append(current_hap)
                        else:
                            b=6

                #---------------------------------------------------------
                #      group extinction
                #---------------------------------------------------------
                total_last_iter_group_lst = [
                    list(reorder_mis_strict_groups_hap_all_iters_dict[iter_ind - d_iters].keys()), \
                         list(reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind - d_iters].keys())]
                total_last_iter_group_lst = [item for sublist in total_last_iter_group_lst for item in sublist]
                total_current_iter_group_lst = [
                    list(reorder_mis_groups_hap_dict.keys()), \
                    list(reorder_mis_strict_group_extension_hap_dict.keys())]
                total_current_iter_group_lst = [item for sublist in total_current_iter_group_lst for item in sublist]

                extincted_group_lst = list(set(total_last_iter_group_lst)
                          - set(total_current_iter_group_lst))

                if len(extincted_group_lst) > 0:
                    extinct_group_all_iters_dict[iter_ind] = extincted_group_lst
                    b = 6

                #----------------------------------------------------------------
                #        rename groups labels - using labels in the last iter
                #----------------------------------------------------------------
                reorder_mis_hap_strict_group_extension_dict = {}
                for gro_ind, hap_lst in reorder_mis_strict_group_extension_hap_dict.items():
                    for hap in hap_lst:
                        reorder_mis_hap_strict_group_extension_dict[hap] = gro_ind

                reorder_mis_strict_groups_hap_all_iters_dict[iter_ind] = \
                    copy.deepcopy(reorder_mis_groups_hap_dict)
                reorder_mis_strict_groups_rnase_all_iters_dict[iter_ind] = \
                    copy.deepcopy(reorder_mis_groups_rnase_dict)
                reorder_mis_hap_strict_groups_all_iters_dict[iter_ind] = \
                    copy.deepcopy(reorder_mis_hap_groups_dict)
                reorder_mis_rnase_strict_groups_all_iters_dict[iter_ind] = \
                    copy.deepcopy(reorder_mis_rnase_groups_dict)
                reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind] = \
                    copy.deepcopy(reorder_mis_strict_group_extension_hap_dict)
                reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind] = \
                    copy.deepcopy(reorder_mis_hap_strict_group_extension_dict)
                group_naming_current_to_new_all_iters_dict[iter_ind] = \
                    copy.deepcopy(group_naming_current_to_new_dict)
                # regroup_current_all_iters_dict[iter_ind] = regroup_current_dict
                current_to_new_group_labeling_all_iters_dict[iter_ind] = \
                    copy.deepcopy(current_to_new_group_labeling_dict)

                #'''
                new_rnases_lst = list(set(RNasesInHaplotypes_all_iters_Dict[iter_ind].keys()) - set(
                    RNasesInHaplotypes_all_iters_Dict[iter_ind - d_iters]))
                if len(new_rnases_lst) > 0:
                    for mutant_rnase_ind in new_rnases_lst:
                        live_enough = True
                        for ii in range(min_duration):
                            if mutant_rnase_ind not in RNasesInHaplotypes_all_iters_Dict[iter_ind + d_iters * ii]:
                                live_enough = False
                                break
                        if live_enough:
                            new_any_rnases_lst.append(mutant_rnase_ind)
                            t_new_any_RNase.append(iter_ind)
                            mutant_hap = RNasesInHaplotypes_all_iters_Dict[iter_ind][mutant_rnase_ind][0]
                            # parent_hap_ = funcs.parent_hap_func(haplotypes_dict, haplotypes_dict_previous, mutant_hap, parents_haps, iter_ind, d_iters, n)
                            parent_hap = funcs.parent_hap_n_iters_before(haplotypes_dict_previous, mutant_hap, parents_haps)
                            # if parent_hap_ != parent_hap:
                            #     b=6
                            mut_gro = np.nan
                            mut_ext_gro = np.nan
                            if mutant_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind]:
                                mut_gro = reorder_mis_hap_strict_groups_all_iters_dict[iter_ind][mutant_hap]
                            if mutant_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind]:
                                mut_ext_gro = reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind][mutant_hap]
                            par_gro = np.nan
                            par_ext_gro = np.nan
                            if parent_hap in haplotypes_dict:
                                parent_RNase = haplotypes_dict[parent_hap].RNasesIndices[0]
                                if parent_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind]:
                                    par_gro = reorder_mis_hap_strict_groups_all_iters_dict[iter_ind][parent_hap]
                                if parent_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind]:
                                    par_ext_gro = reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind][parent_hap]
                            elif parent_hap in haplotypes_dict_previous:
                                parent_RNase = haplotypes_dict_previous[parent_hap].RNasesIndices[0]
                                if parent_hap in reorder_mis_hap_strict_groups_all_iters_dict[iter_ind-d_iters]:
                                    par_gro = reorder_mis_hap_strict_groups_all_iters_dict[iter_ind-d_iters][parent_hap]
                                if parent_hap in reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind-d_iters]:
                                    par_ext_gro = reorder_mis_hap_strict_group_extension_all_iters_dict[iter_ind-d_iters][parent_hap]

                            tmp_dict = {}
                            tmp_dict['orig group'] = par_gro
                            tmp_dict['orig ext group'] = par_ext_gro
                            tmp_dict['orig RNase'] = parent_RNase
                            tmp_dict['mutant RNase'] = mutant_rnase_ind
                            tmp_dict['orig hap'] = parent_hap
                            tmp_dict['mutant hap'] = mutant_hap
                            tmp_dict['mutant group'] = mut_gro
                            tmp_dict['mutant ext group'] = mut_ext_gro
                            tmp_dict['mutant hap birth'] = haplotypes_dict[mutant_hap].BirthIter
                            if iter_ind not in mut_rnases_dict:
                                mut_rnases_dict[iter_ind] = []
                            mut_rnases_dict[iter_ind].append(tmp_dict)

                max_group_ind = np.max([max_group_ind, np.max(list(reorder_mis_strict_groups_hap_all_iters_dict[iter_ind].keys()))])

                #----------------------------------------------------------
                #  save the groupings
                #----------------------------------------------------------
                if iter_ind % 1000 == 0:
                    with open(output_dir + '/mut_rnases_' + str(min_duration*10) + '_iters_dict.pkl', 'wb') as output:
                        pickle.dump(mut_rnases_dict, output, pickle.HIGHEST_PROTOCOL)

                if iter_ind % 1000 == 0:

                    with open(output_dir + '/reorder_mis_strict_groups_hap_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(reorder_mis_strict_groups_hap_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
                    with open(output_dir + '/reorder_mis_strict_groups_rnase_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(reorder_mis_strict_groups_rnase_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
                    with open(output_dir + '/reorder_mis_hap_strict_groups_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(reorder_mis_hap_strict_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
                    with open(output_dir + '/reorder_mis_rnase_strict_groups_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(reorder_mis_rnase_strict_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

                    with open(output_dir + '/reorder_mis_hap_strict_group_extension_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(reorder_mis_hap_strict_group_extension_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
                    with open(output_dir + '/reorder_mis_strict_group_extension_hap_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(reorder_mis_strict_group_extension_hap_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

                    with open(output_dir + '/group_naming_current_to_new_all_iters_dict.pkl', 'wb') as output:
                        pickle.dump(group_naming_current_to_new_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

                    with open(output_dir + '/new_any_rnases_lst.pkl', 'wb') as output:
                        pickle.dump(new_any_rnases_lst, output, pickle.HIGHEST_PROTOCOL)
                    with open(output_dir + '/t_new_any_RNase.pkl', 'wb') as output:
                        pickle.dump(t_new_any_RNase, output, pickle.HIGHEST_PROTOCOL)

            b=6
            with open(output_dir + '/split_one_group_into_two_all_iters_' + str(min_duration*10) + '_iters_dict.pkl', 'wb') as output:
                pickle.dump(split_one_group_into_two_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/extinct_group_all_iters_' + str(min_duration*10) + '_iters_dict.pkl', 'wb') as output:
                pickle.dump(extinct_group_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/mut_rnases_' + str(min_duration*10) + '_iters_dict.pkl', 'wb') as output:
                pickle.dump(mut_rnases_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/new_any_rnases_' + str(min_duration*10) + '_iters_lst.pkl', 'wb') as output:
                pickle.dump(new_any_rnases_lst, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/t_new_any_RNase' + str(min_duration*10) + '_iters_.pkl', 'wb') as output:
                pickle.dump(t_new_any_RNase, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/reorder_mis_strict_groups_hap_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(reorder_mis_strict_groups_hap_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/reorder_mis_strict_groups_rnase_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(reorder_mis_strict_groups_rnase_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/reorder_mis_hap_strict_groups_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(reorder_mis_hap_strict_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/reorder_mis_rnase_strict_groups_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(reorder_mis_rnase_strict_groups_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(output_dir + '/reorder_mis_hap_strict_group_extension_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(reorder_mis_hap_strict_group_extension_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
            with open(output_dir + '/reorder_mis_strict_group_extension_hap_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(reorder_mis_strict_group_extension_hap_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(output_dir + '/group_naming_current_to_new_all_iters_dict.pkl', 'wb') as output:
                pickle.dump(group_naming_current_to_new_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)

# except:
        #     continue