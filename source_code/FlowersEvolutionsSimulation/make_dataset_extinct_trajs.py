import numpy as np
import pickle
import os
import simulationVars as Model
import calc_weighted_num_partners_of_allele as weighted_num_partners
min_duration = 10
d_iters = 10

os.chdir("..")
interaction_thresh = -6
os.chdir("FlowersEvolutionsSimulation_ver13_data")
# directory = "output_data_10haplotypes_k5_18AAs_E_thresh" + str(interaction_thresh) + "_p10-4"  # output_data_10haplotypes_k5_" + str(L) + "AAs_E_thresh-6_p2_2_10-5_mut_2" #test2" #test" #
#
# subdirectories = ['2022-11-12-15-57-44-14'] #list(set(next(os.walk(directory))[1]) -set(['2022-11-12-15-57-45-24']) - set(['2022-11-12-15-57-45-36']))
# subdirectories = list(set(next(os.walk(directory))[1])-set(['2022-11-12-15-57-45-15']))
#
# subdirectories = list(set(next(os.walk(directory))[1]) - set(['2023-04-19-15-41-27-33']) - set(['2023-04-19-15-41-27-36'])                      -set(['2023-03-19-11-56-28-30'])
#     - set(['2023-03-19-11-56-28-25']))

# os.chdir("FlowersEvolutionsSimulation_ver13_data3") # ... data2 -> n_haps = 250, data3 -> n_haps = 750
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4"
subdirectories = list(set(next(os.walk(directory))[1]))

directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_alpha_08_delta_09"
subdirectories = list(set(next(os.walk(directory))[1]) - set(['2023-06-07-10-54-24-41']) \
      -set(['2023-06-07-10-54-24-13']))
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p2times10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_phalf_times10-4"
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_no_dup_no_del"
subdirectories = list(set(next(os.walk(directory))[1]))
subdirectories =  ['2023-06-25-15-01-01-11']
subdirectories = ['2023-06-25-15-01-01-7']
directory = "output_data_10haplotypes_k5_18AAs_E_thresh-6_p10-4_2000haps"
subdirectories = list(set(next(os.walk(directory))[1])
                      - set(['2023-06-28-10-19-50-3']) - set(['2023-06-28-10-19-51-22']))
comp_rnase_became_noncomp_rnase = 0
n_splt_ext = 0
n_cross_SLF_mut_ext = 0
n_neutral_ext = 0
n_cross_SLF_orig_ext = 0
n_non_classified_ext = 0
n_small_ext_groups = []
n_mut_RNases_neutral_dict = {}
n_mut_RNases_SLF_mut_cross_dict = {}
n_mut_RNases_SLF_orig_cross_dict = {}
n_mut_RNases_splt_ext_dict = {}

n_orig_gro_SLF_mut_cross_dict = {}
n_orig_gro_SLF_orig_cross_dict = {}
n_orig_gro_neutral_dict = {}
n_orig_gro_splt_ext_dict = {}
n_ext_group_for_RNase_mut_all_files_dict = {}
p_mut_RNases_SLF_mut_cross_lst = []
p_mut_RNases_SLF_orig_cross_lst = []
p_mut_RNases_neutral_lst = []
p_mut_RNases_splt_ext_lst = []
cum_t_n_gro_dict = {}
for n_gro in range(3,30):
    cum_t_n_gro_dict[n_gro] = 0

cum_t = 0
for ii_dir, subdirectory in enumerate(subdirectories):
    output_dir = directory + "/" + subdirectory
    print(os.path.join(subdirectory))
    if os.path.isfile(output_dir + '/analyze_extinct_group_all_iters_' + str(min_duration * 10) + '_dict.pkl') \
            and os.path.isfile(output_dir + '/mut_rnases_' + str(min_duration * 10) + '_iters_dict.pkl'):
        with open(output_dir + '/RNasesInHaplotypes_all_iters_Dict.pkl', 'rb') as input:
            RNasesInHaplotypes_all_iters_Dict = pickle.load(input)
        with open(output_dir + '/analyze_extinct_group_all_iters_' + str(min_duration * 10) + '_dict.pkl', 'rb') as input:
            analyze_extinct_group_all_iters_dict = pickle.load(input)
        with open(output_dir + '/analyze_split_rnase_cros_slfs_info_' + str(min_duration * 10) + '_iters.pkl', 'rb') as input:
            analyze_split_rnase_cros_slfs_info = pickle.load(input)
        with open(output_dir + '/AA_rnases_dict.pkl', 'rb') as input:
            AA_rnases_dict = pickle.load(input)
        with open(output_dir + '/AA_slfs_dict.pkl', 'rb') as input:
            AA_slfs_dict = pickle.load(input)
        with open(output_dir + '/reorder_mis_strict_groups_hap_all_iters_dict.pkl', 'rb') as input:
            reorder_mis_strict_groups_hap_all_iters_dict = pickle.load(input)
        with open(output_dir + '/reorder_mis_strict_group_extension_hap_all_iters_dict.pkl', 'rb') as input:
            reorder_mis_strict_group_extension_hap_all_iters_dict = pickle.load(input)
        with open(output_dir + '/reorder_mis_rnase_strict_groups_all_iters_dict.pkl', 'rb') as input:
            reorder_mis_rnase_strict_groups_all_iters_dict = pickle.load(input)
        with open(output_dir + '/extinct_group_all_iters_' + str(min_duration * 10) + '_iters_dict.pkl', 'rb') as input:
            extinct_group_all_iters_dict = pickle.load(input)
        with open(output_dir + '/small_ext_groups_all_iters_dict_' + str(min_duration * 10) + '_dict.pkl',
                  'rb') as input:
            small_ext_groups_all_iters_dict = pickle.load(input)

        g_ext_ind = []

        first_iter = list(reorder_mis_strict_groups_hap_all_iters_dict.keys())[0]
        for iter_ind in reorder_mis_strict_groups_hap_all_iters_dict.keys():
            if iter_ind in reorder_mis_strict_group_extension_hap_all_iters_dict:
                gro_arr = \
                    [item for sublist in
                    [reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind], \
                    reorder_mis_strict_groups_hap_all_iters_dict[iter_ind]] for item in
                    sublist]
            gro_arr = np.unique(gro_arr)
            n_gro = len(gro_arr)
            cum_t_n_gro_dict[n_gro] += 10

        n_ext_group_for_RNase_mut_dict = {}
        mut_RNases_lst = []
        ext_gro_not_comp_more_than_one_gro = 0
        analyze_ext_cros_slfs_info = {}
        comp_rnase_became_noncomp_rnase_dict = {}
        for iter_ext_ind in extinct_group_all_iters_dict.keys():
            if iter_ext_ind in analyze_extinct_group_all_iters_dict:
                ext_lst = analyze_extinct_group_all_iters_dict[iter_ext_ind]
                analyze_ext_cros_slfs_info[iter_ext_ind] = []
                if len(ext_lst) > 1:
                    ext_gro_lst = []
                    for ext_dict in ext_lst:
                        ext_gro_lst.append(ext_dict['extinct_group'])
                    if len(np.unique(ext_gro_lst)) != len(ext_gro_lst):
                        ext_gro_not_comp_more_than_one_gro += 1
                for ext_ii, ext_event_ind in enumerate(ext_lst):
                    orig_RNase = ext_event_ind['orig RNase']
                    mut_RNase = ext_event_ind['mut RNase']
                    orig_gro = ext_event_ind['RNase group']
                    ext_gro = ext_event_ind['extinct_group']
                    mut_iter = ext_event_ind['mut iter']


                    #count the number of extincts - mutant RNase at the moment of extinction
                    #--------------------------------------------------
                    #   count RNase_mut at the moment of extinction
                    #--------------------------------------------------

                    with open(output_dir + '/haplotypes_dict_' + str(iter_ext_ind- d_iters) + '.pkl', 'rb') as input:
                        haplotypes_dict = pickle.load(input)
                    with open(output_dir + '/haplotypes_count_dict_' + str(iter_ext_ind- d_iters) + '.pkl', 'rb') as input:
                        haplotypes_count_dict = pickle.load(input)

                    # haps in orig group at the moment of extinction
                    if orig_gro in reorder_mis_strict_groups_hap_all_iters_dict[iter_ext_ind-d_iters]:
                        haps1 = reorder_mis_strict_groups_hap_all_iters_dict[iter_ext_ind-d_iters][orig_gro]
                    else:
                        haps1 = []
                    if orig_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ext_ind-d_iters]:
                        haps2 = reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ext_ind-d_iters][orig_gro]
                    else:
                        haps2 = []

                    haps_orig_gro_iter_ext = \
                        [item for sublist in
                         [haps1, haps2] for item in sublist]

                    n_orig_gro_iter_ext = 0
                    hap_lst = haps_orig_gro_iter_ext
                    for hap in hap_lst:
                        n_orig_gro_iter_ext += haplotypes_count_dict[hap]

                    if ext_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ext_ind - d_iters]:
                        haps1 = reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ext_ind - d_iters][ext_gro]
                    else:
                        haps1 = []
                    if ext_gro in reorder_mis_strict_groups_hap_all_iters_dict[iter_ext_ind - d_iters]:
                        haps2 = reorder_mis_strict_groups_hap_all_iters_dict[iter_ext_ind - d_iters][ext_gro]
                    else:
                        haps2 = []

                    haps_ext_gro_iter_ext = \
                        [item for sublist in
                         [haps1, haps2] for item in
                         sublist]
                    haps_ext_gro = haps1
                    n_ext = 0
                    for hap_ii, hap in enumerate(haps_ext_gro_iter_ext):
                        n_ext += haplotypes_count_dict[hap]
                    n_ext1 = 0
                    for hap_ii, hap in enumerate(haps1):
                        n_ext1 += haplotypes_count_dict[hap]
                    rnase_in_orig_gro_iter_ext = []
                    for hap in haps_orig_gro_iter_ext:
                        rnase_ind = haplotypes_dict[hap].RNasesIndices[0]
                        rnase_in_orig_gro_iter_ext.append(rnase_ind)
                    rnase_in_orig_gro_iter_ext = np.unique(rnase_in_orig_gro_iter_ext)

                    #---------------------------------------------------------------
                    # count mut RNases - not compatible with the extincted group
                    #---------------------------------------------------------------

                    # count the number of mut RNases in the group that is responsinle for the group extinction
                    # those are the RNases that before mutation were compatible with the extincted group
                    # and after the mutation they have lost compatibility

                    RNase_mut_lst = []
                    RNase_orig_lst = []
                    for rnase_ind in rnase_in_orig_gro_iter_ext: #rnases in the group that cause extinction
                        rnase_val = AA_rnases_dict[rnase_ind]
                        n_comp_mut = 0  # initiate
                        is_rnase_comp = np.zeros(len(haps_ext_gro_iter_ext), dtype=bool)
                        for hap_ii, hap in enumerate(haps_ext_gro_iter_ext):  # SLFs from haps in extincted group at extinction iter
                            slf_lst = haplotypes_dict[hap].SLFsIndices  # exterminating
                            for slf_ind in slf_lst:
                                slf_val = AA_slfs_dict[slf_ind]
                                if np.sum(Model.eij[rnase_val, slf_val], axis=0) < interaction_thresh:
                                    is_rnase_comp[hap_ii] = True
                                    n_comp_mut += haplotypes_count_dict[hap]
                                    break
                        if n_ext == 0:
                            b=6
                        if all(is_rnase_comp) or n_comp_mut / n_ext >= 0.75:
                            RNase_orig_lst.append(rnase_ind)
                        elif all(~is_rnase_comp) or n_comp_mut / n_ext <= 0.25:
                            RNase_mut_lst.append(rnase_ind)
                        else:
                            b = 6

                    haps_rnase_mut = []
                    for rnase_ind in RNase_mut_lst:
                        haps_rnase_mut.append(
                            RNasesInHaplotypes_all_iters_Dict[iter_ext_ind - d_iters][rnase_ind])

                    haps_rnase_mut = \
                        [item for sublist in
                         haps_rnase_mut for item in sublist]

                    n_rnase_mut_iter_ext = 0
                    for hap_ii, hap in enumerate(haps_rnase_mut):
                        n_rnase_mut_iter_ext += haplotypes_count_dict[hap]

                    b=6

                    #-----------------------------------------------------------------------------------------
                    # detect split in the group that is responsible to the extinction in addition to extinction
                    #-----------------------------------------------------------------------------------------
                    is_splt = False
                    for iter_splt, splt_lst in analyze_split_rnase_cros_slfs_info.items():
                        if not is_splt:
                            if len(splt_lst) == 1:
                                if orig_RNase == splt_lst[0]['orig_rnase'] \
                                    and orig_gro == splt_lst[0]['orig_gro'] \
                                    and mut_RNase == splt_lst[0]['mut_rnase']:
                                    #extinct + split

                                    analyze_ext_cros_slfs_info[iter_ext_ind].append('ext+split')
                                    analyze_extinct_group_all_iters_dict[iter_ext_ind][ext_ii]['ext_type'] = 'ext+split'
                                    n_splt_ext += 1
                                    is_splt = True
                                    if iter_ext_ind < iter_splt:
                                        if mut_RNase not in n_ext_group_for_RNase_mut_dict:
                                            n_ext_group_for_RNase_mut_dict[mut_RNase] = 0
                                        n_ext_group_for_RNase_mut_dict[mut_RNase] +=1

                                    else:
                                        b=6
                                        break

                                    if n_rnase_mut_iter_ext not in n_mut_RNases_splt_ext_dict:
                                        n_mut_RNases_splt_ext_dict[n_rnase_mut_iter_ext] = 0
                                    n_mut_RNases_splt_ext_dict[n_rnase_mut_iter_ext] += 1
                                    if n_orig_gro_iter_ext not in n_orig_gro_splt_ext_dict:
                                        n_orig_gro_splt_ext_dict[n_orig_gro_iter_ext] = 0
                                    n_orig_gro_splt_ext_dict[n_orig_gro_iter_ext] += 1
                                    if n_orig_gro_iter_ext == 0:
                                        p_mut_RNases_splt_ext_lst.append(np.nan)
                                    else:
                                        p_mut_RNases_splt_ext_lst.append(n_rnase_mut_iter_ext/n_orig_gro_iter_ext)


                            else: #more than one split in an iter
                                b=6
                    b=6

                    if not is_splt:
                        # -----------------------------------------------------------------------------------------
                        # detect cross SLF_mut in the group that is responsible to the extinction in addition to extinction
                        # -----------------------------------------------------------------------------------------
                        interact_mut = False
                        rnge = range(mut_iter, iter_ext_ind, d_iters)
                        is_interact_mut_rnase_iters = np.zeros(len(rnge), dtype=bool) #len(rnge)
                        n_haps_mut_all_iters = np.zeros(len(rnge))
                        n_haps_orig_all_iters = np.zeros(len(rnge))
                        for iter_ii, iter_ind in enumerate(rnge):
                            # iter_ind = mut_iter
                            # iter_ii = 0
                            with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                                haplotypes_dict = pickle.load(input)
                            with open(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl',
                                      'rb') as input:
                                haplotypes_count_dict = pickle.load(input)

                            # haps in the group that cause the extinction in iter 'iter_ind'
                            if orig_gro in reorder_mis_strict_groups_hap_all_iters_dict[iter_ind]:
                                haps1 = reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][orig_gro]
                            else:
                                haps1 = []
                            if orig_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind]:
                                haps2 = reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][orig_gro]
                            else:
                                haps2 = []

                            haps_orig_gro = \
                                [item for sublist in
                                 [haps1, haps2] for item in sublist]

                            rnase_in_orig_gro = []
                            for hap in haps_orig_gro:
                                rnase_ind = haplotypes_dict[hap].RNasesIndices[0]
                                rnase_in_orig_gro.append(rnase_ind)
                            rnase_in_orig_gro = np.unique(rnase_in_orig_gro)

                            # haps in orig group in iter 'iter_ind'
                            if ext_gro in reorder_mis_strict_groups_hap_all_iters_dict[iter_ind]:
                                haps1 = reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][
                                    ext_gro]
                            else:
                                haps1 = []
                            if ext_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[
                                iter_ind]:
                                haps2 = \
                                    reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][
                                        ext_gro]
                            else:
                                haps2 = []

                            haps_ext_gro = \
                                [item for sublist in
                                 [haps1, haps2] for item in sublist]
                            haps_ext_gro = haps1
                            haps_ext_gro = \
                                [item for sublist in
                                [haps1, haps2] for item in sublist]
                            n_ext = 0
                            for hap_ii, hap in enumerate(haps_ext_gro):
                                n_ext += haplotypes_count_dict[hap]

                            n_ext1 = 0
                            for hap_ii, hap in enumerate(haps1):
                                n_ext1 += haplotypes_count_dict[hap]

                            RNase_mut_lst = []
                            RNase_orig_lst = []

                            for rnase_ind in rnase_in_orig_gro:  # rnases in the group that cause extinction
                                rnase_val = AA_rnases_dict[rnase_ind]
                                is_rnase_comp = np.zeros(len(haps_ext_gro), dtype=bool)
                                n_comp_mut = 0
                                for hap_ii, hap in enumerate(haps_ext_gro):  # SLFs from haps in extincted group
                                    slf_lst = haplotypes_dict[hap].SLFsIndices  # exterminating
                                    for slf_ind in slf_lst:
                                        slf_val = AA_slfs_dict[slf_ind]
                                        if np.sum(Model.eij[rnase_val, slf_val],
                                                  axis=0) < interaction_thresh:
                                            is_rnase_comp[hap_ii] = True
                                            n_comp_mut += haplotypes_count_dict[hap]
                                            break
                                if n_ext == 0:
                                    b=6
                                if all(is_rnase_comp) or n_comp_mut/n_ext >=0.75:
                                    RNase_orig_lst.append(rnase_ind)
                                elif all(~is_rnase_comp) or n_comp_mut/n_ext <=0.25:
                                    RNase_mut_lst.append(rnase_ind) #16388
                                else:
                                    b=6

                            haps_rnase_orig = []
                            for rnase_ind in RNase_orig_lst:
                                haps_rnase_orig.append(RNasesInHaplotypes_all_iters_Dict[iter_ind][rnase_ind])

                            haps_rnase_orig = \
                                [item for sublist in
                                 haps_rnase_orig for item in sublist]

                            n_haps = 0
                            for hap in haps_rnase_orig:
                                n_haps += haplotypes_count_dict[hap]

                            n_haps_orig_all_iters[iter_ii] = n_haps

                            haps_rnase_mut = []
                            for rnase_ind in RNase_mut_lst:
                                haps_rnase_mut.append(RNasesInHaplotypes_all_iters_Dict[iter_ind][rnase_ind])

                            haps_rnase_mut = \
                                [item for sublist in
                                 haps_rnase_mut for item in sublist]

                            n_haps = 0
                            for hap in haps_rnase_mut:
                                n_haps += haplotypes_count_dict[hap]

                            n_haps_mut_all_iters[iter_ii] = n_haps

                            if len(RNase_mut_lst) == 0:  # no threat on the ext group, unless it extinct because a different reason
                                if iter_ext_ind in analyze_ext_cros_slfs_info \
                                        or len(analyze_extinct_group_all_iters_dict[iter_ext_ind]) > 1:
                                    b = 6  # make sure that this is not the reason for extintion , and if so dont count it
                                    break
                                else:
                                    b = 6
                            # if len(RNase_orig_lst) == 0:
                            #     b = 6
                            # if len(RNase_orig_lst) == 0:
                            #     b = 6
                            # if mut_RNase in RNase_orig_lst and len(RNase_mut_lst) == 0:
                            #     b = 6
                            # if orig_RNase in RNase_mut_lst and len(RNase_orig_lst) == 0:
                            #     if iter_ext_ind not in comp_rnase_became_noncomp_rnase_dict:
                            #         comp_rnase_became_noncomp_rnase_dict[iter_ext_ind] = []
                            #     comp_rnase_became_noncomp_rnase_dict[iter_ext_ind].append(analyze_extinct_group_all_iters_dict[iter_ext_ind][ext_ii]['extinct_group'])
                            haps_rnase_mut = []
                            for rnase_ind in RNase_mut_lst:
                                #rnase_ind = RNase_mut_lst
                                haps_rnase_mut.append(RNasesInHaplotypes_all_iters_Dict[iter_ind][rnase_ind])

                            haps_rnase_mut = \
                                [item for sublist in
                                 haps_rnase_mut for item in sublist]
                            is_interact_mut_rnase = np.zeros(len(RNase_orig_lst), dtype=bool)
                            for rnase_ii, rnase_ind in enumerate(RNase_orig_lst):
                                # rnase_ind = RNase_orig_lst
                                rnase_val = AA_rnases_dict[rnase_ind]
                                is_interact_mut = np.zeros(len(haps_rnase_mut), dtype=bool)
                                for hap_ii, hap in enumerate(haps_rnase_mut):  # male haps with mut RNase
                                        slf_lst = haplotypes_dict[hap].SLFsIndices  # exterminating
                                        for slf_ind in slf_lst:
                                            slf_val = AA_slfs_dict[slf_ind]
                                            if np.sum(Model.eij[rnase_val, slf_val], axis=0) < interaction_thresh:
                                                # interacting SLF
                                                # extinct + cross mut - mutation of SLF of RNase mut
                                                if haplotypes_dict[hap].RNasesIndices[0] != ext_event_ind[
                                                    'orig RNase']:  # not SC
                                                    is_interact_mut[hap_ii] = True
                                                    break
                                                else:
                                                    b = 6
                                b=6
                                # at least one of the haplotypes is cross compatible
                                is_interact_mut_rnase[rnase_ii] = any(is_interact_mut)
                            b=6
                            # all the RNases are compatible with at least one of the haplotypes
                            is_interact_mut_rnase_iters[iter_ii] = all(is_interact_mut_rnase)

                        # end of for loop
                        interact_mut = \
                            len(is_interact_mut_rnase_iters[is_interact_mut_rnase_iters ==True])/\
                            len(is_interact_mut_rnase_iters) >0.8 or \
                            all(is_interact_mut_rnase_iters[int(len(is_interact_mut_rnase_iters) - 1 \
                            - np.floor(len(is_interact_mut_rnase_iters)*0.2)):]) \
                            == True # 80% of the interactions are cross mut RNases or last 20% of the interactions are cross mut

                        if interact_mut:
                            if mut_RNase not in n_ext_group_for_RNase_mut_dict:
                                n_ext_group_for_RNase_mut_dict[mut_RNase] = 0
                            n_ext_group_for_RNase_mut_dict[mut_RNase] +=1
                            analyze_extinct_group_all_iters_dict[iter_ext_ind][ext_ii]['ext_type'] = 'SLF mut cross'
                            analyze_ext_cros_slfs_info[iter_ext_ind].append('SLF mut cross')
                            n_cross_SLF_mut_ext += 1
                            if n_rnase_mut_iter_ext not in n_mut_RNases_SLF_mut_cross_dict:
                                n_mut_RNases_SLF_mut_cross_dict[n_rnase_mut_iter_ext] = 0
                            n_mut_RNases_SLF_mut_cross_dict[n_rnase_mut_iter_ext] += 1
                            if n_orig_gro_iter_ext not in n_orig_gro_SLF_mut_cross_dict:
                                n_orig_gro_SLF_mut_cross_dict[n_orig_gro_iter_ext] = 0
                            n_orig_gro_SLF_mut_cross_dict[n_orig_gro_iter_ext] += 1
                            if n_orig_gro_iter_ext == 0:
                                p_mut_RNases_SLF_mut_cross_lst.append(np.nan)
                            else:
                                p_mut_RNases_SLF_mut_cross_lst.append(n_rnase_mut_iter_ext / n_orig_gro_iter_ext)

                        if not interact_mut:
                            #if it is not a cross SLFmut before the extinction
                            # ----------------------------------------------------------------------------------------------
                            # detect cross SLF_orig in the group that causes extinction in addition to extinction process
                            # ----------------------------------------------------------------------------------------------
                            interact_orig = False
                            rnge = range(mut_iter, iter_ext_ind, d_iters)
                            is_interact_orig_rnase_iters = np.zeros(len(rnge), dtype=bool) #len(rnge)
                            n_haps_mut_all_iters = np.zeros(len(rnge)) #len(rnge)
                            n_haps_orig_all_iters = np.zeros(len(rnge)) #
                            for iter_ii, iter_ind in enumerate(rnge):
                            # iter_ind = mut_iter
                            # iter_ii = 0
                            # for iter_ind in rnge:
                                with open(output_dir + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'rb') as input:
                                    haplotypes_dict = pickle.load(input)
                                with open(output_dir + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl',
                                          'rb') as input:
                                    haplotypes_count_dict = pickle.load(input)

                                # haps in the group that cause the extinction in iter 'iter_ind'
                                if orig_gro in reorder_mis_strict_groups_hap_all_iters_dict[iter_ind]:
                                    haps1 = reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][orig_gro]
                                else:
                                    haps1 = []
                                if orig_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind]:
                                    haps2 = reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][orig_gro]
                                else:
                                    haps2 = []

                                haps_orig_gro = \
                                    [item for sublist in
                                     [haps1, haps2] for item in sublist]

                                rnase_in_orig_gro = []
                                for hap in haps_orig_gro:
                                    rnase_ind = haplotypes_dict[hap].RNasesIndices[0]
                                    rnase_in_orig_gro.append(rnase_ind)
                                rnase_in_orig_gro = np.unique(rnase_in_orig_gro)

                                # haps in orig group in iter 'iter_ind'
                                if ext_gro in reorder_mis_strict_groups_hap_all_iters_dict[iter_ind]:
                                    haps1 = reorder_mis_strict_groups_hap_all_iters_dict[iter_ind][
                                        ext_gro]
                                else:
                                    haps1 = []
                                if ext_gro in reorder_mis_strict_group_extension_hap_all_iters_dict[
                                    iter_ind]:
                                    haps2 = \
                                    reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ind][
                                        ext_gro]
                                else:
                                    haps2 = []

                                haps_ext_gro = \
                                    [item for sublist in
                                     [haps1, haps2] for item in sublist]
                                # if len(haps1) > 0:
                                #     haps_ext_gro = haps1
                                n_ext = 0
                                for hap_ii, hap in enumerate(haps_ext_gro):
                                    n_ext += haplotypes_count_dict[hap]
                                n_ext1 = 0
                                for hap_ii, hap in enumerate(haps1):
                                    n_ext1 += haplotypes_count_dict[hap]

                                RNase_mut_lst = []
                                RNase_orig_lst = []
                                for rnase_ind in rnase_in_orig_gro:  # rnases in the group that cause extinction
                                    rnase_val = AA_rnases_dict[rnase_ind]
                                    is_rnase_comp = np.zeros(len(haps_ext_gro), dtype=bool)
                                    n_comp_mut = 0
                                    for hap_ii, hap in enumerate(haps_ext_gro):  # SLFs from haps in extincted group
                                        slf_lst = haplotypes_dict[hap].SLFsIndices  # exterminating
                                        for slf_ind in slf_lst:
                                            slf_val = AA_slfs_dict[slf_ind]
                                            if np.sum(Model.eij[rnase_val, slf_val],
                                                      axis=0) < interaction_thresh:
                                                is_rnase_comp[hap_ii] = True
                                                n_comp_mut += haplotypes_count_dict[hap]
                                                break
                                    if n_ext == 0:
                                        b=6
                                    if all(is_rnase_comp) or n_comp_mut / n_ext >= 0.75:
                                        RNase_orig_lst.append(rnase_ind)
                                    elif all(~is_rnase_comp) or n_comp_mut / n_ext <= 0.25:
                                        RNase_mut_lst.append(rnase_ind)
                                    else:
                                        b = 6

                                if len(RNase_mut_lst) == 0: #no threat on the ext group, unless it extinct because a different reason
                                    if iter_ext_ind in analyze_ext_cros_slfs_info\
                                            or len(analyze_extinct_group_all_iters_dict[iter_ext_ind]) > 1:
                                        b=6# make sure that this is not the reason for extintion , and if so dont count it
                                        break
                                    else:
                                        b=6
                                if len(RNase_orig_lst) == 0:
                                    b=6
                                # if mut_RNase in RNase_orig_lst and len(RNase_mut_lst) == 0:
                                #     b=6
                                # if orig_RNase in RNase_mut_lst and len(RNase_orig_lst) == 0:
                                #     b=6
                                #     if iter_ext_ind not in comp_rnase_became_noncomp_rnase_dict:
                                #         comp_rnase_became_noncomp_rnase_dict[iter_ext_ind] = []
                                #     comp_rnase_became_noncomp_rnase_dict[iter_ext_ind].append(
                                #         analyze_extinct_group_all_iters_dict[iter_ext_ind][ext_ii]['extinct_group'])

                                haps_rnase_orig = []
                                for rnase_ind in RNase_orig_lst:
                                    haps_rnase_orig.append(RNasesInHaplotypes_all_iters_Dict[iter_ind][rnase_ind])

                                haps_rnase_orig = \
                                    [item for sublist in
                                     haps_rnase_orig for item in sublist]

                                n_haps = 0
                                for hap in haps_rnase_orig:
                                    n_haps += haplotypes_count_dict[hap]

                                n_haps_orig_all_iters[iter_ii] = n_haps

                                haps_rnase_mut = []
                                for rnase_ind in RNase_mut_lst:
                                    haps_rnase_mut.append(RNasesInHaplotypes_all_iters_Dict[iter_ind][rnase_ind])

                                haps_rnase_mut = \
                                    [item for sublist in
                                     haps_rnase_mut for item in sublist]

                                n_haps = 0
                                for hap in haps_rnase_mut:
                                    n_haps += haplotypes_count_dict[hap]

                                n_haps_mut_all_iters[iter_ii] = n_haps

                                is_interact_orig_rnase= np.zeros(len(RNase_mut_lst), dtype=bool)
                                for rnase_ii, rnase_ind in enumerate(RNase_mut_lst):
                                    rnase_val = AA_rnases_dict[rnase_ind]
                                    # rnase_ind = RNase_mut_lst
                                    is_interact_orig = np.zeros(len(haps_rnase_orig), dtype=bool)
                                    n_orig_rnase_cross = 0
                                    n_orig_rnase_non_cross = 0
                                    for hap_ii, hap in enumerate(haps_rnase_orig): #male haps with orig RNase
                                        slf_lst = haplotypes_dict[hap].SLFsIndices  # exterminating
                                        for slf_ind in slf_lst:
                                            slf_val = AA_slfs_dict[slf_ind]
                                            if np.sum(Model.eij[rnase_val, slf_val], axis=0) < interaction_thresh:
                                                # interacting SLF
                                                # extinct + cross orig - mutation of SLF of RNase orig
                                                if haplotypes_dict[hap].RNasesIndices[0] != ext_event_ind[
                                                    'mut RNase']:  # not SC
                                                    is_interact_orig[hap_ii] = True
                                                    n_orig_rnase_cross += haplotypes_count_dict[hap]
                                                    break
                                                else:
                                                    b = 6
                                        if not is_interact_orig[hap_ii]:
                                            n_orig_rnase_non_cross += haplotypes_count_dict[hap]
                                    b=6
                                    # if any(is_interact_orig):
                                    is_interact_orig_rnase[rnase_ii] = any(is_interact_orig) #at least one of the haplotypes with RNase orig
                                    # is compatible (as a male) with RNase mut
                                    is_interact_orig_rnase_iters[iter_ii] = all(is_interact_orig_rnase) # all the haplotypes in this iter with RNase orig
                                    # are compatible (as male) with RNasemut
                                    # in this iter SLForig was cross mut with RNasemut

                                b=6
                            # end of for loop
                            interact_orig = \
                                len(is_interact_orig_rnase_iters[is_interact_orig_rnase_iters == True]) / \
                                len(is_interact_orig_rnase_iters) > 0.8 or \
                                all(is_interact_orig_rnase_iters[int(len(is_interact_orig_rnase_iters) - 1 \
                                - np.floor(len(is_interact_orig_rnase_iters) * 0.2)):]) \
                                == True  # 80% of the interactions are cross orig RNases or last 20% of the interactions are cross mut

                            if interact_orig:
                                if mut_RNase not in n_ext_group_for_RNase_mut_dict:
                                    n_ext_group_for_RNase_mut_dict[mut_RNase] = 0
                                n_ext_group_for_RNase_mut_dict[mut_RNase] +=1
                                analyze_extinct_group_all_iters_dict[iter_ext_ind][ext_ii]['ext_type'] = 'SLF orig cross'
                                analyze_ext_cros_slfs_info[iter_ext_ind].append('SLF orig cross')
                                n_cross_SLF_orig_ext += 1
                                if n_rnase_mut_iter_ext not in n_mut_RNases_SLF_orig_cross_dict:
                                    n_mut_RNases_SLF_orig_cross_dict[n_rnase_mut_iter_ext] = 0
                                n_mut_RNases_SLF_orig_cross_dict[n_rnase_mut_iter_ext] += 1
                                if n_orig_gro_iter_ext not in n_orig_gro_SLF_orig_cross_dict:
                                    n_orig_gro_SLF_orig_cross_dict[n_orig_gro_iter_ext] = 0
                                n_orig_gro_SLF_orig_cross_dict[n_orig_gro_iter_ext] += 1
                                if n_orig_gro_iter_ext == 0:
                                    p_mut_RNases_SLF_orig_cross_lst.append(np.nan)
                                else:
                                    p_mut_RNases_SLF_orig_cross_lst.append(
                                        n_rnase_mut_iter_ext / n_orig_gro_iter_ext)

                            if not interact_orig:
                                # ----------------------------------------------------------------------------------------------
                                # detect neutral extinction
                                # ----------------------------------------------------------------------------------------------
                                if len(n_haps_mut_all_iters) == 0:
                                    b=6
                                if len(n_haps_orig_all_iters) == 0:
                                    b=6
                                if len(n_haps_mut_all_iters[n_haps_mut_all_iters != 0])/\
                                        len(n_haps_mut_all_iters)<0.75:
                                    if not(iter_ext_ind in analyze_ext_cros_slfs_info \
                                            or len(analyze_extinct_group_all_iters_dict[iter_ext_ind]) > 1):
                                        n_non_classified_ext += 1
                                    continue
                                # if len(n_haps_orig_all_iters[n_haps_orig_all_iters != 0]) /\
                                #     len(n_haps_orig_all_iters)!=1:
                                #     b=6
                                # if not (mut_RNase in RNase_orig_lst or orig_RNase in RNase_mut_lst):
                                analyze_extinct_group_all_iters_dict[iter_ext_ind][ext_ii]['ext_type'] = 'neutral'
                                analyze_ext_cros_slfs_info[iter_ext_ind].append('neutral')
                                n_neutral_ext += 1
                                if mut_RNase not in n_ext_group_for_RNase_mut_dict:
                                    n_ext_group_for_RNase_mut_dict[mut_RNase] = 0
                                n_ext_group_for_RNase_mut_dict[mut_RNase] +=1

                                if n_rnase_mut_iter_ext not in n_mut_RNases_neutral_dict:
                                    n_mut_RNases_neutral_dict[n_rnase_mut_iter_ext] = 0
                                n_mut_RNases_neutral_dict[n_rnase_mut_iter_ext] += 1
                                if n_orig_gro_iter_ext not in n_orig_gro_neutral_dict:
                                    n_orig_gro_neutral_dict[n_orig_gro_iter_ext] = 0
                                n_orig_gro_neutral_dict[n_orig_gro_iter_ext] += 1
                                if n_orig_gro_iter_ext == 0:
                                    p_mut_RNases_neutral_lst.append(np.nan)
                                else:
                                    p_mut_RNases_neutral_lst.append(
                                        n_rnase_mut_iter_ext / n_orig_gro_iter_ext)

                                # else:
                                #     wrong_rnase_classification += 1
                b=6
            else: #not as a result of non-comp RNase-mut
                if iter_ext_ind in small_ext_groups_all_iters_dict:
                    #calc the spontaneus extincted group size(its maximal size)
                    g_lst = small_ext_groups_all_iters_dict[iter_ext_ind]
                    for g in g_lst:
                        # check if at extinct group size is small.
                        # count only those, that at extinct iter group size is small
                        hps_lst = []
                        if g in reorder_mis_strict_groups_hap_all_iters_dict[iter_ext_ind-d_iters]:
                            hps_lst.append(reorder_mis_strict_groups_hap_all_iters_dict[iter_ext_ind-d_iters][g])
                        if g in reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ext_ind-d_iters]:
                            hps_lst.append(reorder_mis_strict_group_extension_hap_all_iters_dict[iter_ext_ind-d_iters][g])
                        hps_lst = [item for sublist in hps_lst for item in sublist]
                        n_hps = 0
                        with open(output_dir + '/haplotypes_count_dict_' + str(iter_ext_ind-d_iters) + '.pkl', 'rb') as input:
                            haplotypes_count_dict = pickle.load(input)
                        for h in hps_lst:
                            n_hps += haplotypes_count_dict[h]

                        # if n_hps > 10:
                        #     break

                        cond = True
                        itr = iter_ext_ind - d_iters
                        n_haps_lst = []
                        while cond:
                            hps_lst = []
                            if g in reorder_mis_strict_groups_hap_all_iters_dict[itr]:
                                hps_lst.append(reorder_mis_strict_groups_hap_all_iters_dict[itr][g])
                            if g in reorder_mis_strict_group_extension_hap_all_iters_dict[itr]:
                                hps_lst.append(reorder_mis_strict_group_extension_hap_all_iters_dict[itr][g])
                            hps_lst = [item for sublist in hps_lst for item in sublist]

                            n_hps = 0
                            with open(output_dir + '/haplotypes_count_dict_' + str(itr) + '.pkl', 'rb') as input:
                                haplotypes_count_dict = pickle.load(input)
                            for h in hps_lst:
                                n_hps += haplotypes_count_dict[h]
                            n_haps_lst.append(n_hps)
                            itr = itr - d_iters
                            if itr in reorder_mis_strict_groups_hap_all_iters_dict:
                                if g in reorder_mis_strict_groups_hap_all_iters_dict[itr] \
                                    or g in reorder_mis_strict_group_extension_hap_all_iters_dict[itr]:
                                    cond = True
                                else:
                                    cond = False
                            else:
                                cond = False
                        if n_haps_lst[0] < 20: #this is the last count of the extincted class (10 iters befor extinction)
                            if iter_ext_ind - first_iter > len(n_haps_lst)*d_iters:
                                n_small_ext_groups.append(np.max(n_haps_lst))
                        else:
                            b=6
        for n_ext_group in n_ext_group_for_RNase_mut_dict.values():
            if n_ext_group not in n_ext_group_for_RNase_mut_all_files_dict:
                n_ext_group_for_RNase_mut_all_files_dict[n_ext_group] = 0
            n_ext_group_for_RNase_mut_all_files_dict[n_ext_group] +=1

        with open(output_dir + '/analyze_ext_cros_slfs_info_' + str(min_duration * 10) + '_iters.pkl', 'wb') as output:
            pickle.dump(analyze_ext_cros_slfs_info, output, pickle.HIGHEST_PROTOCOL)
        with open(output_dir + '/analyze_extinct_group_all_iters_' + str(min_duration * 10) + '_dict.pkl', 'wb') as output:
            pickle.dump(analyze_extinct_group_all_iters_dict, output, pickle.HIGHEST_PROTOCOL)
        b=6
b=6

n_mut_RNases_neutral_lst = []
for n, num_of_times in n_mut_RNases_neutral_dict.items():
    for i in range(num_of_times):
        n_mut_RNases_neutral_lst.append(n)

n_mut_RNases_SLF_mut_cross_lst = []
for n, num_of_times in n_mut_RNases_SLF_mut_cross_dict.items():
    for i in range(num_of_times):
        n_mut_RNases_SLF_mut_cross_lst.append(n)

n_mut_RNases_SLF_orig_cross_lst = []
for n, num_of_times in n_mut_RNases_SLF_orig_cross_dict.items():
    for i in range(num_of_times):
        n_mut_RNases_SLF_orig_cross_lst.append(n)

n_mut_RNases_splt_ext_lst = []
for n, num_of_times in n_mut_RNases_splt_ext_dict.items():
    for i in range(num_of_times):
        n_mut_RNases_splt_ext_lst.append(n)


n_orig_gro_neutral_lst = []
for n, num_of_times in n_orig_gro_neutral_dict.items():
    for i in range(num_of_times):
        n_orig_gro_neutral_lst.append(n)

n_orig_gro_SLF_mut_cross_lst = []
for n, num_of_times in n_orig_gro_SLF_mut_cross_dict.items():
    for i in range(num_of_times):
        n_orig_gro_SLF_mut_cross_lst.append(n)

n_orig_gro_SLF_orig_cross_lst = []
for n, num_of_times in n_orig_gro_SLF_orig_cross_dict.items():
    for i in range(num_of_times):
        n_orig_gro_SLF_orig_cross_lst.append(n)

n_orig_gro_splt_ext_lst = []
for n, num_of_times in n_orig_gro_splt_ext_dict.items():
    for i in range(num_of_times):
        n_orig_gro_splt_ext_lst.append(n)

with open(directory + '/n_orig_gro_neutral_lst.pkl', 'wb') as output:
    pickle.dump(n_orig_gro_neutral_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_orig_gro_splt_ext_lst.pkl', 'wb') as output:
    pickle.dump(n_orig_gro_splt_ext_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_orig_gro_SLF_orig_cross_lst.pkl', 'wb') as output:
    pickle.dump(n_orig_gro_SLF_orig_cross_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_orig_gro_SLF_mut_cross_lst.pkl', 'wb') as output:
    pickle.dump(n_orig_gro_SLF_mut_cross_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_ext_group_for_RNase_mut_all_files_dict.pkl', 'wb') as output:
    pickle.dump(n_ext_group_for_RNase_mut_all_files_dict, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_cross_SLF_mut_ext.pkl', 'wb') as output:
    pickle.dump(n_cross_SLF_mut_ext, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_neutral_ext.pkl', 'wb') as output:
    pickle.dump(n_neutral_ext, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_splt_ext.pkl', 'wb') as output:
    pickle.dump(n_splt_ext, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_cross_SLF_orig_ext.pkl', 'wb') as output:
    pickle.dump(n_cross_SLF_orig_ext, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_non_classified_ext.pkl', 'wb') as output:
    pickle.dump(n_non_classified_ext, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_small_ext_groups.pkl', 'wb') as output:
    pickle.dump(n_small_ext_groups, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_mut_RNases_neutral_dict.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_neutral_dict, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/n_mut_RNases_SLF_mut_cross_dict.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_SLF_mut_cross_dict, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/n_mut_RNases_SLF_orig_cross_dict.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_SLF_orig_cross_dict, output, pickle.HIGHEST_PROTOCOL)
with open(directory + '/n_mut_RNases_splt_ext_dict.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_splt_ext_dict, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_mut_RNases_neutral_lst.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_neutral_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_mut_RNases_SLF_mut_cross_lst.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_SLF_mut_cross_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_mut_RNases_SLF_orig_cross_lst.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_SLF_orig_cross_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/n_mut_RNases_splt_ext_lst.pkl', 'wb') as output:
    pickle.dump(n_mut_RNases_splt_ext_lst, output, pickle.HIGHEST_PROTOCOL)

with open(directory + '/cum_t_n_gro_dict_ext.pkl', 'wb') as output:
    pickle.dump(cum_t_n_gro_dict, output, pickle.HIGHEST_PROTOCOL)