import simulationVars as Model
import simulationLogic as Logic
import calc_weighted_num_partners_of_allele as weighted_num_partners
import mutation_module as Mutation
import allele_duplication_module as Duplicate
import allele_deletion_module as Delete
# import genes_evolution_analysis as Analysis
import pickle
import os
import UNIPROT_calcs
# import time
import numpy as np

def run_simulation(iterationIndex):
    # bioclass_prob_arr = mutate_bioClass.calc_mutation_probs_bio_classes()
    bioclass_uniprot_freqs = UNIPROT_calcs.calc_mutation_probs_bio_classes()
    # print(os.getcwd())
    os.chdir("..")
    # print(os.getcwd())
    os.chdir("FlowersEvolutionsSimulation_ver13_data")
    # print(os.getcwd())
    directory = "output_data_10haplotypes_k5_18AAs_E_thresh-4_p10-4_rerun"
    subdirectories = ['2022-11-12-15-57-45-11'] #['2022-11-12-15-57-45-46']
    subdirectories = list(set(next(os.walk(directory))[1]))
    subdirectory = subdirectories[iterationIndex] #[iterationIndex]
    output_dir = directory + "/" + subdirectory

    with open(path + '/analyze_split_rnase_cros_slfs_info_100_iters.pkl',
              'rb') as input:
        analyze_split_rnase_cros_slfs_info = pickle.load(input)
    split_path = 'second' #'third' #'first'
    if split_path == 'first':
        # first kind of split path (RNase  ->SLFmut -> SLForig)
        splt_iter = 104530
        gro = analyze_split_rnase_cros_slfs_info[splt_iter][0]['orig_gro']
        mut_rnase = analyze_split_rnase_cros_slfs_info[splt_iter][0]['mut_rnase']  # extincted
        orig_rnase = analyze_split_rnase_cros_slfs_info[splt_iter][0]['orig_rnase']
        slf_mut_iter = int(
            np.ceil(analyze_split_rnase_cros_slfs_info[splt_iter][0]['iter_cros_mut_slf_mut_gro'] / 10) * 10) + 10
        last_iter = slf_mut_iter
        rnase_extincted = orig_rnase
        rnase_extincts = mut_rnase
    if split_path == 'second':
        #second kind of split path (RNase -> SLForig ->SLFmut)
        splt_iter = 103210 #93800
        gro = analyze_split_rnase_cros_slfs_info[splt_iter][0]['orig_gro']
        orig_rnase = analyze_split_rnase_cros_slfs_info[splt_iter][0]['orig_rnase']
        mut_rnase = analyze_split_rnase_cros_slfs_info[splt_iter][0]['mut_rnase']  # extincted
        rnase_extincted = mut_rnase
        rnase_extincts = orig_rnase  # extincts
        slf_orig_iter = int(
            np.ceil(analyze_split_rnase_cros_slfs_info[splt_iter][0]['iter_cros_mut_slf_orig_gro'] / 10) * 10) + 10
        last_iter = slf_orig_iter
    if split_path == 'third':
        # third kind of split path (SLForig -> RNase -> SLFmut)
        splt_iter = 101670  # 74600 #58280
        gro = analyze_split_rnase_cros_slfs_info[splt_iter][0]['orig_gro']
        mut_rnase = analyze_split_rnase_cros_slfs_info[splt_iter][0]['mut_rnase']  # extincted
        orig_rnase = analyze_split_rnase_cros_slfs_info[splt_iter][0]['orig_rnase']
        rnase_mut_iter = analyze_split_rnase_cros_slfs_info[splt_iter][0]['rnase_mut_iter']
        last_iter = rnase_mut_iter  # 46250  # 103170 #104400 #101540
        rnase_extincted = mut_rnase
        rnase_extincts = orig_rnase


    # output_dir = directory + "/" + subdirectory + '_splt_gro_86_iter_101670_' + str(iterationIndex)

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
    # '''
    if first_iter == last_iter:
        b=6
        more_than_one_original_ancestor = True
    else:
        b=6
        more_than_one_original_ancestor = False

    dir_rep = '_rerun'
    path2 = output_dir + dir_rep #os.path.join(path, dir_rep)
    os.mkdir(path2)
    print(path2)
    parents_haps = {}
    # more_than_one_original_ancestor = False

    # not_able_to_load = True
    # while not_able_to_load:
    #     try:
    #         open(output_dir + '/haplotypes_dict_' + str(last_iter) + '.pkl', 'rb')
    #         last_iter += 100
    #     except:
    #         not_able_to_load = False
    #         last_iter -= 100

    #last_iter = num_generations
    if more_than_one_original_ancestor:
        with open(output_dir + '/haplotypes_dict_' + str(last_iter) + '.pkl', 'rb') as input:
            haplotypes_dict = pickle.load(input)
        with open(output_dir + '/AA_slfs_dict.pkl', 'rb') as input:
            AA_slfs_dict = pickle.load(input)
        with open(output_dir + '/AA_rnases_dict.pkl', 'rb') as input:
            AA_rnases_dict = pickle.load(input)
        with open(output_dir + '/E_Ri_Fj_dict.pkl', 'rb') as input:
                E_Ri_Fj_dict = pickle.load(input)
        with open(output_dir + '/diploid_count_dict_' + str(last_iter) + '.pkl', 'rb') as input:
            diploid_count_dict = pickle.load(input)
        with open(output_dir + '/haplotypes_count_dict_' + str(last_iter) + '.pkl', 'rb') as input:
            haplotypes_count_dict = pickle.load(input)

        b=6
        # seed_strftime = int(month_int * 1e8 + day_int * 1e6 + hour_int * 1e5 + minute_int * 1e2 + second_int)
        # if Model.use_seed == 1:
        #     random.seed(seed_strftime)

        last_iter = last_iter
        iter_ind = last_iter
        freq_SC_haplotypes = 0

        delta_g_g_dict = {}
        delta_g_r_dict = {}
        cond1_g_g_dict = {}
        cond2_g_g_dict = {}

        d_iters = 500
        first_iter_of_one_ancestor = 0
        freq_SC_haps = 0
        while [more_than_one_original_ancestor or iter_ind < Model.num_of_iters + first_iter_of_one_ancestor] \
                and freq_SC_haps < Model.SC_freq_thresh:

            # -------------------------------------------------------
            #                    duplicate alleles
            # -------------------------------------------------------

            haplotypes_dict, diploid_count_dict, haplotypes_count_dict \
                = Duplicate.duplicate_population(iter_ind, haplotypes_dict, diploid_count_dict, haplotypes_count_dict)

            # -------------------------------------------------------
            #                    delete alleles
            # -------------------------------------------------------

            haplotypes_dict, diploid_count_dict, haplotypes_count_dict \
                = Delete.delete_population(iter_ind, haplotypes_dict, diploid_count_dict,
                                           haplotypes_count_dict)

            # -------------------------------------------------------
            #                    mutate population
            # -------------------------------------------------------

            haplotypes_dict, diploid_count_dict, haplotypes_count_dict, \
            AA_rnases_dict, AA_slfs_dict = \
                Mutation.mutate_population(iter_ind, bioclass_uniprot_freqs, haplotypes_dict,
                                           diploid_count_dict, haplotypes_count_dict, AA_rnases_dict, AA_slfs_dict)

            # -----------------------------------------------------------------
            #                  update interactions
            # -----------------------------------------------------------------

            E_Ri_Fj_dict = Logic.Update_AA_Interaction_dict(E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict,
                                                            haplotypes_count_dict, haplotypes_dict)

            # -------------------------------------------------------
            #   calculate diploids frequencies in the next generation
            # -------------------------------------------------------

            offspring_freq_arr, offspring_dict, delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict \
                = Logic.calc_freq_diploid_next_gen \
                (E_Ri_Fj_dict, diploid_count_dict, \
                 haplotypes_count_dict, haplotypes_dict, \
                 delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict)


            # -------------------------------------------------------
            #                  build next generation
            # -------------------------------------------------------
            # t = time.time()
            diploid_count_dict, haplotypes_count_dict = \
                Logic.build_the_next_generation(offspring_freq_arr, offspring_dict)

            # print("next generation" + str(time.time() - t))

            # --------------------------------------------------------------------------------
            #    update haplotypes_dict
            # --------------------------------------------------------------------------------

            wanted_keys = haplotypes_count_dict.keys()
            bigdict = haplotypes_dict
            haplotypes_dict = dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)

            # --------------------------------------------------------------------------------
            #    update delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict
            # --------------------------------------------------------------------------------

            wanted_keys = diploid_count_dict.keys()
            bigdict = delta_g_g_dict
            delta_g_g_dict = dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)
            bigdict = cond1_g_g_dict
            cond1_g_g_dict = dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)
            bigdict = cond2_g_g_dict
            cond2_g_g_dict = dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)

            g_r = {}
            for dip_keys in diploid_count_dict.keys():
                for hap_keys in haplotypes_count_dict:
                    g_r_key = tuple([item for sublist in [dip_keys, [hap_keys]] for item in sublist])
                    g_r[g_r_key] = []

            wanted_keys = g_r.keys()
            bigdict = delta_g_r_dict
            delta_g_r_dict = dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)
            b = 6

            # ------------------------------------------------------
            #      amount of SC haplotypes
            # ------------------------------------------------------
            num_SC_haps = 0
            for hap_ind, hap in haplotypes_dict.items():
                rnase_ind = hap.RNasesIndices[0]
                rnase_AA = AA_rnases_dict[rnase_ind]
                comp_cond1 = False
                for slf_ind in hap.SLFsIndices:
                    slf_AA = AA_slfs_dict[slf_ind]
                    if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                        E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                    else:
                        E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                    if E < Model.interaction_thresh:
                        # interacting
                        comp_cond1 = True
                        break
                if comp_cond1:
                    hap.IS_SC = 'True'
                    num_SC_haps += haplotypes_count_dict[hap.Haplotype_index_in_pool]

            freq_SC_haps = num_SC_haps / Model.popu_size / 2

            # -------------------------------------------------------------------------
            #     count the number of original ancestor RNases left in the population
            # -------------------------------------------------------------------------
            ancestor_rnases_arr = np.zeros(len(haplotypes_dict), dtype=int)
            for hap_ii, [hap_ind, hap] in enumerate(haplotypes_dict.items()):
                ancestor_rnases_arr[hap_ii] = hap.RNasesAncestralIndices[0]
            if len(np.unique(ancestor_rnases_arr)) == 1:
                more_than_one_original_ancestor = False

            # ------------------------------------------------------------------------
            #   group haplotypes use Maximun Independent Set algorithm
            # ------------------------------------------------------------------------
            if not more_than_one_original_ancestor:
                first_iter_of_one_ancestor = iter_ind
                d_iters = 10
                for current_hap in haplotypes_dict:
                    parents_haps[current_hap] = haplotypes_dict[current_hap].PreviousHaplotype
                if iter_ind % d_iters == 0:
                    with open(path2 + '/parents_haps.pkl', 'wb') as output:
                        pickle.dump(parents_haps, output, pickle.HIGHEST_PROTOCOL)

            # ---------------------------------------------------
            #   write haplotype/Diploid freqs/iters dicts
            # ---------------------------------------------------

            if not more_than_one_original_ancestor:
                first_iter_of_one_ancestor = iter_ind
                d_iters = 10
                for current_hap in haplotypes_dict:
                    parents_haps[current_hap] = haplotypes_dict[current_hap].PreviousHaplotype
                if iter_ind % d_iters == 0:
                    with open(path2 + '/parents_haps.pkl', 'wb') as output:
                        pickle.dump(parents_haps, output, pickle.HIGHEST_PROTOCOL)

            if iter_ind % d_iters == 0:
                print("start generation # %d" % iter_ind)

            if iter_ind % d_iters == 0:
                # print("start generation # %d" % iter_ind)
                Logic.save_files(path2, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

                with open(path2 + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                    pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

                with open(path2 + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                    pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

                with open(path2 + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                    pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)



            iter_ind += 1

        # for hap in haplotypes_iters_dict.keys():
        #     unique_SLFs = np.unique(haplotypes_obj[hap].SLFsIndices)
        #     unique_RNases = np.unique(haplotypes_obj[hap].RNasesIndices)
        #     E_Ri_Fj_dict = {}
        #     for rnase_ind in unique_RNases:
        #         rnase_val = AA_rnases_dict[rnase_ind]
        #         for slf_ind in unique_SLFs:
        #             slf_val = AA_slfs_dict[slf_ind]
        #             E_Ri_Fj_dict[rnase_ind, slf_ind] = np.sum(Model.eij[rnase_val, slf_val], axis=0)
        #
        #     SC = Analyze.is_compatible(E_Ri_Fj_dict, unique_SLFs, unique_RNases)
        #     haplotypes_obj[hap].IS_SC = SC

        Logic.save_files(path2, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

        with open(path2 + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
            pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

        with open(path2 + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
            pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

        with open(path2 + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
            pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)