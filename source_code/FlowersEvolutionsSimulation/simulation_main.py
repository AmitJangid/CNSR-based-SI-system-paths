#import sys
#sys.modules[__name__].__dict__.clear()

import simulationVars as Model
import simulationLogic as Logic
import mutation_module as Mutation
import allele_duplication_module as Duplicate
import allele_deletion_module as Delete
import copy
import analyze_mutations_supportive_functions as Analyze
import sys

# import BLOSUM_calcs
import UNIPROT_calcs

import numpy as np
import pickle
from decimal import Decimal
import time
import os
from numpy import random
from typing import Dict, Any, Tuple, Union

import cProfile, pstats #, StringIO #io

# pr = cProfile.Profile()
# pr.enable()

def run_simulation(iterationIndex):

    #bioclass_prob_arr = mutate_bioClass.calc_mutation_probs_bio_classes()
    # bioclass_blosum_freqs = BLOSUM_calcs.calc_mutation_probs_bio_classes()
    bioclass_uniprot_freqs = UNIPROT_calcs.calc_mutation_probs_bio_classes()
    year = time.strftime("%Y")
    year_int = int(year)
    month = time.strftime("%m")
    month_int = int(month)
    day = time.strftime("%d")
    day_int = int(day)
    hour = time.strftime("%H")
    hour_int = int(hour)
    minute = time.strftime("%M")
    minute_int = int(minute)
    second = time.strftime("%S")
    second_int = int(second)
    millisecond = round(time.monotonic() * 1000)
    timestr = year + "-" + month + "-" + day + "-" + hour + "-" +  minute + "-" + second + "-" + str(iterationIndex)
    seed_strftime = int(month_int * 1e8 + day_int * 1e6 + hour_int * 1e5 + minute_int * 1e2 + second_int)
    if Model.use_seed == 1:
        random.seed(seed_strftime)

    directory = timestr
    print(directory)
    parent_dir = "output_data_10haplotypes"
    # parent_dir = "output_data"
    # os.chdir("..")
    # os.chdir("FlowersEvolutionSimulation_ver9_data")
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    if Model.use_seed == 1:
        np.random.seed(None)
        st0 = np.random.get_state()
        with open(path + '/init_state_RNG.pkl', 'wb') as output:
            pickle.dump(st0, output, pickle.HIGHEST_PROTOCOL)

    log_file = open(path + '/log_file_'
                                + timestr
                                + '.txt', "a")

    log_file.write("Directory: " +directory + "\n"
                   "population size: " + str(Model.popu_size) + "\n" 
                    "num of iteration: " + str(Model.num_of_iters) + "\n"
                    "allele length: " + str(Model.seq_len) + "\n"
                    "initialized num of haplotypes: " + str(Model.num_types_of_haplotypes) + "\n"
                    "initialized num of SLFs in haplotype: " + str(Model.slf_amount_in_haploid) + "\n"
                    "num of AA bio-classes: " + str(Model.alphabet_size) + "\n"
                    "alpha: " + str(Model.alpha) + "\n"
                    "delta: " + str(Model.delta[0]) + "\n"
                    "interaction threshold: " + str(Model.interaction_thresh) + "\n"
                    "p mutation: " +"{:.2E}".format(Decimal(str(Model.p_mutation))) + "\n"
                    "p duplication: " + "{:.2E}".format(Decimal(str(Model.p_duplicate))) + "\n"
                    "p deletion: " + "{:.2E}".format(Decimal(str(Model.p_delete))) + "\n"
                    "5'th class for STOP CODON? - NO \n"  
                    "initialization of bio-classes frequencies: " + Model.initiate_bioclasses_freq + "\n" 
                    "initialized as SC?" + Model.SC_type + "\n"
                    "k(number of female fertilization trials - new model)" + str(Model.num_of_fertilizations_trials) + "\n"   
                    "initialization the population - " + Model.initiation_RNases_costrained_to_log_file + "\n"                                                                                                          
                    "" + Model.initiation_text_to_log_file + "\n" 
                    "RNases duplication is " + Model.enable_dup + "\n"
                    "RNases deletion is " + Model.enable_del + "\n"
                    "seed file in init_state_RNG.pkl, for using seed type: with open('/init_state_RNG.pkl', 'rb') as input:) \n"
                    "and - np.random.set_state(st0) \n"
                    "see: https://stackoverflow.com/questions/32172054/how-can-i-retrieve-the-current-seed-of-numpys-random-number-generator")

    log_file.close()

    #-------------------------------------------------------
    #   initiate mis dicts
    #-------------------------------------------------------
    '''
    mis_groups_hap_all_iters_dict = {}
    mis_groups_rnase_all_iters_dict = {}
    mis_rnase_groups_all_iters_dict = {}
    mis_hap_groups_all_iters_dict = {}
    mis_non_grouped_haps_all_iters_dict = {}
    mis_foreign_comp_groups_rnases_all_iters_dict = {}
    mis_foreign_comp_groups_haps_all_iters_dict = {}
    copy_mis_groups_hap_all_iters_dict = {}
    copy_mis_groups_rnase_all_iters_dict = {}
    copy_mis_rnase_groups_all_iters_dict = {}
    copy_mis_hap_groups_all_iters_dict = {}
    copy_mis_non_grouped_haps_all_iters_dict = {}
    copy_mis_foreign_comp_groups_rnases_all_iters_dict = {}
    copy_mis_foreign_comp_groups_haps_all_iters_dict = {}
    '''
    parents_haps = {}
    # -----------------------------------------------------
    #      Inintiate the population
    #------------------------------------------------------

    if Model.generate_new_population:
        haplotypes_init_dict, all_RNase_ar, all_SLF_ar, \
        SLF_ancestral_Dict, ancestral_SLF_Dict, \
        RNase_ancestral_Dict, ancestral_RNase_Dict\
            = Logic.initiate_flowers_new_alg() #initiate_flowers_no_selection_constrains() #

        E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict = Logic.init_AA_interactions_to_dict(all_RNase_ar, all_SLF_ar)


    with open(path + '/initial_unique_haplotypes_dict.pkl', 'wb') as output:
        pickle.dump(haplotypes_init_dict, output, pickle.HIGHEST_PROTOCOL)

    num_of_init_diploids = int(Model.num_types_of_haplotypes * (Model.num_types_of_haplotypes - 1) / 2)

    # if union probability for all the diploids
    p = 1/num_of_init_diploids * np.ones(num_of_init_diploids)
    # numRNases = len(AA_rnases_dict) #np.size(all_RNase_ar,0)
    # numSLFs = len(AA_slfs_dict) #np.size(all_SLF_ar, 0)
    diploid_count_dict, haplotypes_count_dict, haplotypes_dict, RNases_count_dict, SLFs_count_dict \
        = Logic.build_population_given_prob(haplotypes_init_dict, p)

    # save the initial population using pickle
    # with open(path + '/initial_diploids_obj.pkl', 'wb') as output:
    #     pickle.dump(diploids_obj, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/initial_haplotypes_dict.pkl', 'wb') as output:
        pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/initial_diploid_count_dict.pkl', 'wb') as output:
        pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/initial_haplotypes_count_dict.pkl', 'wb') as output:
        pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/AA_rnases_dict.pkl', 'wb') as output:
        pickle.dump(AA_rnases_dict, output, pickle.HIGHEST_PROTOCOL)
    with open(path + '/AA_slfs_dict.pkl', 'wb') as output:
        pickle.dump(AA_slfs_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/E_Ri_Fj_dict.pkl', 'wb') as output:
        pickle.dump(E_Ri_Fj_dict, output, pickle.HIGHEST_PROTOCOL)

    #------------------------------------------------------------
    #          initializations
    # -----------------------------------------------------------

    haplotypes_SC_count_dict: Dict[Tuple[Any], Union[int, Any]] = {}
    print("start generation # 0")
    iter_ind = 0
    freq_SC_haps = 0
    delta_g_g_dict = {}
    delta_g_r_dict = {}
    cond1_g_g_dict = {}
    cond2_g_g_dict = {}
    # R_g_g_r_dict = {}
    # R_g_h_r_dict = {}

    more_than_one_original_ancestor= True

    d_iters = 500
    first_iter_of_one_ancestor = 0
    while [more_than_one_original_ancestor or iter_ind < Model.num_of_iters + first_iter_of_one_ancestor]\
            and freq_SC_haps < Model.SC_freq_thresh:

        #t = time.time()
        # pr = cProfile.Profile()
        # pr.enable()
        #-------------------------------------------------------
        #                    duplicate alleles
        #-------------------------------------------------------

        haplotypes_dict, diploid_count_dict, haplotypes_count_dict\
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
        #t = time.time()
        E_Ri_Fj_dict = Logic.Update_AA_Interaction_dict\
            (E_Ri_Fj_dict, AA_rnases_dict, AA_slfs_dict, haplotypes_count_dict, haplotypes_dict)

        # -------------------------------------------------------
        #   calculate diploids frequencies in the next generation
        # -------------------------------------------------------
        #t = time.time()
        # pr = cProfile.Profile()
        # pr.enable()
        offspring_freq_arr, offspring_dict, delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict \
            = Logic.calc_freq_diploid_next_gen\
            (E_Ri_Fj_dict, diploid_count_dict, \
            haplotypes_count_dict, haplotypes_dict, \
             delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict)

        # -------------------------------------------------------
        #                  build next generation
        # -------------------------------------------------------
        #t = time.time()
        if np.abs(np.sum(offspring_freq_arr) - 1) > 1e-8:
            print("wrong offspring frequencies")
        diploid_count_dict, haplotypes_count_dict = \
            Logic.build_the_next_generation(offspring_freq_arr, offspring_dict)
        #print("num of haps = " + str(len(haplotypes_count_dict)))
        #print("next generation" + str(time.time() - t))

        # --------------------------------------------------------------------------------
        #    update haplotypes_dict
        # --------------------------------------------------------------------------------

        wanted_keys = haplotypes_count_dict.keys()
        bigdict = haplotypes_dict
        haplotypes_dict = dict((k, bigdict[k]) for k in wanted_keys if k in bigdict)

        #--------------------------------------------------------------------------------
        #    update delta_g_g_dict, delta_g_r_dict, cond1_g_g_dict, cond2_g_g_dict
        #--------------------------------------------------------------------------------

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
        b=6

        # pr.disable()
        # sortby = 'cumulative'
        # ps = pstats.Stats(pr).sort_stats(sortby)  # , stream=s
        # ps.print_stats()

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

        freq_SC_haps = num_SC_haps/Model.popu_size/2

        #-------------------------------------------------------------------------
        #     count the number of original ancestor RNases left in the population
        #-------------------------------------------------------------------------
        ancestor_rnases_arr = np.zeros(len(haplotypes_dict), dtype = int)
        for hap_ii, [hap_ind, hap] in enumerate(haplotypes_dict.items()):
            ancestor_rnases_arr[hap_ii] = hap.RNasesAncestralIndices[0]
        if len(np.unique(ancestor_rnases_arr)) == 1:
            more_than_one_original_ancestor = False

        #------------------------------------------------------------------------
        #   group haplotypes use Maximun Independent Set algorithm
        #------------------------------------------------------------------------
        if not more_than_one_original_ancestor:
            first_iter_of_one_ancestor = iter_ind
            d_iters = 10
            for current_hap in haplotypes_dict:
                parents_haps[current_hap] = haplotypes_dict[current_hap].PreviousHaplotype
            if iter_ind % d_iters == 0:
                with open(path + '/parents_haps.pkl', 'wb') as output:
                    pickle.dump(parents_haps, output, pickle.HIGHEST_PROTOCOL)


            '''
            d_iters = 10
            haps_lst = list(haplotypes_count_dict.keys())
            rnases_arr = np.zeros(len(haps_lst), dtype=int)
            for ii, hap in enumerate(haps_lst):
                rnases_arr[ii] = haplotypes_dict[hap].RNasesIndices[0]
            unique_RNases, counts = np.unique(rnases_arr, return_counts=True)

            mis_groups_hap_dict = {}
            mis_groups_rnase_dict = {}
            mis_non_grouped_haps_dict = {}
            mis_foreign_comp_groups_rnases_dict = {}
            mis_foreign_comp_groups_haps_dict = {}

            # order the haplotypes that are about to be grouped by their amount
            order = np.argsort(list(haplotypes_count_dict.values()))[::-1]  # np.argsort(rnases_count)[::-1]
            haps_arr = np.array(haps_lst)  # /np.sum(list(dt_new_any_RNase_dict.values())))
            haps_arr_sorted = haps_arr[order]

            # compatibility_haps_full_arr = 2 * np.ones((2*Model.popu_size, 2*Model.popu_size))
            # compatibility_haps_types_arr = 2*np.ones((haps_arr_sorted,haps_arr_sorted))
            # add first haplotype to first group (still empty)
            # first check if SI:
            ii = 0
            SC_cond = True
            while SC_cond:
                hap0 = haps_arr_sorted[ii]
                SC_cond = False
                if haplotypes_dict[hap0].IS_SC == 'True':
                    SC_cond = True

                ii += 1
            # if not SC_cond: # SI
            mis_groups_rnase_dict[0] = []
            mis_groups_rnase_dict[0].append(haplotypes_dict[hap0].RNasesIndices[0])
            mis_groups_hap_dict[0] = []
            mis_groups_hap_dict[0].append(hap0)
            for hap_partner in haps_arr_sorted[ii:]:
                # check if hap_partner is SI:
                SC_cond = False
                if haplotypes_dict[hap_partner].IS_SC== 'True':
                    SC_cond = True

                if not SC_cond:  # only SI haplotypes can be labeled in one of the groups
                    comp_check_dict = {}
                    comp_def = np.zeros(len(mis_groups_hap_dict.keys()), dtype=int)

                    partner_rnase_ind = haplotypes_dict[hap_partner].RNasesIndices[0]
                    rnase_AA2 = AA_rnases_dict[partner_rnase_ind]
                    for group_ii, [group_ind, hap_lst] in enumerate(mis_groups_hap_dict.items()):
                        comp_check_dict[group_ind] = np.zeros(len(hap_lst), dtype=int)
                        for hap_ii, hap in enumerate(hap_lst):
                            comp_cond1 = False
                            rnase_ind = haplotypes_dict[hap].RNasesIndices[0]
                            rnase_AA1 = AA_rnases_dict[rnase_ind]
                            for slf_ind in haplotypes_dict[hap_partner].SLFsIndices:
                                slf_AA1 = AA_slfs_dict[slf_ind]
                                if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                    E = np.sum(Model.eij[rnase_AA1, slf_AA1], axis=0)
                                else:
                                    E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                                if E < Model.interaction_thresh:
                                    # interacting
                                    comp_cond1 = True
                                    break
                            comp_cond2 = False
                            for slf_ind in haplotypes_dict[hap].SLFsIndices:
                                slf_AA2 = AA_slfs_dict[slf_ind]
                                if tuple([partner_rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                                    E = np.sum(Model.eij[rnase_AA2, slf_AA2], axis=0)
                                else:
                                    E = E_Ri_Fj_dict[partner_rnase_ind, slf_ind]
                                if E < Model.interaction_thresh:
                                    # interacting
                                    comp_cond2 = True
                                    break
                            if comp_cond1 and comp_cond2:
                                # symetric compatibility
                                # compatibility_haps_types_arr[] = 2
                                comp_check_dict[group_ind][hap_ii] = 2
                            elif not comp_cond1 and not comp_cond2:
                                comp_check_dict[group_ind][hap_ii] = 0
                            else:
                                comp_check_dict[group_ind][hap_ii] = 1

                    for group_ind in comp_check_dict.keys():
                        if all(comp_check_dict[group_ind] == 2):
                            comp_def[group_ind] = 2  # comp

                        elif all(comp_check_dict[group_ind] == 0):
                            comp_def[group_ind] = 0  # non_comp
                        else:
                            comp_def[group_ind] = 1

                    gro_label = np.where(np.array(comp_def) == 0)[0]
                    if all(comp_def == 2):  # new group
                        ind = len(mis_groups_hap_dict)
                        mis_groups_hap_dict[ind] = []
                        mis_groups_hap_dict[ind].append(hap_partner)
                        mis_groups_rnase_dict[ind] = []
                        mis_groups_rnase_dict[ind].append(partner_rnase_ind)
                    elif len(gro_label == 1):  # add to an exsiting already group
                        ind = gro_label[0]
                        mis_groups_hap_dict[ind].append(hap_partner)
                        mis_groups_rnase_dict[ind].append(partner_rnase_ind)
                    elif len(gro_label == 0):
                        b = 6
            mis_rnase_groups_dict = {}
            mis_hap_groups_dict = {}
            for group_ind, rnases_lst in mis_groups_rnase_dict.items():
                for rnase in rnases_lst:
                    mis_rnase_groups_dict[rnase] = group_ind

            for group_ind, haps_lst in mis_groups_hap_dict.items():
                for hap in haps_lst:
                    mis_hap_groups_dict[hap] = group_ind

            for rnase_ind in unique_RNases:
                rnase_AA = AA_rnases_dict[rnase_ind]
                mis_foreign_comp_groups_rnases_dict[rnase_ind] = []
                mis_foreign_comp_groups_haps_dict[rnase_ind] = []
                for group_ind, rnases_lst in mis_groups_rnase_dict.items():
                    if rnase_ind not in np.array(rnases_lst):
                        hap_lst = mis_groups_hap_dict[group_ind]
                        comp_arr = np.zeros(len(hap_lst), dtype=int)
                        for hap_ii, hap in enumerate(hap_lst):
                            comp_cond1 = False
                            for slf_ind in haplotypes_dict[hap].SLFsIndices:
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
                                comp_arr[hap_ii] = 1
                            else:
                                comp_arr[hap_ii] = 0
                                mis_foreign_comp_groups_rnases_dict[rnase_ind].append(group_ind)
                                mis_foreign_comp_groups_haps_dict[rnase_ind].append(hap)
                                b = 6

            non_grouped_haps = list(set(haplotypes_count_dict.keys()) - set(list(mis_hap_groups_dict.keys())))
            for non_grouped_hap in non_grouped_haps:
                mis_non_grouped_haps_dict[non_grouped_hap] = []
                # non_grouped_slfs_arr = np.array(haplotypes_dict[non_grouped_hap].SLFsIndices)
                for rnase_ind, group_ind in mis_rnase_groups_dict.items():
                    rnase_AA = AA_rnases_dict[rnase_ind]
                    comp_cond1 = False
                    for slf_ind in haplotypes_dict[non_grouped_hap].SLFsIndices:
                        slf_AA = AA_slfs_dict[slf_ind]
                        if tuple([rnase_ind, slf_ind]) not in E_Ri_Fj_dict:
                            E = np.sum(Model.eij[rnase_AA, slf_AA], axis=0)
                        else:
                            E = E_Ri_Fj_dict[rnase_ind, slf_ind]
                        if E < Model.interaction_thresh:
                            # interacting
                            comp_cond1 = True
                            break
                    if not comp_cond1:  # non grouped haps do not pollinate a grouped RNase
                        mis_non_grouped_haps_dict[non_grouped_hap].append(group_ind)
            b = 6
            mis_groups_hap_all_iters_dict[iter_ind] = mis_groups_hap_dict
            mis_groups_rnase_all_iters_dict[iter_ind] = mis_groups_rnase_dict
            mis_rnase_groups_all_iters_dict[iter_ind] = mis_rnase_groups_dict
            mis_hap_groups_all_iters_dict[iter_ind] = mis_hap_groups_dict
            mis_non_grouped_haps_all_iters_dict[iter_ind] = mis_non_grouped_haps_dict
            mis_foreign_comp_groups_rnases_all_iters_dict[iter_ind] = mis_foreign_comp_groups_rnases_dict
            mis_foreign_comp_groups_haps_all_iters_dict[iter_ind] = mis_foreign_comp_groups_haps_dict

            if first_iter_with_one_ancestor:
                copy_mis_groups_hap_all_iters_dict[iter_ind] = copy.deepcopy(mis_groups_hap_all_iters_dict[iter_ind])
                copy_mis_groups_rnase_all_iters_dict[iter_ind] = copy.deepcopy(
                    mis_groups_rnase_all_iters_dict[iter_ind])
                copy_mis_rnase_groups_all_iters_dict[iter_ind] = copy.deepcopy(
                    mis_rnase_groups_all_iters_dict[iter_ind])
                copy_mis_hap_groups_all_iters_dict[iter_ind] = copy.deepcopy(mis_hap_groups_all_iters_dict[iter_ind])
                copy_mis_non_grouped_haps_all_iters_dict[iter_ind] = copy.deepcopy(
                    mis_non_grouped_haps_all_iters_dict[iter_ind])
                copy_mis_foreign_comp_groups_rnases_all_iters_dict[iter_ind] = copy.deepcopy(
                    mis_foreign_comp_groups_rnases_all_iters_dict[iter_ind])
                copy_mis_foreign_comp_groups_haps_all_iters_dict[iter_ind] = copy.deepcopy(
                    mis_foreign_comp_groups_haps_all_iters_dict[iter_ind])
            else:

            first_iter_with_one_ancestor = False
            '''

        if iter_ind % d_iters == 0:
            print("start generation # %d" % iter_ind)
            Logic.save_files(path, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

            with open(path + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(path + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
                pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

            with open(path + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
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

    Logic.save_files(path, AA_slfs_dict, AA_rnases_dict, E_Ri_Fj_dict)

    with open(path + '/diploid_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(diploid_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/haplotypes_count_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(haplotypes_count_dict, output, pickle.HIGHEST_PROTOCOL)

    with open(path + '/haplotypes_dict_' + str(iter_ind) + '.pkl', 'wb') as output:
        pickle.dump(haplotypes_dict, output, pickle.HIGHEST_PROTOCOL)

    # [RNasesPartnersDict, SLFsPartnersDict,
    #         RNases_indices_all_iters, SLFs_indices_all_iters,
    #         previous_rnase_all_iters, previous_slf_all_iters,
    #         num_RNases_all_iters, num_SLFs_all_iters,
    #         entropy_RNases_freqs, entropy_SLFs_freqs] \
    #     = Analysis.RNases_SLFs_info_lists(haplotypes_obj, iters_haplotypes_dict,
    #                                       iters_haplotypes_count_dict, all_RNase_ar, all_SLF_ar, E_pool_Ri_Fj)

    # Analysis.write_RNases_SLFs_pkls_all_iters(path, RNases_indices_all_iters, previous_rnase_all_iters,
    #                                           SLFs_indices_all_iters, previous_slf_all_iters,
    #                                           RNasesPartnersDict, SLFsPartnersDict,
    #                                           entropy_RNases_freqs, entropy_SLFs_freqs)
