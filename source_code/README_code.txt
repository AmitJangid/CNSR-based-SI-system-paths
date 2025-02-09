The role of promiscuous molecular recognition in the evolution of RNase-based self-incompatibility
Keren Erez (1), Amit Jangid (1), Ohad Noy Feldheim (2) and Tamar Friedlander (1)

    (1) The Robert H. Smith Institute of Plant Sciences and Genetics in Agriculture
        Faculty of Agriculture, The Hebrew University of Jerusalem,
        P.O. Box 12 Rehovot 7610001, Israel
    (2) The Einstein Institute of Mathematics, Faculty of Natural Sciences,
        The Hebrew University of Jerusalem, Jerusalem 9190401, Israel.
       
    Correspondence: tamar.friedlander@mail.huji.ac.il.

Keren Erez, Amit Jangid: Equal contrubution
###############################################################################################################


This stochastic simulations studies the evolution of a finite population of individuals carrying an S-locus, via rounds of mutation and selection, in search for trajectories of allelic expansion, where crossbreeding between individuals is determined by their allelic content. 


Source Code:

The following is the order of running the python scripts to regenerate the whole datasets. To regenerate the whole datasets for a single run, it takes from 96 hrs to 120 hrs depending on the parameter sets:


1. group_complete_by_mis.py

   or run

   run_group_complete_by_mis_parallel.py -> group_complete_by_mis_parallel.py

2. calc_mis_hap_and_rnase_groups_dict.py

3. main_track_comp_groups.py

   or run

   run_main_track_comp_group_parallel.py -> main_track_comp_group_parallel.py

4. analyze_sequence_events_split.py

5. analyze_sequence_events_extenction.py

6. make_dataset_extinct_trajs.py

7. statistics_splt_num_groups.py

8. statistics_ext_num_groups.py

9. count_tot_time_in_gro.py


###############################################################################################################
