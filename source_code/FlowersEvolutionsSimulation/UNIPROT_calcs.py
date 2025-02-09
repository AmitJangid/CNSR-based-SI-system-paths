import numpy as np
import itertools
import scipy.io
from scipy.io import savemat
import simulationVars as Model

def calc_mutation_probs_bio_classes():
    #----------------------------------------------------------------------------
    #   AA order in the matrices - AA labeling and bio-chemical classification
    #----------------------------------------------------------------------------
    #           C   F   Y   W   M   L   I   V   G   P   A   T   S   N   H   Q   E   D   R   K
    #           4  13  18  17  12  10   9   19  7  14   0  16  15   2   8   5   6   3   1  11
    #Horvitz:   1   0   1   0?  0   0   0   0   0   0   0   1   1   1   1   1   3   3   2   2

    p_dioclass_inds_dict = {}
    p_dioclass_inds_dict[0] = [1,4,5,6,7,8,9,10,3]
    p_dioclass_inds_dict[1] = [0, 2, 11, 12, 13, 14, 15]
    p_dioclass_inds_dict[2] = [18,19]
    p_dioclass_inds_dict[3] = [16,17]

    p_AA1 = [0.0137, 0.0386, 0.0292, 0.0109, 0.0241, 0.0965, 0.0593, 0.0686, 0.0708, 0.0472, 0.0826, 0.0535, 0.066, \
             0.0406, 0.0227, 0.0393, 0.0674, 0.0546, 0.0553, 0.0582]

    p_AA1 = np.array(p_AA1)
    #------------------------------------------------------
    #   reorder q rows and columns by their bio-classes
    #------------------------------------------------------

    new_order = [p_dioclass_inds_dict[0], p_dioclass_inds_dict[1],
                 p_dioclass_inds_dict[2], p_dioclass_inds_dict[3]]

    new_order = [item for sublist in new_order for item in sublist]

    p_bioclass_inds_reordered_dict = {}
    tmp_len = len(p_dioclass_inds_dict[0])
    p_bioclass_inds_reordered_dict[0] = np.arange(0,tmp_len) #[1,4,5,6,7,8,9,10]
    tmp_len2 = tmp_len + len(p_dioclass_inds_dict[1])
    p_bioclass_inds_reordered_dict[1] = np.arange(tmp_len, tmp_len2)#[0, 2, 11, 12, 13, 14, 15]
    tmp_len3 = tmp_len2 + len(p_dioclass_inds_dict[2])
    p_bioclass_inds_reordered_dict[2] = np.arange(tmp_len2, tmp_len3)#[18,19]
    tmp_len4 = tmp_len3 + len(p_dioclass_inds_dict[3])
    p_bioclass_inds_reordered_dict[3] = np.arange(tmp_len3, tmp_len4) #[16,17]

    p_AA1_reordered = p_AA1[new_order]
    p_AA_bioclass = np.zeros((4))
    for ii in range(4):
        tmp_p = 0
        for ind in p_bioclass_inds_reordered_dict[ii]:
            tmp_p += p_AA1_reordered[ind]
        p_AA_bioclass[ii] = tmp_p

    p_AA_bioclass /= np.sum(p_AA_bioclass)
    # check if the p's are summed to 1
    sum_p_AA_bioclass = np.sum(p_AA_bioclass)

    #-----------------------------------
    #    save .mat files
    #-----------------------------------

    scipy.io.savemat('p_AA_bioclass.mat',
                        {'p_AA_bioclass' \
                            :p_AA_bioclass})


    return p_AA_bioclass

