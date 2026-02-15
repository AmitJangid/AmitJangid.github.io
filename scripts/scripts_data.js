const scriptData = {
    "scripts/cnsr_diploid_common/initial_prior.py": `#
#		SLF-RNAse interaction in Non-self recognition self-incompatible (SI) System
#
##########################################################################################

from pylab import genfromtxt
from multiprocessing import Pool#, TimeoutError, Process, Queue, Pipe, Value, Array
import os, math, pdb, time, random as ra, numpy as np, matplotlib.pyplot as plt, sys
import collections, pickle
from datetime import datetime
########################################################################################################
start_time = datetime.now()
print ('Starts|--------|', start_time)



#
#   This is the one we have used for the simulation -----------------------------------------------
#   #########################################################################
#
#   Interaction: One to More
#   In this initialization, each Haplotype will be having some constant number of SLF genes (dinsinct RNAse - 1)

#   First, generate certain number of distinct RNAse

#   Second, generate pools as many as number of distinct RNAse:
#   (Many-to-Many interactions)
#   1. generate one SLF and check its all partner RNAse, then put this SLF in the pools of partner RNAse  \\
#       repeat the same process untill each pool has set of equal number of SLF genes
#   
#   (One-to-One interaction)
#   1. generate one SLF and check its partner RNAse (randomly, check only for one partner), then put this SLF in the pool of partner RNAse \\
#       repeat the same process untill each pool has set of equal number of SLF genes
#   ** In any case, through the SLF if that SLF interact either with all RNAse or none

#   Third, choose SLF from different pools one by one checking it does not interact with same RNAse (where you are placing the SLF) /
#   and choosen SLF also does not have same copy on the same Haplotype
#
#
###############################################################################################
#--------------------------------------------------------------------------------------------------------------------------------------------------------
def initialization_4(input_Values):
    # index, to save the file, total population
    SampleNo, energyThreshold, folder_to_Save, seed_Index = input_Values

    #ra.seed(int(seed_Index))
    np.random.seed(int(seed_Index+1))

    mainPath = folder_to_Save
    pathToSave = mainPath
    #
    #
    #--------------------------------------------------------
    rnase_frequency = np.array([[0.3403, 0.3567, 0.1596, 0.1434],\\
                                [0.3350, 0.3314, 0.1616, 0.1720],\\
                                [0.3046, 0.2887, 0.1907, 0.2161],\\
                                [0.2710, 0.295,  0.2045, 0.2295],\\
                                [0.2209, 0.2946, 0.2290, 0.2555],\\
                                [0.2031, 0.2681, 0.2376, 0.2912],\\
                                [0.1717, 0.2576, 0.2619, 0.3087],\\
                                [0.1297, 0.2392, 0.2985, 0.3326],\\
                                [0.0958, 0.2042, 0.3178, 0.3822]])
    
    rnase_weight_index = int(energyThreshold+10)

    distinct_RNAse = 10                     # number of initial distinct haplotype
    length_RF = 18                          # length of RNAse and F box genes
    weight = [0.5, 0.265, 0.113, 0.122]     # weight of four types of Amino Acid (AA) occurance
    #weight = rnase_frequency[rnase_weight_index] # weight of four types of Amino Acid (AA) occurance in the simulation from previous work
    #weight = weight/sum(weight)
    fourTypesAA = [0,1,2,3]
    
    # energy of interacting pair; one RNAse and another F box gene
    energy_interaction = np.array([[-1,0,0,0],[0,-0.75,-0.25,-0.25],[0,-0.25,1,-1.25],[0,-0.25,-1.25,1]])
    def energyOfIntPair(RNAse, FBoxg):  return sum([energy_interaction[int(x),int(y)] for x, y in zip(RNAse, FBoxg)])

    no_Of_SLF_Each_haplotype = []

    #
    #
    # initialization of the population                                                                                          BLOCK-I--------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------
    #popRNAse, popFbox = [0 for i in range(tot_Population)], [0 for i in range(tot_Population)] # parental population
    tempRNAse = [np.random.choice(fourTypesAA, length_RF, p=weight) for i in range(distinct_RNAse)] # distinct RNAse, need to generate based on similarity
    tempFbox = [[] for i in range(distinct_RNAse)]                                     # distinct F box genes
    SLF_Pools = [[] for i in range(distinct_RNAse)]

    Exact_Size_Of_Pool = 500 if energyThreshold in [-10,-9,-8,-7] else 500 if energyThreshold in [-6, -5] else 500 if energyThreshold in [-4, -3] else 500

    Trials_Limit = 1000000
    trial = 1
    while trial <= Trials_Limit: # finding the pools of fixed size
        #FBoxg = ra.choices(fourTypesAA, weights=weight, k=length_RF)
        FBoxg = np.random.choice(fourTypesAA, length_RF, p=weight)

        interIndex = []
        for RNAse, indexRNAse in zip(tempRNAse, range(distinct_RNAse)):
            if energyOfIntPair(RNAse, FBoxg) < energyThreshold:
                interIndex.append(indexRNAse);
                # break # if do not want shared SLF on different haplotypes

        if len(interIndex) > 0 and len(interIndex) < distinct_RNAse:
            for index_inter in interIndex:
                if len(SLF_Pools[index_inter]) < Exact_Size_Of_Pool: # to make every pool at equal size
                    SLF_Pools[index_inter].append(FBoxg)

        length_Of_Each_SLF_Pool = [len(each_SLF_Pool) for each_SLF_Pool in SLF_Pools]

        if sum(length_Of_Each_SLF_Pool) == Exact_Size_Of_Pool*distinct_RNAse:
            #print (sum(length_Of_Each_SLF_Pool), Exact_Size_Of_Pool*distinct_RNAse)
            #print (length_Of_Each_SLF_Pool, trial)
            break

        trial += 1

    # generating haplotypes by selecting SLF from SLF pools
    length_Of_Each_SLF_Pool = [len(each_SLF_Pool) for each_SLF_Pool in SLF_Pools]
    if length_Of_Each_SLF_Pool == [Exact_Size_Of_Pool for i in range(distinct_RNAse)]:
        print ('pools are ready')

        inter_Pollen = []
        for hap_index in range(distinct_RNAse):
            inter_Each_Pollen = [0 for _ in range(distinct_RNAse)]
            RNASe_Index = [i for i in range(distinct_RNAse) if i != hap_index]
            for RNASe_index in RNASe_Index:
                for attempt in range(Exact_Size_Of_Pool):
                    ra_Selected = ra.randint(0,Exact_Size_Of_Pool-1)
                    selected_SLF = SLF_Pools[RNASe_index][ra_Selected]
                    all_SLF_On_Current_Hap = [tempFbox[hap_index][ij*length_RF:(ij+1)*length_RF] for ij in range(int(len(tempFbox[hap_index])/length_RF))]
                    checkPoint_Of_Identical = 0
                    for each_SLF_On_Current_Hap in all_SLF_On_Current_Hap: # to make sure no two copies are at the same haplotype
                        if np.array_equal(each_SLF_On_Current_Hap, selected_SLF) is True:
                            checkPoint_Of_Identical += 1

                    if energyOfIntPair(tempRNAse[hap_index], selected_SLF) >= energyThreshold and checkPoint_Of_Identical == 0:
                        tempFbox[hap_index] = np.append(tempFbox[hap_index], selected_SLF, 0)
                        inter_Each_Pollen[RNASe_index] = 1 #RNASe_index
                        break

            inter_Pollen.append(np.array(inter_Each_Pollen))

        no_Of_SLF_Each_haplotype = [len(each_Pollen)/length_RF for each_Pollen in tempFbox]

        print (np.array(inter_Pollen))
        print (no_Of_SLF_Each_haplotype)
        print (trial, SampleNo, energyThreshold, '\\n')

        # to make sure every haplotype has fixed number of SLF
        if sum(no_Of_SLF_Each_haplotype) == distinct_RNAse*(distinct_RNAse-1):
            print ('first conditions loop')
            # to make sure every haplotype is interacting with distinct_RNAse-1 hap
            if np.sum(np.array(inter_Pollen), axis=None) == distinct_RNAse*(distinct_RNAse-1):
                print ('second conditions loop')

                np.savetxt(pathToSave+'/slf_initial_{}.dat'.format(SampleNo), np.array(tempFbox), fmt='%i')
                np.savetxt(pathToSave+'/rnase_initial_{}.dat'.format(SampleNo), np.array(tempRNAse), fmt='%i')
                np.savetxt(pathToSave+'/slf_no_initial_{}.dat'.format(SampleNo), np.array(no_Of_SLF_Each_haplotype).T, fmt='%i')

                # with open(pathToSave+'/slf_initial_{}.dat'.format(SampleNo), 'w') as fileWrite:
                #     for each_Pollen in tempFbox:
                #         for data in each_Pollen:  fileWrite.write('%.0f ' % (data))
                #         fileWrite.write('\\n')

    print (SampleNo, 'is done:  Seed-', seed_Index, '\\n', '\\n')


if __name__ == '__main__':


    all_thresholds = [-10, -9.75,  -9.5,  -9.25,  -9,  -8.75,  -8.5,  -8.25,\\
                      -8,  -7.75,  -7.5,  -7.25,  -7,  -6.75,  -6.5,  -6.25,\\
                      -6,  -5.75,  -5.5,  -5.25,  -5,  -4.75,  -4.5,  -4.25,\\
                      -4,  -3.75,  -3.5,  -3.25,  -3,  -2.75,  -2.5,  -2.25,\\
                      -2]
    #all_thresholds = [-10]

    for threshold in all_thresholds:

        no_of_samples = 25

        seeds = [ra.randint(0,100000000) for junk in range(no_of_samples*5)]
        seeds = np.unique(seeds)
        seeds = seeds.astype(int)

        if len(seeds) >= no_of_samples:

            N_Proces = 1
            threshold, folder = threshold, './ethp{}'.format(threshold)

            # Finding indices that has not been saved
            # -----------------------------------------------------------
            def find_SampleNo(Samples):
                Sample_Index = []
                for SampleNo in range(Samples):
                    try:
                        data = np.genfromtxt('./ethp{}/slf_no_initial_{}.dat'.format(threshold,SampleNo))
                        #Sample_Index.append(SampleNo)
                    except:
                        'do THE following'
                        Sample_Index.append(SampleNo)
                return Sample_Index

            Samples = find_SampleNo(no_of_samples)

            Sample_Index    = [i for i in Samples]  # index of sample number
            Threshold       = [threshold for i in Samples] # threshold value
            Folder          = [folder for i in Samples] # path to save the output
            Pool(N_Proces).map(initialization_4, zip(Sample_Index, Threshold, Folder, seeds))






print ('Ends|----|H:M:S|{}'.format(datetime.now() - start_time), '\\n')
###############################################################################################################
###############################################################################################################

`
};