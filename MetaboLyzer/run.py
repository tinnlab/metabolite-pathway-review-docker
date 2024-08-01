#!/usr/bin/env python
# coding: utf-8

#document: https://sites.google.com/a/georgetown.edu/fornace-lab-informatics/home/metabolyzer

import csv,math,sys,pickle,re,random,os
import numpy as np
from scipy.stats import ks_2samp,anderson,scoreatpercentile,poisson,binom,pearsonr,t,norm
import rpy2.robjects as robjects
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from rpy2.robjects.packages import importr
from multiprocessing import Process, Queue, cpu_count
from lxml import html


ctrlRaw = csv.reader(open(sys.argv[1]), delimiter=',')
expRaw = csv.reader(open(sys.argv[2]), delimiter=',')

print("\n                                 _        _             __                          (\__/)  .~    ~. ))")
print("                      /\/\   ___| |_ __ _| |__   ___   / / _   _ _______ _ __       /O O  ./      .'")
print("                     /    \ / _ \ __/ _` | '_ \ / _ \ / / | | | |_  / _ \ '__|     {O__,   \    {")
print("                    / /\/\ \  __/ || (_| | |_) | (_) / /__| |_| |/ /  __/ |          / .  . )    \\")
print("                    \/    \/\___|\__\__,_|_.__/ \___/\____/\__, /___\___|_|          |-| '-' \    } ))")
print("                                                           |___/                    .(   _(   )_.'")
print("                                                                                   '---.~_ _ _&         ")
print("\n---------------------------------------MetaboLyzer (v7.2.0.p3)--------------------------------------\n")


print("*************************************************************************************************")
print("****                                   WARNING: SQUIRRELS                                    ****")
print("****                                     -------------                                       ****")
print("**** > \033[94mOUTPUT FILE NAMES ARE GENERIC AND WILL BE OVERWRITTEN IF NOT CHANGED\033[0m                  ****")
print("**** > \033[94mTO QUIT PRESS CTRL-C\033[0m                                                                  ****")
print("**** > \033[94mPRESSING ENTER IF NOTHING IS TYPED IN THE PROMPTS WILL UTILIZE THE INDICATED DEFAULTS\033[0m ****")
print("**** > \033[94mNORMALIZATION IS ONLY ALLOWED WHEN LOG TRANSFORM IS SELECTED\033[0m                          ****")
print("*************************************************************************************************\n")

print("<----------------------------------------Your terminal window should be at least this wide before proceeding-------------------------------------------->\n")

if len(sys.argv)==4:
    userfilterRaw = csv.reader(open(sys.argv[3]))
    userFilter_List = [set(),set()]
    for thing in userfilterRaw:
        if len(thing)>1:
            userFilter_List[0].add(thing[0])
            userFilter_List[1].add(thing[1])
        else:
            userFilter_List[0].add(thing[0])
    print("User defined filter list loaded!")
else:
    userFilter_List = False

folderName = "output"

zeros_Tr = input(" > Ion presence percentage cutoff [0.70 default]? ")
zeros_Tr = 0.3 if zeros_Tr=='' else 1-float(zeros_Tr)

transf_bool = input(" > Data transformation for statistical analysis: (1)none (2)log (3)inverse hyperbolic sine [2 default]? ")
transf_bool = 2 if transf_bool=='' else int(transf_bool)

if transf_bool == 2:
    norm_bool = input(" > Data normalization for statistical analysis: (1)none (2)standard gaussian (3)KDE max center-scaling [1 default]? ")
    norm_bool = 1 if norm_bool=='' else int(norm_bool)
else:
    norm_bool = 'NA'

outlier_bool = input(" > Remove outliers via 1.5 IQR (not recommended for low sample counts) [Y/n]? ")
outlier_bool = True if ((outlier_bool=="") | (outlier_bool.lower() == "y")) else False

if transf_bool == 1:
    Zexcl_bool = input(" > Exclude zero values in standard statistical analyses [Y/n]? ")
    Zexcl_bool = True if ((Zexcl_bool=="") | (Zexcl_bool.lower() == "y")) else False
elif (transf_bool == 2):
     Zexcl_bool = True
elif transf_bool == 3:
    Zexcl_bool = False
#elif (transf_bool == 4) | (transf_bool == 5):#--------HIDDEN FEATURE!
#    Zexcl_bool = True

sig = input(" > P-value [0.05 default]? ")
sig = 0.05 if sig=='' else float(sig)

bioFilter_bool = input(" > Filter metabolites via HMDB/KEGG/LIPIDMAPS/BioCyc putative identification [Y/n]? ")
bioFilter_bool = True if ((bioFilter_bool=="") | (bioFilter_bool.lower() == "y")) else False

if bioFilter_bool==True:

    species = input("  >> KEGG Species: (1) human (2) mouse (3) rat [1 default]? ")
    if species == '1':
        species = 'hsa'
    elif species == '2':
        species = 'mmu'
    elif species == '3':
        species = 'rno'
    else:
        species = 'hsa'
    
    BioCycDB = input("  >> BioCyc database: (1) HumanCyc (2) MouseCyc (3) MetaCyc [1 default]? ")
    if BioCycDB == '1': BioCycDB = 0
    elif BioCycDB == '2': BioCycDB = 1
    elif BioCycDB == '3': BioCycDB = 2
    else: BioCycDB = 0


    bioFilter_name = input("  >> Datafile pathname for HMDB accession ['HMDB_metabolites.pkl' default]? ")
    bioFilter_name = 'HMDB_metabolites.pkl' if bioFilter_name=='' else bioFilter_name
    pkl_file = open(bioFilter_name, 'rb')
    bioFilter_list = pickle.load(pkl_file, encoding="latin1")
    #print "HMDB total metabolites: %s" % len(bioFilter_list)

    bioFilter_name = input("  >> Datafile pathname for KEGG accession ['KEGG_metabolites.pkl' default]? ")
    bioFilter_name = 'KEGG_metabolites.pkl' if bioFilter_name=='' else bioFilter_name
    pkl_file2 = open(bioFilter_name, 'rb')
    keggMainfile = pickle.load(pkl_file2, encoding="latin1")
    keggFilter_list = keggMainfile[0]
    if species == 'hsa': pathways = keggMainfile[1]
    elif species == 'mmu': pathways = keggMainfile[2]
    elif species == 'rno': pathways = keggMainfile[3]
    #print "KEGG total metabolites: %s" % len(keggFilter_list)

    bioFilter_name = input("  >> Datafile pathname for LIPIDMAPS accession ['LIPIDMAPS_lipids.pkl' default]? ")
    bioFilter_name = 'LIPIDMAPS_lipids.pkl' if bioFilter_name=='' else bioFilter_name
    pkl_file3 = open(bioFilter_name, 'rb')
    lipidFilter_list = pickle.load(pkl_file3, encoding="latin1")
    #print "LIPIDMAPS total metabolites: %s" % len(lipidFilter_list)

    bioFilter_name = input("  >> Datafile pathname for BioCyc accession ['BIOCYC_metabolites.pkl' default]? ")
    bioFilter_name = 'BIOCYC_metabolites.pkl' if bioFilter_name=='' else bioFilter_name
    pkl_file4 = open(bioFilter_name, 'rb')
    biocycFilter_masterList = pickle.load(pkl_file4, encoding="latin1")
    biocycFilter_list = biocycFilter_masterList[BioCycDB]
    #print "BioCyc total metabolites: %s" % len(biocycFilter_list)

    bioMode = input(" >> Molecular species: (1) Positive (2) negative (3) merged [1 default]? ")
    bioMode = 1 if (bioMode=='') else int(bioMode)

    if bioMode==2:
        adductchoose = input("  >> Adducts (pick all that apply): (1) H (2) Cl (3) All [3 default]? ")
        adductchoose = "3" if adductchoose == '' else adductchoose
    elif bioMode==1:
        adductchoose = input("  >> Adducts (pick all that apply): (1) H (2) Na (3) NH4 (4) All [4 default]? ")
        adductchoose = "4" if adductchoose == '' else adductchoose
    elif bioMode==3:
        Nadductchoose = input("  >> Negative mode adducts (pick all that apply): (1) H (2) Cl (3) All [3 default]? ")
        Nadductchoose = "3" if Nadductchoose == '' else Nadductchoose

        Padductchoose = input("  >> Positive mode adducts (pick all that apply): (1) H (2) Na (3) NH4 (4) All [4 default]? ")
        Padductchoose = "4" if Padductchoose == '' else Padductchoose

        adductchoose = [Padductchoose,Nadductchoose]

    else: print("INVALID CHOICE MADE!!!!")

    mw_tol = input("  >> Molecular weight tolerance (ppm) [20 ppm default]? ")
    mw_tol = 20.0 if mw_tol=='' else float(mw_tol)


class bcolors:
    TEAL = '\033[96m'
    PINK = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'


def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack 
    result = [None] * k
    for i in range(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return result

def mannw(ctrlSample,expSample,ranktest,choose="pval"):
    RvectCtrl = robjects.FloatVector(ctrlSample)
    RvectExp = robjects.FloatVector(expSample)

    ranktestVect = ranktest(RvectCtrl,RvectExp)


    p_val = ranktestVect[2][0]
    _W = ranktestVect[0][0]

    #print p_val
    #print ctrlSample
    #print expSample

    if choose == "pval": return p_val
    else: return _W

def kolsmir(ctrlSample,expSample,choose="pval"):
    ksstat,p_val = ks_2samp(ctrlSample,expSample)

    if choose == "pval": return p_val
    else: return ksstat

def welchT(ctrlSample,expSample,ttest):
    RvectCtrl = robjects.FloatVector(ctrlSample)
    RvectExp = robjects.FloatVector(expSample)
    
    forwelch = {'var.equal': False}
    #welch's correction
    correctRvect = ttest(RvectCtrl,RvectExp, **forwelch)
    
    return correctRvect[0][0]

def for_ever(combinedData,combinedIndex,allStats,Nsize,testtype,depth=0,upIndex=0,tempOne="NA"):
    for _i in range(upIndex,len(combinedData)-((Nsize-1)-depth)):
        if depth==0:
            combinedIndex.remove(_i)
            tempOne = [_i]
        else:
            combinedIndex.remove(_i)
            tempOne.append(_i)

        if (depth+1)<Nsize: for_ever(combinedData,combinedIndex,allStats,Nsize,testtype,depth+1,_i+1,tempOne)
        else:
            set1 = [combinedData[indy] for indy in tempOne]
            set2 = [combinedData[indy] for indy in combinedIndex]

            tempStat = mannw(set1,set2,ranktest,"stat") if testtype==2 else welchT(set1,set2,ttest)
            allStats.append(tempStat)

        combinedIndex.append(_i)
        tempOne.remove(_i)



#-------------SMP----------------------------#
def smp_multiVar(data,qu,quCount,mshapiro):
    #robjects.packages.quiet_require("ICSNP")
    r = robjects.r

    for piece in data:
        A_ID = piece[0][0]
        B_ID = piece[0][1]
        AfiltCTRL = piece[1]
        BfiltCTRL = piece[2]
        AfiltEXP = piece[3]
        BfiltEXP = piece[4]

        readTemp = list(AfiltCTRL); readTemp.extend(BfiltCTRL)
        RhotellCTRL = r.matrix(r.c(robjects.FloatVector(readTemp)),ncol=2)

        readTemp = list(AfiltEXP); readTemp.extend(BfiltEXP)
        RhotellEXP = r.matrix(r.c(robjects.FloatVector(readTemp)),ncol=2)

        CTRLnorm = mshapiro(r.t(RhotellCTRL))
        EXPnorm = mshapiro(r.t(RhotellEXP))

        #if (AnormC=="Normal") & (BnormC=="Normal") & (AnormE=="Normal") & (BnormE=="Normal"):
            #print "combined: %s CTRL: %s EXP: %s" % (combCor,CTRLcor,EXPcor)

        if (CTRLnorm[1][0] > 0.10) & (EXPnorm[1][0] > 0.10):
            quCount.put(2)

            readTemp = list(AfiltCTRL); readTemp.extend(BfiltCTRL)
            RhotellCTRL = r.matrix(r.c(robjects.FloatVector(readTemp)),ncol=2)

            readTemp = list(AfiltEXP); readTemp.extend(BfiltEXP)
            RhotellEXP = r.matrix(r.c(robjects.FloatVector(readTemp)),ncol=2)

            RhotellOUT = r.HotellingsT2(RhotellCTRL,RhotellEXP)

            T2p_val = RhotellOUT[1][0]

            if T2p_val < 0.01:
                qu.put([A_ID,B_ID,T2p_val,piece])
        else: quCount.put(1)

def smp_MCresampler(data,sampN,R,qu):
    #robjects.packages.quiet_require("ICSNP")
    r = robjects.r

    for piece in data:
        A_ID = piece[0][0]; B_ID = piece[0][1]
        AfiltCTRL = piece[1]; BfiltCTRL = piece[2]
        AfiltEXP = piece[3]; BfiltEXP = piece[4]

        successRate = 0

        for i in range(0,R):
            random.seed()
            randIndexCTRL = [random.randint(0,len(AfiltCTRL)-1) for temp in range(0,sampN)]
            randIndexEXP = [random.randint(0,len(AfiltEXP)-1) for temp in range(0,sampN)]

            randSampA_C = [AfiltCTRL[index] for index in randIndexCTRL]
            randSampB_C = [BfiltCTRL[index] for index in randIndexCTRL]
            randSampA_E = [AfiltEXP[index] for index in randIndexEXP]
            randSampB_E = [BfiltEXP[index] for index in randIndexEXP]

            readTemp = list(randSampA_C); readTemp.extend(randSampB_C)
            RhotellCTRL = r.matrix(r.c(robjects.FloatVector(readTemp)),ncol=2)

            readTemp = list(randSampA_E); readTemp.extend(randSampB_E)
            RhotellEXP = r.matrix(r.c(robjects.FloatVector(readTemp)),ncol=2)

            RhotellOUT = r.HotellingsT2(RhotellCTRL,RhotellEXP)

            T2p_val = RhotellOUT[1][0]

            if T2p_val < 0.01:
                successRate += 1

        SucPerc = successRate/float(R)
        
        qu.put([A_ID,B_ID,SucPerc])
        #quCount.put(1)

def smp_permuteExact(ctrlSet,expSet,testtype,qu):
    r = robjects.r
    r.options(warn=-1)
    ttest = r('t.test')
    ranktest = r('wilcox.test')

    sub_count = len(ctrlSet)
    for i in range(0,sub_count):
        tempCTRL = [val for val in ctrlSet[i][1:] if val!='NA']
        tempEXP = [val for val in expSet[i][1:] if val!='NA']

        refStat1 = mannw(tempCTRL,tempEXP,ranktest,"stat") if testtype==2 else welchT(tempCTRL,tempEXP,ttest)
        refStat2 = mannw(tempEXP,tempCTRL,ranktest,"stat") if testtype==2 else welchT(tempEXP,tempCTRL,ttest)

        combinedData = list(tempCTRL); combinedData.extend(tempEXP)
        combinedIndex = list(range(0,len(combinedData)))

        allStats = []

        #print "    Datasize: %s %s" % (len(tempCTRL),len(tempEXP))

        for_ever(combinedData,combinedIndex,allStats,len(tempCTRL),testtype)

        statmean = np.mean(allStats)

        perm_pval = 0

        for nummy in allStats:
            if refStat1 > statmean:
                if (nummy >= refStat1) | (nummy <= refStat2): perm_pval += 1
            else: 
                if (nummy <= refStat1) | (nummy >= refStat2): perm_pval += 1

        perm_pval = (perm_pval-1)/float(len(allStats))

        qu.put([ctrlSet[i][0],perm_pval])

def smp_MCpermute(ctrlSet,expSet,Rep,enuff,testtype,qu):
    r = robjects.r
    r.options(warn=-1)
    ranktest = r('wilcox.test')

    sub_count = len(ctrlSet)
    for i in range(0,sub_count):
        tempCTRL = [val for val in ctrlSet[i][1:] if val!='NA']
        tempEXP = [val for val in expSet[i][1:] if val!='NA']

        allStats = []

        if enuff == True:
            refStat1 = kolsmir(tempCTRL,tempEXP,"stat") if testtype==2 else mannw(tempCTRL,tempEXP,ranktest,"stat")
            refStat2 = kolsmir(tempEXP,tempCTRL,"stat") if testtype==2 else mannw(tempEXP,tempCTRL,ranktest,"stat")
        else:
            refStat1 = mannw(tempCTRL,tempEXP,ranktest,"stat") if testtype==2 else welchT(tempCTRL,tempEXP,ttest)
            refStat2 = mannw(tempEXP,tempCTRL,ranktest,"stat") if testtype==2 else welchT(tempEXP,tempCTRL,ttest)

        combinedData = list(tempCTRL); combinedData.extend(tempEXP)
        lenCD = len(combinedData)
        cutL = len(tempCTRL)
        
        for j in range(0,Rep):
            random.seed()
            _remixed = random.sample(combinedData,lenCD)

            if enuff == True: tempStat = kolsmir(_remixed[:cutL],_remixed[cutL:],"stat") if testtype==2 else mannw(_remixed[:cutL],_remixed[cutL:],ranktest,"stat")
            else: tempStat = mannw(_remixed[:cutL],_remixed[cutL:],ranktest,"stat") if testtype==2 else welchT(_remixed[:cutL],_remixed[cutL:],ttest)
            
            allStats.append(tempStat)

        statmean = np.mean(allStats)

        perm_pval = 0

        for nummy in allStats:
            if not ((enuff==True) & (testtype==2)):
                if refStat1 > statmean:
                    if (nummy > refStat1) | (nummy < refStat2): perm_pval += 1
                else: 
                    if (nummy < refStat1) | (nummy > refStat2): perm_pval += 1
            else: #only for kolsmir case
                if (nummy > refStat1): perm_pval += 1

        perm_pval = perm_pval/float(len(allStats))

        qu.put([ctrlSet[i][0],perm_pval])

def smp_dualCorr(the_dataA,AdataN,the_dataB,BdataN,start,finish,analyze_count,qu,quCheck,quWarning):

    summy = 0
    Warnings = []

    for i in range(start,finish):
        unoA = the_dataA[i][1:]
        unoB = the_dataB[i][1:]

        for j in range(i+1, analyze_count):
            dosA = the_dataA[j][1:]
            dosB = the_dataB[j][1:]

    #        TheShitCount = 0
            unoAC = []; dosAC = []
            unoBC = []; dosBC = []

            for k in range(0, AdataN):
                if (unoA[k]!='NA') & (dosA[k]!='NA'):
                    unoAC.append(unoA[k])
                    dosAC.append(dosA[k])

            if len(unoAC) > 2: (pcorA,_pvalA) = pearsonr(unoAC,dosAC)
            else: pcorA = 0

            for k in range(0, BdataN):
                if (unoB[k]!='NA') & (dosB[k]!='NA'):
                    unoBC.append(unoB[k])
                    dosBC.append(dosB[k])

            if len(unoBC) > 2: (pcorB,_pvalB) = pearsonr(unoBC,dosBC)
            else: pcorB = 0

            goodEnough = True if ((len(unoBC)-3) > 0) & ((len(unoAC)-3) > 0) else False

            if (len(unoAC)-2) > 0: 
                tVarA = pcorA*math.sqrt(len(unoAC)-2)/math.sqrt(1-pcorA**2) #conversion to a t distributed variable
                pcorA_pval = 2*(t.cdf(tVarA,len(unoAC)-2)) if tVarA < 0 else 2*(1-t.cdf(tVarA,len(unoAC)-2))
                Za = 0.5*math.log((1+pcorA)/(1-pcorA)) # Fisher's transformation
            else:
                Warnings.append(1)
                pcorA_pval = 1

            if (len(unoBC)-2) > 0: 
                tVarB = pcorB*math.sqrt(len(unoBC)-2)/math.sqrt(1-pcorB**2) #conversion to a t distributed variable
                pcorB_pval = 2*(t.cdf(tVarB,len(unoBC)-2)) if tVarB < 0 else 2*(1-t.cdf(tVarB,len(unoBC)-2))
                Zb = 0.5*math.log((1+pcorB)/(1-pcorB)) # Fisher's transformation
            else:
                Warnings.append(1)
                pcorB_pval = 1

            if goodEnough == True:
                Zscore = (Za - Zb)/math.sqrt((1/float(len(unoAC)-3))+(1/float(len(unoBC)-3)))
                different_pval = 2*norm.cdf(Zscore) if Zscore < 0 else 2*(1-norm.cdf(Zscore))
            else:
                Warnings.append(1)
                different_pval = 1

            #signalDiff = random.random()/1000

                #print "Signal loss: %s %s %s" % (pcorA, pcorB, different_pval)
              #  signalDiff = abs(pcorB) - abs(pcorA)

            signalDiff = abs(pcorB) - abs(pcorA)
            qu.put([signalDiff,i,j,pcorA,pcorA_pval,pcorB,pcorB_pval,different_pval,the_dataA[i][0],the_dataA[j][0]])

            summy += 1
            if summy >= 100:
                quCheck.put(summy)
                summy = 0

            if len(Warnings) > 0:
                quWarning.put(1)
                Warnings = []

    quCheck.put(summy)


#--------------------------------------------#

def IQR(temp_total):
    outlier_count = 0
    total_shift = []
    all_shift = []
    for temp_shift in temp_total:
        filt_shift = [val for val in temp_shift if val!='NA']
        shift_Q1 = scoreatpercentile(filt_shift,25)
        shift_Q3 = scoreatpercentile(filt_shift,75)
        shift_IQR = shift_Q3 - shift_Q1
        LL = shift_Q1 - 1.5 * shift_IQR
        UL = shift_Q3 + 1.5 * shift_IQR

        final_shift = []
        all_temp = []

        for f in range(0,len(temp_shift)):

            if (temp_shift[f] != 'NA'):
                if (temp_shift[f] < LL):
                    outlier_count += 1
                    all_temp.append('NA')
                elif(temp_shift[f] > UL):
                    outlier_count += 1
                    all_temp.append('NA')
                else:
                    final_shift.append(temp_shift[f])
                    all_temp.append(temp_shift[f])

            else:
                all_temp.append('NA')

        total_shift.append(final_shift)
        all_shift.append(all_temp)

    return all_shift, outlier_count


def dataExtractor(CTRLreader, EXPreader, zeros_Tr, transf_bool,Zexcl_bool,norm_bool="NA"):
    print("\nReading in ctrl and exp datasets...")
    CTRLraw_data = []; EXPraw_data = []
    CTRLpart_data = []; EXPpart_data = []
    count = 0
    raw_count = 0

    if norm_bool == 3: CTRLall_data = []; EXPall_data = []
    
    CTRLsample_count = 0; EXPsample_count = 0
    
    bothCount = 0
    
    for CTRLrow in CTRLreader:
        EXProw = next(EXPreader)
        if raw_count >=1: #skip the first row
            CTRLzeros_perc = sum([1 if float(i)==0 else 0 for i in CTRLrow[1:]])/float(len(CTRLrow)-1)
            EXPzeros_perc = sum([1 if float(i)==0 else 0 for i in EXProw[1:]])/float(len(EXProw)-1)
            
            if (CTRLzeros_perc < zeros_Tr) & (EXPzeros_perc < zeros_Tr):
                bothCount += 1
                CTRLraw_data.append(CTRLrow)
                EXPraw_data.append(EXProw)
            elif (CTRLzeros_perc < zeros_Tr) | (EXPzeros_perc < zeros_Tr):
                count += 1
                CTRLpart_data.append(CTRLrow)
                EXPpart_data.append(EXProw)

            if norm_bool == 3: CTRLall_data.append(CTRLrow); EXPall_data.append(EXProw)

        else:
            names = CTRLrow[1:]
            names.extend(EXProw[1:])
            
        raw_count += 1
    
    
    CTRLsample_count = len(CTRLraw_data[0])-1
    EXPsample_count = len(EXPraw_data[0])-1

    print("   %s samples in the ctrl set and %s samples in the exp set" % (CTRLsample_count,EXPsample_count))
    print("   ions with >%s%% presence in both datasets (analyzed via standard statistics): %s" % ((100-zeros_Tr*100), bothCount))
    print("   ions with >%s%% presence in only one dataset (analyzed categorically):        %s" % ((100-zeros_Tr*100),count))
    print("   %s out of %s ions excluded from all analysis due to insufficient data" % (raw_count-(count+bothCount),raw_count))
    
    CTRL_catmatrix = np.zeros((count, CTRLsample_count), float)
    EXP_catmatrix = np.zeros((count, EXPsample_count), float)
    
    #scipy based array creation
    for i in range(0,count):        
        #CTRL_catmatrix[i,:] = [1 if float(abun) > 0 else 0 for abun in CTRLpart_data[i][1:]]
        #EXP_catmatrix[i,:] = [1 if float(abun) > 0 else 0 for abun in EXPpart_data[i][1:]]
        CTRL_catmatrix[i,:] = [float(abun) for abun in CTRLpart_data[i][1:]]
        EXP_catmatrix[i,:] = [float(abun) for abun in EXPpart_data[i][1:]]

    CTRL_cat = CTRL_catmatrix.tolist()
    EXP_cat = EXP_catmatrix.tolist()
    
    [CTRL_cat[j].insert(0,CTRLpart_data[j][0]) for j in range(0,count)] #insert the names back in to the dataframe!
    [EXP_cat[j].insert(0,EXPpart_data[j][0]) for j in range(0,count)]
    
    CTRL_matrix = np.zeros((bothCount, CTRLsample_count), float)
    EXP_matrix = np.zeros((bothCount, EXPsample_count), float)
    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    if transf_bool==2: # log transformation
        CTRLnormParams = []
        EXPnormParams = []
        if norm_bool == 1:
            for i in range(0,bothCount):
                CTRL_matrix[i,:] = [math.log(float(abun)) if float(abun)>0 else -1*np.inf for abun in CTRLraw_data[i][1:]]
                EXP_matrix[i,:] = [math.log(float(abun)) if float(abun)>0 else -1*np.inf for abun in EXPraw_data[i][1:]]
        else:
            for i in range(0,bothCount):
                CTRL_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in CTRLraw_data[i][1:]]
                EXP_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in EXPraw_data[i][1:]]

        if norm_bool == 2:
            print("   Normalizing via standard gaussian center-scaling...")
            for i in range(0,CTRLsample_count):
                CTRLtemp_sample = [math.log(CTRL_matrix[j,i]) for j in range(0,len(CTRL_matrix[:,i])) if (CTRL_matrix[j,i] > 0)]
                CTRLstd_dev = np.std(CTRLtemp_sample)
                CTRLmean = np.mean(CTRLtemp_sample)
                
                CTRL_matrix[:,i] = [(math.log(abun)-CTRLmean)/CTRLstd_dev if abun>0 else -1*np.inf for abun in CTRL_matrix[:,i]]
                CTRLnormParams.append([CTRLmean,CTRLstd_dev])
                
            for i in range(0,EXPsample_count):
                EXPtemp_sample = [math.log(EXP_matrix[j,i]) for j in range(0,len(EXP_matrix[:,i])) if (EXP_matrix[j,i] > 0)]
                EXPstd_dev = np.std(EXPtemp_sample)
                EXPmean = np.mean(EXPtemp_sample)
                
                EXP_matrix[:,i] = [(math.log(abun)-EXPmean)/EXPstd_dev if abun>0 else -1*np.inf for abun in EXP_matrix[:,i]]
                EXPnormParams.append([EXPmean,EXPstd_dev])

        elif norm_bool == 3:
            print("   Normalizing via kernel density estimation maximum center-scaling (KDEMAX)...")
            r = robjects.r
            CTRLall_matrix = np.zeros((raw_count-1, CTRLsample_count), float)
            EXPall_matrix = np.zeros((raw_count-1, EXPsample_count), float)

            for i in range(0,raw_count-1):
                CTRLall_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in CTRLall_data[i][1:]]
                EXPall_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in EXPall_data[i][1:]]

            for i in range(0,CTRLsample_count):
                CTRLtemp_sample = [math.log(CTRLall_matrix[j,i]) for j in range(0,len(CTRLall_matrix[:,i])) if (CTRLall_matrix[j,i] > 0)]
                CTRLr_sample = r.matrix(r.c(robjects.FloatVector(CTRLtemp_sample)), byrow=True, nrow=1)
                KDE_1024c = r.density(CTRLr_sample,n=1024)
                KDE_x = [KDE_1024c.rx2('x')[_g] for _g in range(0,1024)]; KDE_y = [KDE_1024c.rx2('y')[_g] for _g in range(0,1024)]
                _tempCoor = [(KDE_x[_coor],KDE_y[_coor]) for _coor in range(0,1024)]
                CTRL_findCenter = sorted(_tempCoor, key=lambda psig: psig[1], reverse=True)
                CTRLmean = CTRL_findCenter[0][0]

                #CTRLtemp_sampleTEST = [math.log(CTRL_matrix[j,i]) for j in range(0,len(CTRL_matrix[:,i])) if (CTRL_matrix[j,i] > 0)]
                #CTRLstd_devTEST = np.std(CTRLtemp_sampleTEST)
                #CTRLmeanTEST = np.mean(CTRLtemp_sampleTEST)

                _bigger = [_Temp for _Temp in CTRLtemp_sample if _Temp > CTRLmean]
                CTRLstd_dev = math.sqrt(float(sum([(_Temp-CTRLmean)**2 for _Temp in _bigger]))/len(_bigger))
                
                #print "mean %s %s" % (CTRLmeanTEST, CTRLmean)
                #print "sd %s %s" % (CTRLstd_devTEST,CTRLstd_dev)
                #print CTRLtemp_sampleTEST
                #print CTRLtemp_sample
                #print ""

                CTRL_matrix[:,i] = [(math.log(abun)-CTRLmean)/CTRLstd_dev if abun>0 else -1*np.inf for abun in CTRL_matrix[:,i]]
                CTRLnormParams.append([CTRLmean,CTRLstd_dev])
               # print CTRLtemp_sample
               # print [abun for abun in CTRL_matrix[:,i] if abun!=-1*np.inf]
                
            for i in range(0,EXPsample_count):
                EXPtemp_sample = [math.log(EXPall_matrix[j,i]) for j in range(0,len(EXPall_matrix[:,i])) if (EXPall_matrix[j,i] > 0)]
                EXPr_sample = r.matrix(r.c(robjects.FloatVector(EXPtemp_sample)), byrow=True, nrow=1)
                KDE_1024e = r.density(EXPr_sample,n=1024)
                KDE_x = [KDE_1024e.rx2('x')[_g] for _g in range(0,1024)]; KDE_y = [KDE_1024e.rx2('y')[_g] for _g in range(0,1024)]
                _tempCoor = [(KDE_x[_coor],KDE_y[_coor]) for _coor in range(0,1024)]
                EXP_findCenter = sorted(_tempCoor, key=lambda psig: psig[1], reverse=True)
                EXPmean = EXP_findCenter[0][0]

                _bigger = [_Temp for _Temp in EXPtemp_sample if _Temp > EXPmean]
                EXPstd_dev = math.sqrt(float(sum([(_Temp-EXPmean)**2 for _Temp in _bigger]))/len(_bigger))
                
                EXP_matrix[:,i] = [(math.log(abun)-EXPmean)/EXPstd_dev if abun>0 else -1*np.inf for abun in EXP_matrix[:,i]]
                EXPnormParams.append([EXPmean,EXPstd_dev])

        else: #no normalization
            CTRLnormParams = [[0,1] for _i in range(0,CTRLsample_count)]
            EXPnormParams = [[0,1] for _i in range(0,EXPsample_count)]


        CTRL_data = CTRL_matrix.tolist()
        EXP_data = EXP_matrix.tolist()
        
        #CTRL_temp = [[val for val in row if val != float("-inf")] for row in CTRL_data] #remove zeros
        #EXP_temp = [[val for val in row if val != float("-inf")] for row in EXP_data] #remove zeros
        CTRL_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in CTRL_data] #remove zeros
        EXP_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in EXP_data] #remove zeros
        
        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp


        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]
            
        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names, CTRLnormParams, EXPnormParams
    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    elif transf_bool==3: #inverse hyperbolic sine transformation
        
        for i in range(0,bothCount):
            CTRL_matrix[i,:] = [math.log(float(abun)+math.sqrt(float(abun)**2+1)) for abun in CTRLraw_data[i][1:]]
            EXP_matrix[i,:] = [math.log(float(abun)+math.sqrt(float(abun)**2+1)) for abun in EXPraw_data[i][1:]]

        CTRL_temp = CTRL_matrix.tolist()
        EXP_temp = EXP_matrix.tolist()

        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp
        
        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]

        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names

    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    elif (Zexcl_bool) & (transf_bool==1): #no transform but exclude zeros
        for i in range(0,bothCount):
            CTRL_matrix[i,:] = [float(abun) for abun in CTRLraw_data[i][1:]]
            EXP_matrix[i,:] = [float(abun) for abun in EXPraw_data[i][1:]]

        CTRL_data = CTRL_matrix.tolist()
        EXP_data = EXP_matrix.tolist()
        
        #CTRL_temp = [[val for val in row if val != 0] for row in CTRL_data] #remove zeros
        #EXP_temp = [[val for val in row if val != 0] for row in EXP_data] #remove zeros
        CTRL_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in CTRL_data] #remove zeros
        EXP_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in EXP_data] #remove zeros

        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp
        
        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]

        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names
    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    elif (Zexcl_bool==False) & (transf_bool==1): #no transform and no zeros removal
        for i in range(0,bothCount):
            CTRL_matrix[i,:] = [float(abun) for abun in CTRLraw_data[i][1:]]
            EXP_matrix[i,:] = [float(abun) for abun in EXPraw_data[i][1:]]

        CTRL_temp = CTRL_matrix.tolist()
        EXP_temp = EXP_matrix.tolist()

        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp
        
        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]
        
        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names
    
def dataMergedExtractor(CTRLreader, EXPreader, zeros_Tr, transf_bool,Zexcl_bool,norm_bool="NA"):
    print("\nReading in ctrl and exp datasets (merged modes)...")
    CTRLraw_data = []; EXPraw_data = []
    CTRLpart_data = []; EXPpart_data = []
    count = 0
    raw_count = 0

    modeList = []
    modeList_part = []

    if norm_bool == 3: CTRLall_data = []; EXPall_data = []; modeListAll = []
    
    CTRLsample_count = 0; EXPsample_count = 0
    
    bothCount = 0
    
    for CTRLrow_tmp in CTRLreader:
        chemMode = 2 if CTRLrow_tmp[-1]=="Neg" else 1

        EXProw_tmp = next(EXPreader)

        CTRLrow = CTRLrow_tmp[:-1]
        EXProw = EXProw_tmp[:-1]

        if raw_count >=1: #skip the first row
            CTRLzeros_perc = sum([1 if float(i)==0 else 0 for i in CTRLrow[1:]])/float(len(CTRLrow)-1)
            EXPzeros_perc = sum([1 if float(i)==0 else 0 for i in EXProw[1:]])/float(len(EXProw)-1)
            
            if (CTRLzeros_perc < zeros_Tr) & (EXPzeros_perc < zeros_Tr):
                bothCount += 1
                CTRLraw_data.append(CTRLrow)
                EXPraw_data.append(EXProw)
                modeList.append(chemMode)

            elif (CTRLzeros_perc < zeros_Tr) | (EXPzeros_perc < zeros_Tr):
                count += 1
                CTRLpart_data.append(CTRLrow)
                EXPpart_data.append(EXProw)
                modeList_part.append(chemMode)

            if norm_bool == 3: CTRLall_data.append(CTRLrow); EXPall_data.append(EXProw); modeListAll.append(chemMode)

        else:
            names = CTRLrow_tmp[1:]
            names.extend(EXProw_tmp[1:])
            
        raw_count += 1
    
    
    CTRLsample_count = len(CTRLraw_data[0])-1
    EXPsample_count = len(EXPraw_data[0])-1

    print("   %s samples in the ctrl set and %s samples in the exp set" % (CTRLsample_count,EXPsample_count))
    print("   ions with >%s%% presence in both datasets (analyzed via standard statistics): %s" % ((100-zeros_Tr*100), bothCount))
    print("   ions with >%s%% presence in only one dataset (analyzed categorically):        %s" % ((100-zeros_Tr*100),count))
    print("   %s out of %s ions excluded from all analysis due to insufficient data" % (raw_count-(count+bothCount),raw_count))
    
    CTRL_catmatrix = np.zeros((count, CTRLsample_count), float)
    EXP_catmatrix = np.zeros((count, EXPsample_count), float)
    
    #scipy based array creation
    for i in range(0,count):        
        #CTRL_catmatrix[i,:] = [1 if float(abun) > 0 else 0 for abun in CTRLpart_data[i][1:]]
        #EXP_catmatrix[i,:] = [1 if float(abun) > 0 else 0 for abun in EXPpart_data[i][1:]]
        CTRL_catmatrix[i,:] = [float(abun) for abun in CTRLpart_data[i][1:]]
        EXP_catmatrix[i,:] = [float(abun) for abun in EXPpart_data[i][1:]]

    CTRL_cat = CTRL_catmatrix.tolist()
    EXP_cat = EXP_catmatrix.tolist()
    
    [CTRL_cat[j].insert(0,CTRLpart_data[j][0]) for j in range(0,count)] #insert the names back in to the dataframe!
    [EXP_cat[j].insert(0,EXPpart_data[j][0]) for j in range(0,count)]
    
    CTRL_matrix = np.zeros((bothCount, CTRLsample_count), float)
    EXP_matrix = np.zeros((bothCount, EXPsample_count), float)
    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    if transf_bool==2: # log transformation
        CTRLnormParams = []
        EXPnormParams = []
        if norm_bool == 1:
            for i in range(0,bothCount):
                CTRL_matrix[i,:] = [math.log(float(abun)) if float(abun)>0 else -1*np.inf for abun in CTRLraw_data[i][1:]]
                EXP_matrix[i,:] = [math.log(float(abun)) if float(abun)>0 else -1*np.inf for abun in EXPraw_data[i][1:]]
        else:
            for i in range(0,bothCount):
                CTRL_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in CTRLraw_data[i][1:]]
                EXP_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in EXPraw_data[i][1:]]

        if norm_bool == 2:
            print("   Normalizing via standard gaussian center-scaling...")
            for i in range(0,CTRLsample_count):
                CTRLtemp_sampleN = [math.log(CTRL_matrix[j,i]) for j in range(0,len(CTRL_matrix[:,i])) if ((CTRL_matrix[j,i] > 0) & (modeList[j]==2))]
                CTRLtemp_sampleP = [math.log(CTRL_matrix[j,i]) for j in range(0,len(CTRL_matrix[:,i])) if ((CTRL_matrix[j,i] > 0) & (modeList[j]==1))]
                
                CTRLstd_dev = [np.std(CTRLtemp_sampleP), np.std(CTRLtemp_sampleN)]
                CTRLmean = [np.mean(CTRLtemp_sampleP), np.mean(CTRLtemp_sampleN)]
                
                CTRL_matrix[:,i] = [(math.log(CTRL_matrix[j,i])-CTRLmean[modeList[j]-1])/CTRLstd_dev[modeList[j]-1] if (CTRL_matrix[j,i] > 0) else -1*np.inf for j in range(0,len(CTRL_matrix[:,i]))]
                CTRLnormParams.append([CTRLmean,CTRLstd_dev])
                
            for i in range(0,EXPsample_count):
                EXPtemp_sampleN = [math.log(EXP_matrix[j,i]) for j in range(0,len(EXP_matrix[:,i])) if ((EXP_matrix[j,i] > 0) & (modeList[j]==2))]
                EXPtemp_sampleP = [math.log(EXP_matrix[j,i]) for j in range(0,len(EXP_matrix[:,i])) if ((EXP_matrix[j,i] > 0) & (modeList[j]==1))]
                
                EXPstd_dev = [np.std(EXPtemp_sampleP), np.std(EXPtemp_sampleN)]
                EXPmean = [np.mean(EXPtemp_sampleP), np.mean(EXPtemp_sampleN)]
                
                EXP_matrix[:,i] = [(math.log(EXP_matrix[j,i])-EXPmean[modeList[j]-1])/EXPstd_dev[modeList[j]-1] if (EXP_matrix[j,i] > 0) else -1*np.inf for j in range(0,len(EXP_matrix[:,i]))]
                EXPnormParams.append([EXPmean,EXPstd_dev])

        elif norm_bool == 3:
            print("   Normalizing via kernel density estimation maximum center-scaling (KDEMAX)...")
            r = robjects.r
            CTRLall_matrix = np.zeros((raw_count-1, CTRLsample_count), float)
            EXPall_matrix = np.zeros((raw_count-1, EXPsample_count), float)

            for i in range(0,raw_count-1):
                CTRLall_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in CTRLall_data[i][1:]]
                EXPall_matrix[i,:] = [float(abun) if float(abun)>0 else -1*np.inf for abun in EXPall_data[i][1:]]

            for i in range(0,CTRLsample_count):
                CTRLtemp_sampleN = [math.log(CTRLall_matrix[j,i]) for j in range(0,len(CTRLall_matrix[:,i])) if ((CTRLall_matrix[j,i] > 0) & (modeListAll[j]==2))]
                CTRLr_sampleN = r.matrix(r.c(robjects.FloatVector(CTRLtemp_sampleN)), byrow=True, nrow=1)
                KDE_1024N = r.density(CTRLr_sampleN,n=1024)
                KDE_xN = [KDE_1024N.rx2('x')[_g] for _g in range(0,1024)]; KDE_yN = [KDE_1024N.rx2('y')[_g] for _g in range(0,1024)]
                _tempCoorN = [(KDE_xN[_coor],KDE_yN[_coor]) for _coor in range(0,1024)]
                CTRL_findCenterN = sorted(_tempCoorN, key=lambda psig: psig[1], reverse=True)
                CTRLmeanN = CTRL_findCenterN[0][0]
                _biggerN = [_Temp for _Temp in CTRLtemp_sampleN if _Temp > CTRLmeanN]
                CTRLstd_devN = math.sqrt(float(sum([(_Temp-CTRLmeanN)**2 for _Temp in _biggerN]))/len(_biggerN))

                CTRLtemp_sampleP = [math.log(CTRLall_matrix[j,i]) for j in range(0,len(CTRLall_matrix[:,i])) if ((CTRLall_matrix[j,i] > 0) & (modeListAll[j]==1))]
                CTRLr_sampleP = r.matrix(r.c(robjects.FloatVector(CTRLtemp_sampleP)), byrow=True, nrow=1)
                KDE_1024P = r.density(CTRLr_sampleP,n=1024)
                KDE_xP = [KDE_1024P.rx2('x')[_g] for _g in range(0,1024)]; KDE_yP = [KDE_1024P.rx2('y')[_g] for _g in range(0,1024)]
                _tempCoorP = [(KDE_xP[_coor],KDE_yP[_coor]) for _coor in range(0,1024)]
                CTRL_findCenterP = sorted(_tempCoorP, key=lambda psig: psig[1], reverse=True)
                CTRLmeanP = CTRL_findCenterP[0][0]
                _biggerP = [_Temp for _Temp in CTRLtemp_sampleP if _Temp > CTRLmeanP]
                CTRLstd_devP = math.sqrt(float(sum([(_Temp-CTRLmeanP)**2 for _Temp in _biggerP]))/len(_biggerP))
                
                CTRLstd_dev = [CTRLstd_devP, CTRLstd_devN]
                CTRLmean = [CTRLmeanP, CTRLmeanN]

                CTRL_matrix[:,i] = [(math.log(CTRL_matrix[j,i])-CTRLmean[modeList[j]-1])/CTRLstd_dev[modeList[j]-1] if (CTRL_matrix[j,i] > 0) else -1*np.inf for j in range(0,len(CTRL_matrix[:,i]))]
                CTRLnormParams.append([CTRLmean,CTRLstd_dev])
               # print CTRLtemp_sample
               # print [abun for abun in CTRL_matrix[:,i] if abun!=-1*np.inf]
                
            for i in range(0,EXPsample_count):
                EXPtemp_sampleN = [math.log(EXPall_matrix[j,i]) for j in range(0,len(EXPall_matrix[:,i])) if ((EXPall_matrix[j,i] > 0) & (modeListAll[j]==2))]
                EXPr_sampleN = r.matrix(r.c(robjects.FloatVector(EXPtemp_sampleN)), byrow=True, nrow=1)
                KDE_1024N = r.density(EXPr_sampleN,n=1024)
                KDE_xN = [KDE_1024N.rx2('x')[_g] for _g in range(0,1024)]; KDE_yN = [KDE_1024N.rx2('y')[_g] for _g in range(0,1024)]
                _tempCoorN = [(KDE_xN[_coor],KDE_yN[_coor]) for _coor in range(0,1024)]
                EXP_findCenterN = sorted(_tempCoorN, key=lambda psig: psig[1], reverse=True)
                EXPmeanN = EXP_findCenterN[0][0]
                _biggerN = [_Temp for _Temp in EXPtemp_sampleN if _Temp > EXPmeanN]
                EXPstd_devN = math.sqrt(float(sum([(_Temp-EXPmeanN)**2 for _Temp in _biggerN]))/len(_biggerN))

                EXPtemp_sampleP = [math.log(EXPall_matrix[j,i]) for j in range(0,len(EXPall_matrix[:,i])) if ((EXPall_matrix[j,i] > 0) & (modeListAll[j]==1))]
                EXPr_sampleP = r.matrix(r.c(robjects.FloatVector(EXPtemp_sampleP)), byrow=True, nrow=1)
                KDE_1024P = r.density(EXPr_sampleP,n=1024)
                KDE_xP = [KDE_1024P.rx2('x')[_g] for _g in range(0,1024)]; KDE_yP = [KDE_1024P.rx2('y')[_g] for _g in range(0,1024)]
                _tempCoorP = [(KDE_xP[_coor],KDE_yP[_coor]) for _coor in range(0,1024)]
                EXP_findCenterP = sorted(_tempCoorP, key=lambda psig: psig[1], reverse=True)
                EXPmeanP = EXP_findCenterP[0][0]
                _biggerP = [_Temp for _Temp in EXPtemp_sampleP if _Temp > EXPmeanP]
                EXPstd_devP = math.sqrt(float(sum([(_Temp-EXPmeanP)**2 for _Temp in _biggerP]))/len(_biggerP))
                
                EXPstd_dev = [EXPstd_devP, EXPstd_devN]
                EXPmean = [EXPmeanP, EXPmeanN]
                
                EXP_matrix[:,i] = [(math.log(EXP_matrix[j,i])-EXPmean[modeList[j]-1])/EXPstd_dev[modeList[j]-1] if (EXP_matrix[j,i] > 0) else -1*np.inf for j in range(0,len(EXP_matrix[:,i]))]
                EXPnormParams.append([EXPmean,EXPstd_dev])

        else: #no normalization
            CTRLnormParams = [[[0,0],[1,1]] for _i in range(0,CTRLsample_count)]
            EXPnormParams = [[[0,0],[1,1]] for _i in range(0,EXPsample_count)]


        CTRL_data = CTRL_matrix.tolist()
        EXP_data = EXP_matrix.tolist()
        
        #CTRL_temp = [[val for val in row if val != float("-inf")] for row in CTRL_data] #remove zeros
        #EXP_temp = [[val for val in row if val != float("-inf")] for row in EXP_data] #remove zeros
        CTRL_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in CTRL_data] #remove zeros
        EXP_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in EXP_data] #remove zeros
        
        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp


        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]
            
        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names, CTRLnormParams, EXPnormParams, modeList, modeList_part
    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    elif transf_bool==3: #inverse hyperbolic sine transformation
        
        for i in range(0,bothCount):
            CTRL_matrix[i,:] = [math.log(float(abun)+math.sqrt(float(abun)**2+1)) for abun in CTRLraw_data[i][1:]]
            EXP_matrix[i,:] = [math.log(float(abun)+math.sqrt(float(abun)**2+1)) for abun in EXPraw_data[i][1:]]

        CTRL_temp = CTRL_matrix.tolist()
        EXP_temp = EXP_matrix.tolist()

        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp
        
        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]

        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names, modeList, modeList_part

    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    elif (Zexcl_bool) & (transf_bool==1): #no transform but exclude zeros
        for i in range(0,bothCount):
            CTRL_matrix[i,:] = [float(abun) for abun in CTRLraw_data[i][1:]]
            EXP_matrix[i,:] = [float(abun) for abun in EXPraw_data[i][1:]]

        CTRL_data = CTRL_matrix.tolist()
        EXP_data = EXP_matrix.tolist()
        
        #CTRL_temp = [[val for val in row if val != 0] for row in CTRL_data] #remove zeros
        #EXP_temp = [[val for val in row if val != 0] for row in EXP_data] #remove zeros
        CTRL_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in CTRL_data] #remove zeros
        EXP_temp = [[val if val != float("-inf") else 'NA' for val in row] for row in EXP_data] #remove zeros

        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp
        
        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]

        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names, modeList, modeList_part
    #--------------------------------------------------------------------------------------------------------------------------------------------------------#
    elif (Zexcl_bool==False) & (transf_bool==1): #no transform and no zeros removal
        for i in range(0,bothCount):
            CTRL_matrix[i,:] = [float(abun) for abun in CTRLraw_data[i][1:]]
            EXP_matrix[i,:] = [float(abun) for abun in EXPraw_data[i][1:]]

        CTRL_temp = CTRL_matrix.tolist()
        EXP_temp = EXP_matrix.tolist()

        if outlier_bool==True:
            CTRL_data,ctrlOUTLIER = IQR(CTRL_temp)
            EXP_data,expOUTLIER = IQR(EXP_temp)
            print("   %s outliers detected in the ctrl set and %s in the exp set via 1.5 IQR" % (ctrlOUTLIER,expOUTLIER))
        else:
            CTRL_data = CTRL_temp
            EXP_data = EXP_temp
        
        [CTRL_data[j].insert(0,CTRLraw_data[j][0]) for j in range(0,bothCount)] #insert the names back in to the dataframe!
        [EXP_data[j].insert(0,EXPraw_data[j][0]) for j in range(0,bothCount)]
        
        return CTRL_data, EXP_data, CTRL_cat, EXP_cat, bothCount, count, CTRLraw_data, EXPraw_data, names, modeList, modeList_part

def normality(sample):
    Astat,critVals,sig = anderson(sample)
    if Astat > critVals[1]: #for the 10% critical value
        return "Non-Normal"
    else:
        return "Normal"

def bothT(ctrlSample,expSample,ttest):
    
    RvectCtrl = robjects.FloatVector(ctrlSample)
    RvectExp = robjects.FloatVector(expSample)
    
    forstd = {'var.equal': True}
    #equivariance
    equiRvect = ttest(RvectCtrl,RvectExp, **forstd)
    
    forwelch = {'var.equal': False}
    #welch's correction
    correctRvect = ttest(RvectCtrl,RvectExp, **forwelch)
    
    return equiRvect[2][0],correctRvect[2][0]

def fishtest(ctrlCAT,expCAT,r,fe_test):
    ctrlNZcount = sum(ctrlCAT)
    expNZcount = sum(expCAT)
    
    ctrlZcount = len(ctrlCAT)-ctrlNZcount
    expZcount = len(expCAT)-expNZcount
    
    feRvect = fe_test(r.matrix(r.c(ctrlNZcount,expNZcount,ctrlZcount,expZcount), nr = 2))
    
    return feRvect[0][0],float(ctrlNZcount)/float(len(ctrlCAT)),float(expNZcount)/float(len(expCAT))

def confidenceInt(ctrlSample,expSample):
    ctrlN = len(ctrlSample)
    expN = len(expSample)
    ctrlMean = np.mean(ctrlSample)
    expMean = np.mean(expSample)

    ctrlSD = np.std(ctrlSample)
    expSD = np.std(expSample)

    ctrl_CI = [ctrlMean-1.96*ctrlSD/math.sqrt(ctrlN),ctrlMean+1.96*ctrlSD/math.sqrt(ctrlN)]
    exp_CI = [expMean-1.96*expSD/math.sqrt(expN),expMean+1.96*expSD/math.sqrt(expN)]

    if ((ctrlMean > exp_CI[0]) & (ctrlMean < exp_CI[1])) | ((expMean > ctrl_CI[0]) & (expMean < ctrl_CI[1])):
        status = "Insignificant"
    else:
        status = "Significant"

    return status,ctrl_CI,exp_CI,ctrlSD,expSD

def bioIdentifier(KEGGref_pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,bioMode,adductchoose,mw_tol, pathAgg,BpathAgg, pathDisp = False):

    mammalian_count = 0
    HMDBputative_names = []
    
    total_count = 0
    KEGG_count = 0
    BIOCYC_count = 0

    KEGGputative_names = []
    KEGGindv_paths = []
    KEGG_pathways = set()
    StraightKEGGID = []

    LIPIDMAPSputative_names = []

    BIOCYCputative_names = []
    BIOCYCindv_paths = []
    BIOCYC_pathways  =set()
    StraightBIOCYCID = []
    if BioCycDB==2:
        BIOCYCindv_tax = dict()
        BIOCYC_tax = set()


    ion_id = re.split("_",ID[0])
    massR = float(ion_id[0])
    pw_id = re.compile('\d+')

    if bioMode==2: #negative mode
        #for negative mode, we only concentrate on Cl and H for now
        massR_perturbations = []
        if (adductchoose=="3"): massR_perturbations = [massR+1.00728,massR-34.9927]
        if re.search('1',adductchoose) != None: massR_perturbations.append(massR+1.00728)
        if re.search('2',adductchoose) != None: massR_perturbations.append(massR-34.9927)

        #print massR_perturbations
    else: #positive mode
        massR_perturbations = []
        #for positive mode, we only concentrate on H, Na, and NH4 for now
        if (adductchoose=="4"): massR_perturbations = [massR-1.00728,massR-22.9892,massR-18.0338]
        if re.search('1',adductchoose) != None: massR_perturbations.append(massR-1.00728)
        if re.search('2',adductchoose) != None: massR_perturbations.append(massR-22.9892)
        if re.search('3',adductchoose) != None: massR_perturbations.append(massR-18.0338)


    for perb in massR_perturbations:

        #for HMDB
        for met in bioFilter_list:
            temp_mass = met[1]
            check_ppm = (abs(temp_mass-perb)/perb)*1000000

            if (check_ppm < mw_tol):
                if met[2] == "Mammalian Metabolite":
                    mammalian_count += 1
                    HMDBputative_names.append(met[3])
                total_count += 1

        #for KEGG        
        for cpd in keggFilter_list:
            temp_mass = cpd[3]
            check_ppm = (abs(temp_mass-perb)/perb)*1000000

            if (check_ppm < mw_tol):
                temp_pathlist = set()
                [temp_pathlist.add(re.search(pw_id,thing[0]).group(0)) for thing in cpd[2]]
                hum_set = temp_pathlist & KEGGref_pathways
                if len(hum_set) > 0:
                    KEGGputative_names.append(cpd[0]+": "+cpd[1])
                    StraightKEGGID.append(cpd[0])
                    KEGGindv_paths.append(hum_set)
                    [KEGG_pathways.add(thingies) for thingies in cpd[2] if re.search(pw_id,thingies[0]).group(0) in hum_set]
                KEGG_count += 1

        #for lipidmaps
        for LMID in lipidFilter_list:
            temp_mass = LMID[1]
            check_ppm = (abs(temp_mass-perb)/perb)*1000000

            if (check_ppm < mw_tol): LIPIDMAPSputative_names.append([LMID[0],LMID[2],LMID[3],LMID[4]])

        #for BioCyc        
        for uid in biocycFilter_list:
            temp_mass = uid[2]
            check_ppm = (abs(temp_mass-perb)/perb)*1000000

            if (check_ppm < mw_tol):
                BIOCYCputative_names.append(uid[0]+": "+uid[1])
                StraightBIOCYCID.append(uid[0])

                if len(uid[3]) > 0:
                    temp_pathlist = set()
                    [temp_pathlist.add(thing[0]) for thing in uid[3]]

                    BIOCYCindv_paths.append(temp_pathlist)
                    #print uid[3]
                    if BioCycDB!=2: [BIOCYC_pathways.add(thingies) for thingies in uid[3]]
                    else:
                        [BIOCYC_pathways.add((thingies[0],thingies[1])) for thingies in uid[3]]

                        [[BIOCYC_tax.add(str(thing[0])+"\t: "+str(thing[1])+" > "+str(thing[2])) for thing in pathy[2]] for pathy in uid[3]]

                        for thingies in uid[3]:
                            if thingies[0] not in BIOCYCindv_tax:
                                temp_taxlist = []
                                for thing in thingies[2]:
                                    temp_taxlist.append(thing[0])
                                BIOCYCindv_tax[thingies[0]] = temp_taxlist
                else:
                    BIOCYCindv_paths.append([])

                BIOCYC_count += 1



    if pathDisp != False:
        #HMDB output
        print((bcolors.TEAL+"      HMDB (mammalian/total): %s/%s"+bcolors.ENDC) % (len(HMDBputative_names),total_count))
        [sys.stdout.write("       "+nam+"\n") for nam in HMDBputative_names]

        #KEGG output
        print((bcolors.GREEN+"      KEGG (species/total):   %s/%s"+bcolors.ENDC) % (len(KEGGputative_names),KEGG_count))
        #[sys.stdout.write("       "+nam+"\n") for nam in KEGGputative_names]
        [sys.stdout.write("       "+KEGGputative_names[ind]+" "+str(list(KEGGindv_paths[ind]))+"\n") for ind in range(len(KEGGputative_names))]

        if len(KEGGputative_names) != 0:
            print(bcolors.GREEN+"      KEGG pathways: "+bcolors.ENDC)
            [sys.stdout.write("       "+pathinfo[0]+" "+pathinfo[1]+"\n") for pathinfo in KEGG_pathways]

        #LIPIDMAPS output
        if len(LIPIDMAPSputative_names) != 0:
            print(bcolors.YELLOW+"      LIPIDMAPS: "+bcolors.ENDC)
            [sys.stdout.write("       "+nam[0]+": "+nam[1]+"\n") for nam in LIPIDMAPSputative_names]

        #BioCyc output
        print((bcolors.PINK+"      BioCyc:   %s"+bcolors.ENDC) % BIOCYC_count)
        [sys.stdout.write("       "+html.fromstring(BIOCYCputative_names[ind]).text_content()+" "+str(list(BIOCYCindv_paths[ind]))+"\n") for ind in range(len(BIOCYCputative_names))]
        if len(BIOCYCputative_names) != 0:
            print((bcolors.PINK+"      BioCyc pathways: "+bcolors.ENDC))
            if BioCycDB != 2: [sys.stdout.write("       "+pathinfo[0]+": "+html.fromstring(pathinfo[1]).text_content()+"\n") for pathinfo in BIOCYC_pathways]
            else:
                [sys.stdout.write("       "+pathinfo[0]+": "+html.fromstring(pathinfo[1]).text_content()+" "+str(list(BIOCYCindv_tax[pathinfo[0]]))+"\n") for pathinfo in BIOCYC_pathways]
                print((bcolors.PINK+"      MetaCyc based NCBI taxonomy identification:"+bcolors.ENDC))
                [sys.stdout.write("       "+taxinfo+"\n") for taxinfo in BIOCYC_tax]


            
    for pathinfo in KEGG_pathways:
        pathID = (pathinfo[0],pathinfo[1])
        if pathID in pathAgg:
            pathAgg[pathID].append(ID[0])
        else:
            pathAgg[pathID] = [ID[0]]

    for pathinfo in BIOCYC_pathways:
        pathID = (pathinfo[0],pathinfo[1])
        if pathID in BpathAgg:
            BpathAgg[pathID].append(ID[0])
        else:
            BpathAgg[pathID] = [ID[0]]


    HMDBfinal_status = str(mammalian_count)+"/"+str(total_count)
    KEGGfinal_status = str(len(KEGGputative_names))+"/"+str(KEGG_count)


    return HMDBfinal_status,HMDBputative_names, KEGGfinal_status,KEGGputative_names,StraightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID

def MDS(DataMat,totSamp,dist,cmdscale,isoMDS,type):
    params1 = {'eig':True, 'k':3}

    d = dist(DataMat) #calculate the euclidean distances

    if type==1:
        fit = cmdscale(d,**params1)
    else:
        fit = isoMDS(d, k=3)

    projections = fit.rx2('points')

    pc1 = [projections[i] for i in range(0,totSamp)]
    pc2 = [projections[i] for i in range(totSamp,totSamp*2)]
    pc3 = [projections[i] for i in range(totSamp*2,totSamp*3)]

    return pc1,pc2,pc3

def kPCA(DataMat,totSamp,kern,slot,r):
    kType = input("  >> Utilize the (1)Laplacian kernel (2)nth degree polynomial kernel (3)ANOVA kernel (4)Gaussian kernel [2 default]? ")
    kType = 2 if ((kType=='') | (kType == 2)) else int(kType)

    if kType == 1:
        zigma = input("     ~> Sigma value [0.2 default]? ")
        zigma = 0.2 if (zigma=='') else float(zigma)
        kpc = kern.kpca(DataMat,kernel="laplacedot",kpar=r.list(sigma=zigma),features=3)

    if kType == 2: #polynomial kernel
        nth = input("     ~> nth Degree for polynomial kernel [3 default]? ")
        nth = 3 if (nth=='') else int(nth)

        settings = {"degree":nth, "scale":1}
        kpc = kern.kpca(DataMat,kernel="polydot",kpar=r.list(**settings),features=3)
    if kType == 3:
        #kpc = kern.kpca(dataQframe,kernel="besseldot",kpar=r.list(sigma=1, order=1,degree=1),features=3)
        kpc = kern.kpca(DataMat,kernel="anovadot",kpar=r.list(sigma=1,degree=1),features=3)
    #if kType == 4:
    #    kpc = kern.kpca(dataQframe,kernel="splinedot",kpar=r.list(),features=3)
    if kType == 4:
        zigma = input("     ~> Sigma value [0.05 default]? ")
        zigma = 0.05 if (zigma=='') else float(zigma)
        kpc = kern.kpca(DataMat,kernel="rbfdot",kpar=r.list(sigma=zigma),features=3)

    projections = slot(kpc,"rotated")
    
    pc1 = [projections[i] for i in range(0,totSamp)]
    pc2 = [projections[i] for i in range(totSamp,totSamp*2)]
    pc3 = [projections[i] for i in range(totSamp*2,totSamp*3)]

    return pc1,pc2,pc3

def ICA(DataMat,totSamp,fastICA):
    params = {'alg.typ':"parallel","fun":"exp","alpha":1,"method":"C","row.norm":False,"maxit":200,"tol":0.0001,"verbose":False}

    ICAresult = fastICA(DataMat,3,**params)

    projections = ICAresult.rx2("S")

    pc1 = [projections[i] for i in range(0,totSamp)]
    pc2 = [projections[i] for i in range(totSamp,totSamp*2)]
    pc3 = [projections[i] for i in range(totSamp*2,totSamp*3)]

    return pc1,pc2,pc3

def heatMapMaker(forGFX,gfx_nu,totSamp,datatype,ctrlN,expN,r):

    robjects.packages.quiet_require("gplots")

    colscheme = input('  >> Color scheme: (1)green-black-red (2)blue-black-orange [1 default]? ')


    #r.library("gplots")
    colscheme = r.greenred(75) if (colscheme=="1") | (colscheme=="") else r.colorpanel(75, "blue","black","orange")
    #r.library("RColorBrewer")
    #hmcols=r.colorRampPalette(r.c("blue4", "cyan", "white", "yellow", "red4"))

    Rread = []

    #Rprint = []
    #[Rprint.append([(j-np.mean(forGFX[i][datatype]))/np.std(forGFX[i][datatype]) for j in forGFX[i][datatype]]) for i in range(1,gfx_nu)]
    #for i in range(1,gfx_nu):
    #    print "-------------------"
    #    print forGFX[i][0]
    #    print forGFX[i][datatype]
    #    print Rprint[i-1]
    if (transf_bool==1) | (transf_bool==3): #you do range control for things with zeros in them, or if it's not log transformed,
        colNorm_bool = True
    else:
        colNorm_bool = input("  >> Balance color (1)across all features or (2)within each feature [1 default]? ")
        colNorm_bool = False if (colNorm_bool=="") | (colNorm_bool=="1") else True

    if colNorm_bool == True:  # local color balance
        [Rread.extend(
            [(j - np.mean(forGFX[i][datatype])) / np.std(forGFX[i][datatype]) 
            if abs((j - np.mean(forGFX[i][datatype])) / np.std(forGFX[i][datatype]))<1.3 else math.copysign(1.3, (j - np.mean(forGFX[i][datatype])) / np.std(forGFX[i][datatype]))
            for j in forGFX[i][datatype]]
        ) for i in range(1, gfx_nu)]

    else:  # global color balance
        [Rread.extend(
            [(j - np.mean(forGFX[i][datatype])) if abs(j - np.mean(forGFX[i][datatype])) < 1.7 else math.copysign(1.7, (j - np.mean(forGFX[i][datatype]))) for j in forGFX[i][datatype]]
        ) for i in range(1,gfx_nu)]
    #-----------------
    #[Rread.extend([(j-np.mean(forGFX[i][datatype])) for j in forGFX[i][datatype]]) for i in range(1,gfx_nu)]

    #[Rread.extend(forGFX[i][transfPCA]) for i in range(1,gfx_nu)]
    #[Rread.extend([1 if (j-np.mean(forGFX[i][datatype]))/np.std(forGFX[i][datatype])>0 else -1 for j in forGFX[i][datatype]]) for i in range(1,gfx_nu)]

    MatParams = {'ncol':totSamp, 'dimnames':r.list(r.c([forGFX[i][0] for i in range(1,gfx_nu)]),r.c(forGFX[0][1:])),'byrow':True}
    HeatMat =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)

    heatmap = r("heatmap.2")
    #hcols = r('topo.colors')
    r.png(file=folderName+"/sigHeatmap.png",width=3000,height=3000,res=154)
    colorBar = []; colorBar.extend(["#FF0000" for i in range(0,ctrlN)]); colorBar.extend(["#0000FF" for i in range(0,expN)])
    aschar = r('as.character')

    heatParam = {'col':colscheme,'ylab':"Ions",'xlab':"Samples","main":"Heat Map for Significant Ions",
        'margins': r.c(8, 8), 'Colv': r('NA'), 'ColSideColors': aschar(r.c(colorBar)),'scale':'none','key':True,'trace':"none",'keysize':0.7, 'srtCol': 45}

    heatmap(HeatMat, **heatParam)

    print("  Creating PNG file for heatmap to \033[94msigHeatmap.png\033[0m")
    devoff = robjects.r('dev.off')
    devoff()


def multivariateTesting(analyze_count,CTRL_data,EXP_data,ctrlN,expN,sig_list,r):
    robjects.packages.quiet_require("mvnormtest")
    robjects.packages.quiet_require("ICSNP")
    mshapiro = r("mshapiro.test")

    wrHT = csv.writer(open('Hotellings_sig.csv', 'w'), delimiter=',')
    wrHT.writerow(["ion A","Ion B","Hotelling's T2 p-value"])

    T2numtest = 0
    T2sigCount = 0

    newSig_dict = {}

    bin1 = []
    bin2 = []
    bin3 = []
    bin4 = []
    countoff = 0

    sigPiece = []

    print("    %s univariately insignificant full presence ions to be analyzed bivariately" % (analyze_count-len(sig_list)))

    for i in range(0,analyze_count-1):
        AtempCTRL = CTRL_data[i][1:]
        AtempEXP = EXP_data[i][1:]
        A_ID = CTRL_data[i][0]

        for j in range(i+1,analyze_count):
            BtempCTRL = CTRL_data[j][1:]
            BtempEXP = EXP_data[j][1:]
            B_ID = CTRL_data[j][0]

            if (A_ID not in sig_list) & (B_ID not in sig_list):

                AfiltCTRL = []; BfiltCTRL = []; AfiltEXP = []; BfiltEXP = []
                for k in range(0,ctrlN):
                    if (AtempCTRL[k]!= 'NA') & (BtempCTRL[k]!='NA'):
                        AfiltCTRL.append(AtempCTRL[k]); BfiltCTRL.append(BtempCTRL[k])

                for k in range(0,expN):
                    if (AtempEXP[k]!= 'NA') & (BtempEXP[k]!='NA'):
                        AfiltEXP.append(AtempEXP[k]); BfiltEXP.append(BtempEXP[k])

                if countoff%4==0: bin1.append([[A_ID,B_ID],AfiltCTRL,BfiltCTRL,AfiltEXP,BfiltEXP])
                if countoff%4==1: bin2.append([[A_ID,B_ID],AfiltCTRL,BfiltCTRL,AfiltEXP,BfiltEXP])
                if countoff%4==2: bin3.append([[A_ID,B_ID],AfiltCTRL,BfiltCTRL,AfiltEXP,BfiltEXP])
                if countoff%4==3: bin4.append([[A_ID,B_ID],AfiltCTRL,BfiltCTRL,AfiltEXP,BfiltEXP])

                countoff += 1

    qu = Queue()
    quCount = Queue()
    
    p = Process(target=smp_multiVar, args=(bin1,qu,quCount,mshapiro))
    q = Process(target=smp_multiVar, args=(bin2,qu,quCount,mshapiro))
    r = Process(target=smp_multiVar, args=(bin3,qu,quCount,mshapiro))
    s = Process(target=smp_multiVar, args=(bin4,qu,quCount,mshapiro))

    p.start()
    q.start()
    r.start()
    s.start()

    progress = 0

    itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

    while (itsAlive | (not qu.empty()) | (not quCount.empty())):

        if(quCount.empty() != True):
            progress += 1
            temp = quCount.get()
            if temp == 2: T2numtest += 1
            sys.stdout.write("\r    Progress:  %s out of %s (%s%%) ion pairs analyzed        " % (progress,countoff,100*(float(progress)/countoff)))
            sys.stdout.flush()

        if (qu.empty() != True):
            temp1 = qu.get()

            A_ID = temp1[0]
            B_ID = temp1[1]
            T2p_val = temp1[2]

            sigPiece.append(temp1[3])

            T2sigCount += 1
            if A_ID in newSig_dict: newSig_dict[A_ID] += 1
            else: newSig_dict[A_ID] = 1
            if B_ID in newSig_dict: newSig_dict[B_ID] += 1
            else: newSig_dict[B_ID] = 1

            wrHT.writerow([A_ID,B_ID,T2p_val])

        itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

    p.join()
    q.join()
    r.join()
    s.join()

    qu.close()
    quCount.close()


    print("\n    %s normally distributed bivariables identified via multivariate Shapiro Wilks test" % T2numtest)
    print("    %s bivariately significant ion pairs identified from %s univariately insignificant ions" % (T2sigCount,len(list(newSig_dict.keys()))))

    #calculate probability of significant results
    UsefulProb = poisson.cdf(T2sigCount,0.01*T2numtest)

    #wrHT.writerow([])
    #wrHT.writerow(["ion","probability of significance"])
    T2printout = []
    for key in newSig_dict:
        T2ionCount = newSig_dict[key]
        T2ionProb = poisson.cdf(T2ionCount,(1.0/len(list(newSig_dict.keys())))*T2sigCount)
        T2printout.append([key,T2ionCount,1-T2ionProb])
        #wrHT.writerow([key,T2ionProb])
    T2printoutS = sorted(T2printout, key=lambda psig: psig[1], reverse=True)

    print("\n    -------> Top ion frequencies <-------\n    ion\t\t\tcount\tp-value")
    for row in T2printoutS:
        if row[2] < 0.01: print("    %s\t%s\t%s" % (row[0],row[1],row[2]))

    print("")

    print("Results saved to Hotellings_sig.csv")
    print("Complete! Probability of validity: %s%%"% (UsefulProb*100))

    MCveri = input("Verify bivariate results via Monte Carlo resampling [Y/n]? ")
    MCveri = True if ((MCveri=="") | (MCveri.lower() == "y")) else False

    if MCveri==True:
        wrMC = csv.writer(open(folderName+'/Hotellings_MCverify.csv', 'w'), delimiter=',')
        wrMC.writerow(["ion A","Ion B","Monte Carlo resampling pass rate"])

        binA = []; binB = []; binC = []; binD = []
        countoff = 0
        MC_r = input(" > Replicates [10000 default]? ")
        MC_r = 10000 if (MC_r=="") else int(MC_r)

        veriCount = 0

        for sub in sigPiece:
            if countoff%4 == 0: binA.append(sub)
            if countoff%4 == 1: binB.append(sub)
            if countoff%4 == 2: binC.append(sub)
            if countoff%4 == 3: binD.append(sub)

            countoff += 1

        sampN = int(round(float(min(expN,ctrlN))*0.75))

        print("    Subsampling size: %s\n" % sampN)
        print("    ------> Verification results <------")

        qu = Queue()
        
        p = Process(target=smp_MCresampler, args=(binA,sampN,MC_r,qu))
        q = Process(target=smp_MCresampler, args=(binB,sampN,MC_r,qu))
        r = Process(target=smp_MCresampler, args=(binC,sampN,MC_r,qu))
        s = Process(target=smp_MCresampler, args=(binD,sampN,MC_r,qu))

        p.start()
        q.start()
        r.start()
        s.start()

        itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

        while (itsAlive | (not qu.empty())):

            if(qu.empty() != True):
                temp = qu.get()
                wrMC.writerow([temp[0],temp[1],temp[2]])
                if temp[2] < 0.90: 
                    if (len(temp[0]) + len(temp[1])) > 27: print("    %s <<>> %s :\t %s" % (temp[0],temp[1],temp[2]))
                    else: print("    %s <<>> %s :\t\t %s" % (temp[0],temp[1],temp[2]))
                else:
                    if (len(temp[0]) + len(temp[1])) > 27: print("    %s <<>> %s :\t %s\t <-- VERIFIED" % (temp[0],temp[1],temp[2]))
                    else: print("    %s <<>> %s :\t\t %s <-- VERIFIED" % (temp[0],temp[1],temp[2]))
                    veriCount += 1

            itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

        p.join()
        q.join()
        r.join()
        s.join()

        qu.close()

        print("Results saved to Hotellings_MC.csv")
        print("Complete! %s verified significant bivariate ion pairs" % veriCount)

def permute(analyze_count,CTRL_data,EXP_data,ctrlN,expN,testtype,enuff):

    permSigCount = 0

    if ctrlN <= 10:

        if testtype==2: print("ctrl sample size is 10 or less! Proceeding with exact permuted Mann-Whitney U test...")
        else: print("ctrl sample size is 10 or less! Proceeding with exact permuted Welch's T-test...")
        
        Cbin1 = []; Cbin2 = []; Cbin3 = []; Cbin4 = []
        Ebin1 = []; Ebin2 = []; Ebin3 = []; Ebin4 = []

        countoff = 0
        permute_pvals = {}

        for i in range(0,analyze_count):
            if countoff%4 == 0: Cbin1.append(CTRL_data[i]); Ebin1.append(EXP_data[i])
            if countoff%4 == 1: Cbin2.append(CTRL_data[i]); Ebin2.append(EXP_data[i])
            if countoff%4 == 2: Cbin3.append(CTRL_data[i]); Ebin3.append(EXP_data[i])
            if countoff%4 == 3: Cbin4.append(CTRL_data[i]); Ebin4.append(EXP_data[i])

            countoff += 1

        qu = Queue()
        
        p = Process(target=smp_permuteExact, args=(Cbin1,Ebin1,testtype,qu))
        q = Process(target=smp_permuteExact, args=(Cbin2,Ebin2,testtype,qu))
        r = Process(target=smp_permuteExact, args=(Cbin3,Ebin3,testtype,qu))
        s = Process(target=smp_permuteExact, args=(Cbin4,Ebin4,testtype,qu))

        p.start()
        q.start()
        r.start()
        s.start()

        progress = 0

        itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

        while (itsAlive | (not qu.empty())):


            if (qu.empty() != True):
                progress += 1
                sys.stdout.write("\r    Progress:  %s out of %s (%s%%) ions analyzed        " % (progress,countoff,100*(float(progress)/countoff)))
                sys.stdout.flush()
                temp1 = qu.get()
                #print "New ion complete! %s" % (temp1)

                ionName = temp1[0]
                ion_permPval = temp1[1]
                permute_pvals[ionName] = ion_permPval
                if ion_permPval < sig: permSigCount += 1


            itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

        p.join()
        q.join()
        r.join()
        s.join()

        qu.close()

        print("\nComplete! %s significant ions found\n" % permSigCount)

        return permute_pvals


    else:
        if (enuff==True) & (testtype==2): print("ctrl sample size exceeds 10! Proceeding with Monte Carlo permuted K-S test...")
        elif ((enuff==True) & (testtype==1)) | ((enuff==False) & (testtype==2)): print("ctrl sample size exceeds 10! Proceeding with Monte Carlo permuted Mann-Whitney U test...")
        elif (enuff==False) & (testtype==1): print("ctrl sample size exceeds 10! Proceeding with Monte Carlo permuted Welch's T-test...")

        MC_r = input(" > Replicates [10000 default]? ")
        MC_r = 10000 if (MC_r=="") else int(MC_r)

        Cbin1 = []; Cbin2 = []; Cbin3 = []; Cbin4 = []
        Ebin1 = []; Ebin2 = []; Ebin3 = []; Ebin4 = []

        countoff = 0
        permute_pvals = {}

        for i in range(0,analyze_count):
            if countoff%4 == 0: Cbin1.append(CTRL_data[i]); Ebin1.append(EXP_data[i])
            if countoff%4 == 1: Cbin2.append(CTRL_data[i]); Ebin2.append(EXP_data[i])
            if countoff%4 == 2: Cbin3.append(CTRL_data[i]); Ebin3.append(EXP_data[i])
            if countoff%4 == 3: Cbin4.append(CTRL_data[i]); Ebin4.append(EXP_data[i])

            countoff += 1

        qu = Queue()
        
        p = Process(target=smp_MCpermute, args=(Cbin1,Ebin1,MC_r,enuff,testtype,qu))
        q = Process(target=smp_MCpermute, args=(Cbin2,Ebin2,MC_r,enuff,testtype,qu))
        r = Process(target=smp_MCpermute, args=(Cbin3,Ebin3,MC_r,enuff,testtype,qu))
        s = Process(target=smp_MCpermute, args=(Cbin4,Ebin4,MC_r,enuff,testtype,qu))

        p.start()
        q.start()
        r.start()
        s.start()

        progress = 0

        itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

        while (itsAlive | (not qu.empty())):


            if (qu.empty() != True):
                progress += 1
                sys.stdout.write("\r    Progress:  %s out of %s (%s%%) ions analyzed        " % (progress,countoff,100*(float(progress)/countoff)))
                sys.stdout.flush()
                temp1 = qu.get()
                #print "New ion complete! %s" % (temp1)

                ionName = temp1[0]
                ion_permPval = temp1[1]
                permute_pvals[ionName] = ion_permPval
                if ion_permPval < sig: permSigCount += 1


            itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

        p.join()
        q.join()
        r.join()
        s.join()

        qu.close()

        print("\nComplete! %s significant ions found" % permSigCount)

        return permute_pvals

def correlationHeatMapMaker(the_data,names,dataN,analyze_count,filename, dendroInfo, corType):
    robjects.packages.quiet_require("gplots")
    r = robjects.r
    heatmap = r("heatmap.2")

    Rread = []
    [Rread.extend([thing if thing!='NA' else robjects.NA_Real for thing in the_data[i][1:]]) for i in range(0,analyze_count)]

    MatParams = {'ncol':dataN, 'dimnames':r.list(r.c([the_data[i][0] for i in range(0,analyze_count)]),r.c(names)),'byrow':True}
    DataMat =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)
    DataMat = r.t(DataMat) #I know this is lazy but it's just easier!

    print("   * Calculating correlation matrix")
    U = r.cor(DataMat, use="pairwise.complete.obs")
    

    #colscheme = raw_input('  >> Color scheme: (1)green-black-red (2)blue-black-orange [1 default]? ')
    #colscheme = r.greenred(75) if (colscheme=="1") | (colscheme=="") else r.colorpanel(75, "blue","black","orange")

    #r.library("RColorBrewer")
    #hmcols=r.colorRampPalette(r.c("blue4", "cyan", "white", "yellow", "red4"))
    #colscheme = hmcols(75)

    #colscheme = r.bluered(75)
    #colscheme = r.colorpanel(75, "green","black","red")

    #hcols = r('topo.colors')
    print("   * Creating PNG file of heatmap to \033[94m%s\033[0m" % filename)
    if analyze_count < 1000:
        r.png(file=folderName+"/"+filename,width=8000,height=8000,res=120)
    else:
        r.png(file=folderName+"/"+filename,width=12000,height=12000,res=120)
    if corType == 1:
        colscheme = input('  >> Color scheme: (1)green-black-red (2)teal-black-orange [2 default]? ')
        colscheme = r.colorpanel(75, "#0BB5FF","black","#F44A00")  if (colscheme=="2") | (colscheme=="") else r.colorpanel(75, "#14D9C2","black","#EB263D")
        heatParam = {'Rowv':dendroInfo,'Colv':dendroInfo,'col':colscheme,'ylab':"Ions",'xlab':"Ions","main":"Correlation Heat Map (Pearson's)",
            'margins':r.c(10,10),'symm':True,'symbreaks':True,'scale':'none','key':True,'symkey':True,'trace':"none",'keysize':0.7}

        heatmap(U, **heatParam)

    elif corType == 2:
        dissimilarity = ((r.abs(U)).ro*-1).ro+1
        colscheme = r.colorpanel(75, "#660000","white")
        heatParam = {'Rowv':dendroInfo,'Colv':dendroInfo,'col':colscheme,'ylab':"Ions",'xlab':"Ions","main":"Dissimilarity Heat Map",
            'margins':r.c(10,10),'symm':True,'scale':'none','key':True,'trace':"none",'keysize':0.7}
    
        heatmap(dissimilarity, **heatParam)


    devoff = robjects.r('dev.off')
    devoff()

def pureCorMap(theMat,names,analyze_count,filename, dendroInfo, corType, title):
    robjects.packages.quiet_require("gplots")
    r = robjects.r
    heatmap = r("heatmap.2")

    Rread = []
    [Rread.extend(theMat[i,:].tolist()) for i in range(0,analyze_count)]
    MatParams = {'ncol':analyze_count, 'dimnames':r.list(names,names),'byrow':True}
    Ufiltered =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)
    

    print("   * Creating PNG file of heatmap to \033[94m%s\033[0m" % filename)
    if analyze_count < 1000:
        r.png(file=folderName+"/"+filename,width=8000,height=8000,res=120)
    else:
        r.png(file=folderName+"/"+filename,width=12000,height=12000,res=120)
    if corType == 1:
        colscheme = input('  >> Color scheme: (1)green-black-red (2)teal-black-orange [2 default]? ')
        colscheme = r.colorpanel(75, "#0BB5FF","black","#F44A00")  if (colscheme=="2") | (colscheme=="") else r.colorpanel(75, "#14D9C2","black","#EB263D")
        heatParam = {'Rowv':dendroInfo,'Colv':dendroInfo,'col':colscheme,'ylab':"Ions",'xlab':"Ions","main":title,
            'margins':r.c(10,10),'symm':True,'symbreaks':True,'scale':'none','key':True,'symkey':True,'trace':"none",'keysize':0.7}

        heatmap(Ufiltered, **heatParam)

    elif corType == 2:
        dissimilarity = ((r.abs(Ufiltered)).ro*-1).ro+1
        colscheme = r.colorpanel(75, "#660000","white")
        heatParam = {'Rowv':dendroInfo,'Colv':dendroInfo,'col':colscheme,'ylab':"Ions",'xlab':"Ions","main":title,
            'margins':r.c(10,10),'symm':True,'scale':'none','key':True,'trace':"none",'keysize':0.7}
    
        heatmap(dissimilarity, **heatParam)


    devoff = robjects.r('dev.off')
    devoff()

def dendroCreator(the_data,dataN,analyze_count,names):
    sys.stdout.write("  Calculating hierarchical clustering...")
    sys.stdout.flush()
    r = robjects.r
    hdist = robjects.r('as.dist')
    asdendro = robjects.r('as.dendrogram')
    devoff = robjects.r('dev.off')

    Rread = []
    [Rread.extend([thing if thing!='NA' else robjects.NA_Real for thing in the_data[i][1:]]) for i in range(0,analyze_count)]

    MatParams = {'ncol':dataN, 'dimnames':r.list(r.c([the_data[i][0] for i in range(0,analyze_count)]),r.c(names)),'byrow':True}
    DataMat =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)
    DataMat = r.t(DataMat) #I know this is lazy but it's just easier!

    U = r.cor(DataMat, use="pairwise.complete.obs")
    #print U
    #dissimilarity = ((U.ro*-1).ro+1).ro/2 #Rpy vector basic operations
    dissimilarity = ((r.abs(U)).ro*-1).ro+1

    distanceR = hdist(dissimilarity)
    #print distanceR
    hclusty = r.hclust(distanceR)

    sys.stdout.write("complete!\n")
    sys.stdout.flush()

    print("  Creating PNG file of hierarchical clustering dendrogram to \033[94mionclust.png\033[0m")
    if analyze_count < 1000: r.png(file=folderName+"/ionclust.png",width=10500,height=1700,res=84)
    else: r.png(file=folderName+"/ionclust.png",width=15500,height=1700,res=84)
    r.plot(hclusty,main="Hierarchical Clustering of Ions")

    devoff()

    #convert to dendrogram
    dendroInfo = asdendro(hclusty)

    return dendroInfo

def dualCorr_parallelized(the_dataA,AdataN,the_dataB,BdataN,analyze_count,sig,corcut):
    size = analyze_count
    SMP = cpu_count()
    corMatrix = np.zeros((analyze_count,analyze_count))
    sigDiffCor = np.zeros((analyze_count,analyze_count))
    ApcorMat = np.eye(analyze_count)
    BpcorMat = np.eye(analyze_count)
    qu = Queue()
    quCount = Queue()
    quWarning = Queue()
    totaltocalc = ((size**2)-size)/2
    progress = 0

    lossList = []
    gainList = []
    Asig = 0
    Bsig = 0
    sigGain = 0
    sigLoss = 0

    uhOh = 0
    
    if SMP <= 6:
        print("       %s logical processors detected! Spawning 4 slave processes" % SMP)
        e1 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/4)))/2)) #1/4
        e2 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/2)))/2)) #2/4
        e3 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((3*size*(size-1))/4)))/2)) #3/4

        p = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,0,e1,analyze_count,qu,quCount,quWarning))
        q = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e1,e2,analyze_count,qu,quCount,quWarning))
        r = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e2,e3,analyze_count,qu,quCount,quWarning))
        s = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e3,size,analyze_count,qu,quCount,quWarning))
        
    if (SMP > 6) & (SMP < 11):
        print("       %s logical processors detected! Spawning 6 slave processes" % SMP)
        e1 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/6)))/2))
        e2 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/3)))/2))
        e3 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/2)))/2))
        e4 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((2*size*(size-1))/3)))/2))
        e5 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((5*size*(size-1))/6)))/2))

        p = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,0,e1,analyze_count,qu,quCount,quWarning))
        q = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e1,e2,analyze_count,qu,quCount,quWarning))
        r = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e2,e3,analyze_count,qu,quCount,quWarning))
        s = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e3,e4,analyze_count,qu,quCount,quWarning))
        t = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e4,e5,analyze_count,qu,quCount,quWarning))
        u = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e5,size,analyze_count,qu,quCount,quWarning))

    if SMP >= 12:
        print("       %s logical processors detected! Spawning 10 slave processes" % SMP)
        e1 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/10)))/2))
        e2 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/5)))/2))
        e3 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((3*size*(size-1))/10)))/2))
        e4 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((2*size*(size-1))/5)))/2))
        e5 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((size*(size-1))/2)))/2))
        e6 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((3*size*(size-1))/5)))/2))
        e7 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((7*size*(size-1))/10)))/2))
        e8 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((4*size*(size-1))/5)))/2))
        e9 = int(math.floor(((2*size-1)-math.sqrt(((2*size-1)**2)-4*((9*size*(size-1))/10)))/2))


        p = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,0,e1,analyze_count,qu,quCount,quWarning))
        q = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e1,e2,analyze_count,qu,quCount,quWarning))
        r = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e2,e3,analyze_count,qu,quCount,quWarning))
        s = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e3,e4,analyze_count,qu,quCount,quWarning))
        t = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e4,e5,analyze_count,qu,quCount,quWarning))
        u = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e5,e6,analyze_count,qu,quCount,quWarning))
        v = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e6,e7,analyze_count,qu,quCount,quWarning))
        w = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e7,e8,analyze_count,qu,quCount,quWarning))
        x = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e8,e9,analyze_count,qu,quCount,quWarning))
        y = Process(target=smp_dualCorr,args=(the_dataA,AdataN,the_dataB,BdataN,e9,size,analyze_count,qu,quCount,quWarning))


    
    p.start()
    q.start()
    r.start()
    s.start()

    if SMP > 6:
        t.start()
        u.start()
    if SMP >= 12:
        v.start()
        w.start()
        x.start()
        y.start()

    itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()

    if SMP > 6:
        itsAlive = itsAlive | t.is_alive() | u.is_alive()

    if SMP >= 12:
        itsAlive = itsAlive | v.is_alive() | w.is_alive() | x.is_alive() | y.is_alive()
    
    while (itsAlive | (not qu.empty()) | (not quCount.empty()) | (not quWarning.empty())):
        if(qu.empty() != True):

            
            temp = qu.get()
            _pcorA = temp[3]
            _pcorA_pval = temp[4]
            _pcorB = temp[5]
            _pcorB_pval = temp[6]
            _diff_pval = temp[7]
            _Xname = temp[8]
            _Yname = temp[9]

            corMatrix[temp[1],temp[2]] = temp[0]
            corMatrix[temp[2],temp[1]] = temp[0]


            #plain correlation
            random.seed()
            randomator = random.random()/1000
            if _pcorA_pval < sig: ApcorMat[temp[1],temp[2]] = _pcorA; ApcorMat[temp[2],temp[1]] = _pcorA
            else: ApcorMat[temp[1],temp[2]] = randomator; ApcorMat[temp[2],temp[1]] = randomator
            if _pcorB_pval < sig: BpcorMat[temp[1],temp[2]] = _pcorB; BpcorMat[temp[2],temp[1]] = _pcorB
            else: BpcorMat[temp[1],temp[2]] = randomator; BpcorMat[temp[2],temp[1]] = randomator

            if _pcorA_pval < sig: Asig += 1
            if _pcorB_pval < sig: Bsig += 1
            
            atLeastOneSig = (_pcorA_pval < sig) | (_pcorB_pval < sig)
            crossTheLine = ((_pcorA > 0) & (_pcorB < 0)) | ((_pcorA < 0) & (_pcorB > 0))
            sigDiffCor[temp[1],temp[2]] = randomator; sigDiffCor[temp[2],temp[1]] = randomator

            if ((_diff_pval < sig) & (atLeastOneSig) & ((abs(_pcorA - _pcorB) > corcut))):

                if ((_pcorA_pval > sig) & (_pcorB_pval < sig)):
                    sigGain += 1
                    gainList.append([_Xname,_Yname,_pcorA,_pcorA_pval,_pcorB,_pcorB_pval,_diff_pval])
                    sigDiffCor[temp[1],temp[2]] = temp[0]; sigDiffCor[temp[2],temp[1]] = temp[0]
                if ((_pcorA_pval < sig) & (_pcorB_pval > sig)):
                    sigLoss += 1
                    lossList.append([_Xname,_Yname,_pcorA,_pcorA_pval,_pcorB,_pcorB_pval,_diff_pval])
                    sigDiffCor[temp[1],temp[2]] = temp[0]; sigDiffCor[temp[2],temp[1]] = temp[0]

                bothSig = ((_pcorA_pval < sig) & (_pcorB_pval < sig))
                if ((bothSig) & (not crossTheLine) & (abs(_pcorA) > abs(_pcorB))):
                    sigLoss += 1
                    lossList.append([_Xname,_Yname,_pcorA,_pcorA_pval,_pcorB,_pcorB_pval,_diff_pval])
                    sigDiffCor[temp[1],temp[2]] = temp[0]; sigDiffCor[temp[2],temp[1]] = temp[0]
                if  ((bothSig) & (not crossTheLine) & (abs(_pcorA) < abs(_pcorB))):
                    sigGain += 1
                    gainList.append([_Xname,_Yname,_pcorA,_pcorA_pval,_pcorB,_pcorB_pval,_diff_pval])
                    sigDiffCor[temp[1],temp[2]] = temp[0]; sigDiffCor[temp[2],temp[1]] = temp[0]

        if (quWarning.empty() != True):
            ohshit = quWarning.get()
            uhOh += 1

        if(quCount.empty() != True):
            temp2 = quCount.get()
            progress += temp2
            sys.stdout.write("\r       Progress: %s out of %s (%s%%)" % (progress,totaltocalc,100*(float(progress)/totaltocalc)))
            sys.stdout.flush()
            
        itsAlive = p.is_alive() | q.is_alive() | r.is_alive() | s.is_alive()
        if SMP > 6:
                itsAlive = itsAlive | t.is_alive() | u.is_alive()
        if SMP >= 12:
            itsAlive = itsAlive | v.is_alive() | w.is_alive() | x.is_alive() | y.is_alive()

    p.join()
    q.join()
    r.join()
    s.join()
    
    if SMP > 6:
        t.join()
        u.join()

    if SMP >= 12:
        v.join()
        w.join()
        x.join()
        y.join()
    
    qu.close()
    quCount.close()
    quWarning.close()

    if uhOh == 0:
        print("")
    else:
        print("\n       \033[91mWarning! %s incomplete calculations due to insufficient sample pairs! Strongly advise rerunning after increasing ion presence threshold and disabling outlier filtering!\033[0m" % uhOh)

    return corMatrix,sigDiffCor,ApcorMat,BpcorMat,Asig,Bsig,sigGain,sigLoss,gainList,lossList

def dualCorrelator(the_dataA,AdataN,the_dataB,BdataN,analyze_count,sig,dendroInfo):
    #robjects.packages.quiet_require("gplots")
    r = robjects.r
    #hdist = robjects.r('as.dist')
    #asdendro = robjects.r('as.dendrogram')

    corrWr_gains = csv.writer(open(folderName+'/correlation_SigGains.csv', 'w'), delimiter=',')
    corrWr_losses = csv.writer(open(folderName+'/correlation_SigLosses.csv', 'w'), delimiter=',')
    header = ["ion X","ion Y","ctrl set correlation", "ctrl correlation p-value","exp set correlation","exp correlation p-value", "differential correlation p-value"]
    corrWr_gains.writerow(header); corrWr_losses.writerow(header)


    print("   * Calculating differential correlation statistics (this can take a while)...")
    corcut = input("  >> Minimum correlation difference cutoff [0.1 default]? ")
    corcut = 0.1 if corcut == "" else float(corcut)
    baseFrame, sigDiffCor, AcorMat, BcorMat, Asig, Bsig, sigGain, sigLoss, gainList, lossList = dualCorr_parallelized(the_dataA,AdataN,the_dataB,BdataN,analyze_count,sig,corcut)

    corrWr_gains.writerows(gainList)
    corrWr_losses.writerows(lossList)


    #Rread = []
    #[Rread.extend(baseFrame[i,:].tolist()) for i in range(0,analyze_count)]
    RionNames = r.c([the_dataA[i][0] for i in range(0,analyze_count)])
    #MatParams = {'ncol':analyze_count, 'dimnames':r.list(RionNames,RionNames),'byrow':True}
    #CorMat =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)
    
    #print "   * making filtered ctrl correlation heatmap"
    #pureCorMap(AcorMat,RionNames,analyze_count,"ctrlCor_onlySig.png",dendroInfo,1)

    #print "   * making filtered exp correlation heatmap"
    #pureCorMap(BcorMat,RionNames,analyze_count,"expCor_onlySig.png",dendroInfo,1)

    print("     Significant correlations: ctrl - %s   exp - %s" % (Asig,Bsig))
    print("     Statistically significant correlation losses: %s (read out to \033[94mcorrelation_SigLosses.csv\033[0m)" % sigLoss)
    print("     Statistically significant correlation gains: %s (read out to \033[94mcorrelation_SigGains.csv\033[0m)" % sigGain)

    #clustChoice = raw_input("  >> Cluster differential correlation heatmap (1)internally (2)previous choice (ctrl or exp) [1 default]? ")
    #clustChoice = True if ((clustChoice=="") | (clustChoice=="1")) else dendroInfo

    #if clustChoice == True:
    #    Rread = []
    #    [Rread.extend(baseFrame[i,:].tolist()) for i in range(0,analyze_count)]
    #    MatParams = {'ncol':analyze_count, 'dimnames':r.list(RionNames,RionNames),'byrow':True}
    #    UforDist =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)
    #    Udist = r.dist(UforDist)
    #    internalClust = r.hclust(Udist)
    #    clustChoice = asdendro(internalClust)

    print("   * Creating full heatmap with internal clustering")
    pureCorMap(baseFrame,RionNames,analyze_count,"diffCor_internalClust.png",True,1,"Absolute Correlation Differential Heat Map")

    print("   * Creating full heatmap with user defined clustering")
    pureCorMap(baseFrame,RionNames,analyze_count,"diffCor_commonClust.png",dendroInfo,1,"Absolute Correlation Differential Heat Map")

    print("   * Creating statistically filtered heatmap with user defined clustering")
    pureCorMap(sigDiffCor,RionNames,analyze_count,"diffCor_onlySig.png",dendroInfo,1,"Statistically significant Absolute Correlation Differential Heat Map")


    
#--------------------------main script---------------------------#

if bioMode != 3:
    if transf_bool != 2:
        CTRL_data, EXP_data, CTRL_cat, EXP_cat, analyze_count, fisher_count, cSPRAW, eSPRAW, names = dataExtractor(ctrlRaw,expRaw,zeros_Tr,transf_bool,Zexcl_bool)
    else:
        CTRL_data, EXP_data, CTRL_cat, EXP_cat, analyze_count, fisher_count, cSPRAW, eSPRAW, names,CTRLlognormParams,EXPlognormParams = dataExtractor(ctrlRaw,expRaw,zeros_Tr,transf_bool,Zexcl_bool,norm_bool)
else:
    if transf_bool != 2:
        CTRL_data, EXP_data, CTRL_cat, EXP_cat, analyze_count, fisher_count, cSPRAW, eSPRAW, names, modeList, modeList_part = dataMergedExtractor(ctrlRaw,expRaw,zeros_Tr,transf_bool,Zexcl_bool)
    else:
        CTRL_data, EXP_data, CTRL_cat, EXP_cat, analyze_count, fisher_count, cSPRAW, eSPRAW, names,CTRLlognormParams,EXPlognormParams, modeList, modeList_part = dataMergedExtractor(ctrlRaw,expRaw,zeros_Tr,transf_bool,Zexcl_bool,norm_bool)

print("Complete!")

ctrlN = len(cSPRAW[0])-1
expN= len(eSPRAW[0])-1

#print "\nCorrelation analysis of data..."
#dualCorrelator(CTRL_data,ctrlN,EXP_data,expN,analyze_count,sig)
#print "Complete!"


if transf_bool != 2:
    if bioMode != 3:
        CTRLlognormParams = [[0,1] for _i in range(0,ctrlN)]
        EXPlognormParams = [[0,1] for _i in range(0,expN)]
    else:
        CTRLlognormParams = [[[0,0],[1,1]] for _i in range(0,ctrlN)]
        EXPlognormParams = [[[0,0],[1,1]] for _i in range(0,expN)]

#mainly for Hotelling's testing filtering
sig_list = set()
CONTsig_list = []

if min(ctrlN,expN) < 15:
    enoughBool = False
    if min(ctrlN,expN) > 5: print("\nSample size below 15! Only a subset of statistical tests will be conducted!")
    else: print("\n\033[91mSample size is exceptionaly low! All calculated statistical results are questionable! Highly suggest exiting program and going for ice cream instead!\033[0m")
    testtype = input(" > Test filter: (1)Welch's T-test (2)Mann-Whitney U test (3)either [2 default]? ")
    testtype = 2 if testtype=='' else int(testtype)
else:
    enoughBool = True
    print("\nSample size 15 or above! Full gamut of statistical tests will be conducted!")
    testtype = input(" > Test filter: (1)Mann-Whitney U test (2)K-S test (3)either [2 default]? ")
    testtype = 2 if testtype=='' else int(testtype)




statResult = [[CTRL_data[j][0]] for j in range(0,analyze_count)] #feed in the names
gfx_nu = 1

forGFX = [["ID"]]; forGFX[0].extend(names)
#print forGFX[0]
#forGFX[0].extend(["ctrl_samp"+str(i) for i in range(1,ctrlN+1)]); 
#forGFX[0].extend(["exp_samp"+str(i) for i in range(1,expN+1)]); 
forAna = [["ID"]]; forAna[0].extend(names)

if transf_bool!=1:
    tforAna = [["ID"]]; tforAna[0].extend(names)
    wrlog = csv.writer(open(folderName+'/transFilteredData.csv', 'w'), delimiter=',')

futureExcludeLog = csv.writer(open(folderName+'/MT_SigList.csv', 'w'), delimiter=',')

if bioMode != 3:
    if enoughBool ==True:
        if testtype==1: statResult.insert(0,["ion","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Ctrl Anderson-Darling Test","Exp Anderson-Darling Test","Student's T-Test","Welch's T-Test","Mann-Whitney U Test","K-S Test","Permuted Mann-Whitney U Test","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
        if testtype==2: statResult.insert(0,["ion","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Ctrl Anderson-Darling Test","Exp Anderson-Darling Test","Student's T-Test","Welch's T-Test","Mann-Whitney U Test","K-S Test","Permuted K-S Test","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
    else:
        if testtype==1: statResult.insert(0,["ion","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Student's T-Test","Welch's T-Test","Mann-Whitney U Test","Permuted Welch's Test","95%% CI hypothesis testing","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
        if testtype==2: statResult.insert(0,["ion","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Student's T-Test","Welch's T-Test","Mann-Whitney U Test","Permuted Mann-Whitney U Test","95%% CI hypothesis testing","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
else:

    if enoughBool ==True:
        if testtype==1: statResult.insert(0,["ion","Mode","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Ctrl Anderson-Darling Test","Exp Anderson-Darling Test","Student's T-Test","Welch's T-Test","Mann-Whitney U Test","K-S Test","Permuted Mann-Whitney U Test","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
        if testtype==2: statResult.insert(0,["ion","Mode","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Ctrl Anderson-Darling Test","Exp Anderson-Darling Test","Student's T-Test","Welch's T-Test","Mann-Whitney U Test","K-S Test","Permuted K-S Test","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
    else:
        if testtype==1: statResult.insert(0,["ion","Mode","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Student's T-Test","Welch's T-Test","Mode","Mann-Whitney U Test","Permuted Welch's Test","95%% CI hypothesis testing","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])
        if testtype==2: statResult.insert(0,["ion","Mode","HMDB identification","KEGG identification", "LIPIDMAPS identification","BioCyc identification",
            "Student's T-Test","Welch's T-Test","Mann-Whitney U Test","Permuted Mann-Whitney U Test","95%% CI hypothesis testing","Log2 fold change","Ctrl mean","Ctrl CI","","Ctrl sd","Exp mean","Exp CI","","Exp sd"])

r = robjects.r
r.options(warn=-1)
ttest = r('t.test')
ranktest = r('wilcox.test')

#permutation testing
perm_bool = input(" > Conduct Permutation testing for selected test filter [y/N]? ")
perm_bool = False if ((perm_bool=="") | (perm_bool.lower() == "n")) else True
if perm_bool == True:
    print("\nCommencing permutation testing...")
    perm_dict = permute(analyze_count,CTRL_data,EXP_data,ctrlN,expN,testtype,enoughBool)


if bioFilter_bool == True:
    KEGGIDwr = csv.writer(open(folderName+'/KEGG_CID.csv', 'w'), delimiter=',')
        

#conduct standard battery of statistical tests
print("Commencing battery of standard statistical tests for the %s ions with >%s%% presence in both sets..." % (analyze_count,100*(1-zeros_Tr)))

NonSig_volcano_data = [[],[]]
Sig_volcano_data = [[],[]]

#Pathway aggregation counting
pathAgg = {}
crapAgg = {}
BpathAgg = {}
BcrapAgg = {}

#go through each ion and do statistical tests and stuff
for i in range(0,analyze_count):

    tempCTRL = [val for val in CTRL_data[i][1:] if val!='NA']
    tempEXP = [val for val in EXP_data[i][1:] if val!='NA']
    
    CTRLmean = np.mean(tempCTRL)
    EXPmean = np.mean(tempEXP)

    HMDB_info = 'NA'
    KEGG_info = 'NA'
    LIPIDMAPS_info = 'NA'
    BIOCYC_info = 'NA'

    if bioMode == 3:
        chemMode = modeList[i]
        if chemMode == 2: chemString = "ESI-"
        else: chemString = "ESI+"

    if enoughBool == True: # for big sample cases
        CTRL_ADtest = normality(tempCTRL) 
        EXP_ADtest = normality(tempEXP)
        kstest_pval = kolsmir(tempCTRL,tempEXP)
    else:
        kstest_pval = 1

    CI_stat,ctrlCI,expCI,ctrlSD,expSD = confidenceInt(tempCTRL,tempEXP)

    ttest_pval,welch_pval = bothT(tempCTRL,tempEXP,ttest)
    #print "%s" % statResult[i+1][0]
    ranktest_pval = mannw(tempCTRL,tempEXP,ranktest)

    if enoughBool == True:
        if testtype==2: sigVal = kstest_pval 
        elif testtype==1: sigVal = ranktest_pval
        else: sigVal = min(ranktest_pval,kstest_pval)
    else:
        sigVal= ranktest_pval if testtype==2 else welch_pval
        if testtype==2: sigVal = ranktest_pval 
        elif testtype==1: sigVal = welch_pval
        else: sigVal = min(ranktest_pval,welch_pval)

    ID = statResult[i+1]
    if userFilter_List==False:
        if perm_bool == False: goGo = True if sigVal < sig else False
        else: goGo = True if perm_dict[ID[0]] < sig else False
    else:
        if perm_bool == False:
            goGo = True if ((ID[0] not in userFilter_List[0]) & (sigVal < sig)) else False
        else: 
            goGo = True if ((ID[0] not in userFilter_List[0]) & (perm_dict[ID[0]] < sig)) else False


    if goGo==False:
        if (transf_bool==2): fld_chg = math.log(math.exp(EXPmean-CTRLmean),2)
        else: fld_chg = math.log(EXPmean/CTRLmean,2)

        #if fld_chg > 1: print sigVal; print "----------------- means: %s %s" % (CTRLmean,EXPmean)

        NonSig_volcano_data[0].append(fld_chg)
        NonSig_volcano_data[1].append(-1*math.log(sigVal,10))

        if bioFilter_bool == True:
            #ID = statResult[i+1]
            if bioMode != 3: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID = bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,bioMode,adductchoose,mw_tol,crapAgg,BcrapAgg)
            else: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID = bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,chemMode,adductchoose[chemMode-1],mw_tol,crapAgg,BcrapAgg)

            HMDB_info = HMDB_ident + ''.join([",["+name+"]" for name in HMDBputative_names])
            KEGG_info = KEGG_ident + ''.join([",["+name+"]" for name in KEGGputative_names])
            LIPIDMAPS_info = str(len(LIPIDMAPSputative_names)) + ''.join([",{"+name[0]+": "+name[1]+" < "+name[3]+" < "+name[2]+"}" for name in LIPIDMAPSputative_names])
            BIOCYC_info = str(len(BIOCYCputative_names)) + ''.join([",["+name+"]" for name in BIOCYCputative_names])



    elif goGo==True:
        sig_list.add(ID[0])
        CONTsig_list.append(ID[0])

        if bioMode != 3:
            if bioFilter_bool == True: print("\n    Ion %s is significant!" % ID[0])
        else: 
            if bioFilter_bool == True: print("\n    Ion %s is significant! %s" % (ID[0],chemString))

        if (transf_bool==2): fld_chg = math.log(math.exp(EXPmean-CTRLmean),2)
        else: fld_chg = math.log(EXPmean/CTRLmean,2)

        Sig_volcano_data[0].append(fld_chg)
        Sig_volcano_data[1].append(-1*math.log(sigVal,10))

        #output for PCA without logging
        forGFX.append([CTRL_data[i][0]])
        cTEMP = [float(j) for j in cSPRAW[i][1:]]
        eTEMP = [float(j) for j in eSPRAW[i][1:]]
        forGFX[gfx_nu].append(cTEMP)
        forGFX[gfx_nu][1].extend(eTEMP)

        if transf_bool!=1:
            #potential processing for PCA with log-type transforming and filling in zeros
            allmean = (len(tempCTRL)*CTRLmean+len(tempEXP)*EXPmean)/(len(tempCTRL)+len(tempEXP))

            if bioMode != 3:
                ctTEMPnorm = [(math.log(float(cSPRAW[i][j+1]))-CTRLlognormParams[j][0])/CTRLlognormParams[j][1] if float(cSPRAW[i][j+1])!=0 else allmean for j in range(0,ctrlN)]#lognorm
                etTEMPnorm = [(math.log(float(eSPRAW[i][j+1]))-EXPlognormParams[j][0])/EXPlognormParams[j][1] if float(eSPRAW[i][j+1])!=0 else allmean for j in range(0,expN)]
            else:
                ctTEMPnorm = [(math.log(float(cSPRAW[i][j+1]))-CTRLlognormParams[j][0][chemMode-1])/CTRLlognormParams[j][1][chemMode-1] if float(cSPRAW[i][j+1])!=0 else allmean for j in range(0,ctrlN)]#lognorm
                etTEMPnorm = [(math.log(float(eSPRAW[i][j+1]))-EXPlognormParams[j][0][chemMode-1])/EXPlognormParams[j][1][chemMode-1] if float(eSPRAW[i][j+1])!=0 else allmean for j in range(0,expN)]

            forGFX[gfx_nu].append(ctTEMPnorm)
            forGFX[gfx_nu][2].extend(etTEMPnorm)

            ctTEMPsinh = [math.log(float(j)+math.sqrt(float(j)**2+1)) for j in cSPRAW[i][1:]]#IHC
            etTEMPsinh = [math.log(float(j)+math.sqrt(float(j)**2+1)) for j in eSPRAW[i][1:]]
            forGFX[gfx_nu].append(ctTEMPsinh)
            forGFX[gfx_nu][3].extend(etTEMPsinh)



        #---output for CSV logged filtered data with outlier IDing and for plain filtered data---
        if transf_bool!=1: tforAna.append([CTRL_data[i][0]]); cATEMP = []

        forAna.append([CTRL_data[i][0]]); cNTEMP = []

        lognormIndex = 0
        for j in cSPRAW[i][1:]:
            lnum = float(j)
            if lnum!=0:
                if transf_bool==2:
                    if bioMode != 3: lnnum = (math.log(lnum)-CTRLlognormParams[lognormIndex][0])/CTRLlognormParams[lognormIndex][1]
                    else:  lnnum = (math.log(lnum)-CTRLlognormParams[lognormIndex][0][chemMode-1])/CTRLlognormParams[lognormIndex][1][chemMode-1]
                elif transf_bool==3: lnnum = math.log(lnum+math.sqrt(lnum**2+1))
                else: lnnum = lnum

                if lnnum in tempCTRL:
                    if transf_bool!=1: cATEMP.append(lnnum)
                    cNTEMP.append(j)
                else:
                    if transf_bool!=1: cATEMP.append('>>'+str(lnnum)+"<<")
                    cNTEMP.append('>>'+str(j)+'<<')
            else:
                if (transf_bool==2): cATEMP.append('NA')
                elif transf_bool==3:
                    if lnum in tempCTRL: cATEMP.append('0')
                    else: cATEMP.append('>>0<<')
                
                if (Zexcl_bool==True) | (transf_bool==2): cNTEMP.append('<<0>>')
                else:
                    if lnum in tempCTRL: cNTEMP.append('0')
                    else: cNTEMP.append('>>0<<')
            lognormIndex += 1

        lognormIndex = 0
        eATEMP = []; eNTEMP = []
        for j in eSPRAW[i][1:]:
            lnum = float(j)
            if lnum!=0:
                if transf_bool==2:
                    if bioMode != 3:
                        lnnum = (math.log(lnum)-EXPlognormParams[lognormIndex][0])/EXPlognormParams[lognormIndex][1]
                    else:
                        lnnum = (math.log(lnum)-EXPlognormParams[lognormIndex][0][chemMode-1])/EXPlognormParams[lognormIndex][1][chemMode-1]
                elif transf_bool==3: lnnum = math.log(lnum+math.sqrt(lnum**2+1))
                else: lnnum = lnum

                if lnnum in tempEXP:
                    if transf_bool!=1: eATEMP.append(lnnum)
                    eNTEMP.append(j)
                else:
                    if transf_bool!=1: eATEMP.append('>>'+str(lnnum)+"<<")
                    eNTEMP.append('>>'+str(j)+'<<')
            else:
                if (transf_bool==2): eATEMP.append('NA')
                elif transf_bool==3:
                    if lnum in tempEXP: eATEMP.append('0')
                    else: eATEMP.append('>>0<<')
                
                if (Zexcl_bool==True) | (transf_bool==2): eNTEMP.append('<<0>>')
                else:
                    if lnum in tempEXP: eNTEMP.append('0')
                    else: eNTEMP.append('>>0<<')
            lognormIndex += 1

        if transf_bool!=1: tforAna[gfx_nu].extend(cATEMP); tforAna[gfx_nu].extend(eATEMP)
        forAna[gfx_nu].extend(cNTEMP); forAna[gfx_nu].extend(eNTEMP)
        #--------------

        if bioFilter_bool == True:
            #ID = statResult[i+1]
            if bioMode != 3: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID= bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,bioMode,adductchoose,mw_tol, pathAgg,BpathAgg, True)
            else:  HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID= bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,chemMode,adductchoose[chemMode-1],mw_tol, pathAgg,BpathAgg, True)

            if len(straightKEGGID)!=0: straightKEGGID.insert(0,statResult[i+1][0]); KEGGIDwr.writerow(straightKEGGID)

            HMDB_info = HMDB_ident + ''.join([",["+name+"]" for name in HMDBputative_names])
            KEGG_info = KEGG_ident + ''.join([",["+name+"]" for name in KEGGputative_names])
            LIPIDMAPS_info = str(len(LIPIDMAPSputative_names)) + ''.join([",{"+name[0]+": "+name[1]+" < "+name[3]+" < "+name[2]+"}" for name in LIPIDMAPSputative_names])
            BIOCYC_info = str(len(BIOCYCputative_names)) + ''.join([",["+name+"]" for name in BIOCYCputative_names])

        gfx_nu += 1


    if perm_bool == True:
        perm_pval = perm_dict[ID[0]]
    else:
        perm_pval = 'NA'

    if enoughBool == True:
        if bioMode == 3: statResult[i+1].extend([chemString,HMDB_info,KEGG_info,LIPIDMAPS_info,BIOCYC_info,CTRL_ADtest,EXP_ADtest,ttest_pval,welch_pval,ranktest_pval,kstest_pval,perm_pval,fld_chg,CTRLmean,ctrlCI[0],ctrlCI[1],ctrlSD,EXPmean,expCI[0],expCI[1],expSD])
        else: statResult[i+1].extend([HMDB_info,KEGG_info,LIPIDMAPS_info,BIOCYC_info,CTRL_ADtest,EXP_ADtest,ttest_pval,welch_pval,ranktest_pval,kstest_pval,perm_pval,fld_chg,CTRLmean,ctrlCI[0],ctrlCI[1],ctrlSD,EXPmean,expCI[0],expCI[1],expSD])
    else:
        if bioMode == 3: statResult[i+1].extend([chemString,HMDB_info,KEGG_info,LIPIDMAPS_info,BIOCYC_info,ttest_pval,welch_pval,ranktest_pval,perm_pval,CI_stat,fld_chg,CTRLmean,ctrlCI[0],ctrlCI[1],ctrlSD,EXPmean,expCI[0],expCI[1],expSD])
        else: statResult[i+1].extend([HMDB_info,KEGG_info,LIPIDMAPS_info,BIOCYC_info,ttest_pval,welch_pval,ranktest_pval,perm_pval,CI_stat,fld_chg,CTRLmean,ctrlCI[0],ctrlCI[1],ctrlSD,EXPmean,expCI[0],expCI[1],expSD])

print("Complete! %s statistically significant complete presence ions identified" % (gfx_nu-1))
pathSig_P = gfx_nu-1

#volcano plot output
print("Outputting volcano plot (close window to continue)...")
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(NonSig_volcano_data[0],NonSig_volcano_data[1],s=20,marker='o',c="0.75")
ax1.scatter(Sig_volcano_data[0],Sig_volcano_data[1],s=20,marker='o',c='r')
plt.xlabel(r'log2 fold change', fontsize=10)
plt.ylabel(r'-log10 p-value', fontsize=10)
plt.title('p-value vs. fold change')
plt.ylim([0,max(max(NonSig_volcano_data[1]),max(Sig_volcano_data[1]))])
plt.grid(True)
plt.show()


print("\nConducting Fisher's exact test of categorical presence for the %s ions with >%s%% presence in only one set..." % (fisher_count,100*(1-zeros_Tr)))

catResult = [[CTRL_cat[j][0]] for j in range(0,fisher_count)]
if bioMode != 3: catResult.insert(0,["ion","HMDB identification","KEGG identification","LIPIDMAPS identification","BioCyc identification","Ctrl presence %","Exp presence %","Fisher's exact test",""])
else: catResult.insert(0,["ion","Mode","HMDB identification","KEGG identification","LIPIDMAPS identification","BioCyc identification","Ctrl presence %","Exp presence %","Fisher's exact test",""])
catResult[0].extend(names)

fisherexact_test = robjects.r('fisher.test')

CATsig = 0
CATsig_list = []

for i in range(0,fisher_count):
    tempCTRL = [1 if abun > 0 else 0 for abun in CTRL_cat[i][1:]]
    tempEXP = [1 if abun > 0 else 0 for abun in EXP_cat[i][1:]]

    HMDB_info = 'NA'
    KEGG_info = 'NA'
    LIPIDMAPS_info = 'NA'
    BIOCYC_info = 'NA'


    if bioMode == 3:
        chemMode = modeList_part[i]
        if chemMode == 2: chemString = "ESI-"
        else: chemString = "ESI+"


    fe_pval,CTRLpres,EXPpres = fishtest(tempCTRL,tempEXP,r,fisherexact_test)

    ID = catResult[i+1]
    if userFilter_List==False:
        goGo = True if fe_pval < sig else False
    else:
        goGo = True if ((ID[0] not in userFilter_List[1]) & (fe_pval < sig)) else False

    if (goGo==True):
        CATsig += 1
        CATsig_list.append(ID[0])

        if bioMode != 3:
            if bioFilter_bool == True: print("\n    Ion %s is significant!" % ID[0])
        else: 
            if bioFilter_bool == True: print("\n    Ion %s is significant! %s" % (ID[0],chemString))

        if bioFilter_bool == True:
            if bioMode != 3: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID = bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,bioMode,adductchoose,mw_tol,pathAgg,BpathAgg, True)
            else: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID = bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,chemMode,adductchoose[chemMode-1],mw_tol,pathAgg,BpathAgg, True)

            HMDB_info = HMDB_ident + ''.join([",["+name+"]" for name in HMDBputative_names])
            KEGG_info = KEGG_ident + ''.join([",["+name+"]" for name in KEGGputative_names])
            LIPIDMAPS_info = str(len(LIPIDMAPSputative_names)) + ''.join([",{"+name[0]+": "+name[1]+" < "+name[3]+" < "+name[2]+"}" for name in LIPIDMAPSputative_names])
            BIOCYC_info = str(len(BIOCYCputative_names)) + ''.join([",["+name+"]" for name in BIOCYCputative_names])

            if len(straightKEGGID)!=0: straightKEGGID.insert(0,catResult[i+1][0]); KEGGIDwr.writerow(straightKEGGID)

    elif bioFilter_bool == True:
        ID = catResult[i+1]
        if bioMode != 3: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID = bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,bioMode,adductchoose,mw_tol,crapAgg,BcrapAgg)
        else: HMDB_ident,HMDBputative_names, KEGG_ident,KEGGputative_names,straightKEGGID, LIPIDMAPSputative_names, BIOCYCputative_names,StraightBIOCYCID = bioIdentifier(pathways,ID,bioFilter_list,keggFilter_list,lipidFilter_list,biocycFilter_list,chemMode,adductchoose[chemMode-1],mw_tol,crapAgg,BcrapAgg)

        HMDB_info = HMDB_ident + ''.join([",["+name+"]" for name in HMDBputative_names])
        KEGG_info = KEGG_ident + ''.join([",["+name+"]" for name in KEGGputative_names])
        LIPIDMAPS_info = str(len(LIPIDMAPSputative_names)) + ''.join([",{"+name[0]+": "+name[1]+" < "+name[3]+" < "+name[2]+"}" for name in LIPIDMAPSputative_names])
        BIOCYC_info = str(len(BIOCYCputative_names)) + ''.join([",["+name+"]" for name in BIOCYCputative_names])
    
    if bioMode == 3: catResult[i+1].extend([chemString,HMDB_info,KEGG_info,LIPIDMAPS_info,BIOCYC_info,CTRLpres,EXPpres,fe_pval,""])
    else: catResult[i+1].extend([HMDB_info,KEGG_info,LIPIDMAPS_info,BIOCYC_info,CTRLpres,EXPpres,fe_pval,""])
    catResult[i+1].extend(CTRL_cat[i][1:])
    catResult[i+1].extend(EXP_cat[i][1:])

print("Complete! %s statistically significant partial presence ions identified" % CATsig)
pathSig_P += CATsig
pathSig_P = float(pathSig_P)/(analyze_count+fisher_count)


def gbar(ax, x, y, width=0.7, bottom=0):
   X = [[.7, .7],[.9,.9]]
   for left,top in zip(x, y):
        right = left+width
        ax.imshow(X, interpolation='bicubic', cmap=cm.PuBu_r,
                 extent=(left, right, bottom, top), alpha=1)

if bioFilter_bool==True:
    print("Outputting KEGG metabolic pathway hits (close window to continue)...")
    pathCount = [[wut[0]+": "+wut[1],len(pathAgg[wut]), wut] for wut in pathAgg]
    pathCountSort = sorted(pathCount, key=lambda psig: psig[1], reverse=True)

    ind = []; num = 1; ybar=[]; group_labels = []
    ybarSig = []
    for pathies in pathCountSort[1:]:
        if pathies[1]>2:
            num += 1
            ind.append(num)
            ybar.append(pathies[1])
            group_labels.append(pathies[0])
            if pathies[2] in crapAgg: unsigCount = len(crapAgg[pathies[2]])
            else: unsigCount = 0

            lambda_val = (unsigCount + pathies[1]) * pathSig_P
            #pathsigPval = 1-poisson.cdf(pathies[1]-1,lambda_val)
            Binompathsig = 1-binom.cdf(pathies[1]-1,unsigCount + pathies[1], pathSig_P)
            #ybarSig.append(-1*math.log10(pathsigPval))
            ybarSig.append(-1*math.log10(Binompathsig))
            #else: ybarSig.append(-1*math.log(0.0001))

            #print "%s sig: %s   unsig: %s    lambda: %s   poisson p-val: %s binomial p-val: %s" % (pathies[0], pathies[1], unsigCount, lambda_val, pathsigPval, Binompathsig)


    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.bar(ind,ybar,facecolor='#777777', align='center')
    ax2.set_ylabel("Putative Metabolite Count")
    ax2.set_xlabel("Pathway")
    ax2.set_title("Putative Metabolite Hits by KEGG Pathway")
    ax2.set_xticks(ind)
    ax2.set_xticklabels(group_labels, fontsize=9)
    fig2.autofmt_xdate()
    plt.show()

    #reformat the data for proper CSV output
    AggOutput = [["KEGG Pathway ID","Description","ions"]]
    for ag in pathAgg:
        templine = [ag[0],ag[1]]
        templine.extend(pathAgg[ag])
        AggOutput.append(templine)

    #print AggOutput
    wrp = csv.writer(open(folderName+'/KEGGpathlist.csv', 'w'), delimiter=",")

    print("Outputting KEGG metabolic pathways significance (close window to continue)...")

    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111)
    ax5.bar(ind,ybarSig,facecolor='#0066cc', align='center')
    #gbar(ax5, ind, ybarSig)
    ax5.set_title("Putative Pathway Enrichment Significance")
    ax5.set_xlabel("Pathway")
    ax5.set_ylabel("P-Value (-1*log)")
    l = plt.axhline(y=-1*math.log10(sig),linewidth=2, color="#FF6600")
    #ax5.grid()
    ax5.set_xticks(ind)
    ax5.set_xticklabels(group_labels, fontsize=9)
    fig5.autofmt_xdate()
    ax5.set_aspect('auto')
    plt.show()

    #--------------------------------
    
    print("Outputting BioCyc metabolic pathway hits (close window to continue)...")
    BpathCount = [[wut[0]+": "+wut[1],len(BpathAgg[wut]), wut] for wut in BpathAgg]
    BpathCountSort = sorted(BpathCount, key=lambda psig: psig[1], reverse=True)

    ind = []; num = 1; ybar=[]; Bgroup_labels = []
    ybarSig = []
    for pathies in BpathCountSort[1:]:
        if pathies[1]>2:
            num += 1
            ind.append(num)
            ybar.append(pathies[1])
            Bgroup_labels.append(html.fromstring(pathies[0]).text_content())
            if pathies[2] in BcrapAgg: unsigCount = len(BcrapAgg[pathies[2]])
            else: unsigCount = 0

            lambda_val = (unsigCount + pathies[1]) * pathSig_P
            #pathsigPval = 1-poisson.cdf(pathies[1]-1,lambda_val)
            Binompathsig = 1-binom.cdf(pathies[1]-1,unsigCount + pathies[1], pathSig_P)
            #ybarSig.append(-1*math.log10(pathsigPval))
            ybarSig.append(-1*math.log10(Binompathsig))

    fig6 = plt.figure()
    ax6 = fig6.add_subplot(111)
    ax6.bar(ind,ybar,facecolor='#777777', align='center')
    ax6.set_ylabel("Putative Metabolite Count")
    ax6.set_xlabel("Pathway")
    ax6.set_title("Putative Metabolite Hits by BioCyc Pathway")
    ax6.set_xticks(ind)
    ax6.set_xticklabels(Bgroup_labels, fontsize=9)
    fig6.autofmt_xdate()
    plt.show()

    #reformat the data for proper CSV output
    BAggOutput = [["BioCyc Pathway ID","Description","ions"]]
    for ag in BpathAgg:
        templine = [ag[0],ag[1]]
        templine.extend(BpathAgg[ag])
        BAggOutput.append(templine)

    #print AggOutput
    wrq = csv.writer(open(folderName+'/BIOCYCpathlist.csv', 'w'), delimiter=",")

    print("Outputting BioCyc metabolic pathways significance (close window to continue)...")

    fig7 = plt.figure()
    ax7 = fig7.add_subplot(111)
    ax7.bar(ind,ybarSig,facecolor='#0066cc', align='center')
    #gbar(ax7, ind, ybarSig)
    ax7.set_title("Putative Pathway Enrichment Significance")
    ax7.set_xlabel("Pathway")
    ax7.set_ylabel("P-Value (-1*log)")
    l = plt.axhline(y=-1*math.log10(sig),linewidth=2, color="#FF6600")
    #ax7.grid()
    ax7.set_xticks(ind)
    ax7.set_xticklabels(Bgroup_labels, fontsize=9)
    fig7.autofmt_xdate()
    ax7.set_aspect('auto')
    plt.show()





wrg = csv.writer(open(folderName+'/bothPresent_stats.csv', 'w'), delimiter=',')
wrf = csv.writer(open(folderName+'/onePresent_stats.csv', 'w'), delimiter=',')
wrh = csv.writer(open(folderName+'/rawFilteredData.csv', 'w'), delimiter=',')


wrg.writerows(statResult)
wrf.writerows(catResult)


#selective bivariate Hotelling's T-squared statistical analysis
if (enoughBool == True):
    multBool = input(" > Conduct Hotelling's T2-based bivariate testing [y/N]? ")
    multBool = False if ((multBool=="") | (multBool.lower()=="n")) else True
    if multBool == True:
        print("\nCommencing selective bivariate statistical analysis...")
        multivariateTesting(analyze_count,CTRL_data,EXP_data,ctrlN,expN,sig_list,r)

print("\nAll analyses complete! Results written out to CSV formatted files \033[94mbothPresent_stats.csv\033[0m and \033[94monePresent_stats.csv\033[0m\n")


repeatR=True
PCAtime=0
while((repeatR==True) & (PCAtime!=5)):
    PCAtime = input(" > Proceed with (1)standard PCA (2)kernel PCA (3)MDS (4)ICA for significant ions or (5)skip entirely [1 default]? ")
    PCAtime = 1 if (PCAtime=="") else int(PCAtime)

    totSamp = ctrlN + expN

    if PCAtime != 5:
        if transf_bool!=1:
            if norm_bool==3: transfPCA = input("  >> Data visualization: (1)raw (2)KDEMAX transnormalized (3)inverse hyperbolic sine transformed [2 default]? ")
            elif norm_bool==2: transfPCA = input("  >> Data visualization: (1)raw (2)std-gaussian transnormalized (3)inverse hyperbolic sine transformed [2 default]? ")
            elif (norm_bool=="NA") | (norm_bool==1): transfPCA = input("  >> Data visualization: (1)raw (2)log transformed (3)inverse hyperbolic sine transformed [%s default]? " % transf_bool)
            transfPCA = transf_bool if (transfPCA=="") else int(transfPCA)
        else: transfPCA = transf_bool


        Rread = []
        [Rread.extend(forGFX[i][transfPCA]) for i in range(1,gfx_nu)]


        #MatParams = {'nrow':totSamp,'col.names':r.c(forGFX[0][1:])}

        #MatParams = {'nrow':totSamp, 'dimnames':r.list(r.c(forGFX[0][1:]),r.c([forGFX[i][0] for i in range(1,gfx_nu)]))}
        DataMat =  r.matrix(r.c(robjects.FloatVector(Rread)),nrow=totSamp)
        #DataMat =  r.matrix(r.c(robjects.FloatVector(Rread)),**MatParams)


        if PCAtime == 1:
            print("\nConducting PCA via singular value decomposition...")

            pca = r.prcomp

            pca_params = {'retx': True, 'center': True, 'scale.': True}

            mydata_pca = pca(DataMat, **pca_params)

            scores = mydata_pca.rx2('x')

            pc1 = [scores[i] for i in range(0,totSamp)]
            pc2 = [scores[i] for i in range(totSamp,totSamp*2)]
            pc3 = [scores[i] for i in range(totSamp*2,totSamp*3)]

        elif PCAtime==2:

            kern = importr('kernlab')
            slot = r.slot

            pc1,pc2,pc3 = kPCA(DataMat,totSamp,kern,slot,r)

        elif PCAtime==3:
            print("\nConducting MDS...")
            mType = input("  >> Utilize (1)euclidean metric (2)Kruskal's non-metric [1 default]? ")
            mType = 1 if ((mType=='') | (mType == 1)) else int(mType)

            dist = r.dist
            cmdscale = r.cmdscale
            
            MASS = importr('MASS')
            isoMDS = MASS.isoMDS

            pc1,pc2,pc3 = MDS(DataMat,totSamp,dist,cmdscale,isoMDS,mType)

            #r.plot(pc1,pc2,xlab="one",ylab="two")

        elif PCAtime==4:
            print("\nConducting ICA...")

            ICAlib = importr('fastICA')
            fastICA = ICAlib.fastICA

            pc1,pc2,pc3 = ICA(DataMat,totSamp,fastICA)


        print("Outputting 3D projections plot (close window when finished)...") #matplotlib stuff here
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        ax.scatter(pc1[0:ctrlN], pc2[0:ctrlN], pc3[0:ctrlN], c='r', marker='o',s=60)
        ax.scatter(pc1[ctrlN:], pc2[ctrlN:], pc3[ctrlN:], c='b', marker='^',s=60)

        if (PCAtime!=3) & (PCAtime!=4):
            ax.set_xlabel('Principal Component 1')
            ax.set_ylabel('Principal Component 2')
            ax.set_zlabel('Principal Component 3')
        else:
            ax.set_xlabel('Dimension 1')
            ax.set_ylabel('Dimension 2')
            ax.set_zlabel('Dimension 3')

        plt.show()
        repeatR = input(" > Repeat multidimensional visualization analysis [Y/n]? ")
        repeatR = True if (repeatR == "") | (repeatR.lower()=="y") else False

#R heatmap creation
print("\nGenerating heatmap for significant ions...")
heatMapMaker(forGFX,gfx_nu,totSamp,transf_bool,ctrlN,expN,r)
print("Complete!")

if enoughBool==True:
    goThru = input("\nSample size is sufficient for correlation heatmap creation! Proceed [Y/n]? ")
    goThru = True if (goThru == "") | (goThru.lower()=="y") else False
else:
    if min(ctrlN,expN) >= 5:
        goThru = input("\n\033[91mSample size below recommended threshold for correlation heatmap creation! This may result in errors and/or unreliable results\033[0m\n > Proceed anyway (don't blame me if things break!!) [y/N]? ")
        goThru = False if (goThru == "") | (goThru.lower()=="n") else True
    else:
        goThru = False
        print("Sample size is %s! Don't even think about doing correlation related statistics!" % min(ctrlN,expN))

if goThru == True:
    print("\nCommence correlation based calculations and heatmap generation...")

    clustChoice = input(" > Hierarchical clustering by (1)ctrl (2)exp data [1 default]? ")
    clustChoice = 1 if clustChoice=="" else int(clustChoice)
    if clustChoice == 1:
        dendroInfo = dendroCreator(CTRL_data,ctrlN,analyze_count,names[:ctrlN])
    else:
        dendroInfo = dendroCreator(EXP_data,expN,analyze_count,names[ctrlN:])

    corType = input(" > Heatmap type: (1)raw Pearson's correlation (2)dissimilarity [1 default]? ")
    corType = 1 if corType=="" else int(corType)

    if corType==1: print("  Generating ctrl correlation heatmap"); pngname = "CorrHeatmap_ctrl.png"
    else: print("  Generating ctrl dissimilarity heatmap"); pngname = "dissHeatmap_ctrl.png"
    correlationHeatMapMaker(CTRL_data,names[:ctrlN],ctrlN,analyze_count,pngname, dendroInfo, corType)

    if corType==1: print("  Generating exp correlation heatmap"); pngname = "CorrHeatmap_exp.png"
    else: print("  Generating exp dissimilarity heatmap"); pngname = "dissHeatmap_exp.png"
    correlationHeatMapMaker(EXP_data,names[ctrlN:],expN,analyze_count,pngname, dendroInfo, corType)

    if min(ctrlN,expN) > 6:
        goThru2 = input(" > Proceed with statistical analysis of correlation significance (highly discouraged if sample size is small) [Y/n]? ")
        goThru2 = True if (goThru2 == "") | (goThru2.lower()=="y") else False
    else:
        print("  Correlation significance analysis bypassed due to low sample count")
        goThru2 = False

    if goThru2 == True:
        print("  Conducting differential correlation analysis...")
        dualCorrelator(CTRL_data,ctrlN,EXP_data,expN,analyze_count,sig,dendroInfo)
    print("")

if transf_bool!=1:
    print("Reading out transformed data into \033[94mtransFilteredData.csv\033[0m")
    wrlog.writerows(tforAna)


print("Reading out filtered data for significant ions to \033[94mrawFilteredData.csv\033[0m")

wrh.writerows(forAna)

print("Reading out significant M_T list to \033[94mMT_SigList.csv\033[0m")
futureExcludeLog.writerow(["complete presence significant ions","partial presence significant ions"])
for i in range(0,max(len(CONTsig_list),len(CATsig_list))):
    if ((i+1 < len(CONTsig_list)) & (i+1 < len(CATsig_list))):
        tempRow = [CONTsig_list[i],CATsig_list[i]]
    elif((i+1 < len(CONTsig_list)) & (i+1 > len(CATsig_list))):
        tempRow = [CONTsig_list[i],""]
    elif((i+1 > len(CONTsig_list)) & (i+1 < len(CATsig_list))):
        tempRow = ["",CATsig_list[i]]

    futureExcludeLog.writerow(tempRow)

if bioFilter_bool==True:
    print("Reading out KEGG and BioCyc pathway aggregate list to \033[94mKEGGpathlist.csv and BIOCYCpathlist.csv\033[0m")
    wrp.writerows(AggOutput)
    wrq.writerows(BAggOutput)

print("\n***** Please change filenames before the next run to avoid overwriting results! *****\n")

