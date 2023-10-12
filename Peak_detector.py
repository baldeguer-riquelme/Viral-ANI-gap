import warnings
warnings.filterwarnings( "ignore" )
import numpy as np
import pandas as pd
from scipy.signal import find_peaks,find_peaks_cwt, savgol_filter
from scipy.stats import gaussian_kde, normaltest
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from diptest import diptest
import itertools
import os
import multiprocessing


#Function to detect valleys
def get_valleys(valley,pivot):
    list_valleys=[]
    for i in np.arange(0,len(valley)):
        list_valleys.append(round(pivot['ANI'][valley[i]],1))
    return(list_valleys)


#Function to classify the species
def Group(gap_valleys_smooth, gap_peaks_smooth, pval_dip, pval_norm, average_ANI):
    if isinstance(pval_norm, float) == True and isinstance(pval_dip,float) == True: # Try first with normal/unimodal data
        if gap_valleys_smooth == True and pval_dip < 0.05 and pval_norm < 0.05 and average_ANI < 99.8: # Species that show an ANI gap (valley) at 99.3 - 99.8
            group="Group 1" 
        elif average_ANI >= 99.8: # Too clonal to say
            group="Group 2" 
        elif gap_valleys_smooth == False and gap_peaks_smooth == True: # Species that peak at the gap (99.3 - 99.8)
            group="Group 4"
        elif (gap_valleys_smooth == False and gap_peaks_smooth == False) or (pval_dip >= 0.05) or (pval_norm >= 0.05): # Species with a pattern that is compatible with the existence of the gap
            group="Group 3"
    else: # if normal/unimodal data was not calculated, proceed without these data
        if gap_valleys_smooth == True and average_ANI < 99.8: # Species that show an ANI gap (valley) at 99.3 - 99.8
            group="Group 1" 
        elif average_ANI >= 99.8: # Too clonal to say
            group="Group 2" 
        elif gap_valleys_smooth == False and gap_peaks_smooth == False: # Species with a pattern that is compatible with the existence of the gap
            group="Group 3"
        elif gap_valleys_smooth == False and gap_peaks_smooth == True: # Species that peak at the gap (99.3 - 99.8)
            group="Group 4"
    return(group)


#Function to analyze each species
def analyze_sp(species, table, pdf):
    table_sp=table[table['Species'] == species]
    round_df=table_sp['ANI'].round(1)
    average_ANI=round_df.mean()
    pivot=round_df.value_counts().to_frame()
    pivot.columns=['Freq']
    pivot['ANI']=pivot.index
    pivot=pivot.sort_values(by='ANI', ascending=False)
    pivot=pivot.reset_index()
    range=np.arange(pivot['ANI'].min(),100.1,0.1)
    pivot_smooth=pd.DataFrame()
    pivot_smooth['ANI']=range
    pivot_smooth['Freq_smooth']=gaussian_kde(round_df,bw_method=0.18).evaluate(range) # kernel density. PROBABLY I SHOULD CHANGE "pivot['ANI']" BY RANGE DEFINED AS "range=np.arange(pivot['ANI'].min(),100,0.1)"
    #Hartigan's diptest of unimodality test and normality test
    if len(pivot_smooth['Freq_smooth']) >= 8:
        #Diptest
        range_diptest=np.arange(pivot['ANI'].min(),100.1,0.001) # Create range of ANI values
        df_diptest=pd.DataFrame()
        df_diptest['ANI']=range_diptest # Append values of ANI created above
        df_diptest['Freq_smooth']=gaussian_kde(round_df,bw_method=0.18).evaluate(range_diptest) # kernel density estimation
        df_diptest['Normalized']= df_diptest['Freq_smooth'] / df_diptest['Freq_smooth'].max() * 1000 # Normalize data 0 to 1. Then, multiply * 1000 to get ‰ (that is, per mille)
        arr=[]
        for i in np.arange(0,len(df_diptest['ANI'])):
            arr.extend(itertools.repeat(str(df_diptest['ANI'][i]), int(round(df_diptest['Normalized'][i],0)))) # Repeat each ANI values based on the normalized data obtained above. Bins below 0.05% (0.5 per mille) of frequency will not be considered
        df_diptest_input=pd.DataFrame(arr, columns=['Smooth']) # Array to df that will be the input for diptest
        dip, pval_dip = diptest(df_diptest_input['Smooth'].astype(float)) #diptest using smooth data
        #Normality test
        stat, pval_norm =normaltest(table_sp['ANI']) # D’Agostino and Pearson’s normal test 
    else:
        pval_norm = "Low number of bins"
        pval_dip = "Low number of bins"
    #Find valley
    valley_smooth, valley_smooth_ = find_peaks(-pivot_smooth['Freq_smooth']) #Negative symbol before pivot['Freq'] searchs for valleys, instead of peaks. 
    try:
        valleys_smooth = get_valleys(valley_smooth, pivot_smooth)
    except ValueError:
        valleys_smooth="No valleys detected"
    #Find peak
    peak_smooth, peak_smooth_ = find_peaks(pivot_smooth['Freq_smooth']) #Negative symbol before pivot['Freq'] searchs for valleys, instead of peaks. 
    try:
        peaks_smooth=get_valleys(peak_smooth,pivot_smooth)
    except ValueError:
        peaks_smooth="No peaks detected"
    #Check if there is any value between 99.3 and 99.8
    i, j= 99.3, 99.8
    gap_valleys_smooth=any(i <= x <= j for x in valleys_smooth if isinstance(valleys_smooth,list) == True)
    gap_peaks_smooth=any(i <= x <= j for x in peaks_smooth if isinstance(peaks_smooth,list) == True)
    #Classify species in groups
    group=Group(gap_valleys_smooth, gap_peaks_smooth, pval_dip, pval_norm, average_ANI)
    #Print results    
    print('\t'.join([species, group,str(average_ANI), 
                     str(pval_norm), str(pval_dip),
                     str(valleys_smooth), str(gap_valleys_smooth),
                     str(peaks_smooth), str(gap_peaks_smooth)]))
    #Plot figures
    plt.figure(figsize=(14,10))
    plt.plot(pivot_smooth['ANI'], pivot_smooth['Freq_smooth']*(pivot['Freq'].max() / pivot_smooth['Freq_smooth'].max()), '-', color="sandybrown", label='Smooth')
    plt.bar(pivot['ANI'],pivot['Freq'], width=0.05, color="deepskyblue", edgecolor="grey", label='Original')
    if isinstance(valleys_smooth, list) == True:
        plt.vlines(x=valleys_smooth, ymin=0,ymax=pivot['Freq'].max(),linestyles ="dotted", colors="black",  label="Valleys")
    if isinstance(peaks_smooth, list) == True:
        plt.vlines(x=peaks_smooth, ymin=0,ymax=pivot['Freq'].max(),linestyles ="dashdot",colors="black",label="Peaks")
    plt.xlabel("ANI")
    plt.ylabel("Frequency")
    plt.rc('axes', labelsize=15)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    plt.title('\n'.join([species,group]), fontsize=15, loc="left")
    plt.title('\n'.join([' '.join(["Norm_test",str(pval_norm)]),
                            ' '.join(["Dip_test",str(pval_dip)])]), fontsize=15, loc="right")
    plt.legend(loc="upper left", fontsize=15)
    plt.savefig(pdf, format="pdf", bbox_inches="tight")
    plt.close()


def main():
    table=pd.read_table("ani_table_final_sp_included.txt", header=None,delimiter="\t")
    table.columns=["Seq1","Seq2","ANI","Frag_used","Frag_total","Species"]
    table_filt = table[(table['ANI'] >= 90) & (table['Seq1'] != table['Seq2'])] # Filter hits by ANI. I observed that there are many pair around 92-94% ANI for some species, so I decided to show pairs up to 90%.
    p = PdfPages("Peak_Valley_per_species.pdf")
    print("Species\tGroup\tavg_ANI\tpval_normaltest\tpval_diptest\tValley_smooth\tValley_smooth_in_gap\tPeak_smooth\tPeak_smooth_in_gap\t")
    for sp in table_filt['Species'].unique():
        analyze_sp(sp, table_filt, p)
    p.close()



if __name__ == '__main__':
    main()