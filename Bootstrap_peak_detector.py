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

#Function to detect valleys
def get_valleys(valley,pivot):
    list_valleys=[]
    for i in np.arange(0,len(valley)):
        list_valleys.append(round(pivot['ANI'][valley[i]],1))
    return(list_valleys)

#Function to subset species to get the same number of pairs per species
def rand_subsamp(table_filt):
    df_subset=pd.DataFrame()
    for sp in table_filt['Species'].unique():
        table_sp=table_filt[table_filt['Species'] == sp]
        if len(table_sp) > 150:
            sub=table_sp.sample(n=150)
            df_subset=df_subset.append(sub)
    round_df=df_subset['ANI'].round(1)
    average_ANI=round_df.mean()
    pivot=round_df.value_counts().to_frame()
    pivot.columns=['Freq']
    pivot['ANI']=pivot.index
    pivot=pivot.sort_values(by='ANI', ascending=False)
    pivot=pivot.reset_index()
    range=np.arange(pivot['ANI'].min(),100.1,0.1)
    pivot_smooth=pd.DataFrame()
    pivot_smooth['ANI']=range
    pivot_smooth['Freq_smooth']=gaussian_kde(round_df,bw_method=0.12).evaluate(range) # kernel density. PROBABLY I SHOULD CHANGE "pivot['ANI']" BY RANGE DEFINED AS "range=np.arange(pivot['ANI'].min(),100,0.1)"
    return(pivot,pivot_smooth)

#Function to find peaks and valleys
def find_peaks_valleys(table):
    #Find valley
    valley_smooth, valley_smooth_ = find_peaks(-table['Freq_smooth']) #Negative symbol before pivot['Freq'] searchs for valleys, instead of peaks. 
    if len(valley_smooth) > 0:
        valleys_smooth = get_valleys(valley_smooth, table)
    else:
        valleys_smooth="No valleys detected"
    #Find peak
    peak_smooth, peak_smooth_ = find_peaks(table['Freq_smooth']) #Negative symbol before pivot['Freq'] searchs for valleys, instead of peaks. 
    try:
        peaks_smooth=get_valleys(peak_smooth,table)
    except ValueError:
        peaks_smooth="No peaks detected"
    return(peaks_smooth, valleys_smooth)


#Call rand_subsamp() function and identify peaks and valleys for each subsampling
def bootstrap(iter, table):
    num_iterations=iter # iter is the number of iterations to perform
    list_peaks=[]
    list_valleys=[]
    for i in np.arange(0,num_iterations):
        pivot, pivot_smooth = rand_subsamp(table)
        peaks, valleys = find_peaks_valleys(pivot_smooth)
        list_peaks.extend(peaks)
        list_valleys.extend(valleys)
    return(list_peaks,list_valleys)


def main():
    #Read table
    table=pd.read_table("ani_table_final_sp_included.txt", header=None,delimiter="\t")
    table.columns=["Seq1","Seq2","ANI","Frag_used","Frag_total","Length_Seq1","Length_Seq2","Perc_genome_shared","Species"]

    #Filter by identity and remove self-hits
    table_filt = table[(table['ANI'] >= 95) & (table['Seq1'] != table['Seq2'])] # Filter hits by ANI. I observed that there are many pair around 92-94% ANI for some species, so I decided to show pairs up to 90%.

    #Random subsample and smooth data
    pivot, pivot_smooth = rand_subsamp(table_filt) # Random subsampling to plot
    average_freq=pivot_smooth['Freq_smooth'].mean() # Average frequency of the smoothed data
    list_peaks, list_valleys = bootstrap(1000, table_filt) # Perform bootstrap analysis. Subsample x number of times and find valleys and peaks

    #Plot figures
    p=PdfPages("Bootstrap_plot.pdf")

    figure, axis = plt.subplots(2, 1, figsize=(15, 15))

    axis[0].hist(list_valleys, range=[95, 100], bins=51, color="palegreen", edgecolor="black", label="Valleys")
    axis[0].hist(list_peaks, range=[95, 100], bins=51, color="deepskyblue", edgecolor="black", label="Peaks")
    axis[0].set_xlabel("ANI", fontsize=18)
    axis[0].set_ylabel("Number of bootstraps", fontsize=18)
    axis[0].tick_params(axis='both', labelsize=18)
    axis[0].legend(loc="upper left", fontsize=18, bbox_to_anchor=(1, 1))

    axis[1].plot(pivot_smooth['ANI'], pivot_smooth['Freq_smooth']*((pivot['Freq'].max() / pivot_smooth['Freq_smooth'].max())/2), '-', color="sandybrown", label='Smooth', linewidth=3)
    axis[1].bar(pivot['ANI'],pivot['Freq'], width=0.05, color="deepskyblue", edgecolor="grey", label='Original')
    axis[1].vlines(x=list_valleys, ymin=0,ymax=pivot['Freq'].max(),linestyles ="dotted", colors="palegreen", label="Valleys", linewidth=2)
    axis[1].vlines(x=list_peaks, ymin=0,ymax=pivot['Freq'].max(),linestyles ="dashdot",colors="deepskyblue",label="Peaks", linewidth=2)
    axis[1].hlines(y=average_freq * ((pivot['Freq'].max() / pivot_smooth['Freq_smooth'].max())/2) ,xmin=95,xmax=100,linestyles ="dashed",colors="black",label="Average", linewidth=2)
    axis[1].set_xlabel("ANI", fontsize=18)
    axis[1].set_ylabel("Frequency", fontsize=18)
    axis[1].tick_params(axis='both', labelsize=18)
    axis[1].legend(loc="upper left", fontsize=18, bbox_to_anchor=(1, 1))

    plt.savefig(p, format="pdf", bbox_inches="tight")
    plt.close()

    p.close()


if __name__ == '__main__':
    main()