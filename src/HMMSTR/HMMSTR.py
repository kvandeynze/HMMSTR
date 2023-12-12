#I dont need all of these, remove all unused at the end, many will be encapsulated by other files
#from collections import Counter, defaultdict
import argparse
import os, glob
import gzip
#import numpy as np
import pandas as pd
from time import perf_counter
import multiprocessing as mp
import pickle as pkl
import sys
#custom imports
from profile_HMM.profile_HMM import ProfileHMM
from HMMSTR_utils.HMMSTR_utils import *
from process_read.process_read import Process_Read
from GMM_stats.GMM_stats import GMMStats
from KDE_stats.KDE_stats import KDE_cluster

class LoadFromFile (argparse.Action):
    #for reading from input file instead of command line
    def __call__ (self, parser, namespace, values, option_string = None):
        with values as f:
            # parse arguments in the file and store them in the target namespace
            parser.parse_args(f.read().split(), namespace)

def read_my_file(argv):
    # if there is a '-A' string in argv, replace it, and the following filename
    # with the contents of the file (as strings)
    # you can adapt code from _read_args_from_files
    new_argv = []
    with open(argv,"r") as params:
        for a in params: #read each file
            print(a.rstrip().split(" "))
            for b in a.rstrip().split(" "):
                new_argv.append(b)
    return new_argv
def write_input_command(args):
    with open(args.out + "_run_input.txt","w") as outfile:
        outfile.write(vars(args)["subcommand"] + "\n")
        #check subcommand 
        if vars(args)["subcommand"] == "coordinates":
            #write out coordinates specific arguments
            outfile.write(vars(args)["coords_file"] + "\n")
            outfile.write(vars(args)["chrom_sizes_file"] + "\n")
            outfile.write(vars(args)["ref"] + "\n")
            outfile.write("--input_flank_length"+" "+str(vars(args)["input_flank_length"]) + "\n")
        else:
            #targets
            outfile.write(vars(args)["targets"] + "\n")
        outfile.write(vars(args)["out"] + "\n")
        outfile.write(vars(args)["inFile"] + "\n")

        #TODO supress bootstrapping options if bootstrapping isn't selected
        for key,val in vars(args).items():
            if key in ["subcommand","out","inFile","coords_file","targets","ref","chrom_sizes_file","input_flank_length"]: #these are positional fields, put these where they go first
                continue
            if val in [None, False]: #check if we want this parameter in our input
                continue
            if val == True: #stores true
                outfile.write("--"+key +"\n")
                continue
            outfile.write("--"+key + " "+ str(val)+"\n")
    return
def convert_to_fasta(tsv,out):
    '''
    Function to write "read" entry for a given target for both prefix and suffix
    fasta files to be used as queries in blastn

    Args:
        tsv (Pandas Dataframe): input tsv containing name, prefix, repeat, suffix
        out (str): prefix/path for outputted files

    Returns:
        nothing, writes prefix and suffix fasta files to out_prefix.fa and out_suffix.fa respectively

    '''
    for row in tsv.itertuples():
        curr_out = open(out+"_prefix.fa",'a')
        curr_out.write(">"+row.name+"\n")
        curr_out.write(row.prefix + "\n")
        curr_out.close()
        curr_out = open(out +"_suffix.fa",'a')
        curr_out.write(">"+row.name+"\n")
        curr_out.write(row.suffix + "\n")
        curr_out.close()
    return
def build_all(curr_target, out, background,alphabet,transitions,mismatch_probs, repeat_probs,flanking_size):
    '''
    Function to initialize a ProfileHMM object for a row in input targets tsv

    Parameters
    ----------------------------------------------------------------------------
    curr_target: pandas Series
        Row from input tsv corresponding to a single target
    alphabet: list
        list of characters used to define alphabet of HMM
    background: dictionary
        Dictionary of background probabilities to use as background frequencies
    out: str
        output prefix

    Returns
    ------------------------------------------------------------------------------
    None, writes ProfileHMM in format needed for modified umdhmm-v1.02
    '''
    name = curr_target[0]
    prefix = curr_target[1].rstrip()
    repeat = curr_target[2].rstrip()
    suffix = curr_target[3].rstrip()
    #initialize ProfileHMM object
    curr_hmm = ProfileHMM(name=name,alphabet=alphabet, prefix=prefix, repeat=repeat, suffix=suffix, background=background, out=out, transitions=transitions,mismatch_probs=mismatch_probs, repeat_probs=repeat_probs,flanking_size=flanking_size)
    curr_hmm.build_profile("forward") #builds and writes forward ProfileHMM object
    curr_hmm.build_profile("reverse") #builds and writes reverse ProfileHMM object

    #write out the transition and emission matrices to pickle files for original scaling protocol CAN BE REMOVED IN FINAL
    # with open(out + "_" + name + "_A.pkl", "wb") as outfile:
    #     pkl.dump(curr_hmm.A, outfile)
    # with open(out + "_" + name + "_E.pkl", "wb") as outfile:
    #     pkl.dump(curr_hmm.E, outfile)
    return
def process_read(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,out,targets, build_pre, mode, cutoff, k, w, use_full_seq, flanking_size,output_labelled_seqs):
    '''
    Wrapper function to call all Process_Read methods. This will initialize all operations
    for a given read including (1) alignment (2) target identification and filtering 
    (3) repeat identification with viterbi (4) repeat counting and (5) recording relevant metrics

    Parameters
    ----------------------------------------------------------------------------------------------
    header: str
        header for read from fasta file
    seq: str
        raw read sequence
    hmm_file: str
        path to file containing written profile hmm from build_all
    rev_hmm_file: str
        path to file continaing the reverse strand profile hmm from build_all
    hidden_states: str
        path to hidden states file 
    hidden_states_rev: str
        path to hidden states file for reverse strand
    out: str
        Output prefix for run including full or relative path
    targets: pandas DataFrame
        DataFrame containing input information for each target to search for in reads
    build_pre: str
        Path and prefix for hmm files, used in case targets have already been used and user doesn't want duplicate files
    mode: str
        Sequencing technology used, map-ont for Nanopore, pb for PacBio
    cutoff: int
        Mapq cutoff for mappy aligner to determine read maping location
    
    Returns
    ------------------------------------------------------------------------------------------------
    None, all outputs written to out
    '''
    curr_read = Process_Read(header, seq, cutoff, mode, out, k,w, use_full_seq)
    curr_read.align_mappy() #sets prefix_df and suffix_df fields for read
    # if len(curr_read.suffix_df.name) > 0:
    #     print(curr_read.prefix_df)
    #     print(curr_read.suffix_df) 
    curr_read.assign_targets(targets) #sets target_info for current read, this is a dictionary of keys=target name and values=info, all information needed for the target
    #The check for no assigned targets is in run_vit
    targets_found = curr_read.run_viterbi(hmm_file=hmm_file,rev_hmm_file=rev_hmm_file,hidden_states=hidden_states, rev_states=hidden_states_rev, out=out,build_pre=build_pre, prefix_idx = flanking_size,output_labelled_seqs=output_labelled_seqs)
    return
def ratio_gmm_stats(row, out, out_count_name,plot_hists,max_peaks, filter_outliers, bootstrap, CI_width, resample_size, allele_specific_plots,allele_specif_CIs,filter_quantile): 
    '''
    This is a wrapper function to run stats methods on all targets. This version corresponds to our original stats methods
    which includes normalizing by a perfect match path through the hmm followed by GMM clustering and allele calls

    Args:
        target(pandas series): row of target dataframe
        out(str): output suffix
        out_count_name(str): output suffix for final result tsv

    Returns:
        None, results written out to tsv per target

    '''
    #check existence of count file, may not exist if no reads identified for a given target
    name = row[0]
    out_count_file = out + "_" + name + out_count_name
    #print(out_count_file)
    if os.path.exists(out_count_file) == False:
        print(out_count_file + " does not exist, writing null row for", name,"...")
        # returning null at the beginning of the run may be causing issues with apply, try returning a null row of the right dimensions
        curr_dict = {"name":name}
        for i in range(1,max_peaks+1):
            #curr_mean = "A"+str(i)+":mean"
            curr_median = "A"+str(i)+":median"
            curr_mode = "A"+str(i)+":mode"
            curr_sd = "A" + str(i) + ":SD"
            curr_support = "A" + str(i) + ":supporting_reads"
            #if curr_mean not in curr_dict.keys():
             #   curr_dict[curr_mean] = 0
            if curr_median not in curr_dict.keys():
                curr_dict[curr_median] = 0
            if curr_mode not in curr_dict.keys():
                curr_dict[curr_mode] = 0
            if curr_sd not in curr_dict.keys():
                curr_dict[curr_sd] = 0
            if curr_support not in curr_dict.keys():
                curr_dict[curr_support] = 0
            #curr_dict["A"+str(i)] = -1
            curr_dict["num_supporting_reads"] = 0

            #check if we were suppose to bootstrap for this run so we can return the correct number of columns
            if bootstrap:
                curr_dict["A"+str(i)+":median_CI"] = None
            if allele_specif_CIs:
                curr_dict["A"+ str(i)+":median_CI_allele_specific"] = None
        return pd.Series(curr_dict)
    #initialize GMMStats object
    gmm_stats = GMMStats(target_row=row) #contains all target attributes as well as E and A dictionaries
    final_data = gmm_stats.get_stats(out_count_file, out,plot_hists, filter_outliers, filter_quantile)

    final_data.to_csv(out +"_"+ gmm_stats.name +"_read_info.tsv", index = False, sep="\t") #previously final_out.tsv

    #peak calling
    final_data = final_data[final_data.outlier == False].copy()
    curr_row, final_data['cluster_assignments'] = gmm_stats.call_peaks(final_data, out, max_peaks, plot=plot_hists, save_allele_plots = allele_specific_plots)
    if curr_row is None:
        print("current row doesnt exist")
        return
    
    if bootstrap:
        print("Starting bootstrap...")
        curr_row, cluster_assignments_bootstrap = gmm_stats.bootstrap_gmm(curr_row,final_data[final_data.outlier == False], resample_size, CI_width, max_peaks, out)
    if allele_specif_CIs:
        for assignment in final_data['cluster_assignments'].unique():
            curr_row["A"+ str(assignment+1)+":median_CI_allele_specific"] = gmm_stats.bootstrap_gmm_allele_specific(final_data[final_data.cluster_assignments == assignment],resample_size,CI_width,out)

    #write out cluster assignments to file
    assignments = final_data[['read_id','counts']]
    assignments["Name"] = gmm_stats.name
    assignments['assignment'] = final_data['cluster_assignments']
    
    #check if this is the first entry, if it is, save the header, if its not, append to existing file
    if os.path.exists(out + "read_assignments_gmm.tsv"):
        assignments[['Name','read_id','counts','assignment']].to_csv(out + "read_assignments_gmm.tsv",header=None, sep="\t", index=False, mode="a")
    else:
        assignments[['Name','read_id','counts','assignment']].to_csv(out + "read_assignments_gmm.tsv", sep="\t", index=False, mode="a")

    return curr_row #return genotype

def decide_method(target,out,out_count_name):
    #read in counts for current row
    name = target[0]
    out_count_file = out + "_" + name + out_count_name
    if os.path.exists(out_count_file):
        counts_file= pd.read_csv(out_count_file, sep=" ",header=None)
        counts_file.columns = ["read_id","strand","align_score","neg_log_likelihood","subset_likelihood","repeat_likelihood","repeat_start","repeat_end","align_start", "align_end","counts"]
        counts=counts_file.counts
    else:
        return None
    
    QR_range, IQR_range = get_IQR(counts,0.10)
    # print(IQR_range)
    # print(len(counts[counts < IQR_range[0]]),len(counts[counts > IQR_range[1]]))

    #number of outliers according to IQR
    lower_counts = len(counts[counts <= IQR_range[0]])
    upper_counts = len(counts[counts >= IQR_range[1]])
    #number of outliers according to quartiles
    lower_Q1_counts = len(counts[counts <= QR_range[0]])
    upper_Q3_counts = len(counts[counts <= QR_range[1]])
    #gonna try to play with these metrics to figure out what combination works well
    # print(QR_range[0],IQR_range[0],min(counts))
    # print(QR_range[1],IQR_range[1],max(counts))
    data_width = max(counts) - min(counts)
    QR_width = QR_range[1]-QR_range[0]
    IQR_width = IQR_range[1] - IQR_range[0]

    #3 cases: (1) KDE--throw-outliers (2) KDE (3) GMM
    #Try defining "small quantile range" as QR_width < 3
    #Try defining outliers conservatively as lower_counts or upper_counts > 0, use lower_Q1_counts and upper_Q2_counts for calling more outliers -- would probably filter these if we do use the KDE-throw-outliers since its not as dependent on normality but definitely not for GMM
    #Try defining "large range" proportionally to the QR range since we only care about this in defining when we want to use KDE-throw-outliers vs GMM, ie large + small QR vs larger + large QR (or similar)
    #set cutoff as an arbitrary low number since KDE is very sensitive and will be thrown off by most outliers

    #check if we should use KDE, no outliers found -- need to play with data_width definition here
    #print("QR_width: ",QR_width, " data_width: ",data_width)
    if data_width < 4 or (QR_width <= 3 and data_width < 2*QR_width): #not including checking for outliers here rn, might need to
        # print("using KDE")
        decision = "kde"
    elif QR_width <= 4 and (lower_counts > 0 or upper_counts > 0): #might need to play with this
        # print("using KDE but throwing outliers")
        decision = "kde_throw_outliers"
    else:
        # print("using GMM") #default
        decision = "gmm"
    return decision

#assumes that decision has been added to targets dataframe as an additional column
def call_peaks(row, out, out_count_name, plot_hists, max_peaks, filter_outliers=False, filter_quantile=0.25,bootstrap=False, CI_width=0.95, resample_size=100,allele_specific_plots=False,allele_specif_CIs=False, bandwidth='scott',kernel="gaussian",flanking_like_filter=False):
    '''
    This method is the wrapper function for calling both KDE and GMM classes based on the decision per target, it also insures the outputs are uniform across all methods
    '''
    #make cluster object based on decision
    decision = row[4]

    if decision is not None and "kde" in decision:
        if decision == "kde_throw_outliers":
            filter_outliers = True #override when this option is chosen, will be the same as kde if user overrides
        curr_kde = KDE_cluster(row,out,out_count_name, filter_outliers,flanking_like_filter) #discard outliers will still be included in all of these regardless of if KDE-throw-outliers is chosen, that is the only option where it is the default
        if curr_kde.data is None: #FIXME I am not sure if this will require column names to be returned or if apply doesnt care
            allele_calls = {}
            allele_calls = {"name":curr_kde.name}
            allele_calls['bandwidth'] = -1
            allele_calls["peak_calling_method"] = -1
            for i in range(max_peaks):
                allele_calls["A"+str(i + 1)+":median"] = 0
                allele_calls["A"+str(i + 1)+":mode"] = 0
                allele_calls["A"+str(i + 1)+":supporting_reads"] = 0
                allele_calls["A"+str(i+1)+":SD"] = 0
                if bootstrap:
                    allele_calls["A"+str(i)+":median_CI"] = (0,0)
                if allele_specif_CIs:
                    allele_calls["A"+ str(i+1)+":median_CI_allele_specific"] = (0,0)
            allele_calls_series = pd.Series(allele_calls)
            cols = allele_calls_series.index.tolist()
            cols.sort()
            final = allele_calls_series[cols]
            return final
        #get read info tsv
        curr_kde.data = curr_kde.get_stats(plot_hists, filter_outliers, filter_quantile, flanking_like_filter)
        #curr_kde.data.to_csv(out +"_"+ curr_kde.name +"_final_out.tsv", index = False, sep="\t") # combined with read_assignments
        
        clusters, allele_calls, outliers, flanking_outliers = curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=max_peaks, output_plots =  plot_hists, filter_quantile=filter_quantile)
        assignments = curr_kde.assign_clusters(clusters,outliers, flanking_outliers)

        #plot allele_specific plots if applicable
        if allele_specific_plots:
            for assignment in assignments.cluster_assignments.unique():
                curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=1, output_plots =  True, subset=assignments[assignments.cluster_assignments == assignment], allele_specific=True, allele=assignment, filter_quantile=filter_quantile)
    
        #write out read assignments
        assignments["name"] = curr_kde.name
        assignments["peak_calling_method"] = decision
        
        #only output flanking outlier if applicable
        if flanking_like_filter:
            if os.path.exists(out + "read_assignments.tsv") == False:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        else:
            if os.path.exists(out + "read_assignments.tsv") == False:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        
        #account for differing column number for less than max_k peaks
        #print("unique assignements: ", assignments.assignment.unique())
        #print("max_k: ", max_k)

        if len(assignments[assignments.outlier == False][assignments.flanking_outlier == False].cluster_assignments.unique()) < max_peaks:
            for i in range(max_peaks - len(assignments[assignments.outlier == False][assignments.flanking_outlier == False].cluster_assignments.unique())):
                j = len(assignments[assignments.outlier == False][assignments.flanking_outlier == False].cluster_assignments.unique())+1
                allele_calls["A"+str(i+j)+":median"] = 0
                allele_calls["A"+str(i+j)+":mode"] = 0
                allele_calls["A"+str(i+j)+":supporting_reads"] = 0
                allele_calls["A"+str(i+j)+":SD"] = 0
        allele_calls["name"] = curr_kde.name
        allele_calls_series = pd.Series(allele_calls)
        #bootstrap
        assignments = pd.merge(left=assignments, right=curr_kde.data[["read_id","freq"]], how="left")
        if bootstrap:
            if clusters != -1:
                print("Starting bootstrap...")
                #TODO I don't think I should be subsetting this beforehand since I call clustering from the kde object which has the fields for filtering outliers, change this once original test is done
                median_CIs = curr_kde.bootstrap_KDE(assignments, resample_size, CI_width, max_peaks, out)
                allele_calls_series = pd.concat([allele_calls_series,pd.Series(median_CIs)])
        if allele_specif_CIs: 
            for assignment in assignments.cluster_assignments.unique():
                if clusters == -1:
                    break
                if pd.isna(assignment):
                    #print("hi")
                    continue
                print("Starting allele-specific bootstrap...")
                allele_calls_series["A"+ str(assignment)+":median_CI_allele_specific"] = curr_kde.bootstrap_KDE_allele_specific(assignments[assignments.cluster_assignments == assignment], resample_size, CI_width, out)
        
        #get total number of supporting reads
        allele_calls_series["num_supporting_reads"] = 0
        for k in range(max_peaks):
            allele_calls_series["num_supporting_reads"] = allele_calls_series["num_supporting_reads"] + allele_calls["A"+str(k+1) + ":supporting_reads"]
        allele_calls_series["peak_calling_method"] = decision
        cols = allele_calls_series.index.tolist()
        cols.sort()
        final = allele_calls_series[cols]
        return final #remove this later, we want all of the outputs to be uniform just checking this works for now
        
    else:#GMM is the only other option
        #check existence of count file, may not exist if no reads identified for a given target
        name = row[0]
        out_count_file = out + "_" + name + out_count_name
        #print(out_count_file)
        if os.path.exists(out_count_file) == False:
            print(out_count_file + " does not exist, writing null row for", name,"...")
            # returning null at the beginning of the run may be causing issues with apply, try returning a null row of the right dimensions
            curr_dict = {"name":name}
            curr_dict['bandwidth'] = -1
            curr_dict["peak_calling_method"] = -1
            for i in range(1,max_peaks+1):
                #curr_mean = "A"+str(i)+":mean"
                curr_median = "A"+str(i)+":median"
                curr_mode = "A"+str(i)+":mode"
                curr_sd = "A" + str(i) + ":SD"
                curr_support = "A" + str(i) + ":supporting_reads"
                #if curr_mean not in curr_dict.keys():
                 #   curr_dict[curr_mean] = 0
                if curr_median not in curr_dict.keys():
                    curr_dict[curr_median] = 0
                if curr_mode not in curr_dict.keys():
                    curr_dict[curr_mode] = 0
                if curr_sd not in curr_dict.keys():
                    curr_dict[curr_sd] = 0
                if curr_support not in curr_dict.keys():
                    curr_dict[curr_support] = 0
                if bootstrap:
                    curr_dict["A"+str(i)+":median_CI"] = (0,0)
                if allele_specif_CIs:
                    curr_dict["A"+ str(i+1)+":median_CI_allele_specific"] = (0,0)
                #curr_dict["A"+str(i)] = -1
                curr_dict["num_supporting_reads"] = 0
            
            return pd.Series(curr_dict)
        #initialize GMMStats object
        gmm_stats = GMMStats(target_row=row) #contains all target attributes as well as E and A dictionaries
        final_data = gmm_stats.get_stats(out_count_file, out,plot_hists, filter_outliers, filter_quantile=filter_quantile,flanking_like_filter=flanking_like_filter)
    
        #final_data.to_csv(out +"_"+ gmm_stats.name +"_final_out.tsv", index = False, sep="\t")
    
        #NEW check the range of the data to determine if we should be calling this homozygous
        # if final_data[final_data.outlier == False].counts.max()-final_data[final_data.outlier == False].counts.min() < 2:
        #     #change max peaks to 1, this should be called homozygous and the gmm won't do that
        #     max_peaks = 1
        
        #peak calling
        final_data2 = final_data[final_data.outlier == False][final_data.flanking_outlier == False][final_data.counts != 0].copy()
        curr_row, final_data2['cluster_assignments'] = gmm_stats.call_peaks(final_data2, out, max_peaks, plot=plot_hists, save_allele_plots = allele_specific_plots)
        if curr_row is None:
            print("current row doesnt exist")
            return
        
        if bootstrap:
            print("Starting bootstrap...")
            curr_row = gmm_stats.bootstrap_gmm(curr_row,final_data2, resample_size, CI_width, max_peaks, out)
        if allele_specif_CIs:
            for assignment in final_data2['cluster_assignments'].unique():
                print("Starting allele-specific bootstrap...")
                median_CIs = gmm_stats.bootstrap_gmm_allele_specific(final_data2[final_data2.cluster_assignments == assignment],resample_size,CI_width,out)
                curr_row["A"+ str(assignment+1)+":median_CI_allele_specific"] = median_CIs
    
        #write out cluster assignments to file
        #assignments = pd.merge(left = final_data2[['read_id','counts','cluster_assignments']],right=final_data[["read_id","counts","outlier","flanking_outlier"]], on = ["read_id","counts"], how="right")
        assignments = pd.merge(left = final_data2[['read_id','counts','cluster_assignments']],right=final_data, on = ["read_id","counts"], how="right")
        assignments["name"] = gmm_stats.name
        assignments["peak_calling_method"] = decision
        assignments["cluster_assignments"] = assignments["cluster_assignments"]+1
        #assignments['assignment'] = final_data['cluster_assignments']
        if flanking_like_filter:
            if os.path.exists(out + "read_assignments.tsv") == False:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end','counts', 'freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts','freq','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        else:
            if os.path.exists(out + "read_assignments.tsv") == False:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts','freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','strand', 'align_score','neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood','repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts','freq','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
    #add bandwidth column for kde calls
    curr_row['bandwidth'] = -1
    curr_row["peak_calling_method"] = decision
    return curr_row #return genotype
    
def call_peaks_stranded(row, out, out_count_name, plot_hists, max_peaks, filter_outliers=False, filter_quantile=0.25,bootstrap=False, CI_width=0.95, resample_size=100,allele_specific_plots=False,allele_specif_CIs=False, bandwidth='scott',kernel="gaussian",flanking_like_filter=False, strand=None):
    '''
    This method is the wrapper function for calling both KDE and GMM classes based on the decision per target, it also insures the outputs are uniform across all methods
    This is the stranded version of call_peaks, it does everything exactly the same but separates count data by strand before performing any operations
    '''
    #make cluster object based on decision
    decision = row[4]
    if decision is not None and "kde" in decision:
        if decision == "kde_throw_outliers":
            filter_outliers = True #override when this option is chosen, will be the same as kde if user overrides
        curr_kde = KDE_cluster(row,out,out_count_name, filter_outliers,flanking_like_filter,strand) #discard outliers will still be included in all of these regardless of if KDE-throw-outliers is chosen, that is the only option where it is the default
        if curr_kde.data is None: #FIXME I am not sure if this will require column names to be returned or if apply doesnt care
            allele_calls = {}
            allele_calls = {"name":curr_kde.name}
            allele_calls['bandwidth'] = -1
            allele_calls["peak_calling_method"] = -1
            for i in range(max_peaks):
                allele_calls["A"+str(i + 1)+":median"] = 0
                allele_calls["A"+str(i + 1)+":mode"] = 0
                allele_calls["A"+str(i + 1)+":supporting_reads"] = 0
                allele_calls["A"+str(i+1)+":SD"] = 0
                if bootstrap:
                    allele_calls["A"+str(i)+":median_CI"] = (0,0)
                if allele_specif_CIs:
                    allele_calls["A"+ str(i+1)+":median_CI_allele_specific"] = (0,0)
            allele_calls["strand"] = strand
            allele_calls_series = pd.Series(allele_calls)
            cols = allele_calls_series.index.tolist()
            cols.sort()
            final = allele_calls_series[cols]
            #final_stranded.append(pd.DataFrame(final).T)
            #continue #FIXME test this out
            return final #TODO figure out how to return in a stranded manner
        #get read info tsv
        curr_kde.data = curr_kde.get_stats(plot_hists, filter_outliers, filter_quantile, flanking_like_filter, strand=strand)
        #curr_kde.data.to_csv(out +"_"+ curr_kde.name +"_"+strand+"_final_out.tsv", index = False, sep="\t") #stranded out
        
        clusters, allele_calls, outliers, flanking_outliers = curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=max_peaks, output_plots =  plot_hists, filter_quantile=filter_quantile,strand=strand)
        assignments = curr_kde.assign_clusters(clusters,outliers, flanking_outliers)

        #plot allele_specific plots if applicable
        if allele_specific_plots:
            for assignment in assignments.cluster_assignments.unique():
                curr_kde.call_clusters(kernel=kernel, bandwidth=bandwidth,max_k=1, output_plots =  True, subset=assignments[assignments.cluster_assignments == assignment], allele_specific=True, allele=assignment, filter_quantile=filter_quantile, strand=strand)
    
        #write out read assignments
        assignments["name"] = curr_kde.name
        assignments["peak_calling_method"] = decision
        #only output flanking outlier if applicable
        if flanking_like_filter:
            if os.path.exists(out +"_"+strand+ "_read_assignments.tsv") == False:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out + "_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        else:
            if os.path.exists(out +"_"+strand+ "_read_assignments.tsv") == False:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out + "_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        
        #account for differing column number for less than max_k peaks
        #print("unique assignements: ", assignments.assignment.unique())
        #print("max_k: ", max_k)
        if len(assignments[assignments.outlier == False].cluster_assignments.unique()) < max_peaks:
            for i in range(max_peaks - len(assignments[assignments.outlier == False].cluster_assignments.unique())):
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":median"] = 0
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":mode"] = 0
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":supporting_reads"] = 0
                allele_calls["A"+str(i+len(assignments[assignments.outlier == False].cluster_assignments.unique())+1)+":SD"] = 0
        allele_calls["name"] = curr_kde.name
        allele_calls["strand"] = strand
        allele_calls_series = pd.Series(allele_calls)

        #bootstrap
        assignments = pd.merge(left=assignments, right=curr_kde.data[["read_id","freq"]], how="left")
        if bootstrap:
            if clusters != -1:
                print("Starting bootstrap...")
                #TODO I don't think I should be subsetting this beforehand since I call clustering from the kde object which has the fields for filtering outliers, change this once original test is done
                median_CIs = curr_kde.bootstrap_KDE(assignments, resample_size, CI_width, max_peaks, out)
                allele_calls_series = pd.concat([allele_calls_series,pd.Series(median_CIs)])
        if allele_specif_CIs: 
            for assignment in assignments.cluster_assignments.unique():
                if clusters == -1:
                    break
                if pd.isna(assignment):
                    #print("hi")
                    continue
                #print(assignment)
                #print(assignments[assignments.assignment == assignment])
                print("Starting allele-specific bootstrap...")
                allele_calls_series["A"+ str(assignment)+":median_CI_allele_specific"] = curr_kde.bootstrap_KDE_allele_specific(assignments[assignments.cluster_assignments == assignment], resample_size, CI_width, out)
        
        #get total number of supporting reads
        allele_calls_series["num_supporting_reads"] = 0
        for k in range(max_peaks):
            allele_calls_series["num_supporting_reads"] = allele_calls_series["num_supporting_reads"] + allele_calls["A"+str(k+1) + ":supporting_reads"]
        allele_calls_series["peak_calling_method"] = decision
        cols = allele_calls_series.index.tolist()
        cols.sort()
        final = allele_calls_series[cols]
        return final #TODO figure out how I want to return genotypes in strand specific manner
    
    else:#GMM is the only other option
        #check existence of count file, may not exist if no reads identified for a given target
        name = row[0]
        out_count_file = out + "_" + name + out_count_name
        #print(out_count_file)
        if os.path.exists(out_count_file) == False:
            print(out_count_file + " does not exist, writing null row for", name,"...")
            # returning null at the beginning of the run may be causing issues with apply, try returning a null row of the right dimensions
            curr_dict = {"name":name}
            curr_dict['bandwidth'] = -1
            curr_dict["peak_calling_method"] = -1
            for i in range(1,max_peaks+1):
                #curr_mean = "A"+str(i)+":mean"
                curr_median = "A"+str(i)+":median"
                curr_mode = "A"+str(i)+":mode"
                curr_sd = "A" + str(i) + ":SD"
                curr_support = "A" + str(i) + ":supporting_reads"
                #if curr_mean not in curr_dict.keys():
                #   curr_dict[curr_mean] = 0
                if curr_median not in curr_dict.keys():
                    curr_dict[curr_median] = 0
                if curr_mode not in curr_dict.keys():
                    curr_dict[curr_mode] = 0
                if curr_sd not in curr_dict.keys():
                    curr_dict[curr_sd] = 0
                if curr_support not in curr_dict.keys():
                    curr_dict[curr_support] = 0
                if bootstrap:
                    curr_dict["A"+str(i)+":median_CI"] = (0,0)
                if allele_specif_CIs:
                    curr_dict["A"+ str(i+1)+":median_CI_allele_specific"] = (0,0)
                #curr_dict["A"+str(i)] = -1
                curr_dict["num_supporting_reads"] = 0
            curr_dict["strand"] = strand
            return pd.Series(curr_dict)
        #initialize GMMStats object
        gmm_stats = GMMStats(target_row=row) #contains all target attributes as well as E and A dictionaries
        final_data = gmm_stats.get_stats(out_count_file, out,plot_hists, filter_outliers, filter_quantile=filter_quantile,flanking_like_filter=flanking_like_filter, curr_strand=strand)
    
        final_data.to_csv(out +"_"+ gmm_stats.name +"_"+strand+"_final_out.tsv", index = False, sep="\t")
        
        #peak calling
        final_data2 = final_data[final_data.outlier == False][final_data.flanking_outlier == False][final_data.counts != 0].copy()
        curr_row, final_data2['cluster_assignments'] = gmm_stats.call_peaks(final_data2, out, max_peaks, plot=plot_hists, save_allele_plots = allele_specific_plots, strand=strand)
        if curr_row is None:
            print("current row doesnt exist")
            return
        
        if bootstrap:
            print("Starting bootstrap...")
            curr_row = gmm_stats.bootstrap_gmm(curr_row,final_data2, resample_size, CI_width, max_peaks, out)
        if allele_specif_CIs:
            for assignment in final_data2['cluster_assignments'].unique():
                median_CIs = gmm_stats.bootstrap_gmm_allele_specific(final_data2[final_data2.cluster_assignments == assignment],resample_size,CI_width,out)
                curr_row["A"+ str(assignment+1)+":median_CI_allele_specific"] = median_CIs
    
        #write out cluster assignments to file
        assignments = pd.merge(left = final_data2[['read_id','counts','cluster_assignments']],right=final_data[["read_id","counts","outlier","flanking_outlier"]], on = ["read_id","counts"], how="right")
        assignments["name"] = gmm_stats.name
        assignments["peak_calling_method"] = decision
        assignments["cluster_assignments"] = assignments["cluster_assignments"]+1
        #assignments['assignment'] = final_data['cluster_assignments']
        if flanking_like_filter:
            if os.path.exists(out +"_"+strand+ "_read_assignments.tsv") == False:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out +"_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","flanking_outlier","peak_calling_method"]].to_csv(out +"_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
        else:
            if os.path.exists(out +"_"+strand+ "_read_assignments.tsv") == False:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out +"_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a")
            else:
                assignments[['name','read_id','counts','cluster_assignments',"outlier","peak_calling_method"]].to_csv(out +"_"+strand+ "_read_assignments.tsv", sep="\t", index=False, mode="a", header=None)
    #add bandwidth column for kde calls
    curr_row['bandwidth'] = -1
    curr_row["peak_calling_method"] = decision
    curr_row["strand"] = strand
    return curr_row


def main():
    parser = argparse.ArgumentParser()
    #Required inputs, TODO subcommands for input types, test this
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True
    parser_infile = subparsers.add_parser('targets_tsv')
    # add a required argument
    parser_infile.add_argument(dest ="targets",type=str, help='TSV with name, prefix, repeat, suffix, header required. Recommended at least 100bp with default alignment parameters')
    #coordinates mode parser
    parser_coords = subparsers.add_parser('coordinates')
    parser_coords.add_argument(dest = 'coords_file',type=str,help='Path to input custom bed file with either 4 (chr,start,end,repeat motif) or 5 (chr,start,end,repeat motif, name) columns')
    parser_coords.add_argument(dest='chrom_sizes_file',type=str,help = 'Path to chromosome sizes file')
    parser_coords.add_argument(dest="ref",type=str,help="Path to reference genome to get flanking sequence from")
    parser_coords.add_argument("--input_flank_length",type=int,help="Length of prefix and suffix to get from reference, must be longer than 30bp (Default) or flanking_size parameter (default: %(default)s)",default=200)

    #parser.add_argument(dest ="targets",type=str, help='TSV with name, prefix, repeat, suffix, header required. Recommended at least 100bp with default alignment parameters')

    #same as above but for targets_tsv
    parser.add_argument(dest ="out",type=str, help='Output prefix including directory path')
    parser.add_argument(dest = "inFile", type=str, help= 'Sequence to search and annotate in fasta or fastq, gzipped accepted')


    #coordinate input options ***these inputs are dependent on other inputs and I need to figure out how to deal with that
    #TODO figure out how to make arguments conditional ie how do I make it so you can pass either coordinates or target sheet? do I need a bunch of conditionals to control this or is there a
    #streamlined way to do this?

    #optional, if input, models and accompanying files already produced and use the provided prefix
    parser.add_argument("--background", type=str, help="Pickle of custom background dictionary, must include empty string for deletion character and total must add up to 1", default=None)
    parser.add_argument("--E_probs", type=str, help= 'TSV of custom emissions probability matrix, see format specifications', default=None)
    parser.add_argument("--A_probs", type=str, help= 'TSV of custom transition probabilities, see format specifications', default=None)
    parser.add_argument("--custom_RM", type=str, help= 'TSV of custom repeat match state probability matrix, see format specifications, used for motif mosacism, currently supports single target', default=None)
    parser.add_argument("--hmm_pre", type=str,help="Prefix for files produced by build function, use if running the same targets across multiple input files") #if I clean up the directory after I will get rid of this potentially
    #defaults subject to change
    #parser.add_argument("--output_hist",help="Output supporting read histogram", action='store_true')

    #TODO change all output_hist flags to output_plots and wrap additional plots in this flag, also make plots directory at out prefix
    parser.add_argument("--output_plots",help="Output plots, default outputs supporting read histograms and model of best fit plots per target", action='store_true')
    parser.add_argument("--output_labelled_seqs",help="Output the model path through prefix, repeat, and suffix identified per read as context_labelled.txt per target", action='store_true')

    #parser.add_argument("--exclude_sd", type=int,help="Number of standard deviations from the mean of tallest peaks to exclude neighboring peaks", default=2)
    parser.add_argument("--max_peaks", type=int,help="Maximum number of peaks to calculate AIC and BIC for in peak calling (default: %(default)s)",default=2)
    parser.add_argument("--cpus", type=int,help="Number of cpus to use. If none given, half of available cpus used (default: %(default)s)",default=int(mp.cpu_count()/2))
    parser.add_argument("--flanking_size", type=int,help="Number of basepairs flanking repeat to encode in model, longer recommended for highly repeatitive regions (default: %(default)s)",default=100)
    #alignment parameters
    parser.add_argument("--mode", type=str,help="map-ont (Nanopore), pb (PacBio), or sr (short accurate reads, use for accurate short flanking sequence input) (default: %(default)s)",default="map-ont")
    parser.add_argument("--mapq_cutoff", type=int,help="MapQ cutoff for prefix and suffix alignment (default: %(default)s)",default=30)
    parser.add_argument("--k", type=int,help="Kmer parameter to be passed to mappy",default=None)
    parser.add_argument("--w", type=int,help="Window parameter to be passed to mappy",default=None)
    parser.add_argument("--use_full_read", help="Flag to not subset the read by the alignment, optimal for reads with short prefix and suffix",action='store_true')

    #NEW peakcalling options
    parser.add_argument("--peakcalling_method", type=str, help="If set, this peak calling method will be used instead of chosen by decision process. Options: gmm, kde, or kde_throw_outliers (default: %(default)s)",default="auto")
    parser.add_argument("--bandwidth", type=int, help="Bandwidth to use for KDE (if chosen). Will default to scott method if not provided", default=None)
    parser.add_argument("--kernel", type=str, help="Kernel to use for KDE (if chosen) (default: %(default)s)", default="gaussian")

    #optional boostrapping parameters
    parser.add_argument("--bootstrap", help="Boolean designating to output boostrapped confidence intervals with genotype calls",action='store_true')
    parser.add_argument("--call_width", type=float, help="Decimal percentage designating width of confidence interval for allele calls (default: %(default)s)", default=0.95)
    parser.add_argument("--resample_size", type=int, help="Number of times to resample in bootstrapping (default: %(default)s)", default=500)

    #allele specific stats options
    parser.add_argument("--allele_specific_CIs",  help="Output allele-specific bootsrapped confidence intervals", action='store_true')
    parser.add_argument("--allele_specific_plots", help="Output allele-specific histograms with model of best fit", action='store_true')
    parser.add_argument("--discard_outliers", help="Discard outliers based on quantile", action='store_true')
    parser.add_argument("--filter_quantile", type=float, help="Float designating quantile of count frequency to discard when filtering outliers (default: %(default)s)", default=0.25)

    #post read processing filters/options
    parser.add_argument("--flanking_like_filter", help="If set, filter reads per target that are designated as outliers based on flanking sequence likelihood",action='store_true')
    parser.add_argument("--stranded_report", help="If set, genotypes are called for each strand separately as well as together and strand bias is reported if found",action='store_true')

    #for debugging
    parser.add_argument("--cluster_only", help="Run clustering on output corresponding to required arguments, mostly for debugging past runs", action='store_true')
    parser.add_argument("--save_intermediates", help="Flag designating to save intermediate files such as model inputs, raw count files, and state sequence files", action='store_true')

    #check input to see if file or command line was used
    if len(sys.argv) != 2:
        args = parser.parse_args()
    elif len(sys.argv) == 2:
        if sys.argv[1] == "--help" or sys.argv[1] == "-h":
            args = parser.parse_args()
            sys.exit()
        args = parser.parse_args(read_my_file(sys.argv[1]))

    #TODO print out all options chosen and save parameter file for run
    print("Parameters for this HMMSTR run: ")
    with open(args.out + "_run_parameters.txt","w") as outfile:
        for key,val in vars(args).items():
            print(key,":",val)
            outfile.write(key + ":"+ str(val)+"\n")
    #TODO output file with parameters for this run but output in input format so run can be replicated in the future -- not sure if this should just be an option, outputtin 2 of these seems silly
    write_input_command(args)


    
    if args.output_plots:
        directory = args.out + "_plots"
        if os.path.exists(directory) == False:
            os.mkdir(directory) 
            print("Directory '% s' created" % directory)
        else:
            print("Directory '% s' already exists! Plots will be output to existing directory..." % directory)
    if args.output_labelled_seqs: #FIXME check if directory already exists
        directory = args.out + "_labelled_seqs"
        if os.path.exists(directory) == False:
            os.mkdir(directory) 
            print("Directory '% s' created" % directory)
        else:
            print("Directory '% s' already exists! Labelled sequences will be output to existing directory..." % directory)
    #check input mode and check if inputs compatible for coordinates
    if args.subcommand == 'coordinates':
        #check inputs
        if args.input_flank_length < args.flanking_size:
            print("Prefix and suffix length (input_flank_length,%d) are incompatible with model flanking length (flanking_size, %d), exiting..."%(args.input_flank_length,args.flanking_size))
            return
        #generate infile from coordinates given
        print("Coordinates input selected, getting %dbp flanking sequence from reference at %s..."%(args.input_flank_length,args.ref))
        targets = generate_input_sheet(args.coords_file, args.chrom_sizes_file, args.ref, args.input_flank_length)
        #TODO decide if this should be an option
        print("Inputs generated from coordinates! Saving target inputs to %s..."%(args.out+"_inputs.tsv"))
        targets.to_csv(args.out + "_inputs.tsv", sep="\t",index=False)
    else:
        #read in targets
        print("Target tsv input selected! Reading input from  %s..."%(args.targets))
        targets = pd.read_csv(args.targets, sep="\t")
    #convert custom inputs
    if args.E_probs is not None:
        E_probs = read_model_params(args.E_probs,"emissions")
    else:
        E_probs = None
    if args.A_probs is not None:
        A_probs = read_model_params(args.A_probs,"transitions")
    else:
        A_probs = None
    if args.background is not None:
        background = read_model_params(args.background,"background")
    else:
        background = None
    if args.custom_RM is not None:
        custom_RM = read_model_params(args.custom_RM,"repeat")
        #print(custom_RM)
    else:
        custom_RM = None
    
    if args.bandwidth is None:
        args.bandwidth = 'scott'

    #default alphabet that is compatible with all methods, subject to change if users want to have additional inputs such as methylated bases ect in the future
    alphabet = ['A','T','C','G','']
    if args.cluster_only == False:
        if args.hmm_pre is None:
            pool_start = perf_counter()
            targets.apply(build_all,axis=1,args=(args.out, background, alphabet,A_probs, E_probs, custom_RM, args.flanking_size)) #added flanking size
            pool_end = perf_counter()
            print("All models built! finished .... time was: ", str(pool_end-pool_start))
            build_pre = args.out #need to make sure this is compatible with new build all implentation
        else:
            print("Using previously written input files with prefix: " + args.hmm_pre)
            build_pre = args.hmm_pre
        convert_to_fasta(targets,args.out) #write out all target flanking sequences to fasta files to be used by mappy

        #file prefixes and suffixes to use
        hmm_file = ".hmm" #this is a suffix, does not include individual names****
        rev_hmm_file = "_revcomp" + ".hmm"
        hidden_states = ".hidden_states.txt"
        hidden_states_rev = "_revcomp.hidden_states.txt"

        #start multiprocess of read processing
        pool_start = perf_counter()
        print("Starting multiprocess, cpu's in use: " + str(args.cpus))
        pool = mp.Pool(processes=args.cpus)

        #check input file type and pass to appropriate parser
        #check if compressed
        if args.inFile.endswith('gz'):
            #gzipped file
            if args.inFile.endswith("fasta.gz") or args.inFile.endswith("fa.gz"): #fasta file
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read, args.flanking_size, args.output_labelled_seqs)) for header, seq, bool in read_fasta(gzip.open(args.inFile,'rt'))]
            elif args.inFile.endswith("fastq.gz") or args.inFile.endswith("fq.gz"): #fastq
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read, args.flanking_size,args.output_labelled_seqs)) for header, seq, bool in read_fastq(gzip.open(args.inFile,'rt'))]
        else:
            if args.inFile.endswith('a'): #fasta file
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read,args.flanking_size,args.output_labelled_seqs)) for header, seq, bool in read_fasta(open(args.inFile))]
            elif args.inFile.endswith('q'): #fastq
                [pool.apply_async(process_read, args=(header,seq,hmm_file,rev_hmm_file,hidden_states,hidden_states_rev,args.out,targets, build_pre, args.mode, args.mapq_cutoff, args.k, args.w, args.use_full_read,args.flanking_size,args.output_labelled_seqs)) for header, seq, bool in read_fastq(open(args.inFile))]
        pool.close()
        pool.join()
        pool_end = perf_counter()
        print("All reads processed, the pooled run took: ", pool_end-pool_start)

    #stats runs and allele calls
    if args.hmm_pre is None:
        build_pre = args.out #need to make sure this is compatible with new build all implentation
    else:
        build_pre = args.hmm_pre
    out_count_name = "_counts.txt"

    pool_start = perf_counter()

    #TODO test this integration in place of ratio_gmm_stats
    if args.peakcalling_method == "auto":
        targets["peak_call_method"] = targets.apply(decide_method,args=(args.out, out_count_name), axis=1)
    elif args.peakcalling_method in ["kde","kde_throw_outliers","gmm"]:
        targets["peak_call_method"] = args.peakcalling_method
    else:
        print("Invalid method inputed for peak calling, using auto...")
        targets["peak_call_method"] = targets.apply(decide_method,args=(args.out, out_count_name), axis=1)

    #default output
    if args.stranded_report == False:
        geno_df = targets.apply(call_peaks, args=(args.out, out_count_name, args.output_plots,args.max_peaks,args.discard_outliers,args.filter_quantile,args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs, args.bandwidth,args.kernel,args.flanking_like_filter), axis=1)
        #call original allele call procedure
        #geno_df = targets.apply(ratio_gmm_stats,axis=1,args=(args.out,out_count_name, args.output_hist,args.max_peaks,args.discard_outliers, args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs,args.filter_quantile))
        #geno_df.dropna(axis="rows",how="any", inplace=True)
        if isinstance(geno_df, pd.DataFrame):
            #sort outputs
            geno_df_final = geno_df.apply(sort_outputs, args=(args.max_peaks,args.bootstrap,args.allele_specific_CIs), axis=1)
            geno_df_final.to_csv(args.out + "_genotype_calls.tsv",sep="\t",index=False)
            pool_end = perf_counter()
            print("Genotyping run done! Took: ", pool_end-pool_start)
        else:
            print(geno_df)
            print("Results are not a DataFrame! Something went wrong...")
    else: #stranded output
        geno_forward_df = targets.apply(call_peaks_stranded, args=(args.out, out_count_name, args.output_plots,args.max_peaks,args.discard_outliers,args.filter_quantile,args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs, args.bandwidth,args.kernel,args.flanking_like_filter,"forward"), axis=1)
        geno_reverse_df = targets.apply(call_peaks_stranded, args=(args.out, out_count_name, args.output_plots,args.max_peaks,args.discard_outliers,args.filter_quantile,args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs, args.bandwidth,args.kernel,args.flanking_like_filter, "reverse"), axis=1)
        
        #concat stranded results
        geno_df = pd.concat([geno_forward_df,geno_reverse_df]).sort_values(by="name")
        #call original allele call procedure
        #geno_df = targets.apply(ratio_gmm_stats,axis=1,args=(args.out,out_count_name, args.output_hist,args.max_peaks,args.discard_outliers, args.bootstrap, args.call_width, args.resample_size,args.allele_specific_plots,args.allele_specific_CIs,args.filter_quantile))
        #geno_df.dropna(axis="rows",how="any", inplace=True)
        if isinstance(geno_df, pd.DataFrame):
            geno_df_final = geno_df.apply(sort_outputs, args=(args.max_peaks,args.bootstrap,args.allele_specific_CIs,True), axis=1)
            for i in range(1,args.max_peaks+1):
                geno_df_final["A"+str(i)+":strand_std"] = geno_df_final.groupby("name")["A"+str(i)+":median"].transform(np.std)
            geno_df_final.to_csv(args.out + "_genotype_calls.tsv",sep="\t",index=False)
            #geno_df.to_csv(args.out + "_stranded_genotype_calls.tsv",sep="\t",index=False)
            pool_end = perf_counter()
            print("Genotyping done! Took: ", pool_end-pool_start)
        else:
            print(geno_df)
            print("Results are not a DataFrame! Something went wrong...")
    #clean up intermediates
    #TODO make this an option
    #all files produced for models
    if not args.save_intermediates:
        for file_subset in glob.glob(args.out+"*hidden_states.txt"):
            if os.path.exists(file_subset):
                os.remove(file_subset)
        for file_subset in glob.glob(args.out+"*.hmm"):
            if os.path.exists(file_subset):
                os.remove(file_subset)
        for file_subset in glob.glob(args.out+"*labeled_seqs.txt"):
            if os.path.exists(file_subset):
                os.remove(file_subset)
        #FIXME counts file is required for cluster only
        for file_subset in glob.glob(args.out+"*counts.txt"):
            os.remove(file_subset)
        if os.path.exists(args.out+"_prefix.fa"):
            os.remove(args.out+"_prefix.fa")
        if os.path.exists(args.out+"_suffix.fa"):
            os.remove(args.out+"_suffix.fa")
if __name__ == "__main__":
  main()