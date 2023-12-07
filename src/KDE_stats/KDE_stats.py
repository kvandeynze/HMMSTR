import seaborn as sns
import scipy.stats as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from numpy import array, linspace
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import os
import statistics as stat
from collections import Counter

from HMMSTR_utils.HMMSTR_utils import throw_low_cov, remove_outlier_IQR, remove_flanking_outliers
warnings.filterwarnings('ignore')

class KDE_cluster:
    def __init__(self, target, out, out_count_name, discard_outliers,flanking_like_filter,curr_strand=None):
        self.name = target[0]
        self.out = out
        self.discard_outliers = discard_outliers #this will be used when we goto cluster
        self.flanking_like_filter = flanking_like_filter
        out_count_file = out + "_" + self.name + out_count_name
        if os.path.exists(out_count_file):
            counts_data = pd.read_csv(out_count_file, sep=" ",header=None)
            counts_data.columns = ["read_id","strand","align_score","neg_log_likelihood","subset_likelihood","repeat_likelihood","repeat_start","repeat_end","align_start", "align_end","counts"]
            self.data = counts_data[counts_data.counts != 0] #discard 0 count reads
            if curr_strand is not None:
                self.strand = curr_strand
                self.data = counts_data[counts_data.counts != 0][counts_data.strand == curr_strand] #if this is a stranded run, we only want to keep data from the current strand
            else:
                self.strand = None
            if self.data.empty:
                self.data = None
        else:
            self.data = None
    def get_stats(self, plot_hists, filter_outliers=False, filter_quantile=0.25,flanking_like_filter=False, strand=None):
        '''
        Function to calculate summary stats for a single target

        Args:
            out_count_name(str): output suffix for final result tsv

        Returns:
            counts_data(pandas dataframe): final calculations for given count file

        '''
        #read output file
        counts_data = self.data
        counts_data['freq'] = counts_data.groupby(by='counts')['read_id'].transform('count')
        #ADDED outlier filtering
        outliers = []
        flanking_outliers=[]
        if filter_outliers:
            filtered, outliers = remove_outlier_IQR(counts_data.counts, filter_quantile)
            #counts_data = counts_data[counts_data.counts.isin(filtered)]
        counts_data['outlier'] = np.where(counts_data.counts.isin(outliers),True , False)
        #check for flanking likelihood outliers
        if flanking_like_filter:
            flanking_filtered, flanking_outliers = remove_flanking_outliers(counts_data)
        counts_data['flanking_outlier'] = np.where(counts_data.read_id.isin(flanking_outliers),True , False)
        #testing out including plotting here, will need to add a flag
        if plot_hists:
            self.plot_stats(counts_data, self.out, strand)
        return counts_data
    def plot_stats(self,counts_data,out, strand=None):
        '''
        plot the histograms and median likelihood ratio line plots (hopefully with error bars)

        '''
        #(1) unmodified histogram
        to_plot = counts_data[['counts','freq']].drop_duplicates().sort_values(by="counts")
        if to_plot.shape[0] == 0 or to_plot.shape[1] == 0:
            return
        #(1) unmodified histogram
        sns.reset_orig()
        to_plot.plot.bar(x="counts",y="freq",figsize=(15,5),title= str(self.name), color="red", fontsize=7)
        plt.xlabel("Repeat count")
        plt.ylabel("Number of reads")
        fig = plt.gcf()
        if self.strand is not None:
            fig.savefig(str(self.out + "_plots/"+self.name+ "_"+ self.strand +"_supporting_reads_hist.pdf"))
        else:
            fig.savefig(str(self.out + "_plots/"+self.name+"_supporting_reads_hist.pdf"))
        plt.show()
        plt.clf()
        return #currently does not work and isnt called
    def call_clusters(self,kernel="gaussian", bandwidth = 'scott', max_k = 2, output_plots = False, subset=None, filter_quantile=0.25, allele_specific = False, allele = None, strand=None):
        '''
        This function uses kernal density estimation to resolve alleles. It is ideal when data is normally distributed
        and has a relatively small range or overlapping count distributions

        Parameters
        --------------------------------------------------------------------------------------------------------------
        data: pandas Series
            Series containing count data
        kernel: str
            String designating kernel to use for the kde
        bandwidth: float or str
            Bandwidth to use in kde
        max_k: int
            The maximum number of clusters to keep in outputs, all others are filtered by cluster size
        output_plots: bool
            True if density plot will be saved, False if not

        Returns
        -------------------------------------------------------------------------------------------------------------
        clusters: dict
            Dictionary containing all counts assigned to each cluster
        allele_calls: dict
            Dictionary containing the allele call from each cluster as well as the size and standard deviation within the cluster
        outliers: pandas Series
            Series of counts that are identified as outliers and were dicarded before clustering. Empty if filter_outliers is False

        '''
        #filter outliers if called for
        outliers = pd.Series()
        flanking_outliers = pd.Series()
        if self.discard_outliers:
            filtered, outliers = remove_outlier_IQR(self.data.counts,filter_quantile)
            a = np.array(filtered).reshape(-1,1)
        if self.flanking_like_filter:
            flanking_filtered, flanking_outliers = remove_flanking_outliers(self.data)
            #check if we also filtered outliers so we can do both
            if self.discard_outliers:
                a = np.array(flanking_filtered[flanking_filtered.counts.isin(filtered.unique())].counts).reshape(-1,1) #use reads that are both in both of our filtered sets
            else:
                a = np.array(flanking_filtered.counts).reshape(-1,1)

        else:
            a = np.array(self.data.counts).reshape(-1,1)
        if subset is not None:
            a = np.array(subset.counts).reshape(-1,1) #ADDED 080123, allow a subset of data to be passed for reclustering
        #need to see what the difference between using the fit model and the scores are
        #FIXME threw the error the 0.0 can't be raised to a negative power when computing bandwidth
        try:#initial call
            #print(a.shape)
            kde = KernelDensity( bandwidth=bandwidth,kernel=kernel).fit(a) #bandwith 0.5 works for ALS example

            #ADDED 080423 trying minimum bandwidth
            if kde.bandwidth_ < 0.5:
                kde = KernelDensity( bandwidth=0.5,kernel=kernel).fit(a)
        except:
            if len(a) == 0:
                clusters = -1
                #print("data is empty, returning empty allele call for target: ", self.name)
                allele_calls = {"A"+str(k)+":median":-1 for k in range(1,max_k+1)}
                allele_calls.update({"A"+str(k)+":mode":-1 for k in range(1,max_k+1)})
                allele_calls.update({"A"+str(k)+":SD":-1 for k in range(1,max_k+1)})
                allele_calls.update({"A"+str(k)+":supporting_reads":-1 for k in range(1,max_k+1)})
                allele_calls["bandwidth"] = -1
                return clusters,allele_calls , outliers

            #print("an exception occurred!")
            print("Switching to a constant bandwidth=0.5")
            kde = KernelDensity( bandwidth=0.5,kernel=kernel).fit(a)
        #s = linspace(min(a)-5,max(a)+5)
        #s = linspace(min(a),max(a))
        s = linspace(min(a)-2,max(a)+2)
        e = kde.score_samples(s.reshape(-1,1))
        mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

        #run and change bandwidth till # maxima matches or is less than max_k
        while len(mi)+1 > max_k:
            #print(bandwidth)
            #print(len(mi)+1)
            bandwidth = kde.bandwidth_ + 0.1
            kde = KernelDensity( bandwidth=bandwidth,kernel=kernel).fit(a)
            e = kde.score_samples(s.reshape(-1,1))
            mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

        #save cluster info based on minima and maxima
        all_clusters = {}
        maxima = {}
        sizes = {}
        if len(mi) == 0:
            all_clusters["1"] = a
            maxima["1"] = s[ma]
            sizes["1"] = len(a)
        elif len(mi) == 1:
            #print("hi")
            all_clusters["1"] = a[a < s[mi][0]]
            maxima["1"] = s[ma][0]
            sizes["1"] = len(a[a < s[mi][0]])
            all_clusters["2"] = a[a >= s[mi][0]]
            if len(s[ma]) > 1:
                maxima["2"] = s[ma][1]
            else:
                if isinstance(stat.mode(all_clusters["2"]), np.float64):
                    maxima["2"] = stat.mode(all_clusters["2"]) #added, im not sure why this is an issue
                else:
                    maxima["2"] = stat.mode(all_clusters["2"])[0] #if multiple modes returned, pick the first
            sizes["2"] = len(a[a >= s[mi][0]])
        else: # more than 2 clusters
            all_clusters["1"] = a[a < s[mi][0]]
            maxima["1"] = s[ma][0]
            sizes["1"] = len(a[a < s[mi][0]])
            for i_cluster in range(len(mi)-1):
                all_clusters[str(i_cluster+2)] = a[(a >= s[mi][i_cluster]) * (a <= s[mi][i_cluster+1])]
                maxima[str(i_cluster+2)] = s[ma][i_cluster+1]
                sizes[str(i_cluster+2)] = len(a[(a >= s[mi][i_cluster]) * (a <= s[mi][i_cluster+1])])
            all_clusters[str(len(mi)+1)] = a[a >= s[mi][-1]]
            maxima[str(len(mi)+1)] = s[ma][-1]
            sizes[str(len(mi)+1)] = len(a[a >= s[mi][-1]])

        #filter clusters based on max_k
        sizes_series = pd.Series(sizes)
        #sort sizes
        sizes_series_sorted = sizes_series.sort_values(ascending=False)
        #get max_k clusters by size
        #print(sizes_series_sorted)
        k_clusters = list(sizes_series_sorted.head(max_k).index)
        #print(k_clusters)
        clusters = {k: np.unique(all_clusters[k]) for k in k_clusters}
        #FIXME error was thrown here saying that only arrays of size 1 can be cast to float, I'm not sure what case this is yet
        non_one_maxima = False
        for k in k_clusters:
            #print(type(maxima[k]))
            if isinstance(maxima[k], np.float64) == False and len(maxima[k]) != 1:
                #print(k_clusters)
                non_one_maxima = True
                #print("length of maxima array for cluster ", k, " is not 1, the array is: ",maxima[k])
        if non_one_maxima == False:
            #FIXME need to change cluster names to be sequential, this would also be fixed by limiting clusters called by bandwidth
            #allele_calls = {k+"_KDE_maxima": round(float(maxima[k])) for k in k_clusters}
            allele_calls = {}
            allele_calls.update({"A" +str(k)+":supporting_reads":sizes[k] for k in k_clusters})
            allele_calls.update({"A" +str(k)+":SD":np.std(all_clusters[k]) for k in k_clusters})
            #updated 080123, added stats in case the maxima is weird from clustering outliers or a far-off allele
            allele_calls.update({"A" +str(k)+":median":np.median(all_clusters[k]) for k in k_clusters})
            allele_calls.update({"A" +str(k)+":mode":int(st.mode(np.array(all_clusters[k]),keepdims=True)[0][0]) for k in k_clusters}) #I dont know if this will be the same, for some reason my results from margits data dont line up with finding the mode
            
            allele_calls["bandwidth"] = kde.bandwidth_
        elif len(k_clusters) == 1:
            #this is the case when there is only 1 value detected thus it is the only allele that can be called
            #allele_calls = {"1_KDE_maxima":a[0][0]}
            allele_calls = {}
            allele_calls.update({"A1:supporting_reads":sizes["1"]})
            allele_calls.update({"A1:SD":np.std(all_clusters["1"])})
            allele_calls.update({"A1:median":np.median(all_clusters["1"])})
            allele_calls.update({"A1:mode":int(st.mode(np.array(all_clusters["1"]),keepdims=True)[0][0])})
            allele_calls["bandwidth"] = kde.bandwidth_
        else:
            #allele_calls = {k+"_KDE_maxima":-1 for k in k_clusters}
            allele_calls = {}
            allele_calls["bandwidth"] = -1
            allele_calls.update({"A"+str(k)+":SD":0 for k in k_clusters})
            allele_calls.update({"A"+str(k)+":median":0 for k in k_clusters})
            allele_calls.update({"A"+str(k)+":mode":0 for k in k_clusters})
            return clusters,allele_calls , outliers

        if output_plots:
            if len(mi) == 0:
                #print("first case")
                try:
                    plt.hist(a, density=True, histtype='bar',  label='Data', bins=range(int(min(a)),int(max(a))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left")
                    plt.plot(s[:],np.exp(e[:]),color="black",lw=2,linestyle="-",label="Kernel Density Estimate")
                        
                    #for k in k_clusters:
                    if len(s[ma] != 0):
                        plt.plot(s[ma], np.exp(e[ma]), 'bo',label="KDE Maxima")
                    
                except:
                    print("Encounted an exception in homozygous plotting case (KDE)")
                #plt.plot(s[ma], e[ma], 'go') #plot just the points
                if allele_specific:
                    plt.title(self.name + " KDE: Allele " + allele)
                    ymin, ymax = plt.ylim()
                    plt.ylim(ymin, ymax)
                    plt.vlines(x = allele_calls["A1:mode"],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    plt.vlines(x = allele_calls["A1:median"], ymin=ymin, ymax=ymax,linewidth=2, color='purple',linestyle="dashed",label="Median")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+self.name +"_" +strand+"_allele_"+ allele+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_allele_"+ allele+"_KDE.pdf"))
                else:
                    plt.title(self.name + " KDE")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+self.name +"_" +strand+ "_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_KDE.pdf"))
                plt.show()
                plt.clf()
            elif len(mi) == 1:
                #print("case 2")
                try:
                    plt.hist(a, density=True, histtype='bar',  label='Data', bins=range(int(min(a)),int(max(a))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left")
                    plt.plot(s[:],np.exp(e[:]),color="black",lw=2,linestyle="-",label="Kernel Density Estimate")
                    plt.plot(s[ma], np.exp(e[ma]), 'bo',label="KDE Maxima")
                except:
                    print("An exception has occured")
                plt.plot(s[mi], np.exp(e[mi]), 'ro',label="KDE Minima")
                if allele_specific:
                    plt.title(self.name + " KDE: Allele "+ allele)
                    ymin, ymax = plt.ylim()
                    plt.ylim(ymin, ymax)
                    plt.vlines(x = allele_calls["A1:mode"],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    plt.vlines(x = allele_calls["A1:median"], ymin=ymin, ymax=ymax,linewidth=2, color='purple',linestyle="dashed",label="Median")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+ self.name +"_"+ strand+"_allele_"+ allele+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_allele_"+ allele+"_KDE.pdf"))

                else:
                    plt.title(self.name + " KDE")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+ self.name +"_"+ strand+ "_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_KDE.pdf"))
                plt.show()
                plt.clf()
            #dynamically print clusters
            else:
                try:
                    #print("case3")
                    plt.hist(a, density=True, histtype='bar',  label='Data', bins=range(int(min(a)),int(max(a))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left")
                    plt.plot(s[:],np.exp(e[:]),color="black",lw=2,linestyle="-",label="Kernel Density Estimate")
                    plt.plot(s[ma], np.exp(e[ma]), 'bo',label="KDE Maxima")
                except:
                    print("An exception occured in plotting the KDE")
                plt.plot(s[mi], np.exp(e[mi]), 'ro',label="KDE Minima")
                if allele_specific:
                    plt.title(self.name + " KDE: Allele " +allele)
                    ymin, ymax = plt.ylim()
                    plt.ylim(ymin, ymax)
                    plt.vlines(x = allele_calls["A1:mode"],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    plt.vlines(x = allele_calls["A1:median"], ymin=ymin, ymax=ymax,linewidth=2, color='purple',linestyle="dashed",label="Median")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_"+ strand+"_allele_"+ allele+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_allele_"+ allele+"_KDE.pdf"))

                else:
                    plt.title(self.name + " KDE")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_"+ strand+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_KDE.pdf"))
                plt.show()
                plt.clf()
        return clusters, allele_calls, outliers, flanking_outliers
    def assign_clusters(self, clusters, outliers, flanking_outliers):
        '''
        Function to assign individual reads to clusters determined by KDE

        Parameters
        ------------------------------------------------------------------
        clusters: dict
            Dictionary of counts assigned to each cluster
        
        Returns:
        ----------------------------------------------------------------------
        assignemnt_df: pandas DataFrame
            DataFrame equivalent to input count data with additional columns for cluster
            assignemnts and if the count was an outlier
        '''
        count_assignemnts = {}
        for cluster in clusters:
            for curr in clusters[cluster]:
                count_assignemnts[curr] = cluster
        labels = pd.DataFrame()
        labels['counts'] = count_assignemnts.keys()
        labels['cluster_assignments'] = count_assignemnts.values()
        # TODO figure out if we want this to be consistent with gmm assignment output or if this is fine since we label outliers
        assignment_df = self.data.merge(labels, how = 'left', on = 'counts') #I changed this to not assign a cluster to reads that are outliers
        assignment_df['outlier'] = np.where(assignment_df.counts.isin(outliers),True , False)
        assignment_df['flanking_outlier'] = np.where(assignment_df.read_id.isin(flanking_outliers),True , False)
        return assignment_df

    #add bootstrapping to KDE class
    def bootstrap_KDE(self,data, resample_size, CI_width, max_peaks, out):
        #FIXME need to determine if I should only construct CIs from each allele separately instead of like how I am since if there is a large imbalance in supporting reads, the CIs will include both alleles by chance
        #throw out 1 read coverage, this is something to assume will happen if bootstrapping
        #data['freq']=data.groupby(by='counts')['read_id'].transform('count')
        filtered = throw_low_cov(data)
        call_lists = {}
        for i in range(1,max_peaks+1):
            #call_lists["A"+str(i)+"_mean"]=[]
            call_lists["A"+str(i)+":median"]=[]
            #call_lists["A"+str(i)+"_mode"]=[]

        for i in range(resample_size):
            if resample_size > 100 and (i+1)%10 == 0:
                print("Round " + str(i+1)+ "/" + str(resample_size))
            elif resample_size < 100 and resample_size > 10 and (i+1)%5 == 0:
                print("Round " + str(i+1)+ "/" + str(resample_size))
            elif resample_size < 10:
                print("Round " + str(i+1)+ "/" + str(resample_size))
            #sample up to the depth of the experiment with replacement
            curr_sample = filtered.sample(frac=1,axis=0, replace = True)
            #call gmm with original freqs on this subset
            clusters, allele_calls, outliers, flanking_outliers = self.call_clusters(kernel="gaussian", bandwidth = 'scott', max_k = max_peaks, output_plots = False, subset=curr_sample, filter_quantile=0.25) 
            #add current call to lists
            for j in range(1,max_peaks+1):
                if ("A"+str(j)+":median") not in allele_calls.keys():
                    call_lists["A"+str(j)+":median"].append(0)
                    allele_calls["A"+str(j)+":median"] = 0
                    continue
                #call_lists["A"+str(j)+"_mean"].append(curr_row[0]['H'+str(j)+":mean"])
                call_lists["A"+str(j)+":median"].append(allele_calls["A"+str(j)+":median"])
        #mean_CIs = {}
        median_CIs = {}
        #mode_CIs = {}
        adjustment = float((1.0 - CI_width)/2)
        for j in range(1,max_peaks+1): #get confidence interval for each allele called
            #curr_means = call_lists["A"+str(j)+"_mean"]
            curr_medians = call_lists["A"+str(j)+":median"]
            #curr_modes = call_lists["A"+str(j)+"_mode"]
            #curr_row[0]["A"+str(j)+":mean_CI"] = (np.quantile(curr_means,adjustment),np.quantile(curr_means,1-adjustment))
            median_CIs["A"+str(j)+":median_CI"] = (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
            #curr_row[0]["A"+str(j)+":mode_CI"] = (np.quantile(curr_modes,adjustment),np.quantile(curr_modes,1-adjustment))
        return median_CIs

    def bootstrap_KDE_allele_specific(self,data, resample_size, CI_width, out):
            #throw out 1 read coverage, this is something to assume will happen if bootstrapping
            filtered = throw_low_cov(data)
            call_lists = {}
            #call_lists["mean"]=[]
            call_lists["median"]=[]
            #call_lists["mode"]=[]
    
            for i in range(resample_size):
                if resample_size > 100 and (i+1)%10 == 0:
                    print("Round " + str(i+1)+ "/" + str(resample_size))
                elif resample_size < 100 and resample_size > 10 and (i+1)%5 == 0:
                    print("Round " + str(i+1)+ "/" + str(resample_size))
                elif resample_size < 10:
                    print("Round " + str(i+1)+ "/" + str(resample_size))
                #sample up to the depth of the experiment with replacement
                curr_sample = filtered.sample(frac=1,axis=0, replace = True)
                #call gmm with original freqs on this subset
                clusters, allele_calls, outliers, flanking_outliers = self.call_clusters(kernel="gaussian", bandwidth = 'scott', max_k = 1, output_plots = False, subset=curr_sample, filter_quantile=0.25)  #single allele bootstrapping
                #add current call to lists
                #call_lists["mean"].append(curr_row[0]["A1:mean"])
                call_lists["median"].append(allele_calls["A1:median"])
                #call_lists["mode"].append(curr_row[0]["H1:mode"])
            #mean_CIs = {}
            median_CIs = {}
            #mode_CIs = {}
            adjustment = float((1.0 - CI_width)/2)
            #curr_means = call_lists["mean"]
            curr_medians = call_lists["median"]
            #curr_modes = call_lists["mode"]
            #mean_CI = (np.quantile(curr_means,adjustment),np.quantile(curr_means,1-adjustment))
            median_CI= (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
            #mode_CI = (np.quantile(curr_modes,adjustment),np.quantile(curr_modes,1-adjustment))
            return median_CI