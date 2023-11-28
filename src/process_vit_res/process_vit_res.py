#helper functions to process viterbi outputs
from colorama import Fore, Back, Style
def label_states(vit_out, hidden_states):
    '''
    Function to convert output from testvit.c to a single sequence of labeled hidden states

    Args:
        vit_out (text file): text file output of testvit.c
        hidden_states_in (text file): text file of "." delimited states in correct order to label with

    Returns:
        labeled_seq (list): a list corresponding to the input sequence in terms of hidden_states
        pointers (dict): a dictionary of "pointers" to the beginning of a state set (P,R, or S) so we can jump there later

    '''
    labeled_seq = []
    #pointers to record indeces of the start of each part of our model in the read sequence
    pointers = {"P":False, "R":False, "S":False, "G2":False, "D":[], "I":[]}
    #get likelihoods from viterbi output
    MLE = float(vit_out.split("\n")[0].rstrip().split(" ")[6])
    seq = vit_out.split("\n")[3]
    seq_list = seq.rstrip().split(" ")
    seq_list = [int(i) for i in seq_list]
    #print("seq_list: ", seq_list)
    # print("hidden_states: ", hidden_states)

    #build dictionary to convert from numeric to states
    conversion = {}
    for i,state in enumerate(hidden_states):
        conversion[i+1] = state #viterbi labels are one indexed
    #convert labeled_seq
    for i,obs in enumerate(seq_list):
        labeled_seq.append(conversion[seq_list[i]])
        #check if we need to add a "pointer"
        if  "P" in conversion[seq_list[i]] and pointers["P"] == False:
            pointers["P"] = i
        if "R" in conversion[seq_list[i]] and pointers["R"] == False:
            pointers["R"] = i
        if  "S" in conversion[seq_list[i]] and pointers["S"] == False:
            pointers["S"] = i
        if conversion[seq_list[i]] == "G2" and pointers["G2"] == False:
            pointers["G2"] = i
        if "D" in conversion[seq_list[i]]: #get indeces of all deletions
            pointers["D"].append(i)
        if "I" in conversion[seq_list[i]]: #get indeces of all deletions
            pointers["I"].append(i)
    #print(labeled_seq)
    return labeled_seq, pointers,MLE

# def calc_likelihood(vit_out, pointers, labeled_seq,hidden_states, read):
#     '''
#     Function to calculate likelihoods corresponding to subseqeunces to be used later

#     Args:
#         vit_out (text file): text file output of testvit.c
#         pointers (dict): dictionary of indeces in labeled_seq corresponding to beginning of different sections of the model
#         labeled_seq (str): labeled state seqeunce corresponding to viterbi output
#         hidden_states (text file): text file of "." delimited states in correct order to label with
#         read (str): raw read sequence

#     Returns:
#         final (float): final likelihood corresponding to prefix --> suffix viterbi path
#         final_labels(str): seqeunce labels for prefix --> suffix, "-" delimited
#         repeats (str): sequence including deletion labels corresponding to only the target repeat in the seqeunce
#         context (str): sequence including deletion labels corresponding to the path through prefix-->suffix

#     '''
#     #we can calculate the likelihood using the labeled sequence
#     likelihoods_int = vit_out.split("\n")[4]
#     likelihoods_dec = vit_out.split("\n")[5]
#     #process likelihood outputs, combine leading integers with decimals
#     likelihood_int_list = likelihoods_int.rstrip().split(" ")
#     likelihood_dec_list = likelihoods_dec.rstrip().split(" ")
#     likelihood_int_list = [int(i) for i in likelihood_int_list]
#     likelihood_dec_list = [abs(int(i)) for i in likelihood_dec_list]
#     likelihood_list = [float('.'.join(str(i) for i in x)) for x in list(zip(likelihood_int_list,likelihood_dec_list))]
#     #get the last G1 and S likelihoods
#     Sn = pointers["G2"] - 1
#     Gn = pointers["P"] - 1
#     final = likelihood_list[Sn] - likelihood_list[Gn]
#     final_repeat_like = likelihood_list[pointers["S"]-1] - likelihood_list[pointers["R"]-1] #likelihood of just the repeat
#     final_labels = "-".join(labeled_seq[pointers["P"]:pointers["G2"]])
#     #get context and repeat sequence, can't just use pointers cause there may be deletions so labelled seq and read may be different lengths
#     #FIXME need to adjust for number deletions so we subset less of the read, otherwise indeces are off. not sure if only accounting for prefix or all deletions is necessary
#     print("full labelled seq is: ", labeled_seq)
#     print("curr read is: ",read) 
#     print("count of PD is: ",final_labels.count('PD'))
#     print("sequence between pointers: ", read[pointers["P"]:pointers["G2"]])
#     print("sequence between pointers accounting for deletions: ", read[pointers["P"]:pointers["G2"]+len(pointers['D'])])
#     num_pd = final_labels.count('PD')
#     sub_read = read[pointers["P"]+num_pd:] #need to adjust for deletions but I think it maybe postion dependent
#     #get indeces in labelled with deletions

#     deletions = [i - pointers["P"] for i in pointers["D"]] #convert to be relative to prefix start position
#     print("original deletions: ", pointers['D'])
#     print("new deletion indeces: ", deletions)
#     #for each deletion, add a "-" in the raw sequence (this may not be the most efficient way to do this)
#     context = ""
#     offset = 0
#     print(final_labels)
#     for i in range(len(final_labels.split("-"))):
#         if i in deletions:
#             context= context + "-"
#             offset+=1
#         else:
#             context = context + sub_read[i-offset]

#     #now our read is the same length as the labelled seq so we can index with the pointers
#     R_start = pointers["R"] - pointers["P"]
#     R_end = pointers["S"] - pointers["P"]
#     repeats = context[R_start:R_end] #need to see if this includes last index of repeat
#     return final, final_labels, repeats, context, final_repeat_like
def calc_likelihood(vit_out, pointers, labeled_seq,hidden_states, read,subset_start, subset_end):
    '''
    Function to calculate likelihoods corresponding to subseqeunces to be used later

    Args:
        vit_out (text file): text file output of testvit.c
        pointers (dict): dictionary of indeces in labeled_seq corresponding to beginning of different sections of the model
        labeled_seq (str): labeled state seqeunce corresponding to viterbi output
        hidden_states (text file): text file of "." delimited states in correct order to label with
        read (str): raw read sequence

    Returns:
        final (float): final likelihood corresponding to prefix --> suffix viterbi path
        final_labels(str): seqeunce labels for prefix --> suffix, "-" delimited
        repeats (str): sequence including deletion labels corresponding to only the target repeat in the seqeunce
        context (str): sequence including deletion labels corresponding to the path through prefix-->suffix

    '''
    #we can calculate the likelihood using the labeled sequence
    likelihoods_int = vit_out.split("\n")[4]
    likelihoods_dec = vit_out.split("\n")[5]
    #process likelihood outputs, combine leading integers with decimals
    likelihood_int_list = likelihoods_int.rstrip().split(" ")
    likelihood_dec_list = likelihoods_dec.rstrip().split(" ")
    likelihood_int_list = [int(i) for i in likelihood_int_list]
    likelihood_dec_list = [abs(int(i)) for i in likelihood_dec_list]
    likelihood_list = [float('.'.join(str(i) for i in x)) for x in list(zip(likelihood_int_list,likelihood_dec_list))]
    #get the last G1 and S likelihoods
    Sn = pointers["G2"] - 1
    Gn = pointers["P"] - 1
    final = likelihood_list[Sn] - likelihood_list[Gn]
    final_repeat_like = likelihood_list[pointers["S"]-1] - likelihood_list[pointers["R"]-1] #likelihood of just the repeat

    #ADDED 11/13 -- get repeat coordinates per read for downstream motif analysis
    repeat_end =pointers["S"]-1 + subset_start
    repeat_start = pointers["R"] + subset_start

    final_labels = "-".join(labeled_seq[pointers["P"]:pointers["G2"]])
    #get context and repeat sequence, can't just use pointers cause there may be deletions so labelled seq and read may be different lengths
    #FIXME need to adjust for number deletions so we subset less of the read, otherwise indeces are off. not sure if only accounting for prefix or all deletions is necessary
    # print("count of PD is: ",final_labels.count('PD'))
    # print("sequence between pointers: ", read[pointers["P"]:pointers["G2"]])
    # print("sequence between pointers accounting for deletions: ", read[pointers["P"]:pointers["G2"]+len(pointers['D'])])
    num_pd = final_labels.count('PD')
    sub_read = read[pointers["P"]+num_pd:] #need to adjust for deletions but I think it maybe postion dependent
    #get indeces in labelled with deletions

    deletions = [i - pointers["P"] for i in pointers["D"]] #convert to be relative to prefix start position
    #print("length full labelled seq is: ", len(labeled_seq), " curr read length is: ",len(read),"\n", " number deletions: ", len(pointers['D']), " number of insertions: ", "".join(labeled_seq).count("I"), "\n"," number of RD: ", "".join(labeled_seq).count("RD"), " number of RI: ", "".join(labeled_seq).count("RI"),"\n"," number of SD: ", "".join(labeled_seq).count("SD"), " number of SI: ", "".join(labeled_seq).count("SI"))

    # print("difference in number of insertions and deletions: ", str("".join(labeled_seq).count("RD") - "".join(labeled_seq).count("RI") ))
    # print("length states - length read: ", str(len(labeled_seq) - len(read) ))

    #print("new deletion indeces: ", deletions)
    #for each deletion, add a "-" in the raw sequence (this may not be the most efficient way to do this)
    context = ""
    offset = 0
    #print(final_labels)
    for i in range(len(labeled_seq)):
        if 'G' in labeled_seq[i]:
            continue
        if 'D' in labeled_seq[i]:
            context= context + "-"
            offset+=1
        else:
            context = context + read[i-offset]

    #now our read is the same length as the labelled seq so we can index with the pointers
    R_start = pointers["R"] - pointers["P"]
    R_end = pointers["S"] - pointers["P"]
    repeats = context[R_start:R_end] #need to see if this includes last index of repeat
    return final, final_labels, repeats, context, final_repeat_like, repeat_start, repeat_end
def count_repeats(labeled_seq, pointers,repeat_len,seq):
    '''
    Function to count the number of repeats in a sequence given a sequence labeled by its hidden states

    Args:
        labeled_seq (list): hidden state annotation of query sequence
        pointers (dict): dictionary of "pointers" to the first match to prefix, repeat, and suffix
        repeat_len (int): length of target repeat

    Returns:
        count (int): the number of repeats identified in the query sequence

    '''
    repeat_region = labeled_seq[pointers["R"]:pointers["S"]]
    adjusted_length = len(repeat_region)
    #count insertions to account for
    for state in repeat_region:
        if state[1] == "I":
            adjusted_length -= 1
    #calculate counts from adjusted length, deletions already accounted for in labeled sequence
    count = adjusted_length/repeat_len
    return count

def print_labelled(read_id,strand,sub_labels,context,pointers,out):
    '''
    Function to output color-coded context seqeunce

    Args:
        read_id(str): read id for input sequence
        strand(str): strand of read to output
        sub_labels(str): "-" delimited labels for prefix-->suffix
        context(str): context sequence including deletion labels
        pointers(dict): dictionary of pointers corresponding to indeces of insertions, deletions, and start positions
        out(str): output file suffix
    Returns:
    None, outputs new context sequence string to context file for given target
    '''
    #print("context value: ",context)
    R_start = pointers["R"] - pointers["P"]
    R_end = pointers["S"] - pointers["P"]
    context_list = list(context)
    I = [i - pointers["P"] for i in pointers["I"]]
    D = [i - pointers["P"] for i in pointers["D"]]
    for i in I:
        if i < R_start or i > R_end:
            context_list[i] = '\x1b[5;37;42m' + context_list[i] + Style.RESET_ALL + '\x1b[1;30;40m'
        else:
            context_list[i] = '\x1b[5;37;42m' + context_list[i] +  Style.RESET_ALL + '\x1b[1;37;40m'
    for i in D:
        if i < R_start or i > R_end:
            context_list[i] = '\x1b[5;31;41m'+" " + Style.RESET_ALL + '\x1b[1;30;40m'
        else:
            context_list[i] = '\x1b[5;31;41m' + " "+  Style.RESET_ALL + '\x1b[1;37;40m'
    #replace all inserted bases with color + base to be printed
    file = open(out,"a")
    file.write(read_id + " " + strand + " " + '\x1b[1;30;40m' + "".join(context_list[:R_start]) + Style.RESET_ALL + '\x1b[1;37;40m' + "".join(context_list[R_start:R_end]) + Style.RESET_ALL +'\x1b[1;30;40m'+ "".join(context_list[R_end:])+ Style.RESET_ALL + "\n")
    file.close()
    return
    