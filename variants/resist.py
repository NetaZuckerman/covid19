#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 11:23:40 2023

@author: hagar
"""

import pandas as pd
import os 

script_path = os.path.dirname(os.path.abspath(__file__))

def run(variants, fasta):

    resist = pd.read_csv(script_path+"/resistantce.csv")
    variants = variants[["Sample", "nt substitutions","variant (catalog)"]]
    

    final_df = pd.DataFrame(columns=["Sample", "resistance_mutations"])
    for index, row in variants.iterrows():
        sample = str(row["Sample"])
        found_mut=[] # list of found mutations 
        if not row["variant (catalog)"] == "QC fail":
            temp = pd.DataFrame()
            temp["nt mutation"] = str(row["nt substitutions"]).split(";")  #build dataframe for sample's nt mutations
            temp = temp.merge(resist, how='left', on='nt mutation').dropna() #get all resistant mutations found in sample
            
            
            if len(temp) > 0:
                #is aa mut exist
                sample_aa_df = pd.DataFrame(columns=["nt_start_on_genome", "aa_in_seq"]) 
                suspect_pos = set(temp["nt_start_on_genome"].astype(int)) #positions to check aa
                
                for pos in suspect_pos:
                    codon = fasta[sample][pos-1:pos+2] 
                    aa = str(codon.translate().seq)
                    sample_aa_df.loc[len(sample_aa_df.index)] = [pos, aa] 
                temp = temp.merge(sample_aa_df, how='left', on='nt_start_on_genome')# add aa in seq for each suspected mutation
                temp = temp[temp.mut == temp.aa_in_seq] #only keep row where aa mutation in seq exist
                
                for i, mut in temp.iterrows():
                    found_mut.append(mut["gene"] + ":" + mut["aa mutation"] + "(" + mut["treatment"] + ")")
            
        final_df.loc[len(final_df.index)] = [sample, "; ".join(found_mut)]
    return final_df

if __name__ == "__main__":

    run()