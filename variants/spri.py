#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 04:53:40 2024

@author: hagar
"""
import pandas as pd
from sys import argv
from Bio import SeqIO
from datetime import datetime

def get_sequences(fasta_file):
    sequences = {}
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    for sample, record in fasta.items():
        sequences[sample] = str(record.seq).upper()
        sequences[str(sample)] = sequences.pop(sample)
    return sequences

def generate_spri_df(var_csv_path, fasta):
    '''
    process variants table to match SPRI requerments.

    Parameters
    ----------
    var_csv_path : str
        variants.csv file that the pipeline generates with variants.py.

    Returns
    -------
    SPRI Dataframe

    '''
    spri_df = pd.DataFrame(columns= ['sample_year', 'sample_number', 'executive_unit',
                                     'method' ,'sample_substance', 'pathogen', 'date', 'test', 'result'])
    
    var_df = pd.read_csv(var_csv_path)
    for index, row in var_df.iterrows():
        spri_df.loc[len(spri_df)] = ['', row["Sample"], 'NGS', 'Sequencing', 'Nose + throat culture', 'SARS-CoV-2','',
                                     '%Coverage', row['% coverage']]
        spri_df.loc[len(spri_df)] = ['', row["Sample"], 'NGS', 'Sequencing', 'Nose + throat culture', 'SARS-CoV-2','',
                                     'Variant (catalog)', row['variant (catalog)']]
        spri_df.loc[len(spri_df)] = ['', row["Sample"], 'NGS', 'Sequencing', 'Nose + throat culture', 'SARS-CoV-2','',
                                     'Variant (bodek)', row['variant (bodek)']]
        spri_df.loc[len(spri_df)] = ['', row["Sample"], 'NGS', 'Sequencing', 'Nose + throat culture', 'SARS-CoV-2','',
                                     'Sequence', fasta[str(row["Sample"])]]
    
    return spri_df
        

if __name__ == '__main__':
    
    fasta_path = argv[1]
    output_path = argv[2]
    var_csv_path = 'results/' + datetime.now().strftime('%Y%m%d') + '_variants.csv'
    fasta = get_sequences(fasta_path)
    spri_df = generate_spri_df(var_csv_path, fasta)
    spri_df.to_csv(output_path, index=False)