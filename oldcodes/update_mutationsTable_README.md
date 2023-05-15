# Rules for updaing mutationsTable.xlsx


1. Sheet = Variant. sheet name = variant name.

2. Header is the same as other sheets header. \
![img.png](img.png)
   Position = position in genome \
   Reference = Ref-seq nucleotide in that position \
   Mutation = the variant's nuleotide in that position \
   protein = gene \
   variant = name of mutation. for example: G662S is the substitution mutation from G aa to S aa in the gene's 662 amino acid \
   Mutation type = extragenic/SNP/SNP_silent.... \
   annotation = extra info \
   varname = Giorgio output \
   % of sequences = percentage of seqs of the variant that has that mutation.
    
3. Make sure that gaps of deletions and insertions are written as '-' and not '.'. 
4. Make sure that deletions are represented as one line for each position.
