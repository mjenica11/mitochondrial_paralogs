# Libraries
import pandas as pd

# Quality strings from the first four lines of a non-failing control: SRR7426784.fastq
read1=["9",">","9",">",">","?","<","B","B","B","B","B","B","B","B","B","B","B",
       "B","B","B","?","B","B","B","@","B","@","@","<","B","=",";","@","=","B",
       "B","@","B","A","A","C","C","B","@","D","D","D","B","D","?","F","E","F",
       "A","?","7",")",")",":","8","G","F",">","F","F","D","F","D","D","D","D",
       "D","@","?",";"]

read2=["#","D","D","D","D","B","9","D","D","D","D","D","D","B","D","D","D","D",
       "D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D",
       "D","D","D","D","E","E","E","E","E","F","F","F","F","F","F","H","H","H",
       "J","J","J","J","J","I","J","I","J","H","H","G","H","H","F","F","F","F",
       "F","C","C","C"]

read3=["@","@","@","F","F","F","F","D","F","H","H","G","H","I","J","J","F","D",
       "D","D","D","D","D","D","B","D","D","D","D","D","D","D","A","B","@","B",
       "D","B","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#",
       "#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#",
       "#","#","#","#"]

read4=["C","C","C","F","F","F","F","F","H","H","H","H","H","I","J","I","I","I",
       "G","G","H","F","D","D","D","D","B","B","B","D","D","D","B","B","@","B",
       "D","D","B","9","C","D","D","@","B","D","D","&","#","#","#","#","#","#",
       "#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#","#",
       "#","#","#","#"]

# Convert Phred quality scores (ASCII characters) to probabilities
# Function
read_lst=[read1, read2, read3, read4]

def phred_function(arg0):
    res1=[10**-(ord(x)-33/10.0) for x in arg0]
    res2=['{:0.3e}'.format(x) for x in res1]
    return(res2)

phred1=phred_function(arg0=read_lst[0])
phred2=phred_function(arg0=read_lst[1])
phred3=phred_function(arg0=read_lst[2])
phred4=phred_function(arg0=read_lst[3])

df = pd.DataFrame({'quality_string1':read1, 'probability_score1':phred1,
                   'quality_string2':read2, 'probability_score2':phred2,
                   'quality_string3':read3, 'probability_score3':phred3,
                   'quality_string4':read4, 'probability_score4':phred4})

print(df.head)

df = df.reindex(['quality_string1', 'probability_score1',
                 'quality_string2', 'probability_score2',
                 'quality_string3', 'probability_score3',
                 'quality_string4', 'probability_score4'], axis=1)

print(df.head)

df.to_csv('quality_scores.csv', index=True, encoding='utf-8')
