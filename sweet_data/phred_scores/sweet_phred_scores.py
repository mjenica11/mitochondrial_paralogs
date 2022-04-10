# Libraries
import pandas as pd

# Quality strings from the first 4 lines of SRR7426784.fastq
line1=["#","<","<","B","B","F","F","F","F","F","F","F","F","F","F","F","F",
       "F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F",
       "F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F"]

line2=["#","<","<","B","B","F","F","F","F","B","/","F","B","B","/","F","F",
       "F","/","F","F","<","B","F","/","<","F","/","/","F","F","/","F","F",
       "/","<","F","B","/","B","<","B","/","<","F","F","#","#","#","#","#"]

line3=["#","<","<","B","B","F","F","F","F","F","F","F","F","F","F","F","F",
       "F","F","F","F","F","F","F","F","F","F","F","F","F","B","F","F","F",
       "F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F"]

line4=["#","<","<","/","B","B","F","F","F","F","F","F","F","F","F","F","F",
       "F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F",
       "F","F","F","F","F","F","F","F","F","F","F","F","F","F","B","F","F"]

# Convert Phred quality scores (ASCII characters) to probabilities
#phred1=[10**-(ord(x)-33/10.0) for x in line1]

# Function
#phred_lst=[phred1, phred2, phred3, phred4]
ascii_lst=[line1, line2, line3, line4]

def phred_function(arg0):
    res1=[10**-(ord(x)-33/10.0) for x in arg0]
    res2=['{:0.3e}'.format(x) for x in res1]
    return(res2)

phred1=phred_function(arg0=ascii_lst[0])
phred2=phred_function(arg0=ascii_lst[1])
phred3=phred_function(arg0=ascii_lst[2])
phred4=phred_function(arg0=ascii_lst[3])

df = pd.DataFrame({'quality_string1':line1, 'probability_score1':phred1,
                   'quality_string2':line2, 'probability_score2':phred2,
                   'quality_string3':line3, 'probability_score3':phred3,
                   'quality_string4':line4, 'probability_score4':phred4})

print(df.head)

df = df.reindex(['quality_string1', 'probability_score1',
                 'quality_string2', 'probability_score2',
                 'quality_string3', 'probability_score3',
                 'quality_string4', 'probability_score4'], axis=1)

print(df.head)

df.to_csv('quality_scores.csv', index=True, encoding='utf-8')

