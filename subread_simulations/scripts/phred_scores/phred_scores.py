# Libraries
import pandas as pd

# Quality string
line1=["C","C","C","F","F","F","F","F","H","H","H","H","H","J","J","J","J",
        "J","I","J","J","J","J","J","J","G","I","H","I","I","J","J","J","J",
        "J","J","J","I","J","I","I","I","J","J","J","J","J","J","I","J","J",
        "I","J","G","H","J","J","J","J","J","J","I","J","J","H","F","H","H",
        "H","F","F","F","E","E","F","E","E","E","E","E","E","E","E","F","E",
        "D","D","D","D","D","A","D","D","E","E","D","C","D","E","E"]

line2=["C","C","C","F","F","D","F","F","H","H","H","H","H","J","J","I","J",
        "J","J","J","J","I","I","G","H","E","E","H","I","I","I","J","I","I",
        "J","J","J","B","E","G","I","H","G","J","D","I","H","E","C","C","D",
        "D","D","D","C","A","A","C","D","A","A","D","D","D","D","D","D","B",
        "@","D",">","@",">","D","B","@",">","B","?","@","@","C","C","D","D",
        "B","8","4","9","@","D","B","C","D","D","C","C","C","C","A"]

line3=["C","C","C","F","F","F","F","F","H","H","H","H","H","J","J","J","J",
       "J","J","J","J","J","J","J","J","J","J","J","G","I","I","J","J","J",
       "J","J","J","I","J","J","J","J","J","J","J","J","J","J","J","J","H",
       "J","J","J","I","J","J","I","J","I","J","J","J","J","I","J","J","J",
       "I","B","H","E","F","@","C","F","F","@","C","A","D","E","@","C","D",
       "E","C","C","C","D","C","C","C","C","C","D","C","C","C","@"]

line4=["@","C","C","F","F","F","F","D","H","H","H","H","H","J","J","J","J",
       "I","J","J","J","J","I","J","J","I","I","G","I","F","E","G","G",">",
       "G",">","B","G","G","G","/","=",">","B","B","B","D","?","B","D","C",
       "D",">","A","D","D","D","@","C","C","A","@","C","C","B","?","?","B",
       "C","C","C",">","&","5",">","B","B","B","9","@","9","9","@","D","5",
       "<",">","(","+",":","(",":",">",":","@","@","C","9",">","B"]

# Convert Phred quality scores (ASCII characters) to probabilities
phred1 = [10**-(ord(x)-33/10.0) for x in line1]
phred2 = [10**-(ord(x)-33/10.0) for x in line2]
phred3 = [10**-(ord(x)-33/10.0) for x in line3]
phred4 = [10**-(ord(x)-33/10.0) for x in line4]

# Round probabilities to 3 decimal places in scientific notation
rounded_phred1 = ['{:0.3e}'.format(x) for x in phred1]
rounded_phred2 = ['{:0.3e}'.format(x) for x in phred2]
rounded_phred3 = ['{:0.3e}'.format(x) for x in phred3]
rounded_phred4 = ['{:0.3e}'.format(x) for x in phred4]

# Make a dataframe of the probability assignment and original quality string
df = pd.DataFrame({'quality_string1':line1,'rounded_phred1':rounded_phred1,
                   'quality_string2':line2,'rounded_phred2':rounded_phred2,
                   'quality_string3':line3, 'rounded_phred3':rounded_phred3,
                   'quality_string4':line4, 'rounded_phred4':rounded_phred4})

df = df.reindex(['quality_string1', 'rounded_phred1',
                 'quality_string2', 'rounded_phred2',
                 'quality_string3', 'rounded_phred3',
                 'quality_string4', 'rounded_phred4'], axis=1)
print(df)
df.to_csv('quality_scores.csv', index=True, encoding='utf-8')
