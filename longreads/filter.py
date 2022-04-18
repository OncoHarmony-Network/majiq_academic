import pandas as pd
from pygtftk.gtf_interface import GTF
from gtfparse import read_gtf
import gffpandas.gffpandas as gffpd

tpm_path = "/Users/seongwoohan/Desktop/quant2.sf"
gtf_path = "/Users/seongwoohan/Desktop/gencode.v39.annotation.gtf"
gff3_path = "/Users/seongwoohan/Desktop/gencode.v39.annotation.gff3"

# This extracts the TPM < 1 
df_sf2 = pd.read_csv(tpm_path, sep='\t',engine='python')
df_sf2 = df_sf2[df_sf2['TPM'] < 1]
df_sf2_name = df_sf2['Name']
df_sf2_name = set(df_sf2_name)
#print(df_sf2_name)
a = list(df_sf2_name)

df = pd.read_csv(gff3_path, sep='\t',engine='python')
print(df)
# drop this row if it has this transcript_id
#df_gtf = pd.read_csv(gtf_path, sep='\t', engine='python', header=None, skiprows=5)
# print(df_gtf.head())
# print(df_gtf[8].head())
# print(df_gtf.loc[2][8])
# print()

# df_gtf = df_gtf[df_gtf[8].isin(a)]
# print(df_gtf.tail())

# contain this row if it has this transcript_id
# df_gtf2 = read_gtf(gtf_path)
# df_gtf_transcript = set(df_gtf2['transcript_id'])
# print(df_gtf_transcript)
# print()

# filtered_transcript = df_gtf_transcript.difference(df_sf2_name)
# #filtered_transcript = pd.Series(list(filtered_transcript))
# print(filtered_transcript)

# for s in filtered_annotation:
#     assert s in df_gtf_transcript
#     assert s in df_sf2_name


