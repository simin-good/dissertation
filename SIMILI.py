import pandas as pd


path_complex = r"D:\gene_fusion\R_environment\data\alphafold3_fusion_interaction_residue_level.tsv"
path_fusion = r"D:\gene_fusion\R_environment\data\alphafold3_fusions_residue_level.tsv"

df_complex = pd.read_csv(path_complex,sep='\t')
df_fusion = pd.read_csv(path_fusion,sep='\t')

#先创建一个独一无二的id用于标识每一个融合的情况
df_complex['mergeid'] = df_complex['fusion']+'_'+df_complex['gene']+'_'+df_complex['res'].astype(str)
df_fusion['mergeid'] = df_fusion['fusion']+'_'+df_fusion['gene']+'_'+df_fusion['uniprot_res'].astype(str)

#选出我要组合的列
df_complex= df_complex[['fusion','mergeid','gene','aa','sasa_ch','sasa_cx','rsa_ch','rsa_cx','plddt']]
df_fusion= df_fusion[['mergeid','sasa_ch','sasa_split','rsa_ch','rsa_split','plddt']]

#根据mergeid来组合
df_merge = pd.merge(df_complex,df_fusion,how='left',on='mergeid',suffixes=('_complex', '_fusion'))

#算两个BSA
df_merge['complex_BSA'] = df_merge['sasa_ch_complex']-df_merge['sasa_cx']
df_merge['fusion_BSA'] = df_merge['sasa_split']-df_merge['sasa_ch_fusion']

#输出
pathout1 = r"D:\gene_fusion\R_environment\data\siminLI.xlsx"
pathout2 = r"D:\gene_fusion\R_environment\data\interface_both_af3_new.xlsx"
df_merge.to_excel(pathout1,index=False)
#输出过滤空值后的
df_merge.dropna().to_excel(pathout2,index=False)

#Complex BSA = sasa_ch - sasa_cx
#fusion_BSA = sasa_split - sasa_ch,