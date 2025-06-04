import pandas as pd

df = pd.read_csv('TU_result/TU_coverage_list.csv')
df = df[df['coverage']>10]
df['gene_number'] = df['TU'].apply(lambda x: len(x.split('|')))
gene_id =set()
ge_vs_tu_dict={}
for idx, line in df.iterrows():
    TU_list = line['TU'].split('|')
    for item in TU_list:
        gene_id.add(item)
        if item not in ge_vs_tu_dict:
            ge_vs_tu_dict[item]=line['TU']
        else:
            if len(TU_list) > len(ge_vs_tu_dict[item].split('|')):
                ge_vs_tu_dict[item]=line['TU']
longest_tu =pd.DataFrame(list(ge_vs_tu_dict.values()))
longest_tu.drop_duplicates(inplace=True)
longest_tu.to_csv('longest_tu.csv',index=False)
print(1)
