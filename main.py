import kcounter as kc
from sklearn.feature_extraction import DictVectorizer
import pandas as pd
from collections import defaultdict as dd

clusters_df = pd.read_csv('mmseqs_REPEAT_REGIONS_1_cluster_FILTERED.tsv', sep = '\t', names=['0','1'])

sequences = open('Alongiglumis_CN58138_pseudomolecules.v1.fasta.pass.list.REPEAT_REGIONS.gff3.fas')
seq_dict = {}
n = 0
for line in sequences.readlines():
    n+=1
    if n%10000==0:
        print(n)
    if '>' in line:
        name = line.strip()[1:]
    else:
        seq = line.strip()
        seq_dict[name]=seq
sequences.close()

len_df = 0
file = open('out_treetable_v6.txt')
for line in file.readlines():
    if len(line.strip().split('\t')) > len_df:
        len_df = len(line.strip().split('\t'))
last_column = len_df
file.close()
df = pd.read_csv('out_treetable_v6.txt', sep = '\t', names = [i for i in range(0,last_column)])
df = df.fillna(0)

rev_cols = list(df.columns)
print('rev_cols',rev_cols)
rev_cols.reverse()
for col in rev_cols:
    if len(set(df[col].values)) == 1:
        diff_column = int(col)+1
        # print('diff_column',diff_column)
        break

def make_groups(df, column): #создает словарь, где ключ - клада, значение - список строк, которые в нее входят 
    grps = {}
    for i, row in df.iterrows():
        if row[column] not in grps:
            grps[row[column]] = []
        grps[row[column]].append(i)
    return(grps)

def make_tree(df,groups, column, last_col=last_column):
    node = {}
    for k, v in groups.items():
        for col in (range(column, last_col)): # recording the name of specific clade in which we look for kmers
            if (df.iloc[v,col].unique()[0]) == 0:
                needed_column=col
                break
            elif col==last_col:
                needed_column=col
                break
            elif len(df.iloc[v,col].unique()) == 1:
                continue
            else:
                needed_column = col
                break
        subgroups = make_groups(df.iloc[v], needed_column)
        if len(subgroups) == 1 or int(column) == last_col:
            for key, v in subgroups.items(): #k is taken from previous group index
                node[k] = v
        else:
            add_node = make_tree(df, subgroups, column + 1)
            node[k] = add_node
    return(node)

def make_matrix(group_item):
    rows = group_item[1]
    clade_name = group_item[0]
    df_temp = df.iloc[rows]
    # return(df_temp)
    df_temp['cluster_full_name'] = df_temp[df_temp.columns[0:2]].apply(
    lambda x: ':'.join(x.dropna().astype(str)),
    axis=1)
    clusters_names = (df_temp['cluster_full_name'].tolist())
    cluster_temp = clusters_df[clusters_df.iloc[:,0].isin(clusters_names)]
    cluster_temp_listed = cluster_temp.groupby('0')['1'].apply(list).reset_index(name='1')
    cluster_seq_dict = {}
    for i, row in cluster_temp_listed.iterrows():
        cluster_seq_dict[row[0]] = row[1]
    seqs_in_cluster_list = []
    all_dicts = []
    cluster_names = []
    n = 0
    len_clusters = len(cluster_seq_dict)
    for k, v in cluster_seq_dict.items():
        cluster_kmers_amount_dict = dd(int)
        n+=1
        if n%20==0:
            print(f'{n} of {len_clusters}')
        cluster_names.append(k)
        for i in v:
            seqs_in_cluster_list.append(seq_dict[i])
        for seqq in seqs_in_cluster_list:
            kmers = kc.count_kmers(seqq, 13, canonical_kmers=True)
            for key in kmers.keys():
                cluster_kmers_amount_dict[key]+=1
        
        all_dicts.append(cluster_kmers_amount_dict)

    dictvectorizer = DictVectorizer(sparse=False)
    features = dictvectorizer.fit_transform(all_dicts)
    feature_name=dictvectorizer.get_feature_names_out()
    print(len(feature_name))
    print(features)



grp = make_groups(df, diff_column)
for i in grp.items():
    make_matrix(i)
    o=input()



# all_kmers_dicts_list = []
# seq_names = []
# sequences = open('Alongiglumis_CN58138_pseudomolecules.v1.fasta.pass.list.REPEAT_REGIONS.gff3.fas')
# clusters = open('mmseqs_REPEAT_REGIONS_1_cluster_FILTERED.tsv')
# counter = 0
# for line in sequences.readlines():
#     counter+=1
#     if counter%10000==0:
#         print('counter',counter)
#     if '>' in line:
#         name_seq = line.strip()[1:]
#     else:
#         seq = line.strip()
#         kmers_dict = kc.count_kmers(seq, 13)
#         seq_names.append(name_seq)
#         all_kmers_dicts_list.append(kmers_dict)
# print(len(all_kmers_dicts_list), len(seq_names))

# dictvectorizer = DictVectorizer(sparse=False)
# features = dictvectorizer.fit_transform(all_kmers_dicts_list)
# feature_name=dictvectorizer.get_feature_names_out()
# print(feature_name)