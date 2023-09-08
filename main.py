import kcounter as kc
from sklearn.feature_extraction import DictVectorizer
import pandas as pd
from collections import defaultdict as dd
import sys
import numpy as np
# import matplotlib.pyplot as plt
import json


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
df['cluster_full_name'] = df[df.columns[0:2]].apply(
    lambda x: ':'.join(x.dropna().astype(str)),
    axis=1)
all_clusters_names = (df['cluster_full_name'].tolist())
df = df.drop(columns=['cluster_full_name'])
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
            elif col==last_col-1:
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

def get_leaves(branch):
    global leaves
    if type(branch) is list:
        leaves.extend(branch)
    else:
        for k, v in branch.items():
            get_leaves(v)



def make_matrix(group_item):
    max_kmer = 0
    rows = group_item[1]
    clade_name = group_item[0]
    print(df.shape[0], len(rows))
    # o=input()
    df_temp = df.iloc[rows]
    # return(df_temp)
    df_temp['cluster_full_name'] = df_temp[df_temp.columns[0:2]].apply(
    lambda x: ':'.join(x.dropna().astype(str)),
    axis=1)
    clusters_names = (df_temp['cluster_full_name'].tolist())
    # print('103\n',clusters_names)
    cluster_temp = clusters_df[clusters_df.iloc[:,0].isin(clusters_names)]
    cluster_temp_listed = cluster_temp.groupby('0')['1'].apply(list).reset_index(name='1')
    cluster_seq_dict = {}
    print('cluster_temp', cluster_temp.shape[0])
    print('cluster_temp_listed', cluster_temp_listed.shape[0])
    clusters_sizes = {}
    for i, row in cluster_temp_listed.iterrows():
        clusters_sizes[row[0]] = len(row[1])
        cluster_seq_dict[row[0]] = row[1]
    all_dicts = []
    cluster_names = []
    n = 0
    len_clusters = len(cluster_seq_dict)
    print('cluster_seq_dict len', len(cluster_seq_dict))
    for cl in clusters_names: # поменять словарь откуда кластеры брать на тот, который соответствует порядку кластеров в таблице
        cluster_kmers_amount_dict = dd(int)
        seqs_in_cluster_list = []
        n+=1
        if n%20==0:
            print(f'{n} of {len_clusters}')
        cluster_names.append(cl)
        for i in cluster_seq_dict[cl]:
            seqs_in_cluster_list.append(seq_dict[i])
        for seqq in seqs_in_cluster_list:
            kmers = kc.count_kmers(seqq, 13, canonical_kmers=True)
            for key in kmers.keys():
                cluster_kmers_amount_dict[key]+=1
        # if max_kmer<max(cluster_kmers_amount_dict.values()):
        #     max_kmer = max(cluster_kmers_amount_dict.values())
        all_dicts.append(cluster_kmers_amount_dict)
    # print('120\n',cluster_names)
    # print(clusters_names==all_clusters_names) #список всех кластеров в согласии с индексами строк в таблице клад
    # o=input()
    dictvectorizer = DictVectorizer(sparse=False, dtype=np.int16)
    # print(max_kmer, 'max kmer')
    features = dictvectorizer.fit_transform(all_dicts)
    feature_name=dictvectorizer.get_feature_names_out()
    print('first amount of kmers ',len(feature_name))
    # print(features)

    sum_arr = np.sum(features, axis=0)
    sum_list = sum_arr.tolist()
    columns_to_del = np.where(sum_arr < 4)[0]
    features = np.delete(features, columns_to_del, axis=1)
    # feature_name=feature_name.tolist()
    # columns_to_del=columns_to_del.tolist()
    # del feature_name[columns_to_del]
    feature_name = np.delete(feature_name, columns_to_del)
    sum_arr = np.delete(sum_arr, columns_to_del)
    print(features)
    print('kmers after filtering ', len(feature_name))
    # sum_dict = {feature_name[i]:sum_list[i] for i in range(len(sum_list))}
    # with open("kmers_freq.json", "w") as file:
    #     json.dump(sum_dict, file)
    # file.close()
    # plt.hist(sum_list)
    # plt.xlabel('kmer_freq')
    # plt.ylabel('amount')
    # plt.savefig('hist_kmer_freq.png')

    # max_ind = (sum_arr.argmax())
    # print(max_ind, 'index\n', sum_arr[max_ind], 'max value')
    # cluster_names = {v:i for i,v in enumerate(cluster_names)} # для итерации по строкам массива
    return(features, feature_name, clusters_names, sum_arr, clusters_sizes)

def count_uniq(tree, matrix, feature_names, clusters_names, sum_arr, clusters_sizes):
    for k,v in tree.items():
        get_leaves(v)
        temp_matrix = matrix[leaves]
        clade_sum_arr = np.sum(temp_matrix, axis=0)
        current_clade_size = sum([clusters_sizes[clusters_names[i]] for i in leaves])
        # print('leaves', leaves)
        columns_to_del = np.where(clade_sum_arr < int(0.5*current_clade_size))[0] # фильтрация 60% топовых по чувствительности кмеров
        print('размер клады и threshold',current_clade_size, int(0.5*current_clade_size))
        # удаляем новые отфильтрованные столбцы
        clade_sum_arr = np.delete(clade_sum_arr, columns_to_del)
        temp_sum_arr = np.delete(sum_arr, columns_to_del)
        temp_matrix = np.delete(temp_matrix, columns_to_del, axis=1)

        division_sign = clade_sum_arr/temp_sum_arr
        div_sign_desc = np.argsort(-division_sign)[:100]


        # print(div_sign_desc)
        print(division_sign[div_sign_desc.tolist()]) # отбор значений чувствительности топовых по чувствительности кмеров
        print(clade_sum_arr[div_sign_desc.tolist()])
        print(temp_sum_arr[div_sign_desc.tolist()]) # абсолютные значения суммы данных кмеров во всем датасете
        o = input()

grp = make_groups(df, diff_column-1)
tree = make_tree(df, grp, diff_column)
# print(list(tree.keys())[0])
tree = tree[list(tree.keys())[0]]
for i in grp.items():
    # temp_tree = tree[i[0]]
    features, feature_name, clusters_names, sum_arr, clusters_sizes =  make_matrix(i)
    leaves = []
    count_uniq(tree, features, feature_name, clusters_names, sum_arr, clusters_sizes)
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