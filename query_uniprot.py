import numpy as np
import pandas as pd
from bioservices import UniProt 


# get gene ontologies
def query_uniprot(species_id='SBW25',all_locus_tags='PFLU_5352',search_columns='genes,protein names,comment(FUNCTION),go,go(biological process),go(molecular function)'):
    '''
    species_id: identifier for the species given by uniprot
    all_locus_tags: list of locus tags for query
    search_columns: for the valid column names refer to the website: https://www.uniprot.org/help/uniprotkb_column_names
    '''
    query_search = '\"' + species_id + '\" AND ('
    for locus_tag in all_locus_tags[0:-1]:
        query_search = query_search + '\"' + locus_tag + '\" OR '
    query_search = query_search + '\"' + all_locus_tags[-1] + '\")'
    up = UniProt()
    search_result = up.search(query_search, frmt='tab',columns=search_columns)
    str_up_ALL = search_result.split('\n')
    ls_up = []
    for each_line in str_up_ALL[1:]:
        ls_up.append(each_line.split('\t'))
    df_up = pd.DataFrame(ls_up[0:-1])
    df_up.columns = ['gene','protein','function','go','go_bp','go_mf'] # str_up_ALL[0].split('\t')
    return df_up

def process_query(tags, search_names, locus_tags_keep, C, Csorted):
    # uniprot
    
    # can only simultaneously query up to ~400 tags, hence the split below
    df_up = query_uniprot(all_locus_tags=tags,search_columns=search_names)

    # the gene column of df_up contains both gene name and locus tag, split those up (don't run this twice)
    ls_genes = []
    ls_locus_tags = []
    for ii in range(len(df_up)): 
        split_str = df_up.gene[ii].split()
        if len(split_str) == 2: # single gene name and locus_tag
            ls_genes.append(split_str[0])
            ls_locus_tags.append(split_str[1])
        elif len(split_str) == 1: # no gene name, only locus_tag
            ls_genes.append('N/A')
            ls_locus_tags.append(split_str[0])
        elif len(split_str) == 3: # two gene names and locus_tag
            ls_genes.append(split_str[0]+' '+split_str[1])
            ls_locus_tags.append(split_str[2])
        elif len(split_str) == 4: # three gene names and locus_tag
            ls_genes.append(split_str[0]+' '+split_str[1]+' '+split_str[2])
            ls_locus_tags.append(split_str[3])

    df_up['gene'] = ls_genes
    df_up['locus_tag'] = ls_locus_tags

    df_up

    row_list = []
    for ii,tag in enumerate(tags):
        tag_row = np.where(df_up.locus_tag == tag)[0][0]
        row_list.append(df_up.loc[tag_row])

    df_up_new = pd.DataFrame(row_list)
    df_up_new.reset_index(drop=True, inplace=True)

    # if gene name is n/a, replace it with the protein name
    for ii in range(len(df_up_new)):
        if df_up_new.gene[ii] == 'N/A':
            df_up_new.gene[ii] = df_up_new.protein[ii]

    # get observability rank in percentage from the sampling weights
    tag_inds = []
    for tag in tags:
        ind = locus_tags_keep.index(tag)
        tag_inds.append(ind)

    rank_tag_inds = []
    for ii in tag_inds:
        rank_tag_inds.append(len(C) - list(Csorted).index(C[ii]))     
        rank_tag_inds_per = [np.round((rank_ind-1) / len(C) * 100,2) for rank_ind in rank_tag_inds]
        
    # concatenate the ranks to the df
    df_up_new = pd.concat([df_up_new,pd.DataFrame({'obs_rank':rank_tag_inds_per,'tag_ind':tag_inds})],axis=1)

    return df_up_new
