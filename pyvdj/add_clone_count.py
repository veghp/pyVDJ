def add_clone_count(adata):
    # Number of clones in clonotype
    ## add check of adata.obs['vdj_clonotype']
    clone_count = adata.obs['vdj_clonotype'].value_counts()    
    clone_count_dict = dict(zip(clone_count.index, clone_count))

    adata.obs['vdj_clone_count'] = adata.obs['vdj_clonotype']
    adata.obs['vdj_clone_count'].replace(to_replace=clone_count_dict, inplace=True)

    adata.obs['vdj_clone_count'].loc[adata.obs['vdj_is_clone'] == False] = 1 # not a clone of other cells (unique)
    adata.obs['vdj_clone_count'].loc[adata.obs['vdj_has_vdjdata'] == False] = 0 # no data

    return adata
