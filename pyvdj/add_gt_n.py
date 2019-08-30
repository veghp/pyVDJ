def add_gt_n(adata, n):
    # Clonotypes with not less than n clones
    ## add check of adata.obs['vdj_clonotype']
    clon_bool = adata.obs['vdj_clonotype'][adata.obs['vdj_is_clone']].value_counts()

    clon_bool[clon_bool < n] = 'Small'
    clon_bool[clon_bool != 'Small'] = clon_bool[clon_bool != 'Small'].index.tolist()
    clon_dict = dict(zip(clon_bool.index, clon_bool))

    newcol = 'vdj_clonotype_gt_' + str(n)
    adata.obs[newcol] = adata.obs['vdj_clonotype']
    adata.obs[newcol] = adata.obs[newcol].cat.add_categories(['No_data', 'Small'])

    adata.obs[newcol][-adata.obs['vdj_is_clone']] = 'Small' # keep order of this
    adata.obs[newcol][-adata.obs['vdj_has_vdjdata']] = 'No_data' # and this
    adata.obs[newcol].replace(to_replace=clon_dict, inplace=True)
    adata.obs[newcol].unique()
    
    return adata
