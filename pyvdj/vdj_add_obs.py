def add_has_vdjdata(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_has_vdjdata'] = adata.obs[obs_col].isin(df['barcode_meta'])
      # 'barcode_meta' is in the df and corresponds to 'obs_col' in anndata.
    return adata


def add_clonotype(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_clonotype'] = adata.obs[obs_col]
    # Remove missing ones:
    adata.obs['vdj_clonotype'][-adata.obs['vdj_has_vdjdata']] = None

    tcr_dict = dict(zip(df['barcode_meta'], df['clonotype_meta']))
    adata.obs['vdj_clonotype'].replace(to_replace=tcr_dict, inplace=True)

    return adata


def vdj_add_obs(adata, obs = ['has_vdjdata']):
    adder_functions = {
        "has_vdjdata": add_has_vdjdata,
        "clonotype": add_clonotype}

    for e in obs:
        func = adder_functions[e]
        adata = func(adata)
        
    return adata
