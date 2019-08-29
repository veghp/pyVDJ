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


def add_is_clone(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_is_clone'] = adata.obs[obs_col]
    adata.obs['vdj_is_clone'][-adata.obs['vdj_has_vdjdata']] = False
    clone_dict = dict(zip(df['barcode_meta'], df['is_clone']))
    adata.obs['vdj_is_clone'].replace(to_replace=clone_dict, inplace=True)

    return adata


def add_is_productive(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_is_productive'] = adata.obs[obs_col]
    adata.obs['vdj_is_productive'][-adata.obs['vdj_has_vdjdata']] = None
    prod_dict = dict(zip(df['barcode_meta'], df['productive_all']))
    adata.obs['vdj_is_productive'].replace(to_replace=prod_dict, inplace=True)

    return adata


def add_chains(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    chains = adata.uns['pyvdj']['df']['chain'].unique()
    chains = [x for x in chains if str(x) != 'nan']
    
    chain_nested_dict = dict()
    for c in chains:
        print(c)
        adata.uns['pyvdj']['df']['chain_' + c] = adata.uns['pyvdj']['df']['chain'] == c
        has_chain = adata.uns['pyvdj']['df'].groupby('barcode_meta')['chain_' + c].apply(lambda g: any(g))
        chain_nested_dict[c] = dict(zip(has_chain.index, has_chain))

    for c in chain_nested_dict.keys():
        adata.obs['vdj_chain_' + c] = adata.obs[obs_col]
        adata.obs['vdj_chain_' + c][-adata.obs['vdj_has_vdjdata']] = 'No_data'
        adata.obs['vdj_chain_' + c].replace(to_replace=chain_nested_dict[c], inplace=True)

    return adata


def vdj_add_obs(adata, obs = ['has_vdjdata']):
    adder_functions = {
        'has_vdjdata': add_has_vdjdata,
        'clonotype': add_clonotype,
        'is_clone': add_is_clone,
        'is_productive': add_is_productive,
        'chains': add_chains}

    for e in obs:
        func = adder_functions[e]
        adata = func(adata)
        
    return adata
