def add_has_vdjdata(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']
    adata.obs['vdj_has_vdjdata'] = adata.obs[obs_col].isin(df['barcode_meta'])
      # 'barcode_meta' is in the df and corresponds to 'obs_col' in anndata.
    return adata

#def vdj_add_obs(adata, obs = ['has_vdjdata'], pyvdj = None):
def vdj_add_obs(adata, obs = ['has_vdjdata']):
    # pyvdj is the dictionary returned by vdj_load(), uses adata.uns['pyvdj']
    # if unspecified.
#    if pyvdj == None:
#        pyvdj = adata.uns['pyvdj']

#    obs_col = pyvdj['obs_col']
#    df = pyvdj['df']
    
    adder_functions = {
        "has_vdjdata": add_has_vdjdata
        }
#        "x": 'add_'}

    for e in obs:
        func = adder_functions[e]
        adata = func(adata)
        
    return adata
