def vdj_load(paths, samples, adata = None, add_obs = True):
    # Loads 10x V(D)J sequencing data into a dictionary containing a dataframe.
    # If anndata specified, returns it with an .uns['pyvdj'] slot, else
    # returns the dictionary.
    # add_obs: whether to add some default .obs metadata columns
    import pandas as pd
    from pyVDJ.pyvdj.vdj_add_obs import vdj_add_obs
    cat_df = pd.DataFrame()
    for f in paths:
        df = pd.read_csv(f)
        df['barcode_meta'] = df['barcode'] + "_" + samples[f]

        df['clonotype_meta'] = df['raw_clonotype_id'] + "_" + samples[f]
        df['is_clone'] = ~df['raw_clonotype_id'].isin(['None'])    
        df = df.loc[df['is_cell'] == True] # filter step
    
        cat_df = pd.concat([cat_df, df], ignore_index=True)


    # Productive cells:
    d = {'True': True, 'False': False, 'None': False}
    cat_df['productive'] = cat_df['productive'].map(d)
    is_productive = cat_df.groupby('barcode_meta')['productive'].apply(lambda g: all(g))
    product_dict = dict(zip(is_productive.index, is_productive))
    cat_df['productive_all'] = cat_df['barcode_meta']
    cat_df['productive_all'].replace(to_replace=product_dict, inplace=True)


    vdj_dict = {'df':cat_df, 'samples':samples, 'obs_col':'vdj_obs'}

    
    if adata == None:
        return vdj_dict
    else:
        adata.uns['pyvdj'] = vdj_dict

        if add_obs: # then make a few default metadata columns ('vdj_...')
            adata = vdj_add_obs(adata, obs = ['has_vdjdata'])
        return adata
