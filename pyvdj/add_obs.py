# Copyright 2019 Peter Vegh
#
# This file is part of pyVDJ.
#
# pyVDJ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyVDJ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyVDJ.  If not, see <https://www.gnu.org/licenses/>.


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
    adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_clonotype'] = None

    tcr_dict = dict(zip(df['barcode_meta'], df['clonotype_meta']))
    adata.obs['vdj_clonotype'].replace(to_replace=tcr_dict, inplace=True)
    adata.obs['vdj_clonotype'] = adata.obs['vdj_clonotype'].astype('category')

    return adata


def add_is_clone(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_is_clone'] = adata.obs[obs_col]
    adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_is_clone'] = False
    clone_dict = dict(zip(df['barcode_meta'], df['is_clone']))
    adata.obs['vdj_is_clone'].replace(to_replace=clone_dict, inplace=True)

    return adata


def add_is_productive(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_is_productive'] = adata.obs[obs_col]
    adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_is_productive'] = None
    prod_dict = dict(zip(df['barcode_meta'], df['productive_all']))
    adata.obs['vdj_is_productive'].replace(to_replace=prod_dict, inplace=True)

    return adata


def add_chains(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    chains = adata.uns['pyvdj']['df']['chain'].unique()
    chains = [x for x in chains if str(x) != 'nan']

    chain_nested_dict = dict() # for making one .obs column for each chain type
    for c in chains:
        print(c)
        adata.uns['pyvdj']['df']['chain_' + c] = adata.uns['pyvdj']['df']['chain'] == c
        has_chain = adata.uns['pyvdj']['df'].groupby('barcode_meta')['chain_' + c].apply(lambda g: any(g))
        chain_nested_dict[c] = dict(zip(has_chain.index, has_chain))

    for c in chain_nested_dict.keys():
        adata.obs['vdj_chain_' + c] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_chain_' + c] = 'No_data'
        
        adata.obs['vdj_chain_' + c].replace(to_replace=chain_nested_dict[c], inplace=True)

    return adata


def add_genes(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    constant_genes = adata.uns['pyvdj']['df']['c_gene'].unique()
    constant_genes = [x for x in constant_genes if str(x) != 'nan']

    constant_genes_nested_dict = dict()
    for c in constant_genes:
        print(c)
        adata.uns['pyvdj']['df']['vdj_constant_' + c] = adata.uns['pyvdj']['df']['c_gene'] == c
        has_c_gene = adata.uns['pyvdj']['df'].groupby('barcode_meta')['vdj_constant_' + c].apply(lambda g: any(g))
        constant_genes_nested_dict[c] = dict(zip(has_c_gene.index, has_c_gene))

    for c in constant_genes_nested_dict.keys():
        adata.obs['vdj_constant_' + c] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_constant_' + c] = 'No_data'
        adata.obs['vdj_constant_' + c].replace(to_replace=constant_genes_nested_dict[c], inplace=True)

    return adata


def add_obs(adata, obs = ['has_vdjdata']):
    # obs: which of the below metadata to add?
    adder_functions = {
        'has_vdjdata': add_has_vdjdata,
        'clonotype': add_clonotype,
        'is_clone': add_is_clone,
        'is_productive': add_is_productive,
        'chains': add_chains,
        'genes': add_genes}

    for e in obs:
        func = adder_functions[e]
        adata = func(adata)

    return adata
