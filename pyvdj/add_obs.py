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

    if 'is_clone' in df.columns:
        adata.obs['vdj_is_clone'] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_is_clone'] = False
        clone_dict = dict(zip(df['barcode_meta'], df['is_clone']))
        adata.obs['vdj_is_clone'].replace(to_replace=clone_dict, inplace=True)
    else:
        if 'vdj_clone_count' not in adata.obs.columns:
            adata = add_clone_count(adata)
        adata.obs['vdj_is_clone'] = adata.obs['vdj_clone_count'] > 1

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
    grouped_cells = adata.uns['pyvdj']['df'].groupby('barcode_meta')
    for c in chains:
        print(c)
        adata.uns['pyvdj']['df']['chain_' + c] = adata.uns['pyvdj']['df']['chain'] == c
        has_chain = grouped_cells['chain_' + c].apply(any)
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
    grouped_cells = adata.uns['pyvdj']['df'].groupby('barcode_meta')
    for c in constant_genes:
        print(c)
        adata.uns['pyvdj']['df']['vdj_constant_' + c] = adata.uns['pyvdj']['df']['c_gene'] == c
        has_c_gene = grouped_cells['vdj_constant_' + c].apply(any)
        constant_genes_nested_dict[c] = dict(zip(has_c_gene.index, has_c_gene))

    for c in constant_genes_nested_dict.keys():
        print('Preparing annotation for %s' % c)
        adata.obs['vdj_constant_' + c] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_constant_' + c] = 'No_data'
        adata.obs['vdj_constant_' + c].replace(to_replace=constant_genes_nested_dict[c], inplace=True)

    return adata


def add_v_genes(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    v_genes = adata.uns['pyvdj']['df']['v_gene'].unique()
    v_genes = [x for x in v_genes if str(x) != 'nan']

    v_genes_nested_dict = dict()
    grouped_cells = adata.uns['pyvdj']['df'].groupby('barcode_meta')
    for v in v_genes:
        print(v)
        adata.uns['pyvdj']['df']['vdj_v_' + v] = adata.uns['pyvdj']['df']['v_gene'] == v
        has_v_gene = grouped_cells['vdj_v_' + v].apply(any)
        v_genes_nested_dict[v] = dict(zip(has_v_gene.index, has_v_gene))

    for v in v_genes_nested_dict.keys():
        adata.obs['vdj_v_' + v] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_v_' + v] = 'No_data'
        adata.obs['vdj_v_' + v].replace(to_replace=v_genes_nested_dict[v], inplace=True)

    return adata


def add_j_genes(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    j_genes = adata.uns['pyvdj']['df']['j_gene'].unique()
    j_genes = [x for x in j_genes if str(x) != 'nan']

    j_genes_nested_dict = dict()
    grouped_cells = adata.uns['pyvdj']['df'].groupby('barcode_meta')
    for j in j_genes:
        print(j)
        adata.uns['pyvdj']['df']['vdj_j_' + j] = adata.uns['pyvdj']['df']['j_gene'] == j
        has_j_gene = grouped_cells['vdj_j_' + j].apply(any)
        j_genes_nested_dict[j] = dict(zip(has_j_gene.index, has_j_gene))

    for j in j_genes_nested_dict.keys():
        adata.obs['vdj_j_' + j] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_j_' + j] = 'No_data'
        adata.obs['vdj_j_' + j].replace(to_replace=j_genes_nested_dict[j], inplace=True)

    return adata


def add_clone_count(adata):
    # Number of clones in clonotype
    if 'vdj_clonotype' not in adata.obs.columns:
        adata = add_clonotype(adata)

    clone_count = adata.obs['vdj_clonotype'].value_counts()
    clone_count_dict = dict(zip(clone_count.index, clone_count))

    adata.obs['vdj_clone_count'] = adata.obs['vdj_clonotype']
    adata.obs['vdj_clone_count'].replace(to_replace=clone_count_dict, inplace=True)

    adata.obs.loc[(adata.obs['vdj_has_vdjdata'] == False), 'vdj_clone_count'] = 0
      # no data
    adata.obs['vdj_clone_count'] = adata.obs['vdj_clone_count'].astype(int)

    return adata


def add_obs(adata, obs):
    # obs: list. which of the below metadata to add?
    adder_functions = {
        'has_vdjdata': add_has_vdjdata,  # load_vdj() adds this by default
        'clonotype': add_clonotype,
        'is_clone': add_is_clone,
        'is_productive': add_is_productive,
        'chains': add_chains,
        'genes': add_genes,
        'v_genes': add_v_genes,
        'j_genes': add_j_genes,
        'clone_count': add_clone_count,
        }

    for e in obs:
        func = adder_functions[e]
        adata = func(adata)

    return adata
