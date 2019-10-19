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


import pandas as pd


def make_cdr3set(df):  # for finding public CDR3s
    # df is adata.uns['pyvdj']['df']
    # Creates a set of productive CDR3s for each cell in df
    df = df[df['productive']]
    cdr3sets = df.groupby('barcode_meta')['cdr3'].apply(frozenset)

    cdr3_values = cdr3sets.unique()
    cdr3_codes = dict(zip(range(0, len(cdr3_values) ), cdr3_values))
    cdr3_codes_rev = dict(zip(cdr3_values, range(0, len(cdr3_values) )))

    cdr3set_dict = {
        'cdr3sets': cdr3sets,
        'cdr3_codes': cdr3_codes,
        'cdr3_codes_rev': cdr3_codes_rev}

    return cdr3set_dict


def stats(adata, meta):
    obs_col = adata.uns['pyvdj']['obs_col']
    # This function returns a dictionary of statistics on the VDJ data,
    # grouped by categories in the adata.obs[meta] column. Keys:
    # 'meta' stores the adata.obs columname
    # 'cells' count of cells, and cells with VDJ data per category
    # 'clonotype_counts' number of different clonotypes per category
    # 'clonotype_dist' clone count distribution
    # 'shared_cdr3' dictionary of cdr3 - cell
    stats_dict = {}
    stats_dict['meta'] = meta

    n_cells = adata.obs.groupby(meta).size()
    n_hasvdj = adata.obs.groupby(meta)['vdj_has_vdjdata'].apply(sum)
    stats_dict['cells'] = [n_cells, n_hasvdj]

    n_cl = adata.obs.groupby(meta)['vdj_clonotype'].unique()
    n_cl = n_cl.apply(lambda x: ([z for z in x if str(z) != 'nan']))
    len_cl = [len(y) for y in n_cl]
    dict_cl = dict(zip(adata.obs[meta].unique(), len_cl))
    stats_dict['clonotype_counts'] = dict_cl

    dist_cl = adata.obs.groupby(meta)['vdj_clone_count'].value_counts()
    stats_dict['clonotype_dist'] = dist_cl

    if 'cdr3' not in adata.uns['pyvdj']:
        adata.uns['pyvdj']['cdr3'] = make_cdr3set(adata.uns['pyvdj']['df'])
        print('Adding CDR3 to adata.uns')

    # 'shared_cdr3'
    cdr3set_dict = adata.uns['pyvdj']['cdr3']
    cdr3sets = cdr3set_dict['cdr3sets']
    # Group anndata cells (obs_col) by meta cats
    grouped = adata.obs.groupby(meta)
    # For each meta cat cellist: subset cdr3df to these cells, then
    cdr3_dict = {}
    for name, group in grouped:
        group[obs_col].values
        cdr3_list = cdr3sets[cdr3sets.index.isin(group[obs_col].values)]

        cdr3_dict[name] = cdr3_list.unique()  # keys of cdr3_dict
# and their copynumber (cdr3 clones) within groups -- postponed for later
        # Return simplified sets:
        cdr3_dict[name] = cdr3_list.unique()

    stats_dict['cdr3'] = cdr3_dict

    if 'stats' not in adata.uns['pyvdj']:
        adata.uns['pyvdj']['stats'] = {}
    adata.uns['pyvdj']['stats'][meta] = stats_dict

    return adata
