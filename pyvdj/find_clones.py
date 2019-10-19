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


def make_cdr3_nt_set(df_d):
    # For finding clones within donor. Based on stats.py make_cdr3set()
    df_d = df_d[df_d['productive']]
    cdr3sets = df_d.groupby('clonotype_meta')['cdr3_nt'].apply(frozenset)

    cdr3_values = cdr3sets.unique()
    cdr3_values = [n for n in cdr3_values if len(n) > 1]
    cdr3_codes = dict(zip(range(0, len(cdr3_values) ), cdr3_values))
    cdr3_codes_rev = dict(zip(cdr3_values, range(0, len(cdr3_values) )))

    for i, c in enumerate(cdr3sets):
        if len(c) == 1:
            cdr3sets[i] = cdr3sets.index[i]
              # with only 1 CDR3 seq, we cannot decide whether it's a clone
    cdr3_replace = cdr3sets.replace(to_replace=cdr3_codes_rev, inplace=False)

    cdr3_nt_set_dict = {
        'cdr3sets': cdr3sets,
        'cdr3_codes': cdr3_codes,
        'cdr3_codes_rev': cdr3_codes_rev,
        'cdr3_replace': cdr3_replace}

    return cdr3_nt_set_dict


def add_donor_clonotype(adata):
    obs_col = adata.uns['pyvdj']['obs_col']
    df = adata.uns['pyvdj']['df']

    adata.obs['vdj_donor_clonotype'] = adata.obs[obs_col]
    # Remove missing ones:
    adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), 'vdj_donor_clonotype'] = None

    tcr_dict = dict(zip(df['barcode_meta'], df['donor_clonotype_meta']))
    adata.obs['vdj_donor_clonotype'].replace(to_replace=tcr_dict, inplace=True)
    adata.obs['vdj_donor_clonotype'] = adata.obs['vdj_donor_clonotype'].astype('category')

    return adata


def find_clones(adata, sample_donor):
    # This function returns adata with annotation of clonotypes shared between
    # 10x samples within donor (organism, individual).
    # 'sample_donor': dict of sample:donor
    df = adata.uns['pyvdj']['df']
    df['donor'] = df['sample']
    df['donor'].replace(to_replace=sample_donor, inplace=True)

    df['donor_clonotype_meta'] = df['clonotype_meta']
    for d in df['donor'].unique():
        df_d = df.loc[df['donor'] == d]
        print(d)
        # Calculate CDR3_nt sets:
        cdr3_nt_set_dict = make_cdr3_nt_set(df_d)
        cdr3_replace = cdr3_nt_set_dict['cdr3_replace']
        df['donor_clonotype_meta'].replace(to_replace=cdr3_replace, inplace=True)

    df['donor_clonotype_meta'] = df['donor_clonotype_meta'].astype(str)
    df['donor_clonotype_meta'] = df['donor'] + '_' + df['donor_clonotype_meta']
    adata.uns['pyvdj']['df']['donor_clonotype_meta'] = df['donor_clonotype_meta']

    adata = add_donor_clonotype(adata)

    return adata
