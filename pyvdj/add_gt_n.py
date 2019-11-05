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


def add_gt_n(adata, n):
    # Clonotypes with not less than n clones
    if (not isinstance(n, int) or n < 0):
        raise ValueError('n must be a non-negative integer')

    # These two adata.obs columns will be used:
    if 'vdj_clonotype' not in adata.obs.columns:
        raise Exception('Please add \'clonotype\' with add_obs()')
    if 'vdj_clone_count' not in adata.obs.columns:
        raise Exception('Please add \'clone_count\' with add_obs()')

    clon_bool = adata.obs['vdj_clonotype'].value_counts()
    txt_less = 'Less_than_n'
    clon_bool[clon_bool < n] = txt_less
    clon_bool[clon_bool != txt_less] = clon_bool[clon_bool != txt_less].index.tolist()
    clon_dict = dict(zip(clon_bool.index, clon_bool))

    newcol = 'vdj_clonotype_n_' + str(n)
    adata.obs[newcol] = adata.obs['vdj_clonotype']
    adata.obs[newcol] = adata.obs[newcol].cat.add_categories(['No_data', txt_less])
    adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), newcol] = 'No_data'
    adata.obs[newcol].replace(to_replace=clon_dict, inplace=True)

    return adata
