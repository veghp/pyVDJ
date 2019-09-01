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
    ## add check of adata.obs['vdj_clonotype']
    clon_bool = adata.obs['vdj_clonotype'][adata.obs['vdj_is_clone']].value_counts()

    clon_bool[clon_bool < n] = 'Small'
    clon_bool[clon_bool != 'Small'] = clon_bool[clon_bool != 'Small'].index.tolist()
    clon_dict = dict(zip(clon_bool.index, clon_bool))

    newcol = 'vdj_clonotype_gt_' + str(n)
    adata.obs[newcol] = adata.obs['vdj_clonotype']
    adata.obs[newcol] = adata.obs[newcol].cat.add_categories(['No_data', 'Small'])

    adata.obs.loc[(-adata.obs['vdj_is_clone']), newcol] = 'Small' # keep order of this
    adata.obs.loc[(-adata.obs['vdj_has_vdjdata']), newcol] = 'No_data' # and this
    adata.obs[newcol].replace(to_replace=clon_dict, inplace=True)

    return adata
