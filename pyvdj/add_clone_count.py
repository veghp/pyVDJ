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


def add_clone_count(adata):
    # Number of clones in clonotype
    ## add check of adata.obs['vdj_clonotype']
    clone_count = adata.obs['vdj_clonotype'].value_counts()
    clone_count_dict = dict(zip(clone_count.index, clone_count))

    adata.obs['vdj_clone_count'] = adata.obs['vdj_clonotype']
    adata.obs['vdj_clone_count'].replace(to_replace=clone_count_dict, inplace=True)

    adata.obs.loc[(adata.obs['vdj_is_clone'] == False), 'vdj_clone_count'] = 1
      # not a clone of other cells (unique)
    adata.obs.loc[(adata.obs['vdj_has_vdjdata'] == False), 'vdj_clone_count'] = 0 
      # no data

    return adata
