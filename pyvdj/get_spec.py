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


def get_spec(adata, clonotypes):
    df = adata.uns['pyvdj']['df']
    clonotypes = list(set(clonotypes))

    cdr3_dict = {}
    for c in clonotypes:
        print(c)
        cdr3 = df.loc[df['clonotype_meta'] == c, 'cdr3']
        cdr3 = cdr3.unique()
        cdr3 = cdr3.tolist()
        try:
            cdr3.remove('None')
        except ValueError:
            pass
        cdr3_dict[c] = cdr3

    return cdr3_dict
