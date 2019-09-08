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


try:
    import numpy as np
    from Levenshtein import distance
    from scipy.spatial.distance import pdist, squareform
except ImportError:
    has_pkgs = False
else:
    has_pkgs = True

try:
    import igraph
except ImportError:
    has_igraph = False
else:
    has_igraph = True


def get_cdr3(adata):
    df = adata.uns['pyvdj']['df']
    samples = adata.uns['pyvdj']['samples'].values()

    cdr3_dict = {}
    for s in samples:
      
      cdr3 = df.loc[(df['sample'] == s), 'cdr3']
      cdr3 = cdr3.unique()
      cdr3 = cdr3.tolist()
      try:
          cdr3.remove('None')
      except ValueError:
          pass
    
      cdr3_dict[s] = cdr3

    return cdr3_dict


def get_dist(cdr3_dict, sample):
    """Requires numpy, scipy and Levenshtein."""
    if not has_pkgs:
        raise ImportError("Function requires numpy, scipy and python-Levenshtein packages.")

    cdr3 = cdr3_dict[sample]
    cdr3 = np.array(cdr3).reshape(-1, 1)
    dist = pdist(cdr3, lambda x,y: distance(x[0], y[0]))

    dist = (squareform(dist))

    return dist


def graph_cdr3(dist):
    """Requires igraph."""
    if not has_igraph:
        raise ImportError("Function requires the igraph-python package.")

    dist[dist > 1] = 0 # use Levenshtein distance of 1
    g = igraph.Graph.Adjacency((dist > 0).tolist())
    
    return g
