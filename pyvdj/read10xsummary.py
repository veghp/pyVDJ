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


def read10xsummary(paths, samples = None):
    # This function reads a list of 10x Cell Ranger metrics_summary.csv
    # filepaths into a pandas dataframe.
    # samples is a dictionary of path:sample
    cat_df = pd.DataFrame()
    for f in paths:
        df = pd.read_csv(f, thousands=r',')

        if samples == None:
            df.insert(0, 'Sample', f)
        else:
            df.insert(0, 'Sample', samples[f])

        cat_df = pd.concat([cat_df, df], ignore_index=True)

    return cat_df
