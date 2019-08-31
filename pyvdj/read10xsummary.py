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
