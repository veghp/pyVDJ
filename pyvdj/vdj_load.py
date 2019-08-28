def vdj_load(paths, samples, anndata = None):
    import pandas as pd
    cat_df = pd.DataFrame()
    for f in paths:
        df = pd.read_csv(f)
        df['barcode_meta'] = df['barcode'] + "_" + samples[f]

        cat_df = pd.concat([cat_df, df], ignore_index=True)

    return cat_df
