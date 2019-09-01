# pyVDJ

V(D)J sequencing data analysis

This package adds 10x Genomics V(D)J sequencing data to an [AnnData](http://readthedocs.io/en/stable/) object's `.uns` part, and also makes annotation columns in `.obs`.
This enables plotting various V(D)J properties and handling mRNA and V(D)J sequencing data together.


## Install

`pip3 install pyvdj`


## Usage

    import pyvdj
    adata = pyvdj.load_vdj(paths, samples, adata)
    adata = pyvdj.vdj_add_obs(adata, obs=['is_clone'])

For a detailed tutorial, see HTML LINK.


## Details

The package has 3 main functions that
* read `metrics_summary.csv` files into a pandas dataframe.
* load `filtered_contig_annotations.csv` files into an AnnData object.
* create annotations in the AnnData object.


### Read metrics

The `read10xsummary` function requires a list of paths to `metrics_summary.csv` files, and optionally a dictionary of path:samplename. It returns a dataframe of the metrics.


### Load V(D)J data

The `load_vdj` function loads 10x V(D)J sequencing data (`filtered_contig_annotations.csv` files) into an AnnData object's `.uns['pyvdj']` slot, and returns the object. The `adata.uns['pyvdj']` slot is a dictionary with the following elements:
* `'df'`: a dataframe containing V(D)J data
* `'obs_col'`: the `anndata.obs` columname of matching cellnames.
* `'samples'`: a dictionary of filename:samplename

More entries will be added in the future. If an anndata object is not supplied, the function returns the dictionary.

Arguments:
* `paths`: list of paths to filtered_contig_annotations.csv files.
* `samples`: a dictionary of path:samplename.
* `adata`: the AnnData object.
* `add_obs`: whether to add some default .obs metadata columns.


### Add annotations

The `adata.uns['pyvdj']['df']` is a pandas dataframe of the V(D)J data, with two additional columns that contain unique cell barcode and clonotype labels. These are generated using the user-supplied sample names: `cellbarcode + '_' + samplename` and `clonotype + '_' + samplename`.

These unique cell names are used to match the V(D)J cells to the AnnData `.X` cells, using `adata.obs['vdj_obs']`. The user has to prepare this column using the cell barcodes and the sample names.

The `add_obs` function can add the following annotations:
* `'has_vdjdata'`: does the cell have V(D)J sequencing data?
* `'clonotype'`: add clonotype name
* `'is_clone'`: does it have a clone?
* `'is_productive'`: are all chains productive?
* `'chains'`: adds a boolean column for each chain
* `'genes'`: adds a boolean column for each constant gene

More will be added in the future.

