# pyVDJ

V(D)J sequencing data analysis

This package adds 10x Genomics V(D)J sequencing data to an [AnnData](https://anndata.readthedocs.io) object's `.uns` part, and also makes annotation columns in `.obs`.
This enables plotting various V(D)J properties and handling mRNA (GEX) and V(D)J sequencing data together.


## Install

`pip3 install pyvdj`


## Usage

    import pyvdj
    adata = pyvdj.load_vdj(paths, samples, adata)
    adata = pyvdj.add_obs(adata, obs=['is_clone'])

For a detailed description, see the [tutorial](tutorials/pyVDJ_tutorial.html).


## Details

The package has functions that
* read `metrics_summary.csv` files into a pandas dataframe.
* load `filtered_contig_annotations.csv` files into an AnnData object.
* create various statistics and annotations in the AnnData object.


### Read metrics

The `read10xsummary` function requires a list of paths to `metrics_summary.csv` files, and optionally a dictionary of path:samplename. It returns a dataframe of the metrics.


### Load V(D)J data

The `load_vdj` function loads 10x V(D)J sequencing data (`filtered_contig_annotations.csv` files) into an AnnData object's `.uns['pyvdj']` slot, and returns the object. The `adata.uns['pyvdj']` slot is a dictionary which has the following elements:
* `'df'`: a dataframe containing V(D)J data
* `'obs_col'`: the `anndata.obs` columname of matching cellnames.
* `'samples'`: a dictionary of filename:samplename

If an anndata object is not supplied, the function returns the dictionary.

Arguments:
* `paths`: list of paths to filtered_contig_annotations.csv files.
* `samples`: a dictionary of path:samplename.
* `adata`: the AnnData object.
* `add_obs`: whether to add some default .obs metadata columns.


### Add annotations

The `adata.uns['pyvdj']['df']` is a pandas dataframe of the V(D)J data, with two additional columns that contain unique cell barcode and clonotype labels. These are generated using the user-supplied sample names: `cellbarcode + '_' + samplename` and `clonotype + '_' + samplename`.

These unique cell names are used to match the V(D)J cells to the AnnData `.X` cells, using `adata.obs['vdj_obs']`. The user has to prepare this column using the cell barcodes and the sample names.

The `add_obs` function can add the following annotations:
* `'has_vdjdata'`: does the cell have V(D)J sequencing data?
* `'clonotype'`: add clonotype name
* `'is_clone'`: does it have a clone?
* `'is_productive'`: are all chains productive?
* `'chains'`: adds annotation (True, False, No_data) for each chain
* `'genes'`: adds annotation (True, False, No_data) for each constant gene
* `'v_genes'`: adds annotation (True, False, No_data) for each variable gene
* `'j_genes'`: adds annotation (True, False, No_data) for each joining gene
* `'clone_count'`: adds clone count annotation


### Definitions

* Clone: a cell whose TCR is identical to another cell, within the same individual (donor, organism).*
* Clonotype: a set of all cells with the same TCR in the same individual (donor). A clonotype can have 1 or more cells.**
* Clone count (of a clonotype): number of clones in the clonotype.
* Public TCR (or CDR3) sequence: these are common and occur in multiple (or all) donors.
* Private TCR (or CDR3) sequence: these are unique to one donor.
* Condition-specific TCR (or CDR3) sequence: these occur in donors with a condition (disease, treated etc). These are private (unique) to the condition.

_The above definitions are understood in the context of the sequenced cells._

*As determined by Cell Ranger.
**Note that Cell Ranger v2 does not assign a clonotype id to clonotypes with only 1 clone, but uses ‘None’. Cell Ranger v3 does assign a clonotype id to all cells.


### CDR3 specificity

We can retrieve CDR3 amino acid sequences for given clonotypes using

    pyvdj.get_spec(adata, clonotypes = [clonotype1_sampleA', 'clonotype3_sampleB'])

which returns a dictionary. This can be used to find specificity in CDR3 databases, such as [VDJdb](http://vdjdb.cdr3.net) or [McPAS-TCR](http://friedmanlab.weizmann.ac.il/McPAS-TCR/).


### Clonotype statistics

We can generate and [plot various statistics](tutorials/pyVDJ_tutorial.html) on clonotypes and diversity.

    adata = pyvdj.stats(adata, meta)

This function adds a dictionary of statistics on the VDJ data (`adata.uns['pyvdj']['stats'][meta]`),
grouped by categories in the `adata.obs[meta]` column. Keys:

* `'meta'` stores the adata.obs columname
* `'cells'` count of cells, and cells with VDJ data per category
* `'clonotype_counts'` number of different clonotypes per category
* `'clonotype_dist'` clone count distribution
* `'shared_cdr3'` dictionary of cdr3 - cell


### Public and private CDR3 sequences

We can [find TCR-specificity shared between samples](tutorials/pyVDJ_tutorial.html), donors or any other annotation category.

    adata = pyvdj.find_clones(adata, sample_dict)

This function returns AnnData with clonotype annotation, where clonotypes shared between 10x samples within donor (organism, individual) are combined to have the same clonotype ID.
`'sample_dict'` is a dictionary of sample:donor, matching 10x samples (channels, as specified when the 10x VDJ data was loaded) to donors.


### CDR3-similarity graph

A set of prototype functions build CDR3-similarity graphs using [Levenshtein distances](https://en.wikipedia.org/wiki/Levenshtein_distance). The nodes are the CDR3 sequences, and edges connect nodes with Levenshtein distance of 1.

    cdr3_dict = pyvdj.get_cdr3(adata)  # get CDR3s for each sample
    dist = pyvdj.get_dist(cdr3_dict, sample)  # calculate distances (adjacency matrix)
    g = pyvdj.graph_cdr3(dist)  # returns an igraph graph object.

This requires the python-Levenshtein and the igraph-python packages.


### Dependencies

The package was originally developed for data made with Cell Ranger v2.1.1 (Chemistry: Single Cell V(D)J; V(D)J reference: GRCh38-alts-ensembl) and has been tested to work with Cell Ranger v3.1.0 data, with the following Python (v3.6.9) package versions:

    pandas 0.25.1
    anndata 0.6.21
    scanpy 1.4.3

