# Microarray data analysis #

## Gene expression microarrays ##

For Affymetrix gene expression microarrays:
```
    cd /some/directory/with/microarrays
    raw = read_affy_gene_arrays('.CEL.gz', '/data/csb/pipeline/platforms/affy_hg_u133_plus_2/probes');
    raw.mean = rmabackadj(raw.mean);
    raw.mean = normalize_quantiles(raw.mean);
    gene_expr = pset_summary_median_polish(raw);
```

For Agilent gene expression microarrays:
```
    raw = read_agilent_gene_arrays(find_files('.txt'));
    raw.mean = normalize_quantiles(raw.mean);
    gene_expr = pset_summary_median_polish(raw);
```

Using custom probesets:
```
    >> load /worktmp/pipeline/platforms/affy_ht_hg_u133a/gene_probesets
    >> expr = uarray_expression_rma(data, affy_ht_hg_u133a_gene_probesets)
    expr =
        mean: [29595x191 double]
        meta: [1x1 struct]
```

From the resulting data structure, you can see that the raw probe intensities were transformed into gene expression levels for the ~30000 different genes named in the RefSeq annotations (these gene names include non-coding, hypothetical and pseudogenes). The rows of the gene expression matrix correspond to the known genes of the currently selected organism:
```
    >> organism.Genes
    ans =
                 Name: {29595x1 cell}
      TranscriptCount: [29595x1 double]
          Transcripts: [29595x31 double]
  		     EntrezID: [29595x1 double]
```

If a suitable probeset is not available for the genome build that you have selected, you can always build a new probeset. To do this, first load the microarray's probe information into the pipeline:
```
    >> load /worktmp/pipeline/platforms/affy_ht_hg_u133a/probes
```

Then you can create either gene, transcript or exon probesets by using one of the following commands:
```
    >> probesets = create_gene_probesets(affy_ht_hg_u133a_probes);
    >> probesets = create_transcript_probesets(affy_ht_hg_u133a_probes);
    >> probesets = create_exon_probesets(affy_ht_hg_u133a_probes);
```

The microarray probes will be mapped against the organism's currently selected transcriptome build, and new probesets will be constructed based on the mappings. All three probesets can be used with the uarray\_expression\_rma() and related functions.




## CGH microarrays ##

Given a probe definition data structure for a microarray platform, you can generate CGH probesets against a given genome version in the following way:
```
    >> load /worktmp/pipeline/platforms/agilent_hg_cgh_244a/probes
    >> cgh_probesets = create_cgh_probesets(agilent_hg_cgh_244a_probes)
    cgh_probesets =
            Type: 'Copy number'
        Organism: 'Homo sapiens'
         Version: 'RefSeq 38'
      Chromosome: [235893x1 double]
          Offset: [235893x1 double]
      ProbeCount: [235893x1 double]
          Probes: [235893x3 double]
```

For most commonly used CGH microarray platforms, these CGH probesets have been pre-generated and made available under pipeline/platforms, so you only need to do this step for new and exotic CGH arrays.

Many CGH microarrays have been designed so that each probe has a distinct sequence and is complementary to a distinct genomic position. In such CGH microarrays the "ProbeCount" field will always have a value of 1. Conversely, this field will have values greater than 1 if multiple probes have an identical sequence and interrogate the same genomic position.

With CGH probesets at hand, raw probe level CGH data can be visualized in a genomics viewer software such as IGV by exporting it as a track:
```
    >> tumor = query('gbm agilent cgh 244a hms', 'sample type ~ tumor');
    >> normal = query('gbm agilent cgh 244a hms', 'sample type ~ normal');
    >> [tumor, normal] = paired_samples(tumor, normal, 'Patient');
    >> load /worktmp/pipeline/platforms/agilent_hg_cgh_244a/cgh_probesets;
    >> cgh_logratio_track(realize(tumor), realize(normal), ...
           agilent_hg_cgh_244a_probesets, '~/cgh_probe_tracks.igv');
```

This will generate a file that contains IGV scatterplot tracks for every paired CGH sample in the provided dataset. The zero level of the CGH logratios will be normalized, but no segmentation or quantization will be performed.

If you wish to run segmentation on the CGH samples, you can use the function cgh\_segment\_fast() in the following manner:
```
    >> tumor = query('gbm agilent cgh 244a hms', 'sample type ~ tumor');
    >> normal = query('gbm agilent cgh 244a hms', 'sample type ~ normal');
    >> [tumor, normal] = paired_samples(tumor, normal, 'Patient');
    >> tumor = realize(tumor);
    >> normal = realize(normal);
    >> segments = cgh_segment_fast(tumor, normal, cgh_probesets);
    >> cn_seg_to_track(segments, 'gbm_copy_number.seg')
```

In the above example, we first split the dataset into two query sets, one containing all tumor samples, and one containing all adjacent normal samples. We then use the function paired\_samples() to draw sample pairs from these two groups so that both samples in every pair share the same patient ID. So essentially we generated a list of tumor - adjacent normal sample pairs. We then load the CGH probesets from the disk and run a segmentation algorithm on the raw CGH probe intensity data. Finally we export the calculated segments as a .seg copy number track onto the disk. This track can then be visualized with IGV or some other genomic visualization tool.

You can also tune the segmentation algorithm using a number of optional parameters that are documented under HELP CGH\_SEGMENT\_FAST.

If you prefer to use the circular binary segmentation algorithm by Olshen et al, you can use the segmentation function cgh\_segment\_cbs().

When using genomic visualization tools to display copy number tracks, make sure that you have selected the right genome build (usually hg19). The data will look bad if visualized against the wrong genome build.

You can also search for statistically significant copy number alterations:
```
    >> cna_significance_track(segments, agilent_hg_cgh_244a_probesets, ...
           'significant_cnv.igv')
```

This function will also export an IGV track, but this time the track will consist of p-values transformed as log10(1 / p). In other words, the track will contain high peaks where statistically significant copy number alterations are present.

If you're more interested in the percentage-wise recurrence of copy number alterations, you can use another function:
```
    >> cna_recurrence_track(segments, agilent_hg_cgh_244a_probesets, ...
           'cnv_recurrence.igv')
```







## MicroRNA expression microarrays ##

Given a probe definition data structure for a microarray platform, you can generate microRNA probesets against the currently selected microRNA annotations (see `organism.miRNA`) in the following manner:
```
    >> load /worktmp/pipeline/platforms/agilent_human_mirna_8x15k_v2/probes;
    >> mirna_probesets = create_mirna_probesets(agilent_human_mirna_8x15k_v2_probes)
    mirna_probesets =
             miRNA: {904x1 cell}
        ProbeCount: [904x1 double]
            Probes: [904x17 double]
          Organism: 'Homo sapiens'
              Type: 'miRNA expression'
```

For most commonly used microRNA arrays, these probesets have been pre-generated and made available under pipeline/platforms, so you only need to do this step for new and exotic miRNA arrays.

Given microRNA probesets, you can calculate microRNA expression levels from raw microarray samples in the following way (note that in this example we use the pre-built probesets):
```
    >> load /worktmp/pipeline/platforms/agilent_human_mirna_8x15k_v2/probesets;
    >> raw = realize(query('taylor*mirna raw'));
    >> mirna_expr = uarray_expression_rma(raw, agilent_human_mirna_8x15k_v2_probesets)
    mirna_expr =
        Mean: [904x142 double]
        Meta: [1x1 struct]
```

The raw microRNA array samples in this example were loaded from an existing dataset. If your raw samples have not yet been imported as a dataset, you can do so using `import_uarray_data()` (see HELP IMPORT\_UARRAY\_DATA for more details).

The rows of the `mirna_expr.Mean` matrix represent different mature miRNAs annotated for the currently selected organism. The rows are in the same order as the fields under `organism.miRNA`. Columns represent samples.




## Quality control ##

You can visualize the spatial distribution of probe intensities on the microarray slide by using the function render\_uarray\_intensities():
```
    >> raw = realize(query('TCGA/GBM Affy HT HG U133A raw', 1:10));
    >> render_uarray_intensities(raw, affy_ht_hg_u133a_probes, '~/fig_spatial');
```

Many microarray platforms contain so called control probes whose intensity values are not always kept in the raw sample files. So if your microarray data contains !NaN values, you can choose to either visualize those probes as either black or red points. By default those probes will be shown as black, but you can use the optional argument 'RedMissing' to visualize the missing probes as red pixels.

If a visual inspection of the spatially rendered probe intensities indicates that the microarray samples contain significant spatial artifacts, you can apply spatial detrending in the following fashion:
```
    >> raw = realize(query('TCGA/GBM Affy HT HG U133A raw', 1:10));
    >> load /worktmp/pipeline/platforms/affy_ht_hg_u133a/probes;
    >> raw = spatial_detrend(raw, affy_ht_hg_u133a_probes);
```

You can also visualize the probe intensity histogram for a microarray sample:
```
    >> uarray_intensity_hist(raw, '~/histogram');
```

Or you can compare the probe composition of two microarray platforms:
```
    >> load agilent_human_mirna_8x15k/probes
    >> load agilent_human_mirna_8x15k_v2/probes
    >> compare_probes(agilent_human_mirna_8x15k_v2_probes, ...
           agilent_human_mirna_8x15k_probes);
    Probe comparison results:
    - 7120 probes added
    - 4067 probes removed
    - 7013 probes changed position
```