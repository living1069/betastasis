# Pipeline startup and basics #

On servers where the pipeline is installed, you can start it with the command "pipeline". This command starts a customized Matlab instance with access to functions for the analysis of biological data.

On startup, the pipeline automatically loads into memory some genetic information for Homo sapiens. This species specific information is necessary for running computational analyses. If you're interested in running analyses for some other organism, you can change the currently selected organism by typing:
```
    >> select_organism('mus musculus', '2009')
    Reading organism data into memory...
```

A list of available organisms can be viewed with the command:
```
    >> list_organisms
    List of available organisms:
    - Homo sapiens
    - Mus musculus
```

You should now select Homo sapiens again before proceeding with the tutorial:
```
    >> select_organism('homo sapiens', '2009')
```

You can see genetic information about the currently selected organism by accessing the global variable "organism":
```
    >> organism
    organism =
                 Name: 'Homo sapiens'
              Version: '2009'
          Chromosomes: [1x1 struct]
                Genes: [1x1 struct]
          Transcripts: [1x1 struct]
                Exons: [1x1 struct]
                miRNA: [1x1 struct]
            pre_miRNA: [1x1 struct]
                 SNPs: [1x1 struct]
           Ontologies: [1x1 struct]
```

This global variable is also used by many of the pipeline's internal components. This is why it is important to make sure that you have the correct organism selected when performing data analysis. Having said that, the pipeline does perform some sanity checks and will try to notify the user if the current organism does not seem to match with the data.

In this tutorial we're interested in ovarian cancer, so let's try retrieving some raw microarray data first. Here's an example of loading one of the TCGA ovarian cancer microarray data sets (dataset names are case insensitive):
```
    >> qset = query('TCGA/OV Affy HT HG U133A raw')
    qset =
          Type: 'Microarray probe intensities'
        Sample: [1x1 struct]
      Resource: {527x1 cell}
      Platform: {527x1 cell}
       Patient: [1x1 struct]
          Misc: [1x1 struct]
```

From the command output, you can see that this data set includes 527 microarray samples produced with the Affymetrix HT HG U133A gene expression microarray. The data structure you're seeing on the screen is known as a "query set". Query sets are used for filtering data sets down to the desired components, and for studying the metadata associated with data sets. Query sets notably do not contain the actual data, but instead only act as handles to the data. This is useful, because it allows you to perform preliminary filtering on huge data sets without loading gigabytes of data from the disk.

Now, let's have a look at what the query set contains. The "sample\_id" vector contains a list of IDs for the samples that were hybridized to the microarray in the 527 experiments. For this TCGA data set, the sample IDs match with the original TCGA sample IDs. Here are the first three:
```
    >> qset.Sample.ID(1:3)
    ans =
      'TCGA-13-0807-01'
      'TCGA-13-0891-01'
      'TCGA-13-0912-01'
```

The "meta" substructure contains clinical patient information associated with the microarray samples:
```
    >> qset.Patient
    ans =
           Gender: [527x1 char]
               ID: {527x1 cell}
           Status: {527x1 cell}
     SurvivalTime: [527x1 double]
	          ...
```

The "Misc" substructure contains miscellaneous pieces of clinical information that are not supported by the pipeline data model:
```
    >> qset.Misc
    ans =
          ADDITIONALCHEMOTHERAPY: {527x1 cell}
           ADDITIONALDRUGTHERAPY: {527x1 cell}
        ADDITIONALHORMONETHERAPY: {527x1 cell}
	                 ...
```

These fields can still be used for filtering samples, but the pipeline predicate query language does not support them. Thus, filtering based on these fields must be performed manually.

Now that we know what a query set consists of, we can try filtering our query set by using the predicate query language:
```
    >> female_qset = filter_query(qset, 'gender = female')
    female_qset =
          Type: 'Microarray probe intensities'
      Platform: {509x1 cell}
        Sample: [1x1 struct]
      Resource: {509x1 cell}
       Patient: [1x1 struct]
          Misc: [1x1 struct]
```

Here we narrowed our query set down to those samples that came from female patients. Now it should be noted that this data set is for ovarian cancer, so what's going on with the 18 samples that aren't from female patients? Well, let's find out by using the helpful function meta\_summary():
```
    >> meta_summary(qset)
	...
    Patient.Gender:
    - Female (509 items)
    - N/A (18 items)
	...
```

From the output of the function, you can see that the dataset contains 509 samples from female patients (as we saw earlier), and 18 samples from patients of unknown gender. You can see the same result by querying for all samples with unknown gender:
```
    >> qset = query('TCGA/OV Affy HT HG U133A raw', 'gender = unknown')
    qset =
          Type: 'Microarray probe intensities'
      Platform: {18x1 cell}
        Sample: [1x1 struct]
      Resource: {18x1 cell}
       Patient: [1x1 struct]
          Misc: [1x1 struct]
```

As you can see, the dataset contains 18 samples for which the patient gender is not reported in the clinical data. Nonetheless, it is probably a safe bet that these patients are female as well. We're talking about ovarian cancer, after all :)

In order to gain access to the data a query set refers to, the query set must be "realized" with the following command:
```
    >> data = realize(qset)
    Progress: 100%
    data =
        Mean: [247899x18 double]
        Meta: [1x1 struct]
```

A progress indicator will show up on the screen while the data is being retrieved from the storage. For this query set of 18 samples, the retrieval will happen very fast. Retrieving hundreds of microarray samples from the disk will take a couple of seconds, though. You can see this for yourself by realizing the entire microarray data set:
```
    >> data = realize(query('TCGA/OV Affy HT HG U133A'))
    Progress: 100%
    data =
        Mean: [247899x527 double]
        Meta: [1x1 struct]
```

And there you have it. You can now view the raw probe intensity data in the following fashion:
```
    >> data.Mean(1:3, 1:5)
    ans =
     1.0e+03 *
      0.5315    0.3180    0.3050    0.6150    0.2380
      3.2888    1.6508    1.8145    2.8865    1.3953
      6.3203    3.8602    3.8381    5.7050    3.9452
```

By looking at the realized data structure, you can also see that all of the original metadata from the query set is still available behind the attribute "Meta":
```
    >> data.Meta
    ans =
          Type: 'Microarray probe intensities'
      Platform: {527x1 cell}
        Sample: [1x1 struct]
       Patient: [1x1 struct]
          Misc: [1x1 struct]
```

Now, let's have a look at some gene expression data and see how query sets can be merged together. First we need two query sets of gene expression data, both from the TCGA project. One of the query sets will be for the glioblastoma gene expression profiles, one for the ovarian cancer profiles:
```
    >> qset_gbm = query('TCGA/GBM Affy HT HG U133A gene expression');
    >> qset_ov = query('TCGA/OV Affy HT HG U133A gene expression');
```

Since both query sets refer to gene expression data and were generated for the same organism and using the same annotations, the query sets can be merged together in the following fashion:
```
    >> merged = query_union(qset_gbm, qset_ov);
```

This function can take an arbitrary amount of query sets to merge, but in this example we called it with just two query sets. As the name of the function implies, a union of the two query sets is taken. Strictly speaking, the union is taken over the "Resource" field of the two query sets, since resource names uniquely identify pieces of data in the pipeline.

But let's now have a look at the merged query set:
```
    >> merged
    merged =
                     Misc: [1x1 struct]
                 Organism: 'Homo sapiens'
                  Patient: [1x1 struct]
                 Platform: {922x1 cell}
                 Resource: {922x1 cell}
                   Sample: [1x1 struct]
      SummarizationMethod: {922x1 cell}
                     Type: 'Gene expression'
```

The new query set contains both the 395 glioblastoma profiles and the 527 ovarian cancer profiles, for a total of 922 samples. The merged query set can be realized in the normal fashion, and the data resources will be automatically loaded from the two datasets:
```
    >> data = realize(merged)
    Progress: 100%
    data =
      Mean: [29595x922 double]
      Meta: [1x1 struct]
```

When analyzing the merged data set, you can use metadata to differentiate between tumor types:
```
    >> unique(data.Meta.Patient.TumorType)
    ans =
      '-'
      'Serous Cystadenocarcinoma'
      'Treated primary GBM'
      'Untreated primary (De Nova) GBM'
```

The symbol '-' always stands for "unknown" in the metadata.

After filtering and working on datasets, you can store it in the pipeline data store as a fresh copy:
```
    >> expr = realize(query('TCGA/GBM Affy HT HG U133A gene expression', 1:10));
    >> expr.Mean = quantilenorm(expr.Mean);
    >> create_dataset('quantile normalized example data', expr)
```

You can also remove existing datasets:
```
    >> remove_dataset('quantile normalized example data')
```

Note that all new datasets are fully self-contained and include all of the raw data. This means that if you only want to store the result of a filter query, it is a better idea to simply save() the query set object, since you don't want to make an extra copy of the raw data.

Each dataset is also fully self-contained in terms of metadata. This means that if you wish to share data with another group of researchers, you can simply save the data object, and all metadata will be automatically included:
```
    >> expr = realize(query('TCGA/GBM Affy HT HG U133A gene expression', 1:10));
    >> save ~/my_expr_data.mat expr
    >> clear expr
    >> load ~/my_expr_data.mat
    >> expr.Meta
    ans =
                   Type: 'Gene expression'
               Platform: {10x1 cell}
                 Sample: [1x1 struct]
                Patient: [1x1 struct]
                   Misc: [1x1 struct]
    SummarizationMethod: {10x1 cell}
               Organism: 'Homo sapiens'
                Version: 'RefSeq 38'
```

If the other research group is also using Matlab, they can then immediately use the data, and have access to the full metadata. When sharing gene expression data, it is of course also a good idea to provide the other group with the list of gene names for identifying the rows in your gene expression data matrix.






# DATA IMPORT - GENERAL #

In the previous chapters of this tutorial, we have retrieved all of our data directly from the pipeline data store, in a format that is immediately usable in Matlab. In real life, bioinformatics data is stored in a number of different formats and representations. In order to analyze real life data, the data must first be brought into the pipeline and represented in the unified data formats that the pipeline deals with.

To import microarray data into the pipeline, you can use the following command:
```
    >> import_uarray_data('dataset name', 'Affymetrix Human Exon 1.0 ST')
    Importing GSM526134_YX_Exon1_PCA0001.CEL...
    Importing GSM526135_YX_Exon1_PCA0002.CEL...
    Importing GSM526136_YX_Exon1_PCA0003.CEL...
                   ...
```

This function checks the current working directory for microarray sample files, and reads them into the pipeline, creating a new dataset with the requested name. You also need to specify the microarray platform, since the function cannot yet autodetect all microarray platforms based on the sample file contents. The second argument will likely become optional in the future.

You can get a list of supported microarray platform names from the pipeline's internal ontologies:
```
    >> ontologies.uarray_platforms
```

To import sequencing data into the pipeline, you can use the command:
```
    >> import_seq_reads('dataset name', 'ABI SOLiD V3')
```

Again, the command looks for FASTA, FASTQ and colorspace read files in the current working directory, and creates a new dataset out of them. This function will also move the read files into the pipeline data store, so be careful.

A list of sequencing platform names is again provided in the internal ontologies:
```
    >> ontologies.sequencing_platforms
```

If you call import\_uarray\_data() or import\_seq\_reads() as shown above, no metadata will be associated with the samples. The only populated metadata columns will be Meta.Sample.Filename, Meta.Platform and Meta.Type.

To add metadata to any dataset, you must first create or otherwise acquire an Excel spreadsheet containing the clinical metadata for your sample. Then, once you have an Excel spreadsheet, you can use the function `augment_meta_xls()` to associate metadata with your samples:
```
    >> augment_meta_xls('dataset name', 'spreadsheet.xls')
```

Note that the spreadsheet must be in the old `.xls` format. The newer `.xlsx` is not supported by Matlab. If you have any issues importing metadata from the `.xls` format, you can also opt to use the function `augment_meta_tab()`, which uses tab delimited text files rather than `.xls` files. Excel supports exporting your spreadsheet as a tab delimited text file.

The metadata association process first looks at your spreadsheet and tries to find a column with the name "Filename". If such a column is found, clinical metadata is associated with samples based on the sample filenames stored in Meta.Sample.Filename. The association can also happen through the field Meta.Sample.ID, if a filename-to-sample mapping is provided in the following manner:
```
    >> sample_map = containers.Map;
    >> sample_map('filename') = 'sample ID';   % Construct the filename-sample mapping
    >> ...
    >> augment_meta_xls('dataset name', 'metadata.xls', 'SampleMap', sample_map)
```







# DATA IMPORT - TCGA #

The pipeline supports automated importing of TCGA data, given that your TCGA data is stored in a directory structure similar or identical to that used by the TCGA FTP and HTTPS directories. All you need to do is to enter a directory that contains one or more TCGA archives, and run the function import\_tcga\_uarray\_data(). Here's an example:
```
    >> cd /worktmp/TCGA/gbm/cgcc/mskcc.org/hg-cgh-244a/cna
    >> import_tcga_uarray_data('tcga/gbm agilent cgh 244a mskcc test', ...
           'Agilent HG CGH 244A');
```

The function will automatically scan the directory and any subdirectories for SDRF files, and will then look for raw data files for the samples described in the SDRF files. The function then calls import\_uarray\_data() to import the microarray data into the pipeline. As the final step, TCGA patient metadata is automatically associated with the microarray samples.





# DATA IMPORT - GENE EXPRESSION OMNIBUS (GEO) #

Limited support is provided for importing data from the Gene Expression Omnibus service. Raw microarray samples downloaded from GEO can be automatically imported using the generic import\_uarray\_data(). For metadata, you can download a series matrix file from GEO and then parse it:
```
    >> geo_meta = read_geo_series_matrix('GSE21032_series_matrix.txt')
    geo_meta =
                Title: 'Integrative genomic profiling of human prostate cancer'
      SeriesAccession: 'GSE21032'
       SubmissionDate: 'Mar 23 2010'
       LastUpdateDate: 'Jul 06 2010'
          SampleTitle: {185x1 cell}
      SampleAccession: {185x1 cell}
             SampleID: {185x1 cell}
```

For GEO datasets, you can often use the SampleAccession and SampleID fields to build a mapping from filenames to sample IDs. This mapping can then be given to the function import\_uarray\_data() as the optional parameter 'SampleMap'. This allows you to create a spreadsheet that only contains clinical metadata associated with sample IDs, and the association with filenames then happens automatically via the mapping.