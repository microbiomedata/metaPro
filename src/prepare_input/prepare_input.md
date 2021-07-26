To run the workflow, we need the following data files

For NMDC at PNNL, we collect 

    1. RAW MS/MS files from DMS or EMSL_NEXUS for a FICUS datasets (PNNL internal use only.)
       - Capable of downloading 
         - datasets for a datapackage ids.
         - datasets for a set of datasets ids.
         - msgfplus and masic jobs already processed by DMS(need their respective set of jod numbers as input )
       
    2. Fasta and GFF files from Provided by JGI `cori.nersc.gov`
       - manually secure copy using scp from `cori.nersc.gov` to pnnl server.

    3. a manually created mapper file(excel file) provided by scientist that has information about datasets mapping to fasta files.

Note: Identified FICUS datasets were `Hess`, `Stegen`, `Blanchard`.

For external users of this workflow, please have your datasets ready in a folder `storage/data`.  

Note: 
- it's up to the user whether they want to maintain a directory structure for grouping there datasets or not. Workflow runs independent of how you arrange your input. 
- However, you're expected to put 
  - all raw files in `storage/data`
  - fastas `storage/fastas` (must have .faa & .gff files)
  - pre-defined mappings(excel format) in `storage/mappings`
  - parameters and contaminants in `storage/parameters`

Additionally, you can [mount network-share drives](utility/mount_network_drive.sh) assuming you data directory `storage/` lives there or even in [project's subdirectory.](./)

An Example on a `path/to/your/storage/` input data looks like:

```json
.
├── data
│ └── set_of_Dataset_IDs
│     └── stegen
│         └── 500088
│             └── Froze_Core_2015_N2_50_60_34_QE_26May16_Pippin_16-03-39.raw
├── fastas
│ └── stegen
│     └── 1781_100336
│         ├── Ga0482236_functional_annotation.gff
│         └── Ga0482236_proteins.faa
├── mappings
│   └── EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping_2021-01-25.xlsx
└── parameters
    ├── LTQ-FT_10ppm_2014-08-06.xml
    ├── MSGFPlus_PartTryp_MetOx_20ppmParTol.txt
    ├── MSGFPlus_PartTryp_MetOx_20ppmParTol_ModDefs.txt
    ├── MSGFPlus_Tryp_NoMods_20ppmParTol.txt
    ├── Mass_Correction_Tags.txt
    └── Tryp_Pig_Bov.fasta
```