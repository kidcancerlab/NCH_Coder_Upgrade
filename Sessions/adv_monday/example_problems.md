Examples of problematic datasets:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788
    - The SRA only has one read - missing all cell barcodes and UMI data
    - Processed data is already normalized (not integers)
    - Data is just text files of matrices - waste of space, can't read in with Seurat
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249268
    - Data processed through non-standard pipeline (kallisto)
    - Raw SRA entries has 4 SRA numbers per sample
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222510
    - Data are saved as three matrix files, but each has sample name prefix, must be renamed
    - Gene name or features file missing. Expecting features.tsv.gz. Need to rename genes.tsv.gz to features.tsv.gz
    - The genes.tsv.gz file has the wrong format
>[mvc002@r1pl-hpcf-log01 NCH_Coder_Upgrade]$ zcat input/scRNA/GSM6925133_OX1X/features.tsv.gz | head
>Xkr4
>Gm1992
>Gm37381
>Rp1
>Sox17
>Mrpl15
>Lypla1
>Gm37988
>Tcea1
>Rgs20
>[mvc002@r1pl-hpcf-log01 NCH_Coder_Upgrade]$ zcat /home/gdrobertslab/lab/Counts_2/S0055/filtered_feature_bc_matrix/features.tsv.gz | head
>GRCh38_ENSG00000243485  GRCh38_MIR1302-2HG      Gene Expression
>GRCh38_ENSG00000237613  GRCh38_FAM138A  Gene Expression
>GRCh38_ENSG00000186092  GRCh38_OR4F5    Gene Expression
>GRCh38_ENSG00000238009  GRCh38_AL627309.1       Gene Expression
>GRCh38_ENSG00000239945  GRCh38_AL627309.3       Gene Expression
>GRCh38_ENSG00000239906  GRCh38_AL627309.2       Gene Expression
>GRCh38_ENSG00000241860  GRCh38_AL627309.5       Gene Expression
>GRCh38_ENSG00000241599  GRCh38_AL627309.4       Gene Expression
>GRCh38_ENSG00000286448  GRCh38_AP006222.2       Gene Expression
>GRCh38_ENSG00000236601  GRCh38_AL732372.1       Gene Expression

