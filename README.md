# FIND

Despite advances in annotating and interpreting human genetic variants, existing methods to distinguish pathogenic from neutral variants still inadequately capture the nuanced impact of genetic variants on fitness and disease susceptibility. To address this, we introduced FIND, a method that enhances resolution in differentiating trait-modulating alleles from pathogenic or neutral ones by stratifying genetic variants into refined categories based on selection pressures and derived allele frequency.

- We welcome any discussion, suggestion and potential contribution of new functional prediction scores through github or contact Dr. Mulin Jun Li (mulinli{at}connect.hku.hk).

  

## Download

### FIND V1.0

- **Genome build GRCh37 / hg19**

  All possible SNVs of GRCh37 / hg19 [Download]()

  All possible SNVs of GRCh37 / hg19 incl. all annotations [Download]()

- **Genome build GRCh38 / hg38**

​	All possible SNVs of GRCh38 / hg38 [Download]()

​	All possible SNVs of GRCh38 / hg38 incl. all annotations [Download]()



## Usage

### Score interpretation and prioritization

We have expanded the conventional binary partitioning of all human genome variants into four nuanced categories, utilizing the dimensions of natural selection and derived allele frequency (DAF) spectrum. Therefore, for each SNV, FIND will provide the probability that the variant is predicted to be in four categories, named: FIND_F, FIND_I, FIND_N, FIND_D. The raw scores reported by FIND can be used to determining the categories of variant. Similar to CADD C-scores and it's phred-like scores (**PHRED score**), we recommend to use phred-like scores (FIND_F_PHRED, FIND_I_PHRED, FIND_N_PHRED, FIND_D_PHRED) for the likely causal variant prioritization (for each FIND categories) and even for comparison among different models.



## Building FIND model

### Requirements

- python 3.10
- [TabNet](https://github.com/dreamquark-ai/tabnet)
- scikit-learn
- [VarNote](http://www.mulinlab.org/varnote/index.html)

### Procedures

1. Annotated training dataset by using VarNote or other way.

   Annotation

   ```bash
   ## Examples of VarNote annotations
   java -jar /path to/VarNote-1.2.0.jar AnnotationConfig -I ./config/BigAnno_SNV.config
   ```

2. Feature processing.

   ```bash
   cd ./script/01_DataAnno
   python ./Variants_Annotation_process.py -i [/input file] -o [/output file]
   ```

3. Model training.

   ```bash
   cd ./script/02_ModelTraining
   python FIND_TabNet_train.py --ModelP [/path to training dataset] --MoldeF [training dataset file name] --ModelO [/path to output model file]
   ```



## Copyright

Copyright (c) Mulinlab@Tianjin Medical University 2020-2024. All rights reserved.

