# FIND

Despite advances in annotating and interpreting human genetic variants, existing methods to distinguish pathogenic from neutral variants still inadequately capture the nuanced impact of genetic variants on fitness and disease susceptibility. To address this, we introduced FIND, a method that enhances resolution in differentiating trait-modulating alleles from pathogenic or neutral ones by stratifying genetic variants into refined categories based on selection pressures and derived allele frequency (DAF) . We have expanded the conventional binary partitioning of all human genome variants into four nuanced categories, utilizing the dimensions of natural selection and derived allele frequency spectrum (as shown below), including: 1) Fixed/Nearly Fixed (**F**) category; 2) Intermediate Selection (**I**) category; 3) Neutral (**N**) category; 4) Deleteriousness (**D**) category.

![image-20240526124047898](./figure/image-20240526124047898.png)

- FIND accurately classifies genome-wide pathogenic variants.

- FIND improved resolution for identifying favored allele with trait-modulating effect.

- We welcome any discussion, suggestion and potential contribution of new functional prediction scores through github or contact Dr. Mulin Jun Li (mulinli{at}connect.hku.hk).

  

## Download

### FIND V1.0

- **Genome build GRCh37 / hg19**

  All possible SNVs of GRCh37 / hg19 [Download]() (494G)

  All possible SNVs of GRCh37 / hg19 incl. all annotations [Download]() (1.3T)

- **Genome build GRCh38 / hg38**

<<<<<<< HEAD
  All possible SNVs of GRCh38 / hg38 [Download]() (492G)

  All possible SNVs of GRCh38 / hg38 incl. all annotations [Download]() (1.3T)
=======
  All possible SNVs of GRCh38 / hg38 [Download]()

  All possible SNVs of GRCh38 / hg38 incl. all annotations [Download]()
>>>>>>> 399fac80a32f6c49a5ea2afb35b6f8a710040f2c



## Usage

### Score interpretation and prioritization

We have expanded the conventional binary partitioning of all human genome variants into four nuanced categories, utilizing the dimensions of natural selection and derived allele frequency (DAF) spectrum. Therefore, for each SNV, FIND will provide the probability that the variant is predicted to be in four categories, named: FIND_F, FIND_I, FIND_N, FIND_D. The raw scores reported by FIND can be used to determining the categories of variant. Similar to CADD C-scores and it's phred-like scores (**PHRED score**), we recommend to use phred-like scores (FIND_F_PHRED, FIND_I_PHRED, FIND_N_PHRED, FIND_D_PHRED) for the likely causal variant prioritization (for each FIND categories) and even for comparison among different models.

### Get the FIND scores

- Get all pre-calculated FIND scores in [hg19]() (494G)/[hg38 ]()(492G)

  1. Download

     ```bash
     # hg19
     wget -c web-site(file.gz)
     wget -c web-site(file.gz.vanno)
     wget -c web-site(file.gz.vanno.vi)
     
     # hg38
     wget -c web-site(file.gz)
     wget -c web-site(file.gz.vanno)
     wget -c web-site(file.gz.vanno.vi)
     ```

  2. Get the FIND prediction scores by using [VarNote](http://www.mulinlab.org/varnote/index.html)

     ```bash
     # Modify configuration files (FIND_ScoreAnno_SNV.config)
     ## 1.Path of FIND score file
     ## 2.Path and format of query file
     ## 3.Path of output file
     
     # Get score
     java -jar /path to/VarNote-1.2.0.jar AnnotationConfig -I ./config/FIND_ScoreAnno_SNV.config
     ```

     Please refer to this [document](http://www.mulinlab.org/varnote/documentaiton.html) for more details on downloading and using VarNote.

- Query the pre-calculated FIND scores of the interested variants directly from [Vannoportal](http://www.mulinlab.org/vportal/index.html).

- Recalculated FIND raw scores by the trained FIND model.

  Please refer to the following method for feature annotation.

  ```bash
  cd ./script/02_ModelTraining
  python ./script/FIND_predict.py -i [input file] -o [output file]
  ```

## Building FIND model

### Requirements

- python 3.10
- [TabNet](https://github.com/dreamquark-ai/tabnet)
- scikit-learn
- [VarNote](http://www.mulinlab.org/varnote/index.html)

### Procedures

1. Obtain FIND training data.

   The FIND is based on the differences in fitness among genetic variants from different evolutionary processes in modern human genome. Specifically, FIND distinguishes between genetic variants that derived allele have undergone fixed/nearly fiexed (labelled **F**), intermediate selection (labelled **I**), neutral (labelled **N**), and deleterious/*de novo* (labelled **D**). We constructed a training dataset comprising approximately 2 million variants, ensuring a balanced representation across four distinct categories.

   > The FIND training data set is [FIND_Training_Data.tsv](./FIND/train_dataset), and the corresponding relationship of data labels is F:0, I:1, N:2, D:3.

2. Feature annotation and processing.

   To more comprehensively capture the evolutionary patterns and adaptive selection mechanisms of human variations, especially those located in non-coding regions, we systematically computed and integrated 289 base-wise variant annotation features across the whole genome, with a significant percentage achieving allele-specific resolution.

   - Download annotation file in [hg19]()/[hg38]() (1.3T).

   - Annotated training dataset by using VarNote or other way.

     ```bash
     # Examples of VarNote annotations
     # Modify configuration files (BigAnno_SNV.config)
     ## 1.Path of feature file
     ## 2.Path and format of query file
     ## 3.Path of output file
     
     java -jar /path to/VarNote-1.2.0.jar AnnotationConfig -I ./config/BigAnno_SNV.config
     ```

   - Feature processing.

     ```bash
     cd ./script/01_DataAnno
     python ./Variants_Annotation_process.py -i [input file] -o [output file]
     ```

3. Model training.

   We utilized the [TabNet](https://github.com/dreamquark-ai/tabnet) model to classify four genetic variant categories in our training set.

   ```bash
   cd ./script/02_ModelTraining
   python FIND_TabNet_train.py --ModelP [path to training dataset] --MoldeF [training dataset file name] --ModelO [path to output model file]
   ```

## Copyright

Copyright (c) Mulinlab@Tianjin Medical University 2020-2024. All rights reserved.
