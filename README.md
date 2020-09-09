# SPRINT
Protein-peptide interactions are essential for all cellular processes including DNA repair, replication, gene-expression, and metabolism. As most protein-peptide interactions are uncharacterized, it is cost effective to investigate them computationally as the first step. All existing approaches for predicting protein-peptide binding sites, however, are based on protein structures despite the fact that the structures for most proteins are not yet solved. This article proposes the first machine-learning method called SPRINT to make Sequence-based prediction of Protein-peptide Residue-level Interactions. SPRINT yields a robust and consistent performance for 10-fold cross validations and independent test. The most important feature is evolution-generated sequence profiles. For the test set (1056 binding and non-binding residues), it yields a Matthews' Correlation Coefficient of 0.326 with a sensitivity of 64% and a specificity of 68%. This sequence-based technique shows comparable or more accurate than structure-based methods for peptide-binding site prediction. 

Cite: Taherzadeh, G., Yang, Y., Zhang, T., Liew, A. W. C., & Zhou, Y. (2016). Sequence‐based prediction of protein–peptide binding sites using support vector machine. Journal of computational chemistry, 37(13), 1223-1229.

Instruction:
* Protein-peptide dataset are stored in Data directory. Dataset file contains protein sequences labeled as 1 and 0 for binding and non-binding residues, respectively. Train and test files contain actual binding residues used in this study.
* Run ./SPRINT.py for feature extraction and peptide binding site prediction.
* Find the pre-trained model and protein-peptide complexes in pdb format in Release.

