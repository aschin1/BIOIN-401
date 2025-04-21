# BIOIN-401

Usage:
1. Access the csv file of interest for training OR To create the csv file, call batch_processing.py with main4.py (VADAR) or main5.py (DSSP)
2. ensure the csv file, prepare_dataset.py and split_train_val.py are in the same directory and call prepare_dataset.py (this file can be edited to work with different csv files)
3. call split_train_val.py which will create .pt files for training and evaluation
4. call fine_tune_prot_bert5.py which will use the .pt files for training and validation across the epochs
5. once the fine_tuned_prot_bert model is saved, call evaluate_prot_bert4.py to evaluate its performance across the test set
6. call the final_workflow.py and ensure get_accessions4.py, get_alphafold_data4.py, and get_secondary_struct4.py are accessible


Imports:
- Biopython
- Pytorch
- Sklearn

