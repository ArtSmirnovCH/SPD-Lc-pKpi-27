import numpy as np
import pandas as pd

from typing import Callable
from copy import deepcopy

from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold

import os
import sys
import warnings

#-------------------------------------------------------------------
# Add path

try:
    current_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    current_dir = os.getcwd()
    
parent_dir = os.path.join(current_dir, '..')
sys.path.insert(0, parent_dir)

#-------------------------------------------------------------------
# Add local moduls

from analysis_scripts.config import CFG
#-------------------------------------------------------------------
# Configure env

CFG.seed_all(CFG.seed)
warnings.simplefilter('ignore')

#########################################################################
# CrossVal
#########################################################################

class CrossVal():       # CURRENTLY assumes binary classification
    
    def __init__(
        self,
        x_train: pd.DataFrame,
        y_train: pd.Series,
        nfolds: int = 5,
        stratification: bool = True
    ):
        """
        Initialize the CrossVal object and set up cross-validation splits.
        
        Args:
            x_train (pd.DataFrame): Training feature data (n_samples, n_features)
            y_train (pd.Series): Training target values (n_samples,)
            nfolds (int, optional): Number of folds for cross-validation.
            stratification (bool, optional): Whether to use stratified sampling.
                For classification tasks, maintains class distribution in each fold.
                For regression with >100 unique values, bins target into quantiles.
        """
        self.x_train = x_train
        self.y_train = y_train
        self.nfolds = nfolds
        self.stratification = stratification
        self.models = None

        if stratification:

            if y_train.nunique() > 100:
                n_bins = 10
                y_binned = pd.qcut(y_train, q=n_bins, labels=False, duplicates='drop')
            else:
                y_binned = y_train
            
            skf = StratifiedKFold(n_splits=nfolds, shuffle=True, random_state=CFG.seed)
            split_generator = skf.split(x_train, y_binned)
        
        else:
        
            kf = KFold(n_splits=nfolds, shuffle=True, random_state=CFG.seed)
            split_generator = kf.split(x_train)
            
        self.split_generator = split_generator
        
        
    def fit(self, model: Callable, eval_metric: Callable, train_validation: bool = False):
        """
        Train models using cross-validation and compute out-of-fold predictions.
        """
        self.models = []
        
        oof_pred_proba = np.zeros_like(self.y_train, dtype=np.float64)
    
        for fold, (train_idx, val_idx) in enumerate(self.split_generator):
            
            print(f'--- Fold {fold+1}/{self.nfolds} ---')

            x_train_fold, x_val_fold = self.x_train.iloc[train_idx], self.x_train.iloc[val_idx]
            y_train_fold, y_val_fold = self.y_train.iloc[train_idx], self.y_train.iloc[val_idx]
            
            if model.__class__.__name__ == 'XGBClassifier':
            
                fold_model = deepcopy(model)
            
                if train_validation:
                    fold_model.fit(
                        x_train_fold, 
                        y_train_fold,
                        eval_set=[(x_train_fold, y_train_fold), (x_val_fold, y_val_fold)],
                        verbose=0
                    )
                else:
                    fold_model.fit(x_train_fold, y_train_fold)
            
            elif model.__class__.__name__ == 'LGBMClassifier':
        
                fold_model = deepcopy(model)
            
                if train_validation:
                    fold_model.fit(
                        x_train_fold,
                        y_train_fold,
                        eval_metric='auc',
                        eval_set=[(x_val_fold, y_val_fold)],
                    )
                else:
                    fold_model.fit(x_train_fold, y_train_fold)
            
            elif model.__class__.__name__ == 'CatBoostClassifier':
            
                fold_model = deepcopy(model)
            
                if train_validation:
                    fold_model.fit(
                        x_train_fold,
                        y_train_fold,
                        eval_set=(x_val_fold, y_val_fold),
                        plot=False
                    )
                else:
                    fold_model.fit(x_train_fold, y_train_fold)
                    
            elif model.__class__.__name__ in ['LogisticRegression', 'Pipeline']:
            
                fold_model = deepcopy(model)
                fold_model.fit(x_train_fold, y_train_fold)
            
            elif model.__class__.__name__ == 'GradientBoostedTreesLearner':
                
                fold_model = deepcopy(model)
                
                train_fold_df = x_train_fold.copy()
                train_fold_df['target'] = y_train_fold.copy()
            
                fold_model = fold_model.train(train_fold_df)
                    
            elif model.__class__.__name__ == 'RandomForestLearner':
                
                fold_model = deepcopy(model)
                
                train_fold_df = x_train_fold.copy()
                train_fold_df['target'] = y_train_fold.copy()
            
                fold_model = fold_model.train(train_fold_df)
            
        
            if model.__class__.__name__ in ['GradientBoostedTreesLearner', 'RandomForestLearner']:
                y_pred_proba_train = fold_model.predict(x_train_fold)
                y_pred_proba_val = fold_model.predict(x_val_fold)
            else:
                y_pred_proba_train = fold_model.predict_proba(x_train_fold)[:, 1]
                y_pred_proba_val = fold_model.predict_proba(x_val_fold)[:, 1]   
                
            oof_pred_proba[val_idx] = y_pred_proba_val    
            
            fold_train_score = eval_metric(y_train_fold, y_pred_proba_train)
            fold_val_score = eval_metric(y_val_fold, y_pred_proba_val)
            
            self.models.append(fold_model)    
            
            print(f"Fold {fold+1} {eval_metric.__name__}. OOF Train Score: {fold_train_score:.4f} OOF Val Score: {fold_val_score:.4f}")
        
        return self.models, oof_pred_proba
    
    
    def predict(self, x_df: pd.DataFrame):
        """
        Generate predictions from all trained cross-validation models.
        """
        pass
    
