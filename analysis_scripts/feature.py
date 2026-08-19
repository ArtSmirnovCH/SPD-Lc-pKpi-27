import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import chi2
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import SelectKBest

from typing import Tuple, Optional, List


################################################################################
# FeatureClassification Class
################################################################################

class FeatureClassification:
    
    # -- Features
    #   -- Categorical
    #       -- Nominal
    #           -- Binary
    #           -- Multi
    #       -- Ordinal
    #           -- Binary
    #           -- Multi
    #   -- Numeric
    #       -- Discrete
    #           -- Binary
    #           -- Multi
    #       -- Continual
    
    def __init__(self, df: pd.DataFrame, target_name: str, ordinal_columns: List[str] = None):
        """
        Initialize feature classifier by splitting all columns into detailed feature types.
        
        Args:
            df: Input DataFrame
            target_name: Name of target column (excluded from classification)
            ordinal_columns: List of categorical columns
        """
        columns = df.columns.copy()
        columns = columns.drop(target_name)
    
        cat_features, num_features = self.categorical_numeric_split(df=df, columns=columns)
        # num
        discrete_features, continual_features = self.discrete_continual_split(df=df, numeric_columns=num_features, threshold=0.2)
        discrete_binary_features, discrete_multi_features = self.binary_multi_split(df=df, columns=discrete_features)
        # cat
        ordinal_features, nominal_features = self.categorical_ordinal_nominal_split(df=df, categorical_columns=cat_features, ordinal_cols=ordinal_columns)
        ordinal_binary_features, ordinal_multi_features = self.binary_multi_split(df=df, columns=ordinal_features)
        nominal_binary_features, nominal_multi_features = self.binary_multi_split(df=df, columns=nominal_features)
        
        # cat
        self.cat_features = cat_features
        self.ordinal_features = ordinal_features
        self.ordinal_binary_features = ordinal_binary_features
        self.ordinal_multi_features = ordinal_multi_features
        self.nominal_features = nominal_features
        self.nominal_binary_features = nominal_binary_features
        self.nominal_multi_features = nominal_multi_features
        
        # num
        self.num_features = num_features
        self.discrete_features = discrete_features
        self.continual_features = continual_features
        self.discrete_binary_features = discrete_binary_features
        self.discrete_multi_features = discrete_multi_features
        
        print("=" * 50)
        print('GENERAL FEATURES INFO')
        print("=" * 50)
        for col in df.columns:
            print(f"{col + ',':<30} nunique: {str(df[col].nunique()) + ',':<10} NANs: {df[col].isna().sum()}")
        
        print('\n', "=" * 50)
        print("FEATURE CLASSIFICATION SUMMARY")
        print("=" * 50)
        print(f"\nCATEGORICAL FEATURES ({len(self.cat_features)}): {self.cat_features}")
        print(f"  Ordinal ({len(self.ordinal_features)}): {self.ordinal_features}")
        print(f"    - Binary ({len(self.ordinal_binary_features)}): {self.ordinal_binary_features}")
        print(f"    - Multi ({len(self.ordinal_multi_features)}): {self.ordinal_multi_features}")
        print(f"  Nominal ({len(self.nominal_features)}): {self.nominal_features}")
        print(f"    - Binary ({len(self.nominal_binary_features)}): {self.nominal_binary_features}")
        print(f"    - Multi ({len(self.nominal_multi_features)}): {self.nominal_multi_features}")
        
        print(f"\nNUMERIC FEATURES ({len(self.num_features)}): {self.num_features}")
        print(f"  Discrete ({len(self.discrete_features)}): {self.discrete_features}")
        print(f"    - Binary ({len(self.discrete_binary_features)}): {self.discrete_binary_features}")
        print(f"    - Multi ({len(self.discrete_multi_features)}): {self.discrete_multi_features}")
        print(f"  Continual ({len(self.continual_features)}): {self.continual_features}")
        

    @staticmethod
    def categorical_numeric_split(df: pd.DataFrame, columns: List[str]) -> Tuple[List, List]:
        """
        Split column list into categorical and numeric based on conversion success.
        
        Attempts to convert each column to numeric; if successful, column is numeric,
        otherwise it's categorical.
        """
        cat_features = []
        num_features = []
        
        for col in columns:
            
            try:
                pd.to_numeric(df[col])
                num_features.append(col)
            except ValueError:
                cat_features.append(col)
                
        return cat_features, num_features
        
    
    @staticmethod
    def categorical_ordinal_nominal_split(df: pd.DataFrame, categorical_columns: List[str], ordinal_cols: List[str] = None) -> Tuple[List, List]:
        """
        Split categorical columns into ordinal and nominal features.
        
        Args:
            df: Input DataFrame
            columns: List of categorical column names to split
            ordinal_cols: List of ordinal columns
            
        Returns:
            Tuple of (ordinal_features, nominal_features)
        """
        
        ordinal_features = ordinal_cols if ordinal_cols is not None else []
        nominal_features = [col for col in categorical_columns if col not in ordinal_features]
        
        return ordinal_features, nominal_features
    
    
    @staticmethod
    def discrete_continual_split(df: pd.DataFrame, numeric_columns: List[str], threshold: float = 0.2) -> Tuple[List, List]:
        """
        Split columns into discrete and continual features based on unique value ratio.
        
        Args:
            df: Input DataFrame
            columns: List of column names to check
            threshold: Ratio threshold (unique_count / total_count) to determine discretization
            
        Returns:
            Tuple of (discrete_features, continual_features)
        """
        
        discrete_features = []
        continual_features = []
        
        for col in numeric_columns:
            
            data = df[col].dropna()
            if data.nunique() < int(len(data) * threshold):
                discrete_features.append(col)
            else:
                continual_features.append(col)
        
        return discrete_features, continual_features
        
        
    @staticmethod
    def binary_multi_split(df: pd.DataFrame, columns: List[str]) -> Tuple[List, List]:
        """
        Split features into binary and multi-class based on unique value count.
        
        Columns with exactly 2 unique values are binary; others are multi-class.
        """
        binary_features = []
        multi_features = []
        
        for col in columns:
            
            if df[col].nunique() == 2:
                binary_features.append(col)
            else:
                multi_features.append(col)
                
        return binary_features, multi_features
                
                
################################################################################
# FeatureEngineer Class
################################################################################

class FeatureEngineer:

    def __init__(self):

        self.target_encoding_maps = {}

    
    def target_encoder(
        self,
        df: pd.DataFrame,
        categorical_columns: list[str],
        target_name: str,
        train_df: pd.DataFrame, # after taget encoding
        n_fold: int = 5,
        replace : bool = False) -> pd.DataFrame:
        
        """
        Apply target encoding with out-of-fold strategy for training data and mapping for test data.
        
        This method implements target encoding that prevents target leakage by using out-of-fold
        means for training data. For training data, it performs stratified k-fold cross-validation
        where each fold is encoded using means from other folds. For test data, it uses global
        means calculated from the entire training set. Missing values are filled with the global
        target mean.
        
        The method automatically detects whether input is train or test data by checking for the
        presence of the target column. For test data, it requires that target encoding has been
        previously applied to training data to build the encoding maps.
        
        Args:
            df: Input DataFrame (can be training or test data)
            categorical_columns: List of categorical column names to encode
            train_df: Training DataFrame used to calculate encoding maps (for test data)
            n_fold: Number of folds for out-of-fold encoding (default: 5)
            replace: If True, replaces original columns; if False (default), adds 
                    new columns with '_target_encode' suffix
        
        Returns:
            DataFrame with target-encoded categorical columns
            
        Raises:
            KeyError: When encoding test data without first encoding training data
        """

        df = df.copy()

        # For Null filling
        mean_of_target = train_df[target_name].mean()

        # Prepare for stratification
        if train_df[target_name].nunique() > 100:
            n_bins = 10
            y_binned = pd.qcut(train_df[target_name], q=n_bins, labels=False, duplicates='drop')
        else:
            y_binned = train_df[target_name]

        for col in categorical_columns:

            new_col_name = col + '_target_encode'
            df[new_col_name] = np.nan

            skf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=58)

            if target_name in df.columns:
                # Train data processing - OOF encoding
                
                for tr_ind, val_ind in skf.split(df, y_binned):
                    
                    df_tr, df_val = df.iloc[tr_ind], df.iloc[val_ind]
                    
                    map_encoder = df_tr.groupby(col)[target_name].mean().to_dict()
                    
                    df.loc[df.index[val_ind], new_col_name] = df_val[col].map(map_encoder)
                
                # Fill NaN   
                df[new_col_name].fillna(mean_of_target, inplace = True)

                # Store the global mapping for test data
                global_map = train_df.groupby(col)[target_name].mean().to_dict()
                self.target_encoding_maps[col] = global_map

                if replace:
                    df = df.drop(col, axis=1)
                    
            else:
                # For test data - use mapping from train data
                if col not in self.target_encoding_maps:
                    raise KeyError(f'{new_col_name} must be in train DataFrame. Apply target Encoding for train Dataset first!')
                
                map_encoder = self.target_encoding_maps[col]
                
                if replace:
                    df[col] = df[col].map(map_encoder)
                    df[col].fillna(mean_of_target, inplace=True)
                else:
                    df[new_col_name] = df[col].map(map_encoder)
                    df[new_col_name].fillna(mean_of_target, inplace=True)
                    
        return df  


    @staticmethod
    def decimal_encoder(df: pd.DataFrame, numeric: list[str], replace: bool = False) -> pd.DataFrame:
        """ Encode numeric values into separate columns for each digit position.

        Args:
            df (pd.DataFrame): data to encode
            numeric (list[str]): list of features to encode
            replace (bool, optional): replcae old columns in dataset or add new columns. Defaults to False.

        Returns:
            pd.DataFrame: encoded data
        """
        
        df = df.copy()
        
        for col in numeric:
            
            col_data = df[col].astype(str).str.replace('.', '')
        
            max_length = col_data.str.len().max()
        
            for pos in range(max_length):
                
                # digits = col_data.str[pos].fillna('')
                digits = col_data.str[pos]
                
                df[f'{pos}_decimal_encoder'] = digits
                
            if replace:
                df = df.drop(col, axis=1)
        
        return df


    @staticmethod
    def label_encoder(df: pd.DataFrame, categoric: list[str], replace: bool = False) -> pd.DataFrame:
        """
        Encode categorical columns with integer labels.
        
        This method transforms categorical values into integer representations based on
        alphabetical sorting.
        
        Args:
            df: Input DataFrame
            categoric: List of column names to encode
            replace: If True, replaces original columns; if False (default), creates 
                    new columns with '_label_encode' suffix
        
        Returns:
            DataFrame with encoded categorical columns
        """
        
        df = df.copy()
        
        if replace:
            
            for col in categoric:

                dct = {}
                for val, key in enumerate(sorted(df[col].unique(), reverse=False)):
                    dct[key] = val
                
                df[col] = df[col].map(dct)

        else:
     
            for col in categoric:
                
                dct = {}
                for val, key in enumerate(sorted(df[col].unique(), reverse=False)):
                    dct[key] = val
                
                df[col + '_label_encode'] = df[col].map(dct)

        return df


    @staticmethod
    def count_encoder(df: pd.DataFrame, categorical_columns: list[str], target_name: str, replace: bool = False) -> pd.DataFrame:
        """
        Encode categorical columns by their frequency counts.
        
        This method replaces categorical values with their occurrence counts in the dataset.
        
        Args:
            df: Input DataFrame
            categorical_columns: List of column names to encode
            replace: If True, replaces original columns; if False (default), creates 
                    new columns with '_count_encode' suffix
        
        Returns:
            DataFrame with count-encoded categorical columns
        """
        df = df.copy()
        
        if replace:
            
            for col in categorical_columns:
                map_val_count = df.groupby(by=col)[target_name].count().astype(np.int32).to_dict()
                df[col] = df[col].map(map_val_count)
          
        else:
        
            for col in categorical_columns:
                map_val_count = df.groupby(by=col)[target_name].count().astype(np.int32).to_dict()
                df[col + '_count_encode'] = df[col].map(map_val_count)

        return df


    @staticmethod
    def mean_encoder(df: pd.DataFrame, categorical_columns: list[str], target_name: str, replace: bool = False) -> pd.DataFrame:
        """
        Encode categorical columns by the mean of the target variable.
        
        This method replaces categorical values with the average target value for each category.
        
        Args:
            df: Input DataFrame
            categorical_columns: List of column names to encode
            replace: If True, replaces original columns (dropping them); if False (default), 
                    creates new columns with '_mean_encode' suffix while keeping originals
        
        Returns:
            DataFrame with mean-encoded categorical columns
        """
        df = df.copy()
        
        if replace:
            
            for col in categorical_columns:
                map_val_count = df.groupby(by=col)[target_name].mean().astype(np.int32).to_dict()
                df[col + '_mean_encode'] = df[col].map(map_val_count)
                df = df.drop(col, axis=1)
          
        else:
        
            for col in categorical_columns:
                map_val_count = df.groupby(by=col)[target_name].mean().astype(np.int32).to_dict()
                df[col + '_mean_encode'] = df[col].map(map_val_count)

        return df


    @staticmethod
    def one_hot_encoder(df: pd.DataFrame, categoric: list[str], drop_first: bool = True, replace: bool = False, suffix: str = '_onehot_encode') -> pd.DataFrame:
        """
        One-Hot Encoding for categorical variables
        
        Args:
            categoric : list[str]
                List of categorical column names to encode
            drop_first : bool
                Whether to drop the first category to avoid multicollinearity
            replace : bool
                If True, replace original columns with encoded ones
                If False, keep original columns and add encoded ones
        """
        
        df = df.copy()
        
        # Input validation
        missing_cols = [col for col in categoric if col not in df.columns]
        if missing_cols:
            raise KeyError(f'Columns not found in DataFrame: {missing_cols}')
        
        if not replace:

            encoded_df = pd.get_dummies(
                data=df[categoric],
                prefix=categoric,
                prefix_sep='_',
                drop_first=drop_first,
                dtype=np.int32
            )
            
            # Add suffix to encoded columns
            encoded_cols = [col for col in encoded_df.columns if any(col.startswith(f"{cat}_") for cat in categoric)]
            rename_dict = {col: f"{col + suffix}" for col in encoded_cols}
            encoded_df = encoded_df.rename(columns=rename_dict)
            
            df = pd.concat([df, encoded_df], axis=1)
        
        else:
            
            # Replace original columns
            encoded_df = pd.get_dummies(
                data=df,
                prefix=categoric,
                prefix_sep='_',
                columns=categoric,
                drop_first=drop_first,
                dtype=np.int32
            )
            
            # Add suffix to encoded columns
            encoded_cols = [col for col in encoded_df.columns if any(col.startswith(f"{cat}_") for cat in categoric)]
            rename_dict = {col: f"{col + suffix}" for col in encoded_cols}
            df = encoded_df.rename(columns=rename_dict)

        return df      
    

def feature_stat_importance(x: pd.DataFrame, y: pd.Series, method: callable, select_k_best: int = None): 
    
    k = 'all' if select_k_best is None else select_k_best
    selector = SelectKBest(score_func=method, k=k)
    selector.fit(x, y)

    all_features = x.columns
    all_scores = selector.scores_

    res_df = pd.DataFrame({
        'feature': all_features,
        'score': all_scores,
        'selected': selector.get_support()
    }).sort_values(by='score', ascending=False).reset_index(drop=True)
    
    
    fig, ax = plt.subplots(figsize=(8, max(4, res_df.loc[res_df.selected].shape[0] * 0.3)))
                
    sns.barplot(
        data=res_df.loc[res_df.selected],
        x='score',
        y='feature',
        width=0.4,
        orient='h',
        alpha=0.8,
        ax=ax
    )
    
    ax.set_title(f'Feature Importance - {method.__name__} Scores', fontsize=14, fontweight='bold', pad=10)
    ax.set_xlabel(f'{method.__name__} Score', fontweight='bold')
    ax.set_ylabel('Feature', fontweight='bold')
    ax.grid(True, alpha=0.8, linestyle='--', axis='x')
    
    plt.tight_layout()
    
    return res_df

    
def features_correlation(
    df: pd.DataFrame,
    columns: List[str],
    target: str = None,
    corr_with_target: bool = False,
    target_values: List[str] = None,
    min_corr: float = None,
    figsize: Tuple = (6, 6)
) -> plt.Figure:
    """
    Analyze and visualize feature correlations.

    Args:
        df : pd.DataFrame
            Input dataframe containing the features and optionally target column
        columns : List[str]
            List of feature column names to analyze for correlations
        target : str
            Target column name used for filtering the data.
        corr_with_target : bool
            If calculate correlation with target too.
        target_values : Optional[List[str]]
            Values in target column to filter by. If provided, only rows where target column
            matches these values will be used for correlation analysis
        min_corr : Optional[float] 
            Minimum absolute correlation threshold to show.
        figsize : Tuple(int, int)
            figure size
    
    Returns:
        high_corr_df : pd.DataFrame
            DataFrame with form variable_1 x variable_2 -> correlation.
    """

    if target_values is None and target is None:
        df = df[columns].copy()
    elif target_values is None and target is not None:
        df = df[columns + [target]].copy()
    else:
        df = df.loc[df[target].isin(target_values)][columns + [target]].copy()

    corr = df[columns].corr().copy()
    
    print(f'{len(columns)=}')
    print(f'{np.linalg.matrix_rank(corr)=}')
    print(f'{np.linalg.det(corr)=}')

    # List of high correlated variables
    high_corr_pairs = []

    if min_corr is not None:
        
        mask = np.abs(corr) <= min_corr
        
        for i in range(len(corr.columns)):
            for j in range(i+1, len(corr.columns)): 
                corr_value = corr.iloc[i, j]
                if not np.isnan(corr_value) and abs(corr_value) > min_corr:
                    high_corr_pairs.append((
                        corr.columns[i], 
                        corr.columns[j],
                        corr_value
                    ))
    else:
        min_corr = 0.
        mask = None
    
    if high_corr_pairs:
        high_corr_df = pd.DataFrame(high_corr_pairs, columns=['corr_var_1', 'corr_var_2', 'corr'])
    else:
        high_corr_df = pd.DataFrame()
    
    print(f"\nHighly correlated pairs (|corr| > {min_corr}):\n")
    for col1, col2, corr_value in sorted(high_corr_pairs, key=lambda x: abs(x[-1]), reverse=True):
        print(f'{col1:<15} x {col2:<15}: {corr_value:.4f}.')
        if corr_with_target:
            print(f'Target x {col1:<15}: {df[columns + [target]].corr()[col1][target]:.4f}')
            print(f'Target x {col2:<15}: {df[columns + [target]].corr()[col2][target]:.4f}\n')
        

    # Visualization
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        data=corr,
        mask=mask,
        cmap='RdBu_r',
        annot=False,
        fmt='.3f',
        annot_kws={'fontsize': 10},
        linewidths=1,
        linecolor='gray',
        cbar_kws={
            'label': 'Variable Pair Correlation',
            'format': '%.2f'
        },
        center=0,
        square=True,
        ax=ax
    )

    if min_corr is None:
        if target_values is not None:
            ax.set_title(f'Correlation plot \n For target in {target_values}', fontweight='bold', pad=20)
        else:
            ax.set_title(f'Correlation plot', fontweight='bold', pad=20)
    else:
        if target_values is not None:
            ax.set_title(f'Correlation plot \n For target in {target_values} \n\n (abs(corr) > {min_corr})', fontweight='bold', pad=20)
        else:
            ax.set_title(f'Correlation plot \n\n (abs(corr) > {min_corr})', fontweight='bold', pad=20)
        
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # plt.savefig('../plots/all_corr_sig.pdf', bbox_inches='tight')
    plt.tight_layout()
    
    return high_corr_df