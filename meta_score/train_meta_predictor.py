#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression

from mlxtend.feature_selection import ColumnSelector

from joblib import dump


# ## Parameters

# In[ ]:


INPUT='/home/tbrittoborges/Baltica/nanopore_benchmark/results/SJ_annotated.csv'
TEST_SIZE=0.2
RANDOM=1982
NANOPORE_CUTOFF=0.95

np.random.seed(RANDOM)
# ## Setup

# In[ ]:



# In[ ]:


data = pd.read_csv(INPUT)


# In[ ]:


data.head()


# In[ ]:


data.columns = data.columns.str.lower()


# In[ ]:


data = data.drop(columns=['is_novel', 'gene_name', 'transcript_name', 'class_code', 'exon_number', 'coordinates'])


# In[ ]:


data = data.fillna(0)


# In[ ]:


data.head()


# In[ ]:


y = data['orthogonal_na']
y = np.where(y > 0.95, 1, 0) # classification problem
x = data.drop(columns='orthogonal_na')


# ## Data split and preprocess

# In[ ]:


X_train, X_test, y_train, y_test = train_test_split(x, y, stratify=y, test_size=TEST_SIZE, random_state=RANDOM)


# In[ ]:


X_train.info()


# In[ ]:


np.unique(y_train, return_counts=True)[1]


# In[ ]:


np.unique(y_test, return_counts=True)


# ## Pipeline contruction

# ColumnSelector API http://rasbt.github.io/mlxtend/user_guide/feature_selection/ColumnSelector/#api

# In[ ]:


import itertools


# In[ ]:


col_selector = list(itertools.chain.from_iterable([itertools.combinations(range(4), i + 1) for i in range(4)]))


# In[ ]:


col_selector


# ## GBC 

# In[ ]:


pipe = Pipeline(
    [('select', ColumnSelector()),
     ('xboost', GradientBoostingClassifier(random_state=RANDOM))])


# In[ ]:


param_grid = {
    "select__cols": col_selector,
    "xboost__learning_rate": [1, 0.5, 0.25, 0.1, 0.05, 0.01],
    "xboost__max_depth": [1, 2, 4, 8],
    "xboost__subsample": [0.5, 0.8, 1.0],
    "xboost__min_samples_split": [0.2, 0.4, 0.8],
    "xboost__n_estimators": [1, 2, 4, 8, 16, 32, 64, 100, 200]}
print('Starting GCB grid search')
grid = GridSearchCV(pipe, param_grid, cv=10, n_jobs=-1, scoring='roc_auc')


# In[ ]:


grid.fit(X_train, y_train)
print('Best parameters:', grid.best_params_)
print('Best performance:', grid.best_score_)


# In[ ]:


df_grid_results_1 = pd.DataFrame(grid.cv_results_)


# In[ ]:


df_grid_results_1  = df_grid_results_1.sort_values('rank_test_score')


# In[ ]:


df_grid_results_1


# In[ ]:

print('Starting GCB grid search')
df_grid_results_1.to_csv('sklearn_result.csv')


# In[ ]:


df_grid_results_1.shape


# In[ ]:


df_grid_results_1


# In[ ]:


df_grid_results_1.loc[
    :,['param_select__cols', 'mean_test_score', 'rank_test_score']]


# In[ ]:


df_grid_results_1.loc[
    :,['param_select__cols']].value_counts()


# In[ ]:


df_grid_results_1.loc[:, "param_select__cols"] = df_grid_results_1.loc[:, "param_select__cols"].astype(str)


# In[ ]:


df_grid_results_1.query("param_select__cols != '(0, 1, 2, 3)'").loc[
    :, ["param_select__cols", "mean_test_score", "rank_test_score"]]


# In[ ]:


dump(grid.best_estimator_, 'gcb_meta_best_estimator.joblib')
dump(grid, 'gcb_grid_object.joblib')


# ## Logistic regression

# In[ ]:


pipe2 = Pipeline(
    [('select', ColumnSelector()),
     ('log_reg', LogisticRegression(random_state=RANDOM))])


# In[ ]:
print("Starting LR grid search")

param_grid2= {
    "select__cols": col_selector,
    "log_reg__C": np.logspace(-3, 4, 5), 
    "log_reg__penalty":["l1", "l2"] }


# In[ ]:


grid2 = GridSearchCV(pipe2, param_grid2, cv=10, n_jobs=-1, scoring='roc_auc')


# In[ ]:

print("Saving LR grid search")
grid2.fit(X_train, y_train)
print('Best parameters:', grid2.best_params_)
print('Best performance:', grid2.best_score_)


# In[ ]:


df_grid2_results = pd.DataFrame(grid2.cv_results_)


# In[ ]:


df_grid2_results  = df_grid2_results.sort_values('rank_test_score')


# In[ ]:


df_grid2_results.groupby('param_select__cols').first()


# In[ ]:




