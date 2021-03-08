#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 01:08:48 2021

@author: yusuf
"""

from ipynb.fs.full.CovidReceptorExpressionAnalysis_Functions import *

from io import StringIO
import sys

c1_patients = pd.read_csv('c1.txt',header=None)[0].tolist()
c2_patients = pd.read_csv('allpat.txt',header=None)[0].tolist()

for i in range(len(c1_patients)):
    c1_patients[i] = ('TARGET-30-' + c1_patients[i])
for i in range(len(c2_patients)):
    c2_patients[i] = ('TARGET-30-' + c2_patients [i])

c1_samples = clinical_sample_trimmed.loc[clinical_sample_trimmed['#Patient Identifier'].isin(c1_patients)]['Sample Identifier']
c2_samples = clinical_sample_trimmed.loc[clinical_sample_trimmed['#Patient Identifier'].isin(c2_patients)]['Sample Identifier']


old_stdout = sys.stdout
new_stdout = StringIO()
sys.stdout = new_stdout

compareEFS([c1_samples,c2_samples], 'EFS stratified by GSEA rank clusters', 'Cluster', 'Remaining')

output = new_stdout.getvalue()[-1]
