"""
Created on 2017-07-18 by Edwin F. Juarez using the CCAL library created by Kwat Medetgul-Ernar Pablo Tamayo.

This module will grab a .gct file and a cls file to perform differential expression analysis.
"""
# import os
import sys
# sys.path.append(os.getcwd()+"/ccalnoir")
# from sys import path
# path.append('./ccalnoir')
# print(sys.path)

# import os
# out = open('stdout.txt', 'w')
# out.write(str(os.listdir('.')))

# import zipfile
# zip_ref = zipfile.ZipFile('ccal.zip', 'r')
# zip_ref.extractall('.')
# zip_ref.close()
import ccalnoir.ccalnoir as ccal
from ccalnoir.ccalnoir.mathematics.information import information_coefficient
import pandas as pd
import numpy as np
import scipy
# import math
import seaborn as sns
import matplotlib.pyplot as plt

def custom_pearson(x, y):
    return scipy.stats.pearsonr(x, y)[0]

TOP = 10

arg_n = len(sys.argv)
if arg_n == 1:
    err_out = open('stderr.txt', 'w')
    err_out.write("No files were provided. This module needs a GCT and a CLS file to work.")
    err_out.close()
    sys.exit("Error message: No files were provided. This module needs a GCT and a CLS file to work.")
elif arg_n == 2:
    err_out = open('stderr.txt', 'w')
    err_out.write("Only one file was provided (called = {}). This module needs a GCT and a CLS file to work.".format(sys.argv[0]))
    sys.exit("Only one file was provided (called = {}). This module needs a GCT and a CLS file to work.".format(sys.argv[0]))
elif arg_n == 3:
    gct_name = sys.argv[1]
    cls_name = sys.argv[2]
    function_to_call = custom_pearson
elif arg_n == 4:
    dispatcher = {
        "Pearson Correlation": custom_pearson,
        "PC": custom_pearson,
        "pc": custom_pearson,
        "correlation": custom_pearson,
        "Correlation": custom_pearson,
        "corr": custom_pearson,
        "Corr": custom_pearson,
        "Information Coefficient": information_coefficient,
        "IC": information_coefficient,
        "ic": information_coefficient,
    }
    gct_name = sys.argv[1]
    cls_name = sys.argv[2]
    try:

        function_to_call = dispatcher[sys.argv[3]]
        print('Using '+str(function_to_call)+' as the metric for similarity.')
    except KeyError:
        raise ValueError('This function is not supported at the moment, only Pearson Correlation and '
                         'Information Coefficient are supported at the moment.')
else:
    err_out = open('stderr.txt', 'w')
    err_out.write("Too many inputs. This module needs a GCT and a CLS file to work, "
                  "plus an optional input choosing between Pearson Correlation or Information Coefficient.")
    sys.exit("Too many inputs. This module needs a GCT and a CLS file to work, "
             "plus an optional input choosing between Pearson Correlation or Information Coefficient.")


out = open('stdout.txt', 'w')
# gct_name = "all_aml_test.preprocessed.gct"
df = pd.read_csv(gct_name, sep='\t', skiprows=2)
# cls_name = "all_aml_test.cls"
f = open(cls_name)
f.readline()
labels = np.asarray(f.readline().strip('\n').split(' '), dtype=str)[1:]
idx = np.asarray(f.readline().strip('\n').split(' '), dtype=float)
to_target = pd.Series(data=idx, index=list(df)[2:])

target, features, results = ccal.computational_cancer_biology.association.compute_association(
    target=to_target, features=df.iloc[:, 2:], function=function_to_call)

# input(results)


# import pickle
# pickle.dump((target, features, results, df), file=open("temp.p", 'wb'))
#
# import pickle
# (target, features, results, df) = pickle.load(open("temp.p", 'rb'))
# cls_name = "all_aml_test.cls"
# f = open(cls_name)
# f.readline()
# labels = np.asarray(f.readline().strip('\n').split(' '), dtype=str)[1:]

indexes = results.index.values
df = df.reindex(list(indexes))
features = features.reindex(list(indexes))
features = features.rename(df['Name'])


def make_label(row, labels=np.array([0, 1])):
    if row['Score'] > 0:
        idx = labels[1]
    else:
        idx = labels[0]
    return idx

out_df = pd.DataFrame()
out_df['Name'] = df.iloc[np.r_[0:TOP, -TOP:0], :]['Name']
out_df['Description'] = df.iloc[np.r_[0:TOP, -TOP:0], :]['Description']
out_df['Score'] = results.iloc[np.r_[0:TOP, -TOP:0], :]['score']
out_df['Differentially Expressed In'] = out_df.apply(make_label, args=(labels,), axis=1)

# out_df['Marker-of'] = ['ALL' if np.sign(out_df['Score']) > 0 else "AML"]
out = open('features.txt', 'w')
out.write(features.to_csv())
out.close()
out = open('results.txt', 'w')
out.write(results.to_csv())
out.close()
out = open('target.txt', 'w')
out.write(target.to_csv())
out.close()

out = open("scores.txt", 'w')
out.write(out_df.to_csv())
out.close()
# print(df.head(n=TOP))
# print(features.iloc[np.r_[0:TOP, -TOP:0], :])
sns.heatmap(features.iloc[np.r_[0:TOP, -TOP:0], :], cmap='coolwarm')
plt.yticks(rotation=0)
plt.xticks(rotation=90)
plt.savefig('heatmap.png', dpi=300)

plt.clf()
sns.barplot(y='Score', x='Name', data=out_df, hue='Differentially Expressed In')
plt.title('Differential Expression Analysis')
# my_cmap = sns.color_palette(as_cmap=True)
# plt.ylabel('Similarity Metric\nNegative val. diff. ex. in '+labels[0]+' Positive val. diff. ex. ib'+labels[1])
plt.xticks(rotation=45, horizontalalignment="right")
# [t.set_color(i) for (i, t) in zip(cmap(out_df['Score']), plt.gca().xaxis.get_ticklabels())]
plt.axvline(x=TOP-0.5, linestyle='--', color='gray')
plt.axhline(y=0, linestyle='-', color='k')
plt.ylabel('Similarity Metric')
plt.xlabel('Gene Name')
# plt.tight_layout()
plt.savefig('scores.png', dpi=300, bbox_inches="tight")
