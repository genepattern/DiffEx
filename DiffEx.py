"""
Created on 2017-07-18 by Edwin F. Juarez using the CCAL library created by Kwat Medetgul-Ernar Pablo Tamayo.

This module will grab a .gct file and a .cls file to perform differential expression analysis.
"""

import os
import sys
tasklib_path = os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(tasklib_path + "/ccalnoir")
import matplotlib as mpl
mpl.use('Agg')
import ccalnoir as ccal
from ccalnoir.mathematics.information import information_coefficient
import pandas as pd
import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt


def custom_pearson(x, y):
    return scipy.stats.pearsonr(x, y)[0]

# Error handling:
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
    TOP = 10
    function_to_call = custom_pearson
elif arg_n == 4:
    gct_name = sys.argv[1]
    cls_name = sys.argv[2]
    TOP = int(sys.argv[3])
    function_to_call = custom_pearson
elif arg_n == 5:
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
    TOP = int(sys.argv[3])
    try:
        function_to_call = dispatcher[sys.argv[4]]
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

df = pd.read_csv(gct_name, sep='\t', skiprows=2)
f = open(cls_name)
f.readline()
labels = np.asarray(f.readline().strip('\n').split(' '), dtype=str)[1:]
idx = np.asarray(f.readline().strip('\n').split(' '), dtype=float)
to_target = pd.Series(data=idx, index=list(df)[2:])

out.write("Successfully read file {} and {}, will use {} as the similarity metric.\n"
          .format(gct_name, cls_name, function_to_call.__name__))

target, features, results = ccal.computational_cancer_biology.association.compute_association(
    target=to_target, features=df.iloc[:, 2:], function=function_to_call)
out.write("All analyses complete.")
out.close()

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

# TODO: Make these outputs optional.

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
sns.heatmap(features.iloc[np.r_[0:TOP, -TOP:0], :], cmap='coolwarm')
plt.yticks(rotation=0)
plt.xticks(rotation=90)
plt.title('Top {} differentially expressed genes per phenotype'.format(TOP))
plt.ylabel('Genes')
plt.xlabel('Sample')
plt.savefig('heatmap.png', dpi=300, bbox_inches="tight")

plt.clf()
sns.barplot(y='Score', x='Name', data=out_df, hue='Differentially Expressed In')
plt.title('Differential Expression Analysis')
plt.xticks(rotation=45, horizontalalignment="right")
plt.axvline(x=TOP-0.5, linestyle='--', color='gray')
plt.axhline(y=0, linestyle='-', color='k')
plt.ylabel('Similarity Metric')
plt.xlabel('Gene Name')
plt.savefig('scores.png', dpi=300, bbox_inches="tight")
