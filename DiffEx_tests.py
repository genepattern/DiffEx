from subprocess import call
import pandas as pd
from os.path import isfile

# TODO: if the user wants to select a different version of Python.
log = open('testlog.txt', 'w')
error_count = 0
# Call the module
call("Python DiffEx.py test_dataset.gct test_dataset.cls 5", shell=True)
# Check if the outputs were created:
file_list = ['target.txt', 'stdout.txt', 'scores.txt', 'scores.png', 'results.txt', 'heatmap.png', 'features.txt']
file_counter = 0
for file in file_list:
    if isfile(file):
        file_counter += 1
    else:
        print('Test failed! Missing file ', file)
        print('Test failed! Missing file ', file, file=log)
        error_count = 1
print('{} out of {} files found.'.format(file_counter, len(file_list)))
print('{} out of {} files found.'.format(file_counter, len(file_list)), file=log)

df = pd.read_csv("scores.txt", sep=',')
df.rename(columns={'Unnamed: 0': 'gene_index'}, inplace=True)  # Renaming a column
df.sort_index(inplace=True)
right_PC_order = pd.read_csv('scores_PC.txt', sep=',')
right_PC_order.rename(columns={'Unnamed: 0': 'gene_index'}, inplace=True)  # Renaming a column
right_PC_order.sort_index(inplace=True)
# Check if the ranked genes are in the correct order.
if df['Name'].equals(right_PC_order['Name']):
    print('Pearson Correlation orders genes correctly.')
    print('Pearson Correlation orders genes correctly.', file=log)
else:
    print('Test failed! DiffEx is not ordering the genes correctly using Pearson Correlation.')
    print('Test failed! DiffEx is not ordering the genes correctly using Pearson Correlation.', file=log)
    error_count += 1

# Check if the ranked genes are in the correct order.
call("Python DiffEx.py test_dataset.gct test_dataset.cls 5 IC", shell=True)
df = pd.read_csv("scores.txt", sep=',')
df.rename(columns={'Unnamed: 0': 'gene_index'}, inplace=True)  # Renaming a column

right_IC_order = pd.read_csv('scores_IC.txt', sep=',')
right_IC_order.rename(columns={'Unnamed: 0': 'gene_index'}, inplace=True)  # Renaming a column
right_IC_order.sort_index(inplace=True)
if df['Name'].equals(right_IC_order['Name']):
    print('Information Correlation orders genes correctly.')
    print('Information Correlation orders genes correctly.', file=log)
else:
    print('Test failed! DiffEx is not ordering the genes correctly using Information Coefficient.')
    print('Test failed! DiffEx is not ordering the genes correctly using Information Coefficient.', file=log)
    error_count += 1
if error_count == 0:
    print('Yay! DiffEx is working as expected')
    print('Yay! DiffEx is working as expected', file=log)
log.close()
