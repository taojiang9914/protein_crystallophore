from sqlite3 import connect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,colors
from matplotlib.gridspec import GridSpec
from numpy import array, zeros
from math import ceil
from compiler.ast import flatten

#===================================
#Don't touch this part
def plot_scores(complementary_list,proteins,complexes_pairs):
    n_pairs = len(complementary_list)
    complementary_list = array(complementary_list)
    ind = np.arange(n_pairs)
    plt.rc('axes', axisbelow=True)
    p1 = plt.barh(ind, complementary_list[:,0], color=(0.8, 0.8, 0.8), edgecolor='black')
    p2 = plt.barh(ind, complementary_list[:,1], left = complementary_list[:,0], color=(0.4, 0.4, 0.4), edgecolor='black', hatch='.')
    p3 = plt.barh(ind, complementary_list[:,2], left = complementary_list[:,0]+complementary_list[:,1], color=(0.8, 0.8, 0.8), edgecolor='black', hatch='.')
    ylabels = proteins
    plt.yticks(ind, ylabels,rotation=45)
    plt.xlabel('Number of Hits')
    plt.grid(axis='x')
    plt.legend((p1[0],p2[0],p3[0]),('Native only','both','TbXo4 only'))
    plt.show()
#===================================

conn = connect('crystallophore.db')
c = conn.cursor()
complementary_list = []
native_only = []
Complex_only = []
both = []
proteins = ['pb9','ProteinaseK']
complexes_pairs = [['native','Xo4']
                   ]

for protein in proteins:
    for pair in complexes_pairs:
        pair = [protein+"_"+cmplx for cmplx in pair]
        sample_list = []
        both, p1, p2 = 0, 0,0
        for cmplx in pair:
            # Extract data from the database
            # Select the columns 'plate' and 'well' from the table 'sample' where the protein name is protein_cmplx 
            # In this script the protein names are 'pb9_native', 'pb9_Xo4', 'ProteinaseK_native' and 'ProteinaseK_Xo4'
            data = c.execute("SELECT plate, well FROM sample WHERE protein='"+cmplx+"'")
            sample_list.append(data.fetchall())
        for s in sample_list[0]:
            if s in sample_list[1]:
                both += 1
            else:
                p1 += 1
        for s in sample_list[1]:
            if not s in sample_list[0]:
                p2 += 1
        complementary_list.append([p1, both, p2])
    
plot_scores(complementary_list, proteins, complexes_pairs)
