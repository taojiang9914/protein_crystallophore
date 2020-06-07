from sqlite3 import connect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,colors
from matplotlib.gridspec import GridSpec

#===================================
#Don't touch this part
def plot_scores(total_score_cond,total_score_crys,kit_names,conditions,protein,complexes,no_mask,score):
    ncols = 3
    nrows = 2
    if conditions != ['ph']:
        #with Multi
        #vmax = len(conditions)+2
        #no Multi
        vmax = len(conditions)+1
    else:
        vmax = 8

    cmaplist = [[1., 1., 1.], [0., .4, 1.], [0., .8, 1.],[1., .8, 0.], [1., .4, 0.], [1., 0., 0.], [0.6, 0.2, 0.9]]
    if conditions == ['ph']:
            cmaplist = [[1., 1., 1.], [0., .4, 1.], [0., .8, 1.],[0., 1., 0.],[1., .8, 0.], [1., .4, 0.], [1., 0., 0.], [0.6, 0.2, 0.9]]
    if conditions != ['ph']:
        #with Multi
        #cmap = colors.ListedColormap(cmaplist[:len(conditions)+2])
        #no Multi
        cmap = colors.ListedColormap(cmaplist[:len(conditions)+1])
    else:
        cmap = colors.ListedColormap(cmaplist)
    # create the plots
    if no_mask:
        fig = plt.figure(figsize=(10,6))
    else:
        fig = plt.figure(figsize=(10,8))
    axes = [fig.add_axes([0,0.9,1,0])]
    for r in range(nrows):
        for c in range(ncols):
            axes.append(fig.add_axes([c*0.28+0.1, (1-r)*0.4+0.15, 0.25, 0.30]))

    # add some data
    if no_mask:
        total_score_cond[:] = 0
    for i,ax in enumerate(axes[1:]):
        im = ax.pcolor(total_score_cond[i],vmin=0, vmax=vmax, cmap=cmap,edgecolors='k', linewidths=1)
    if not no_mask:
        if conditions != ['ph']:
            cbar = fig.colorbar(im, ax=axes,orientation='horizontal',aspect=65, shrink=0.6,ticks=range(1,len(conditions)+2))
        else:
            cbar = fig.colorbar(im, ax=axes,orientation='horizontal',aspect=65, shrink=0.6,ticks=range(1,8))
        cbar.ax.get_xaxis().set_ticks([])

        if conditions != ['ph']:
            cond_list = ['None']
            cond_list.extend(conditions)
            #cond_list.append('Multi')
        else:
            cond_list = ['None']
            cond_list.extend(['4.0-5.0', '5.0-6.0', '6.0-7.0', '7.0-8.0', '8.0-9.0', '9.0-10.0','>10.0'])
        for j, cond in enumerate(cond_list):
            cond = cond.replace('"','')
            if len(cond) > 10:
                cond = cond[:10]
            cbar.ax.text((2 * j + 1) / (2*float(len(cond_list))), -1.0, cond, ha='center', va='center')



    # remove the x and y ticks
    axes[0].axis('off')
    axes[0].set_title(protein+' score: '+str(score) ,fontdict={'fontsize':20})
    for i,ax in enumerate(axes[1:]):
        ax.set_title(kit_names[i])
        ax.set_xlim([0,12])
        ax.set_xticks([k+0.5 for k in range(0,12)])
        ax.set_xticklabels([str(k) for k in range(1,13)])
        ax.set_yticks([v+0.5 for v in row_dict.values()])
        ax.set_yticklabels(row_dict.keys())
        ax.xaxis.set_ticks_position('none') 
        ax.yaxis.set_ticks_position('none') 

    for i,ax in enumerate(axes[1:]):
        for j in range(12):
            for k in range(8):
                if total_score_crys[i,k,j,0] == 1:
                    ax.fill([j,j+1,j+1,j],(k,k,k+1,k+1),hatch='o',fill=False)
                if total_score_crys[i,k,j,1] == 1:
                    ax.fill([j,j+1,j+1,j],(k,k,k+1,k+1),hatch='\\\\',fill=False)
                if total_score_crys[i,k,j,2] == 1:
                    ax.fill([j,j+1,j+1,j],(k,k,k+1,k+1),hatch='//',fill=False)

    ps = 'o: '+ complexes[0]+'    \: '+complexes[1]
    plt.text(-10.8, -3, ps)
    nametoken = protein
    for i,c in enumerate(complexes):
        nametoken += '_'
        nametoken += complexes[i][:-4]


    plt.show()
#===================================

proteins = ['HEWL','ProteinaseK','TdTau']
complexes = ['TbXo4_15j','LuXo4_15j']
#conditions = ['ph']
conditions = []
#With condition mask or not
no_mask = True
#no_mask = False
score = 'all'
#score = 6

conn = connect('sengilberge.db')
c = conn.cursor()
# Find out the maximum number of drops
# The plates have 8 rows 'A,B,...,H' and 12 columns
# Creating a 6 by 8 by 12 matrix to store the score sum
row_dict = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7}
kit_names = ['Wizard_I+II_rigaku',
            'Salt-Grid_hampton',
            'JCSG_MD',
            'PEGs-I_qiagen',
            'PACT_MD',
            'Classics-Suite_qiagen']
six_kits = [kit.lower() for kit in kit_names]

total_score_cond = np.zeros([6,8,12])
if not no_mask:
    if conditions != ['ph']:
        data = []
        for cond in conditions:
            # Getting data from the detailed_plate table with given chemical composition 
            data.append([d for d in c.execute("SELECT name, well, "+cond+" FROM detailed_plate").fetchall() if d[-1] != '0' ])
        for i,_data in enumerate(data):
            for d in _data:
                plate, well = d[0], d[1]
                ind_plate = kit_names.index(plate)
                ind_well = [row_dict[d[1][0]], int(d[1][1:])-1]
                if total_score_cond[ind_plate, ind_well[0],ind_well[1]] == 0:
                    total_score_cond[ind_plate, ind_well[0],ind_well[1]] = i+1
                else:
                    total_score_cond[ind_plate, ind_well[0],ind_well[1]] = len(conditions)+1
    else:
        # Getting pH value from the detailed_plate table for all wells with pH information
        data = [d for d in c.execute("SELECT name, well, ph FROM detailed_plate").fetchall() if d[-1] != '0' ]
        for d in data:
            plate, well = d[0], d[1]
            ind_plate = kit_names.index(plate)
            ind_well = [row_dict[d[1][0]], int(d[1][1:])-1]
            if float(d[2])<4.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 0
            elif float(d[2])<5.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 1
            elif float(d[2])<6.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 2
            elif float(d[2])<7.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 3
            elif float(d[2])<8.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 4
            elif float(d[2])<9.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 5
            elif float(d[2])<10.0:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 6
            else:
                total_score_cond[ind_plate, ind_well[0],ind_well[1]] = 7

for protein in proteins:
    total_score_crys = np.zeros([6,8,12,3])
    for i,cmplx in enumerate(complexes):
        _score = np.zeros([6,8,12])
        # Getting the plate and well data in which certain protein complex combintion has crystallized with certain score from the sample table
        if score == 'all':
            data = c.execute("SELECT protein, plate, well FROM sample WHERE protein='"+protein+"_"+cmplx+"'")
        elif score == 6:
            data = c.execute("SELECT protein, plate, well FROM sample WHERE protein='"+protein+"_"+cmplx+"' AND score='6'")
        elif score == 5:
            data = c.execute("SELECT protein, plate, well FROM sample WHERE protein='"+protein+"_"+cmplx+"' AND score='5'")
        elif score == 4:
            data = c.execute("SELECT protein, plate, well FROM sample WHERE protein='"+protein+"_"+cmplx+"' AND score='4'")
        native_samples = []
        for d in data:
            plate, well = d[1], d[2]
            ind_plate = kit_names.index(plate)
            ind_well = [row_dict[d[2][0]], int(d[2][1:])-1]
            total_score_crys[ind_plate, ind_well[0],ind_well[1],i] = 1
            _score[ind_plate, ind_well[0],ind_well[1]] = 1

    plot_scores(total_score_cond, total_score_crys, kit_names, conditions,protein,complexes,no_mask,score)
