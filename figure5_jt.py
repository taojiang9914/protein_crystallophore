from sqlite3 import connect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm,colors
from matplotlib.gridspec import GridSpec
from numpy import array, zeros
from math import ceil

#===================================
#Normally you don't need to touch this part
def plot_scores(total_score_crys,conditions, condition_list,kit_names,protein,complexes,name_token, no_mask, show_plot):
    score = total_score_crys.reshape([6,-1,3]).sum(axis=1)
    score = [len(s) for s in condition_list]
    ncols = 3
   
    if conditions == ['ph']:
        full_cmaplist = [[0.9,0.9,0.9], [0., .4, 1.], [0., .8, 1.],[1., .8, 0.], [1., .4, 0.], [1., 0., 0.], [0.6, 0.2, 0.9],[0.5,0.5,0.5],[1,1,1]]
        cmaplist = full_cmaplist[:]
    elif conditions == ['plate']:
        full_cmaplist = [[0., .4, 1.], [0., .8, 1.],[1., .8, 0.], [1., .4, 0.], [1., 0., 0.], [0.6, 0.2, 0.9],[1,1,1]]
        cmaplist = full_cmaplist[:]
    else:
        full_cmaplist = [[0.9,0.9,0.9], [0., .4, 1.], [0., .8, 1.],[1., .8, 0.], [1., .4, 0.], [1., 0., 0.]]
        cmaplist = full_cmaplist[:len(conditions)+1]
        cmaplist.append([0.6, 0.2, 0.9])
        cmaplist.append([1,1,1])
    cmap = colors.ListedColormap(cmaplist[:])
    if conditions == ['ph']:
        vmax = 9
    elif conditions == ['plate']:
        vmax = 7
    else:
        vmax = len(conditions)+3

    fig = plt.figure(figsize=(10,7))
    fig.patch.set_facecolor('white')
    axes = []
    vspace = 0
    vspace_t = 0
    cell_v_size=0.036
    title_v_size=0.018

    #making cells for the r-th plate
    grid = [int(ceil(max(score)/10.)),10]
    vspace += cell_v_size*grid[0] + 0.04
    vspace_t += title_v_size*grid[0] + 0.04
    for c in range(ncols):
        axes.append(fig.add_axes([c*0.28+0.1, 0.9-vspace, 0.25, cell_v_size*grid[0]]))
        if conditions == ['ph']:
            tmp_score = zeros(int(grid[0]*grid[1]))+8
        elif conditions == ['plate']:
            tmp_score = zeros(int(grid[0]*grid[1]))+7
        else:
            tmp_score = zeros(int(grid[0]*grid[1]))+len(conditions)+2
        if no_mask:
            tmp_score[:int(score[c])] = 0
        else:
            tmp_score[:int(score[c])] = [entry for entry in condition_list[c]]
        tmp_score = tmp_score.reshape(grid)
        im = axes[-1].pcolor(tmp_score,edgecolors='w', linewidths=1,vmin=0,vmax=vmax,cmap=cmap)
    axes[0].axis('off')
    for i,ax in enumerate(axes):
        plt.setp(ax.spines.values(), visible=False)
        ax.xaxis.set_visible(False)
        ax.tick_params(left=False, labelleft=False, right=False, labelright=False)
    titles = [complexes[0][0].upper()+complexes[0][1:]+' unique','Common','Tb-'+complexes[1][0].upper()+complexes[1][1:]+' unique']
    for i in range(3):
        axes[i].set_title(titles[i],fontsize=15,loc='center')

    if not no_mask:
        #Add chemical condition on the SOI
        cbaxes = fig.add_axes([0.28, 0.8-vspace, 0.5, 0.03]) 
        if conditions == ['ph']:
            lgd_txt = ['None','4.0-5.0','5.0-6.0','6.0-7.0','7.0-8.0','8.0-9.0','9.0-10.0','>10.0']
            plt.text(-0.1, -0.2-vspace,'pH:')
        elif conditions == ['plate']:
            lgd_txt = [kit_names.keys()[kit_names.values().index(i)] for i in range(6)]
            plt.text(-0.1, 0.37-vspace,'plate:')
        else:
            lgd_txt = ['None']
            for c in conditions:
                lgd_txt.append(c.replace('"',''))
            lgd_txt.append('multi')
        cbar = fig.colorbar(im, cax = cbaxes,ticks = np.arange(0.5,vmax,1),orientation='horizontal')
        cbar.ax.set_xticklabels(lgd_txt,rotation=-70)
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(size=0)

        if conditions == ['plate']:
            xloc = cbar.ax.get_xticks()
            hit_counter = np.zeros([6])
            for cond in condition_list:
                for c in cond:
                    hit_counter[c] += 1
            for i in range(6):
                plt.text(xloc[i]-0.01, 0.77-vspace,str(int(hit_counter[i])),color='white')
    title=protein+' '+complexes[0][0].upper()+complexes[0][1:]+' vs Tb-'+complexes[1][0].upper()+complexes[1][1:]
    fig.suptitle(title,fontsize=16,y=0.95)
    if not no_mask:
        if conditions[0] == 'ph':
            plt.savefig(path_prefix+'/soi_'+protein+'_'+name_token+'_overlay_ph.png',dpi=300)
        elif conditions[0] == 'plate':
            plt.savefig(path_prefix+'/soi_'+protein+'_'+name_token+'_overlay_plate.png',dpi=300)
    else:
        plt.savefig(path_prefix+'/soi_'+protein+'_'+name_token+'.png',dpi=300)
        plt.savefig(path_prefix+'/soi_'+protein+'_'+name_token+'.eps',dpi=300)
    if show_plot:
        plt.show()
#===================================


no_mask = False
show_plot = True

all_conditions = [['plate'],['ph']]

proteins = ['pb9','ProteinaseK', 'TdTau','HEWL','GRHPR','Protease1']

complexes = ['native','Xo4']
name_token  = complexes[0]+'_'+ complexes[1]

path_prefix = 'figures'

conn = connect('sengilberge.db')
c = conn.cursor()
# Find out the maximum number of drops
# The plates have 8 rows 'A,B,...,H' and 12 columns
# Creating a 6 by 8 by 12 array to store the score sum
row_dict = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7}
kit_names = {'Wizard_I+II_rigaku':0,
             'Salt-Grid_hampton':1,
             'JCSG_MD':2,
             'PEGs-I_qiagen':3,
             'PACT_MD':4,
             'Classics-Suite_qiagen':5}
six_kits = [kit.lower() for kit in kit_names]

total_score_crys = np.zeros([6,8,12])


for conditions in all_conditions:
    sample_list = []
    for i,protein in enumerate(proteins):
        for j,cmplx in enumerate(complexes):
            protein_name  = protein + "_" + cmplx
            if conditions != ['plate']:
                #Use table join operation to obtain protein crystallization condition table 
                command = '''SELECT sample.plate, sample.well'''
                for cond in conditions:
                    command += ''',detailed_plate."''' + cond + '''"'''
                command += '''\nFROM sample
                              INNER JOIN detailed_plate
                              ON sample.plate = detailed_plate.name
                              AND sample.well = detailed_plate.well
                              WHERE sample.protein ="'''+protein_name+'''"'''
            else:
                command = '''SELECT plate, well FROM sample WHERE protein="'''+protein_name+'''"'''
            data = c.execute (command)
            sample_list.append(data.fetchall())
            if i != 0: 
                if j == 0:
                    sample_list[0].extend(sample_list[-1])
                elif j == 1:
                    sample_list[1].extend(sample_list[-1])
                sample_list.pop()
    sample_list[0] = sorted(set(sample_list[0]))
    sample_list[1] = sorted(set(sample_list[1]))
    condition_list = []
    native_only = []
    complex_only = []
    both = []
    
    for s in sample_list[0]:
        conds = []
        #Assign each condition a number which will be used as a color code when plotting
        if conditions == ['plate']:
            cond_id = kit_names[s[0]]
        elif conditions == ['ph']:
            if float(s[2])<4.0:
                cond_id = 0
            elif float(s[2])<5.0:
                cond_id = 1
            elif float(s[2])<6.0:
                cond_id = 2
            elif float(s[2])<7.0:
                cond_id = 3
            elif float(s[2])<8.0:
                cond_id = 4
            elif float(s[2])<9.0:
                cond_id = 5
            elif float(s[2])<10.0:
                cond_id = 6
            else:
                cond_id = 7
        else:
            non_zeros = []
            for i,rec in enumerate(s[2:]):
                if rec != '0': non_zeros.append(i)
            if len(non_zeros) > 1:
                cond_id = len(conditions)+1
            elif len(non_zeros) == 1:
                cond_id = non_zeros[0]+1
            elif len(non_zeros) == 0:
                cond_id = 0
         
        #Check the synergetic effect
        if s in sample_list[1]:
            both.append(cond_id)
        else:
            native_only.append(cond_id)
    
    for s in sample_list[1]:
        conds = []
        #Assign each condition a number which will be used as a color code when plotting
        if conditions == ['plate']:
            cond_id = kit_names[s[0]]
        elif conditions == ['ph']:
            if float(s[2])<4.0:
                cond_id = 0
            elif float(s[2])<5.0:
                cond_id = 1
            elif float(s[2])<6.0:
                cond_id = 2
            elif float(s[2])<7.0:
                cond_id = 3
            elif float(s[2])<8.0:
                cond_id = 4
            elif float(s[2])<9.0:
                cond_id = 5
            elif float(s[2])<10.0:
                cond_id = 6
            else:
                cond_id = 7
        else:
            non_zeros = []
            for i,rec in enumerate(s[2:]):
                if rec != '0': non_zeros.append(i)
            if len(non_zeros) > 1:
                cond_id = len(conditions)+1
            elif len(non_zeros) == 1:
                cond_id = non_zeros[0]+1
            elif len(non_zeros) == 0:
                cond_id = 0
            
        #Check the synergetic effect
        if not s in sample_list[0]:
            complex_only.append(cond_id)
    condition_list.append(sorted(native_only))
    condition_list.append(sorted(both))
    condition_list.append(sorted(complex_only))
    
    plot_scores(total_score_crys, conditions, condition_list,kit_names, '6proteins',complexes,name_token, no_mask, show_plot)
