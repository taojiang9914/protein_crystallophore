import csv
import sqlite3
import os
import glob

kit_names = ['Wizard_I+II_rigaku',
            'Salt-Grid_hampton',
            'JCSG_MD',
            'PEGs-I_qiagen',
            'PACT_MD',
            'Classics-Suite_qiagen']
#conn1 = sqlite3.connect('publication.db')
#c1 = conn1.cursor()
##c1.execute('''CREATE TABLE sample (id INTEGER PRIMARY KEY AUTOINCREMENT, protein text, plate text, well text)''')
##conn1.commit() 
##ind_sample = 1
#tmp = c1.execute('''SELECT * FROM sample ORDER BY id DESC''')
#ind_sample = int(next(tmp)[0])+1
#conn2 = sqlite3.connect('sengilberge.db')
#c2 = conn2.cursor()
#
##proteins = ['pb9','ProteinaseK']
##complexes_pairs = [['native','Xo4']
##                   ]
##proteins = ['HEWL','ProteinaseK','TdTau']
##complexes = ['TbXo4_15j','LuXo4_15j']
#
#proteins = ['TdTau','HEWL','GRHPR','Protease1']
#complexes = ['native','Xo4']
#for protein in proteins:
#    for cmplx in complexes:
#        sample_list = []
#        protein_cmplx = protein+"_"+cmplx
#        data = c2.execute("SELECT plate, well FROM sample WHERE protein='"+protein_cmplx+"'")
#        sample_list.extend(data.fetchall())
#        
#        for s in sample_list:
#            print(s[0])
#            val_str = str(ind_sample)
#            val_str += ", '"+protein_cmplx+"', '"+s[0]+"', '"+s[1]+"'"
#            str_total ="INSERT INTO sample (id, protein, plate, well) VALUES ("
#            str_total += val_str
#            str_total += ")"
#            c1.execute(str_total)
#            conn1.commit()
#            ind_sample += 1

conn1 = sqlite3.connect('publication.db')
c1 = conn1.cursor()
#c1.execute('''CREATE TABLE detailed_plate (id INTEGER PRIMARY KEY AUTOINCREMENT, name text, well text, ph text)''')
ind_sample = 1
conn2 = sqlite3.connect('sengilberge.db')
c2 = conn2.cursor()
data = c2.execute('''SELECT name, well, ph FROM detailed_plate''').fetchall()
c2 = conn2.cursor()
for i,d in enumerate(data):
    print(i,d)
    val_str = str(i+1)+",'"+d[0]+"','"+d[1]+"','"+d[2]+"')"
    print('''INSERT INTO detailed_plate (id, name, well, ph) VALUES('''+val_str)
    c1.execute('''INSERT INTO detailed_plate (id, name, well, ph) VALUES('''+val_str)
conn1.commit() 
