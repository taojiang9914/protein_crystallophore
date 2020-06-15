import csv
import sqlite3
import os
import glob


path='.'
sample_files = glob.glob("samples/*15j.csv")
print(sample_files)
test_run = False


kit_names = ['Wizard_I+II_rigaku',
            'Salt-Grid_hampton',
            'JCSG_MD',
            'PEGs-I_qiagen',
            'PACT_MD',
            'Classics-Suite_qiagen']
conn = sqlite3.connect('example_crystallophore.db')
c = conn.cursor()
c.execute('''CREATE TABLE sample (id INTEGER PRIMARY KEY AUTOINCREMENT, protein text, plate text, well text)''')
conn.commit() 

ind_sample = 1

for sample_file in sample_files:
    print('starting '+sample_file)
    samples = csv.reader(open("/".join((path,sample_file)),'rb'), delimiter=';')#, quotechar='|')
    next(samples)
    if '/' in sample_file:
        protein = sample_file.split('/')[-1].split('.')[0]
    else:
        protein = sample_file.split('.')[0]
    for sample in samples:
        val_str = str(ind_sample)

        if len(sample) > 0:
            score = sample[1]
            #if 1:
            if score == '6':
                if sample[5][1] == '0':
                    well = sample[5][0]+sample[5][2]
                else:
                    well = sample[5][:3]
                if sample[10] not in kit_names:
                    raise('Kit name not in the list, please verify!!')
                else:
                    plate = sample[10]
#                try:
#                    int(sample[1])
#                except ValueError:
#                    if test_run:
#                        print(sample[1])
#                        pass
#                        #raise('strange score, please verify!!')
#                    else:
#                        pass

                 # The columns that correspond to the data can vary depending on the provider of the data file, verify that before.
                val_str += ", '"+protein+"', '"+sample[10]+"', '"
                val_str += well
                val_str += "'"

                str_total ="INSERT INTO sample (id, protein, plate, well) VALUES ("
                str_total += val_str
                str_total += ")"
                print str_total
                ind_sample += 1
                if not test_run:
                    c.execute(str_total)
                    conn.commit() 
    print('=======>',sample_file,' file done')
conn.close() 
