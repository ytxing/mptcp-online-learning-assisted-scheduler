#!/usr/bin/python3
import time
import os
import sys
import decodejson

def get_scheduler_csv(dir):

    for (root, subdirs, files) in os.walk(dir):
        subdirs.sort()
        files.sort()
        data = []
        csvFile = []
        column_num = 0

        for file in files:

            if file.find("json.txt") != -1:
                column_num += 1
                # print(file)
                with open( os.path.join(root,file), 'r' ) as f:
                    raw = f.readlines()[:-2]
                raw = [line.strip().split(',') for line in raw]
                
                this_data = [item[3] for item in raw]
                data.append(this_data)

        for item in data:
            if len(csvFile) == 0:
                csvFile = [[x] for x in item]
                continue
            for i in range(len(item)):
                csvFile[i].append(item[i])

        for i in range(len(csvFile)):
            csvFile[i] = ','.join(csvFile[i])
            
        if (len(csvFile) > 0):
            with open( os.path.join(root,root.split('/')[-1]+'.csv'), 'w' ) as f:
                header = ''
                for num in range(column_num):
                    header += root.split('/')[-1]+'-'+str(num)
                    if num < column_num - 1:
                        header += ','
                    else:
                        header += '\n'
                f.write(header)
                f.write('\n'.join(csvFile))

def get_all_csv(dir: str):
    header = []
    final_csv_file = []

    for (root, subdirs, files) in os.walk(dir):
        subdirs.sort()
        files.sort()
        for file in files:
            if file.find("final.csv") != -1:
                print('DELETE OLD CSV FILE')
                return

            if file.find(".csv") != -1:
                # print(file)
                with open( os.path.join(root,file), 'r' ) as f:
                    lines = f.readlines()
                for i in range(len(lines)):
                    lines[i] = lines[i].strip()
                if final_csv_file == []:
                    final_csv_file = lines
                    final_csv_file[0] = 'intervals,{}'.format(final_csv_file[0])
                    for i in range(1, len(final_csv_file)):
                        final_csv_file[i] = '{},{}'.format(i,final_csv_file[i])
                    continue
                for i in range(len(lines)):
                    final_csv_file[i] = final_csv_file[i]+','+lines[i]

    if (len(final_csv_file) > 0):
        name = dir.strip('/').split('/')[-1]
        with open( os.path.join(dir,'{}-final.csv'.format(name.replace('/', '-'))), 'w' ) as f:
            f.write('\n'.join(final_csv_file))

def get_trunk_csv(dir: str):
    header = []
    final_csv_file = []

    for (root, subdirs, files) in os.walk(dir):
        subdirs.sort()
        files.sort()
        for file in files:
            if file.find("final_trunk.csv") != -1:
                print('DELETE OLD CSV FILE')
                return

            if file.find("sum_bytes") != -1:
                # print(file)
                with open( os.path.join(root,file), 'r' ) as f:
                    lines = f.readlines()
                for i in range(len(lines)):
                    lines[i] = lines[i].strip()
                
                final_csv_file.append('{},{},{}'.format(root.split('/')[-1], lines[0].split(',')[3], lines[1].split(',')[2]))

    if (len(final_csv_file) > 0):
        name = dir.strip('/').split('/')[-1]
        with open( os.path.join(dir,'{}-trunk-final.csv'.format(name.replace('/', '-'))), 'w' ) as f:
            f.write('\n'.join(final_csv_file))
# def calssify(csv_file, key):
#     with open(csv_file) as f:
#             reader = csv.DictReader(f)
#             key_fileds = [r for r in reader.fieldnames if r.find(key) != -1]
#             for key_name in key_fileds:
#                 for line 
            
            

if __name__ == '__main__':
    decodejson.get_dir_json(sys.argv[1])
    get_scheduler_csv(sys.argv[1])
    get_trunk_csv(sys.argv[1])