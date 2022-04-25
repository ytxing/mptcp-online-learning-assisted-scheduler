#!/usr/bin/python3
import json
import os
import sys

def decode_json(json_file, file_name):
    with open(json_file, 'r' ) as f:
        json_data_lines = f.readlines()

    if len(json_data_lines) == 0:
        print('EMPTY!! {}'.format(json_file))
        return
    if os.path.isfile('{}.txt'.format(json_file)):
        os.system("rm {}.txt".format(json_file))
    
    count_lb = 0
    count_rb = 0
    json_datas = []
    one_json = ''
    for i in range(len(json_data_lines)):
        line = json_data_lines[i]
        one_json += line
        if line == "{\n":
            count_lb += 1
        if line == "}\n":
            count_rb += 1
            json_datas.append(one_json)
            one_json = ''

    for i in range(len(json_datas)):
        json_data = json_datas[i]
        data = json.loads(json_data)

        if 'error' in data:
            print('ERROR!! {}'.format(json_file))
            continue

        bytes = data['start']['test_start']['bytes']
        duration = data['start']['test_start']['duration']

        intervals = data['intervals']
        for interval in intervals:
            stream = interval['streams'][0]
            # print('{:.2f},{:.2f},{:.2f}'.format(stream['start'], stream['end'], stream['bits_per_second']/(1024*1024)))
            with open(json_file+'.txt'+str(i), 'a' ) as f:
                f.write('interval,{:.2f},{:.2f},{:.2f}\n'.format(stream['start'], stream['end'], stream['bits_per_second']/(1024*1024)))
        sum_rcv = data['end']['sum_received']
        sum_snt = data['end']['sum_sent']
        with open(json_file+'.txt'+str(i), 'a' ) as f:
                f.write('sum_rcv_thp{},{:.2f},{:.2f},{:.2f}\n'.format(file_name, sum_rcv['start'], sum_rcv['end'], sum_rcv['bits_per_second']/(1024*1024)))
                f.write('sum_snt_thp{},{:.2f},{:.2f},{:.2f}\n'.format(file_name, sum_snt['start'], sum_snt['end'], sum_snt['bits_per_second']/(1024*1024)))
                if bytes != '0':
                    f.write('bytes{},0.00,0.00,{}\n'.format(file_name, bytes))
                    f.write('sum_rcv_bytes{},{},{},{}\n'.format(file_name, sum_rcv['start'], sum_rcv['end'], sum_rcv['bytes']))
                    f.write('sum_snt_bytes{},{},{},{}\n'.format(file_name, sum_snt['start'], sum_snt['end'], sum_snt['bytes']))


# dir='/home/a/mptcp-mptcp_v0.95/log'
def get_dir_json(dir):
    for (root, subdirs, files) in os.walk(dir):
        # print(root)
        flag = 0
        subdirs.sort()
        files.sort()
        for file in files:
            if file.find(".json") != -1 and file.find(".txt") == -1:
                # print(file)
                flag = 1
                file_path = os.path.join(root,file)
                decode_json(file_path, file)

        if flag == 1:
            os.system('cat {}/*.txt* | grep sum_rcv_thp > {}/sum.txt'.format(root, root))
            throughputs = []
            with open('{}/sum.txt'.format(root), 'r' ) as f:
                lines = f.readlines()
                for line in lines:
                    throughputs.append(float(line.split(',')[-1]))
            with open('{}/sum.txt'.format(root), 'a' ) as f:
                if len(throughputs) == 0:
                    f.write('avg,0.00,0.00,0.00')
                else:
                    f.write('avg,0.00,0.00,{:.2f}'.format(sum(throughputs)/len(throughputs)))

            os.system('cat {}/*.txt* | grep bytes > {}/sum_bytes.txt'.format(root, root))

if __name__ == '__main__':
    get_dir_json(sys.argv[1])
    # new_decode_json(sys.argv[1], "test")
