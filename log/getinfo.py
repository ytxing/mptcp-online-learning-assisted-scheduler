#!/usr/bin/python3
import time
import os
import sys

def get_info_dir(dir):
    init_pattern = "ol_init"
    meta_pattern = "meta_tp"

    for (root, subdirs, files) in os.walk(dir):
        print(root)
        if root.find('details') != -1:
            continue

        for file in files:                
            tps = []
            tp_type = []
            if file.find("kernlog.log") != -1:
                os.system("mkdir -p {}/{}-details".format(root, file))
                count = 0
                with open( os.path.join(root,file), 'r' ) as f:
                    for line in f.readlines():
                        count += 1
                        if line.find(init_pattern) != -1:
                            print(line)
                            tps.append(line.replace("\n", "")[-16:])
                            if line.find(meta_pattern) != -1:
                                tp_type.append("meta")
                            else:
                                tp_type.append("sub")
                        if count > 999:
                            break
                
                for i in range(len(tps)):
                    os.system("cat {} | grep \"{}\" > {}/{}-details/{}-{}.log".format(os.path.join(root,file), tps[i], root, file, tp_type[i], tps[i]))
                    print("cat {} | grep \"{}\" > {}/{}-details/{}-{}.log".format(os.path.join(root,file), tps[i], root, file, tp_type[i], tps[i]))

if __name__ == '__main__':
    get_info_dir(sys.argv[1])