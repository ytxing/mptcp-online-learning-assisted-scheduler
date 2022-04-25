#!/usr/bin/python3
import time, os, sys


timer = int(sys.argv[1])
for i in range(timer):
    time.sleep(1)
    print("count down: {}/{}".format(i, timer), end='\r')

while True:
    timer = 600
    do_flag = 1
    for i in range(timer):
        time.sleep(1)
        print("count down: {}/{}".format(i, timer), end='\r')

    result = os.popen("ps aux | grep runexps")
    lines = result.readlines()
    for x in lines:
        if x.find("python4 runexps.py") != -1: 
            print(x.strip())
            do_flag = 0
    
    if do_flag == 1:
        break

print("\nrunexps!")
os.system("python3 runexps.py")
# os.system("ifstat")