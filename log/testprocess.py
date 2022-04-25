#!/usr/bin/python3
import time, multiprocessing

#假定这是你的银行存款
timeout = 0

def iperf(t):
    for i in range(t):
        print('{}, iperfing'.format(i), end='\r')
        time.sleep(1)

def timer(t, iperf_t: multiprocessing.Process):
    for i in range(t):
        print('{}, timing'.format(i), end='\r')
        time.sleep(1)
        if iperf_t.is_alive() == False:
            return
    if iperf_t.is_alive():
        iperf_t.terminate()
        print('{}, KILL'.format(i))

def main(t):
    iperf_t = multiprocessing.Process(target=iperf, args=(t+4,))
    
    iperf_t.start()
    timer(t+3, iperf_t)
    iperf_t.join()
    

main(5)