#!/usr/bin/python3
import sys, runexps
import os


if __name__ == '__main__':
    left_pid = runexps.ps_check(runexps.ns3_left_name)
    right_pid = runexps.ps_check(runexps.ns3_right_name)
    os.system('echo a | sudo kill {}'.format(left_pid))
    os.system('echo a | sudo kill {}'.format(right_pid))
    os.system('ps aux | grep nmb')
    runexps.stop_iperfs()