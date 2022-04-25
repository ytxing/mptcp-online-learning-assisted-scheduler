#!/usr/bin/python3
import sys, runexps



if __name__ == '__main__':
    left_pid = runexps.ps_check(runexps.ns3_left_name)
    right_pid = runexps.ps_check(runexps.ns3_right_name)
    runexps.cmd_local('echo a | sudo kill {}'.format(left_pid))
    runexps.cmd_local('echo a | sudo kill {}'.format(right_pid))
    runexps.cmd_local('ps aux | grep nmb')