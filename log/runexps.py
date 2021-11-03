#!/usr/bin/python3
import time
import os
import decodejson
import getinfo
import mycsv
import sys
import multiprocessing

ns3_path_to_waf = '/home/x/ns-allinone-3.33/ns-3.33/'
ns3_left_name = 'nmb1'
ns3_right_name = 'nmb2'
total = 0
now_second = 0
total_round = 0
total_second = 0

server_name = 'a@192.168.56.4'
server_name_root = 'root@192.168.56.4'
client_name = 'a@192.168.57.4'
source_ip = '10.0.0.1'
host_log_path = '/home/a/mylog/'

local_log_path = '/home/x/mptcp-mptcp_v0.95/log/'
local_dir = ''

def cmd_to_host(host_name, cmd):
    print('ssh {} \"{}\"'.format(host_name, cmd))
    os.system('ssh {} \"{}\"'.format(host_name, cmd))
    with open('{}/command.txt'.format(local_dir), 'a' ) as f:
        f.write('ssh {} \"{}\"\n'.format(host_name, cmd))

def cmd_to_host_root(host_name, cmd):
    print('echo a | ssh {} \"{}\"'.format(host_name, cmd))
    os.system('echo a | ssh {} \"{}\"'.format(host_name, cmd))
    with open('{}/command.txt'.format(local_dir), 'a' ) as f:
        f.write('echo a | ssh {} \"{}\"\n'.format(host_name, cmd))

def cmd_to_host_sudo(host_name, cmd):
    print('echo a | ssh -tt {} \"sudo {}\"'.format(host_name, cmd))
    os.system('echo a | ssh -tt {} \"sudo {}\"'.format(host_name, cmd))
    with open('{}/command.txt'.format(local_dir), 'a' ) as f:
        f.write('echo a | ssh -tt {} \"sudo {}\"\n'.format(host_name, cmd))

def cmd_local(cmd):
    print(cmd)
    result = os.system(cmd)
    print("INFO result: {}".format(result))
    return result

def ps_check(name):
    result = os.popen('ps -a | grep {}'.format(name)).readlines()
    if len(result) == 0:
        return -1
    else:
        return result[0].split()[0]

def collect_data():
    cmd_local('scp -r {}:{} {}'.format(server_name, host_dir, local_log_path))
    cmd_local('scp -r {}:{} {}'.format(client_name, host_dir, local_log_path))

def timer(timeout, iperf_t: multiprocessing.Process):
    global now_second, total_second
    now_tmp = now_second
    for i in range(int(timeout) + 20):
        print('Timing: {}/{}+{} [{}s/{}s({:.2f}h): all]'.format(i, timeout, 20, now_second, total_second, float(total_second / 3600)), end='\r')
        time.sleep(1)
        now_second += 1
        if iperf_t.is_alive() == False:
            print('\niperf finished!')
            return False
    if iperf_t.is_alive():
        iperf_t.terminate()
        print('\nKILL iperf')
        now_second = now_tmp
        return True

def restart_vm():
    cmd_local('vboxmanage controlvm left poweroff')
    time.sleep(5)
    cmd_local('vboxmanage startvm left')
    time.sleep(120)
    result = cmd_local('ssh a@192.168.56.4 \"iperf3 -s -D\"')
    print(result)
    while result != 0:
        time.sleep(30)
        result = cmd_local('ssh a@192.168.56.4 \"iperf3 -s -D\"')
        print(result)
    

def go_iperf(scheduler, round, duration, window, ns3_left_argv, ns3_right_argv, host_dir):
    cmd_to_host_sudo(server_name, 'sysctl net.mptcp.mptcp_scheduler={}'.format(scheduler))
    cmd_local('echo a | sudo -S echo > {}ns3_info_left.log 2>&1 &'.format(local_log_path))
    cmd_local('echo a | sudo -S echo > {}ns3_info_right.log 2>&1 &'.format(local_log_path))
    sub_dir = '{}/{}'.format(host_dir, scheduler)
    sub_dir = sub_dir  + '-' + window + duration + '-' + ns3_right_argv.replace(' ', '')
    cmd_to_host(server_name, 'mkdir -p {}'.format(sub_dir))
    cmd_to_host(client_name, 'mkdir -p {}'.format(sub_dir))
    for i in range(round):
        global total, total_round
        total_round += 1
        print('====================================')
        print('\033[1;31;0m {}/{}: {} ({}/{}: all) \033[0m'.format(i+1, round, scheduler, total_round, total))
        print('====================================')
        os.chdir(ns3_path_to_waf)
        cmd_to_host_root(server_name_root, 'echo > /var/log/kern.log')
        cmd_local('echo a | sudo -S ./waf --run \'{} {}\' >> {}ns3_info_left.log 2>&1 &'.format(ns3_left_name, ns3_left_argv, local_log_path))
        cmd_local('echo a | sudo -S ./waf --run \'{} {}\' >> {}ns3_info_right.log 2>&1 &'.format(ns3_right_name, ns3_right_argv, local_log_path))
        for j in range(999):
            left_pid = ps_check(ns3_left_name)
            right_pid = ps_check(ns3_right_name)
            if left_pid == -1 or right_pid == -1:
                time.sleep(0.1)
                continue
            break
        if left_pid == -1 or right_pid == -1:
            print('ns3 scripts GO WRONG!!'.format(ns3_left_name, left_pid))
            if left_pid != -1:
                cmd_local('echo a | sudo kill {}'.format(left_pid))
            if right_pid != -1:
                cmd_local('echo a | sudo kill {}'.format(right_pid))
            i -= 1
            continue
        print('PID {}：{}'.format(ns3_left_name, left_pid))
        print('PID {}：{}'.format(ns3_right_name, right_pid))
        iperf_t = multiprocessing.Process(target=cmd_to_host, args=(client_name, 'iperf3 -c 10.0.0.1 -i 1 -t {} -w {} -J -R --logfile {}/{}iperf.json'.format(duration, window, sub_dir, i),))
        # cmd_to_host(client_name, 'iperf3 -c 10.0.0.1 -i 1 -t {} -w {} -J -R --logfile {}/{}iperf.json'.format(duration, window, sub_dir, i))
        iperf_t.start()
        timeout = timer(duration, iperf_t)
        iperf_t.join()
        if timeout:
            restart_vm()
            cmd_local('echo a | sudo kill {}'.format(left_pid))
            cmd_local('echo a | sudo kill {}'.format(right_pid))
            i -= 1
            continue
        cmd_to_host(server_name, 'cat /var/log/kern.log | grep ytxing >> {}/{}kernlog.log'.format(sub_dir, i))
        cmd_local('echo a | sudo kill {}'.format(left_pid))
        cmd_local('echo a | sudo kill {}'.format(right_pid))
        collect_data()

if __name__ == '__main__':
    window = '30k'
    duration = '200'
    note = sys.argv[1]
    ns3_left_argv = '-loss=0'
    ns3_right_argv = '-owd=30ms'
# ns3_right_argv = '-changetime=200 -change=1'
    experiments = [ 
                    # ['default',     8, '200', '30k',    ns3_left_argv, '-owd=0ms'], \
                    ['default',     3, '200', '100k',   ns3_left_argv, '-owd=0ms'], \
                    ['redundant',   8, '200', '30k',    ns3_left_argv, '-owd=0ms'], \
                    ['redundant',   8, '200', '100k',   ns3_left_argv, '-owd=0ms'], \
                    ['roundrobin',  8, '200', '30k',    ns3_left_argv, '-owd=0ms'], \
                    ['roundrobin',  8, '200', '100k',   ns3_left_argv, '-owd=0ms'], \
                    ['ol',          8, '200', '30k',    ns3_left_argv, '-owd=0ms'], \
                    ['ol',          8, '200', '100k',   ns3_left_argv, '-owd=0ms'], \
                    # =======================================
                    ['default',     8, '200', '30k',    ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['default',     8, '200', '100k',   ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['redundant',   8, '200', '30k',    ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['redundant',   8, '200', '100k',   ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['roundrobin',  8, '200', '30k',    ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['roundrobin',  8, '200', '100k',   ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['ol',          8, '200', '30k',    ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    ['ol',          8, '200', '100k',   ns3_left_argv,  '-owd=20ms -loss=0.2'], \
                    # ========================================
                    ['default',     8, '200', '30k',    ns3_left_argv,    '-loss=0.2'], \
                    ['default',     8, '200', '100k',   ns3_left_argv,    '-loss=0.2'], \
                    ['redundant',   8, '200', '30k',    ns3_left_argv,    '-loss=0.2'], \
                    ['redundant',   8, '200', '100k',   ns3_left_argv,    '-loss=0.2'], \
                    ['roundrobin',  8, '200', '30k',    ns3_left_argv,    '-loss=0.2'], \
                    ['roundrobin',  8, '200', '100k',   ns3_left_argv,    '-loss=0.2'], \
                    ['ol',          8, '200', '30k',    ns3_left_argv,    '-loss=0.2'], \
                    ['ol',          8, '200', '100k',   ns3_left_argv,    '-loss=0.2'], \
                    # ========================================
                    ['default',     8, '200', '30k',    '-loss=0.2',    '-loss=0.2'], \
                    ['default',     8, '200', '100k',   '-loss=0.2',    '-loss=0.2'], \
                    ['redundant',   8, '200', '30k',    '-loss=0.2',    '-loss=0.2'], \
                    ['redundant',   8, '200', '100k',   '-loss=0.2',    '-loss=0.2'], \
                    ['roundrobin',  8, '200', '30k',    '-loss=0.2',    '-loss=0.2'], \
                    ['roundrobin',  8, '200', '100k',   '-loss=0.2',    '-loss=0.2'], \
                    ['ol',          8, '200', '30k',    '-loss=0.2',    '-loss=0.2'], \
                    ['ol',          8, '200', '100k',   '-loss=0.2',    '-loss=0.2'], \
                    # ========================================
                    ['default',     8, '200', '30k',    ns3_left_argv, '-owd=20ms'], \
                    ['default',     8, '200', '100k',   ns3_left_argv, '-owd=20ms'], \
                    ['redundant',   8, '200', '30k',    ns3_left_argv, '-owd=20ms'], \
                    ['redundant',   8, '200', '100k',   ns3_left_argv, '-owd=20ms'], \
                    ['roundrobin',  8, '200', '30k',    ns3_left_argv, '-owd=20ms'], \
                    ['roundrobin',  8, '200', '100k',   ns3_left_argv, '-owd=20ms'], \
                    ['ol',          8, '200', '30k',    ns3_left_argv, '-owd=20ms'], \
                    ['ol',          8, '200', '100k',   ns3_left_argv, '-owd=20ms'], \
                    # ========================================
                    ['default',     8, '200', '30k',    ns3_left_argv, '-owd=40ms'], \
                    ['default',     8, '200', '100k',   ns3_left_argv, '-owd=40ms'], \
                    ['redundant',   8, '200', '30k',    ns3_left_argv, '-owd=40ms'], \
                    ['redundant',   8, '200', '100k',   ns3_left_argv, '-owd=40ms'], \
                    ['roundrobin',  8, '200', '30k',    ns3_left_argv, '-owd=40ms'], \
                    ['roundrobin',  8, '200', '100k',   ns3_left_argv, '-owd=40ms'], \
                    ['ol',          8, '200', '30k',    ns3_left_argv, '-owd=40ms'], \
                    ['ol',          8, '200', '100k',   ns3_left_argv, '-owd=40ms'], \
                    # =========================================
                    ]
    # experiments = [['ol', 1]]
    
    for x in experiments:
        total = (total + x[1]) 
        total_second = total_second + int(x[2]) * x[1]

    now = time.strftime("%Y%m%d-%H%M%S", time.localtime())
    host_dir = host_log_path + now + note
    local_dir = local_log_path + now + note
    print(local_dir)
    cmd_local('mkdir -p {}'.format(local_dir))
    cmd_to_host(client_name, 'mkdir -p {}'.format(host_dir))
    for experiment in experiments:
        go_iperf(experiment[0], experiment[1], experiment[2], experiment[3], experiment[4], experiment[5], host_dir)

    cmd_local('find . -name "*" -type f -size 0c | xargs -n 1 rm -f')
    cmd_local('find -type d -empty | xargs -n 1 rm -rf')

    decodejson.get_dir_json(local_dir)
    getinfo.get_info_dir(local_dir)
    mycsv.get_scheduler_csv(local_dir)
    mycsv.get_all_csv(local_dir)


