#!/usr/bin/python3
import time
import os
import multiprocessing
import expscript

ns3_path_to_waf = '/home/x/ns-allinone-3.33/ns-3.33/'
ns3_left_name = 'nmb1'
ns3_right_name = 'nmb2'
total_exp = 0
total = 0
now_second = 0
total_round = 0
total_second = 0
crash_cnt = 0

server_name = 'a@192.168.56.4'
server_name_root = 'root@192.168.56.4'
client_name = 'a@192.168.57.4'
source_ip = '10.0.0.1'
host_log_path = '/home/a/mylog/'

local_log_path = '/home/x/mptcp-mptcp_v0.95/log/'
local_dir = local_log_path

def mycmd(cmd: str):
    result = os.system(cmd)
    with open('{}/command.txt'.format(local_dir), 'a' ) as f:
        f.write(cmd+'\n')
    return result

def myprint(*args, **kwargs):
   print('[YTXING] ', end='')
   return print(*args, **kwargs)

def cmd_to_host(host_name, cmd):
    myprint('ssh {} \"{}\"'.format(host_name, cmd))
    mycmd('ssh {} \"{}\"'.format(host_name, cmd))
    # with open('{}/command.txt'.format(local_dir), 'a' ) as f:
    #     f.write('ssh {} \"{}\"\n'.format(host_name, cmd))

def cmd_to_host_root(host_name, cmd):
    myprint('echo a | ssh {} \"{}\"'.format(host_name, cmd))
    mycmd('echo a | ssh {} \"{}\"'.format(host_name, cmd))
    # with open('{}/command.txt'.format(local_dir), 'a' ) as f:
    #     f.write('echo a | ssh {} \"{}\"\n'.format(host_name, cmd))

def cmd_to_host_sudo(host_name, cmd):
    myprint('echo a | ssh -tt {} \"sudo {}\"'.format(host_name, cmd))
    mycmd('echo a | ssh -tt {} \"sudo {}\"'.format(host_name, cmd))
    # with open('{}/command.txt'.format(local_dir), 'a' ) as f:
    #     f.write('echo a | ssh -tt {} \"sudo {}\"\n'.format(host_name, cmd))

def cmd_local(cmd):
    myprint(cmd)
    result = mycmd(cmd)
    myprint("INFO result: {}".format(result))
    return result

def ps_check(name):
    result = os.popen('ps -a | grep {}'.format(name)).readlines()
    if len(result) == 0:
        return -1
    else:
        return result[0].split()[0]

def collect_data(local_dir, sub_dir):
    cmd_local('scp -r {}:{} {}'.format(server_name, sub_dir, local_dir))
    cmd_local('scp -r {}:{} {}'.format(client_name, sub_dir, local_dir))

def timer(timeout, iperf_t: multiprocessing.Process):
    global now_second, total_second, crash_cnt
    now_tmp = now_second
    more_time = 20
    for i in range(int(timeout) + more_time):
        myprint('Timing: {}/{}+{}] [crashes:{}]'.format(i, timeout, more_time, crash_cnt), end='\r')
        time.sleep(1)
        now_second += 1
        if iperf_t.is_alive() == False:
            myprint('\niperf finished!')
            return False
    if iperf_t.is_alive():
        iperf_t.terminate()
        myprint('\nKILL iperf')
        now_second = now_tmp
        return True

def restart_server():
    cmd_local('vboxmanage controlvm left poweroff')
    time.sleep(5)
    cmd_local('vboxmanage startvm left')
    time.sleep(120)
    result = cmd_local('ssh a@192.168.56.4 \"iperf3 -s -D\"')
    myprint(result)
    count = 0
    while result != 0:
        time.sleep(30)
        result = cmd_local('ssh a@192.168.56.4 \"iperf3 -s -D\"')
        myprint(result)
        count += 1
        if count > 3:
            restart_server()
    
def restart_iperfs():
    iperfs_pid = -1
    flag = 0
    results = os.popen('ssh {} \"ps aux | grep iperf3\"'.format(server_name))
    for line in results:
        if line.find('iperf3 -s -D') != -1:
            myprint(line.split()[1])
            iperfs_pid = line.split()[1]
    
    cmd_to_host_sudo(server_name, 'sudo kill {}'.format(iperfs_pid))
    cmd_to_host(server_name, 'iperf3 -s -D')
    results = os.popen('ssh {} \"ps aux | grep iperf3\"'.format(server_name))
    for line in results:
        if line.find('iperf3 -s -D') != -1:
            myprint('Restart iperf3 -s'.format(line.split()[1]))
            flag = 1
            time.sleep(2)

    results = os.popen('ssh {} \"ps aux | grep iperf3\"'.format(client_name))
    for line in results:
        if line.find('iperf3 -c 10.0.0.1') != -1:
            myprint(line.split()[1])
            iperfc_pid = line.split()[1]
            cmd_to_host_sudo(client_name, 'sudo kill {}'.format(iperfc_pid))

    return flag

def stop_iperfs():
    iperfs_pid = -1
    flag = 0
    results = os.popen('ssh {} \"ps aux | grep iperf3\"'.format(server_name))
    for line in results:
        if line.find('iperf3 -s -D') != -1:
            myprint(line.split()[1])
            iperfs_pid = line.split()[1]
    
    cmd_to_host_sudo(server_name, 'sudo kill {}'.format(iperfs_pid))

    results = os.popen('ssh {} \"ps aux | grep iperf3\"'.format(client_name))
    for line in results:
        if line.find('iperf3 -c 10.0.0.1') != -1:
            myprint(line.split()[1])
            iperfc_pid = line.split()[1]
            cmd_to_host_sudo(client_name, 'sudo kill {}'.format(iperfc_pid))

    return flag

def go_iperf(cc, scheduler, round, duration, trunk_size, window, ns3_left_argv, ns3_right_argv, host_dir, local_dir):
    global total, total_round, crash_cnt

    cmd_to_host(server_name, "cat /sys/module/mptcp_ol/parameters/DEBUG_USE_NEW_EPOCH")
    cmd_to_host(server_name, "cat /sys/module/mptcp_ol/parameters/DEBUG_USE_GAMMA_TUNING")
    cmd_to_host_sudo(server_name, 'sysctl net.mptcp.mptcp_scheduler={}'.format(scheduler))
    cmd_local('echo a | sudo -S echo > {}ns3_info_left.log 2>&1 &'.format(local_log_path))
    cmd_local('echo a | sudo -S echo > {}ns3_info_right.log 2>&1 &'.format(local_log_path))
    if trunk_size == '':
        sub_dir = window + '-' + str(duration) + 's-L' + ns3_left_argv.replace(' ', '')+ '-R' + ns3_right_argv.replace(' ', '')
        use_trunk = False
    else:
        sub_dir = window + '-' + str(trunk_size) + '-L' + ns3_left_argv.replace(' ', '')+ '-R' + ns3_right_argv.replace(' ', '')
        use_trunk = True
    sub_dir = sub_dir + '-' + scheduler + '-' + cc 
    
    sub_dir = '{}/{}'.format(host_dir, sub_dir)
    
    cmd_to_host(server_name, 'mkdir -p {}'.format(sub_dir))
    cmd_to_host(client_name, 'mkdir -p {}'.format(sub_dir))

    i = 0
    while i < round:    
        cmd_to_host_sudo(server_name, '/home/a/setup-network.sh')
        cmd_to_host_sudo(client_name, '/home/a/setup-network.sh')
        time.sleep(0.3)

        cmd_to_host_sudo(server_name, 'sysctl net.mptcp.mptcp_scheduler={}'.format(scheduler))
        if cc == 'bbr':
            cmd_to_host_sudo(server_name, 'sysctl net.ipv4.tcp_congestion_control={}'.format(cc))
            cmd_to_host_sudo(server_name, 'sysctl net.core.default_qdisc=fq')
        else:
            cmd_to_host_sudo(server_name, 'sysctl net.ipv4.tcp_congestion_control={}'.format(cc))
            cmd_to_host_sudo(server_name, 'sysctl net.core.default_qdisc=pfifo_fast')

        left_pid = ps_check(ns3_left_name)
        right_pid = ps_check(ns3_right_name)
        if left_pid != -1:
            cmd_local('echo a | sudo kill {}'.format(left_pid))
        if right_pid != -1:
            cmd_local('echo a | sudo kill {}'.format(right_pid))
        count = 0
        while restart_iperfs() != 1:
            count += 1
            time.sleep(3)
            if count > 10:
                restart_server()
                crash_cnt += 1
                continue

        total_round += 1
        myprint('====================================')
        myprint('\033[1;31;0m {}/{}: {} ({}/{} + {}: all) \033[0m'.format(i+1, round, scheduler, total_round, total, crash_cnt))
        myprint('====================================')
        os.chdir(ns3_path_to_waf)
        cmd_to_host_root(server_name_root, 'echo > /var/log/kern.log')
        cmd_local('echo a | sudo -S ./waf --run \'{} {}\' >> {}ns3_info_left.log 2>&1 &'.format(ns3_left_name, ns3_left_argv, local_log_path))
        cmd_local('echo a | sudo -S ./waf --run \'{} {}\' >> {}ns3_info_right.log 2>&1 &'.format(ns3_right_name, ns3_right_argv, local_log_path))
        for j in range(200):
            left_pid = ps_check(ns3_left_name)
            right_pid = ps_check(ns3_right_name)
            if left_pid == -1 or right_pid == -1:
                time.sleep(0.1)
                continue
            break
        if left_pid == -1 or right_pid == -1:
            myprint('ns3 scripts GO WRONG!!'.format(ns3_left_name, left_pid))
            if left_pid != -1:
                cmd_local('echo a | sudo kill {}'.format(left_pid))
            if right_pid != -1:
                cmd_local('echo a | sudo kill {}'.format(right_pid))
            continue
        myprint('PID {}: {}'.format(ns3_left_name, left_pid))
        myprint('PID {}: {}'.format(ns3_right_name, right_pid))
        cmd_to_host(client_name, 'touch {}/{}iperf.json'.format(sub_dir, i))
        cmd_to_host(client_name, 'echo > {}/{}iperf.json'.format(sub_dir, i))
        if use_trunk:
            iperf_t = multiprocessing.Process(target=cmd_to_host, args=(client_name, 'iperf3 -c 10.0.0.1 -i 1 -n {} -w {} -J -R --logfile {}/{}iperf.json'.format(trunk_size, window, sub_dir, i),))
        else:
            iperf_t = multiprocessing.Process(target=cmd_to_host, args=(client_name, 'iperf3 -c 10.0.0.1 -i 1 -t {} -w {} -J -R --logfile {}/{}iperf.json'.format(duration, window, sub_dir, i),))
        iperf_t.start()
        wrong = timer(duration, iperf_t)
        iperf_t.join()
        if wrong:
            restart_server()
            crash_cnt += 1
        cmd_to_host(server_name, 'cat /var/log/kern.log | grep ytxing >> {}/{}kernlog.log'.format(sub_dir, i))
        cmd_local('echo a | sudo kill {}'.format(left_pid))
        cmd_local('echo a | sudo kill {}'.format(right_pid))
        collect_data(local_dir, sub_dir)
        i += 1

def run_exp(experiments, note):
    global total, total_second, total_round, total_exp, local_dir
    total_exp += 1
    total = 0
    total_second = 0
    total_round = 0
    note = '-' + note
    # ns3_right_argv = '-changetime=200 -change_owd=1'
    # experiments = [['ol', 1]]
    for x in experiments:
        total = (total + x[2]) 
        total_second = total_second + int(x[3]) * x[2]


    
    # cmd_to_host_sudo(server_name, '/home/a/setup-network.sh')

    now = time.strftime("%Y%m%d-%H%M%S", time.localtime())
    host_dir = host_log_path + now + note
    local_dir = local_log_path + now + note
    os.system('mkdir -p {}'.format(local_dir))
    os.system('touch {}/command.txt'.format(local_dir))
    os.system('chmod 777 {}/command.txt'.format(local_dir))
    cmd_to_host(client_name, 'mkdir -p {}'.format(host_dir))
    for experiment in experiments:
        # cc, scheduler, round, duration, trunk_size, window, ns3_left_argv, ns3_right_argv, host_dir, local_dir
        go_iperf(experiment[0], experiment[1], experiment[2], experiment[3], experiment[4], experiment[5], experiment[6], experiment[7], host_dir, local_dir)
        os.system('echo \"{}\" >> {}/done_exps.txt'.format(experiment, local_dir))

    cmd_local('find . -name "*" -type f -size 0c | xargs -n 1 rm -f')
    cmd_local('find -type d -empty | xargs -n 1 rm -rf')

    # decodejson.get_dir_json(local_dir)
    # getinfo.get_info_dir(local_dir)
    # mycsv.get_scheduler_csv(local_dir)
    # mycsv.get_all_csv(local_dir)

def main():
    global ns3_left_name, ns3_right_name

    ns3_left_name = 'p2pnmb1'
    ns3_right_name = 'p2pnmb2'
    note = "WSP-gimme_20211230_103655_1125"
    experiments = expscript.gimme_20211230_103655_1125()
    run_exp(experiments, note)

if __name__ == '__main__':
    main()