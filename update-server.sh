#!/bin/sh
scp mptcp_ol.c a@192.168.56.4:/home/a/mptcp-mptcp_v0.95/mptcp-olsched
ssh a@192.168.56.4 -p 22  "ls -la /home/a/mptcp-mptcp_v0.95/mptcp-olsched | grep mptcp_ol.c"
ssh a@192.168.56.4 -p 22  "echo a | sudo -S ./remake.sh"