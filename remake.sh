#!/bin/sh
make
rmmod mptcp_ol
insmod mptcp_ol.ko
sudo sysctl net.mptcp.mptcp_scheduler=ol
sudo sysctl net.mptcp.mptcp_path_manager=fullmesh