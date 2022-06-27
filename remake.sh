#!/bin/sh
make clean
sleep 0.5
make
sleep 0.5
rmmod mptcp_ol
sleep 0.5
insmod mptcp_ol.ko