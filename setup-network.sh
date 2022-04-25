#!/bin/sh
ssh a@192.168.56.4 -p 22  "echo a | sudo -S ./setup-network.sh"