#!/bin/bash

mv MoS2_ex.dat plot/.
mv MoS2_sp.dat plot/.

cd plot
sleep 1
python3 conductivity.py
# python3 absorption.py
sleep 1
feh last.png
