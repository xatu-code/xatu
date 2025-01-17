#!/bin/sh

rm libxatu.a log*

make build
echo '--------------------------------------'
echo '|              FINISHED BUILD         |'
echo '--------------------------------------'
sleep 1
make xatu
echo '--------------------------------------'
echo '|              FINISHED XATU          |'
echo '--------------------------------------'

echo 'Done!'
