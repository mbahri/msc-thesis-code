#!/bin/bash

while read -r line; do 
    HOSTNAME="$(echo $line | cut -d ',' -f 1)"
    MU="$(echo $line | cut -d ',' -f 2)"
    LAMBDA="$(echo $line | cut -d ',' -f 3)"
    ID="$(echo $line | cut -d ',' -f 4)"

    echo "Running mu $MU and lambda $LAMBDA with id $ID on host $HOSTNAME"
    ssh -oStrictHostKeyChecking=no "$HOSTNAME" "/homes/mb2215/bitbucket/Thesis/bilinear-rpca/hypersearch/wrapper.sh $MU $LAMBDA $ID" &
done < hosts
