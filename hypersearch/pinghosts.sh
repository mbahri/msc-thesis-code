#!/bin/bash

while read -r line; do 
    HOSTNAME="$(echo $line | cut -d ',' -f 1)"
    

    ping "$HOSTNAME" -c 1
done < hosts
