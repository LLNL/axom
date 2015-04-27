#!/bin/bash

for file in $(find -E . -type f -regex ".*\.(hpp|cpp)" -follow -print0 | xargs -0); do
    ~/Codes/uncrustify/src/uncrustify -c uncrustify.cfg --no-backup $file
    echo $file
done

    
#echo $FILES
