#!/bin/bash
#script to clean files from a mercury simulation
# version 1.0

function clean {
    rm *.out
    rm *.dmp
    rm *.tmp
    
    # To delete the stderr and stdout of a bash scheduler of the server. 
    # The last "." is very important, in order to avoid suppression of the submission script itself.
    rm *.sh.* 
}

function clean-plots {
    rm *.aei
    rm *.clo
}

#echo "are you sure to erase all files from the simulation? (y/n)"
#read confirm
#if [ $confirm = "y" -o $confirm = "o" ]
#then
    #echo "deleting files..."
    #clean
    #echo "done"
#else
    #echo "nothing was suppressed."
#fi


echo "deleting files..."
clean
clean-plots
echo "done"


