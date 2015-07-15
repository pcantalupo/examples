#!/bin/bash

printf "#### FROM BEGINNING ####\n"
VAR=/path/to/filename.txt
printf "Absolute path: $VAR\n"

# shortest ${variable#pattern}
printf "\nShortest (i.e. relative path):\n"
printf "     \${VAR#*/} -> ${VAR#*/}\n"   # path/to/filename.txt

# longest ${variable##pattern}
printf "\nLongest (i.e. basename):\n"
printf "     \${VAR##*/} -> ${VAR##*/}\n" # filename.txt

printf "\n\n"


printf "#### FROM END ####\n"
VAR=foobar.qux.fa
printf "Filename: $VAR\n"

# shortest ${variable#pattern}
printf "\nShortest (i.e. remove suffix):\n"
printf "     \${VAR%%.*} -> ${VAR%.*}\n"    # foobar.qux

# longest ${variable##pattern}
printf "\nLongest (i.e. keep prefix):\n"
printf "     \${VAR%%%%.*} -> ${VAR%%.*}\n"  # foobar




