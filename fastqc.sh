# run fastqc
# Find all .gz files and run fastqc

find . -name "*.gz" | xargs -n 1 fastqc &
