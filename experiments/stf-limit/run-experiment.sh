time ./../../multiedsm -m 4096m -s ../maw/eds/Y.eds -p ../maw/patterns/filtered_hg37.maws -t 1 > results.txt
python test.py ../maw/patterns/filtered_hg37.maws ../maw/eds/Y.eds 10 >> results.txt
