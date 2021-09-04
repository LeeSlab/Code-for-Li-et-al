
## generate .cool file
cooler cload pairs -c1 2 -p1 3 -c2 5 -p2 6 --temp-dir ./tmp/ sample_100000_abs.change.bed sample_allValidPairs sample_allValidPairs.cool
cooler balance sample_allValidPairs.cool

## generate .cool file
## then you can analysis the cool files using coolpuppy.py from https://github.com/Phlya/coolpuppy
