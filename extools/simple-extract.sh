#!/bin/bash

export LANG=C
date
./extractor -i $1 -d X -c 400000 | sort -t $'\t' -k 1 | gzip > ex.output.gz
date
gzcat ex.output.gz | ./mr_stripe_rule_reduce -p | sort -t $'\t' -k 1 | ./mr_stripe_rule_reduce | gzip > phrase-table.gz

