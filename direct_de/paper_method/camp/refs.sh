#!/bin/bash
cd "$(dirname "$0")/.."
for n in eucl500 scat500 eucl2000 scat2000; do nice -n 19 python3 camp/refs.py $n; done
