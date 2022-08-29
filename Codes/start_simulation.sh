#!/bin/bash
for i in {1..150}
do
   ./ising
   ./bootstrap
   python ising_data_analysis.py
done
