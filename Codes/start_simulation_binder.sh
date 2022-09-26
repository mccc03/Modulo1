#!/bin/bash
g++ -Wall -o ising_binder ising_binder.cpp
for i in {1..150}
do
   ./ising_binder
   python bootstrap_binder.py
done
