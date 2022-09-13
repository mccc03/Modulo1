#!/bin/bash
g++ -Wall -o ising ising.cpp
for i in {1..150}
do
   ./ising
   python bootstrap.py
done
