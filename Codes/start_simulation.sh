#!/bin/zsh
g++ -Wall -o ising ising.cpp
for i in {1..160}
do
   ./ising
   python bootstrap.py
done
