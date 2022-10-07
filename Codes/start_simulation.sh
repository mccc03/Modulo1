#!/bin/zsh
g++ -Wall -o ising ising.cpp
for i in {1..100}
do
   ./ising
   python bootstrap.py
done
