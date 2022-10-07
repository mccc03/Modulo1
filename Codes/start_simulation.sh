#!/bin/zsh
g++ -Wall -o ising ising.cpp
for i in {1..200}
do
   ./ising
   python bootstrap.py
done
