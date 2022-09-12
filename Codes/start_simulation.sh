#!/bin/bash
g++ -Wall -o ising ising.cpp
g++ -Wall -o bootstrap bootstrap.cpp
for i in {1..150}
do
   ./ising
   ./bootstrap
done
