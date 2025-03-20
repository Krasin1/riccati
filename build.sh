#!/bin/sh

cmake -S . -B ./build 
cmake --build ./build --parallel $(nproc)
cp build/solver ./
