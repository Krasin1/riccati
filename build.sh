#!/bin/sh

cmake -S src/ -B ./.build
cmake --build ./.build --parallel $(nproc)
cp .build/solver ./
ln -sf ../.build/compile_commands.json src/
