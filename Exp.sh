#!/bin/sh
for var in Karate ;do
    julia -O3 Main.jl $var 1
    julia -O3 Main.jl $var 2
    julia -O3 Main.jl $var 3
done
for var in Karate ;do
    julia -O3 TestTime.jl $var 1
    julia -O3 TestTime.jl $var 2
    julia -O3 TestTime.jl $var 3
done
