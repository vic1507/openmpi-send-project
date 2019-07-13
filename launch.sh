#!/bin/bash
mpicc -o sand-cellular-automata finalProjectScalingVersion.c -lm $(pkg-config allegro-5 allegro_font-5 allegro_primitives-5 --libs --cflags) | mpirun -np 4 ./sand-cellular-automata

