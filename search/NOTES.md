# Černý conjecture: search in progress

Work in progress; no result is claimed.

- `src/rt.h`: reset threshold of a DFA with n <= 28 states and k <= 8 letters, by
  breadth-first search over subsets.
- `src/test_cerny.c`: checks `rt.h` on the Černý automata C_n, n <= 20 (rt = (n-1)^2).
- `src/sa.c`: simulated annealing over k-letter, n-state automata maximising the reset
  threshold; prints any automaton with rt > (n-1)^2 as COUNTEREXAMPLE.
  Usage: `sa n k iterations seed [T0]`.
- `runs/`: logs of runs at (n, k) = (13, 3), (14, 2), (16, 2), (18, 2). No counterexample;
  best values 83 / 144, 169 / 169 (the Černý automaton), 127 / 225, 149 / 289.

Build: `cc -O2 -o sa src/sa.c -lm`.
