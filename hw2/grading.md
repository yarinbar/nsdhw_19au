# HW2 Grading Standard

HW2: https://github.com/yungyuc/nsdhw_19au/tree/master/hw2

1. Your submission PR needs to create a sub-directory using
exactly your GitHub username as your working directory (nsdhw_19au/hw2/$username).
2. The PR needs to be created against the hw2 branch. (Not `master` branch!)
3. If duplication is found, you may deduct **all points** for the duplicated part of
the later submitter.

**You need to submit your PR link on e3.**

## Grading Detail

There are 5 points in this homework.

Working directory hierarchy may look like:
```
- nsdhw_19au
  - hw2
    - username
      - q1
        - Makefile
        - xxx.cc
        - (any other you need)
      - q2
        - Makefile
        - xxx.cc
        - yyy.py
        - (any other you need)
```

Please make sure your PR is compatible with Unix.

### Q1 has 2 points

- Implementation gets 1 point. Do not modify `main()`, TA will `diff` the output to the correct answer.
- Makefile gets 1 point:
  * When a source file changes (you can touch it), `make` needs to pick it up
    and rebuild.
  * `make run` needs to produce the correct terminal output, without
    crashing.
  * `make clean` needs to remove all the built and intermediate files.
  * To get full point one has to do all 3.

TA will read the code and also check:
```
cd q1
make
make run > output
diff outpout TA_answer
make clean
```

### Q2 for 3 points

- Implementation gets 1 point. 
- Makefile gets 1 point.
  * When a source file changes (you can touch it), `make` needs to pick it up
    and rebuild.
  * `make test` runs the Python-based tests.
  * `make clean` removes all the built and intermediate files.
  * To get full point one has to do all 3.
- Unit-testing gets 1 point.
  * a test function for zero-length 2-vector (invalid input).
  * a test function for zero angle.
  * a test function for right angle (90-deg).
  * a test function for one other angle.
  * To get full point, do at least 2 of the above.

TA will read the code and also check:
```
cd q2

make test
make clean
``
