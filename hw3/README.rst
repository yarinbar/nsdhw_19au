======================
Homework assignment #3
======================

Your submission PR needs to create a sub-directory using exactly your GitHub
username as your working directory (``nsdhw_19au/hw3/$username``).  **You also
need to submit your PR link on E3.** The PR needs to be created against the
``hw3`` branch.  (Not ``master``!)

Working directory hierarchy may look like::

  - (nsdhw_19au)
    - (hw3)
      - username (<- your workind directory)
        - Makefile
        - xxx.{cc|cpp|cxx}
        - test_xxx.py (for unit testing)
        - performance.txt
        - (any other files you need)

Everyone should write his/her own code.  It is OK to discuss, but there should
not be duplicated code.  If duplication is found, **all points** for the
duplicated part of the latter submitter may be deducted.

Please make sure your PR is compatible with the latest NSD AMI on AWS (Ubuntu
18.04 LTS).  You are not required to use the AMI for doing the homework, but
the grader will use the AMI.  If your code doesn't build or properly run, you
**may not get any point**.

Question
========

Develop your own matrix-matrix multiplication code, measure the runtime, and
compare with that of BLAS ``DGEMM`` subroutine.  The matrix size should be
larger than or equal to :math:`1000\times1000`.

Grading guideline
=================

There are totally 5 points in this homework assignment.

1. Develop a two-dimensional matrix class (you may reuse the code in the course
   note) in C++.  (2 points)

   1. Implement a free function performing naive matrix-matrix multiplication.
   2. Implement another free frunction that uses ``DGEMM`` for the
      multiplication.
   3. The ``validate.py`` script needs to pass.
2. Make Python unit tests for the C++ matrix class and multiplication
   functions.  (1 point)
3. Make sure the answers from the two multiplication functions match.  (1
   point)
4. Time the two multiplication functions and compare the performance.  (1
   point)

On the latest NSD AMI, the grader will run the following commands:

.. code-block:: bash

  cd nsdhw_19au/hw3
  cd username
  ../validate.py
  make clean
  cp -f performance.txt performance.txt.orig
  make bench
  diff performance.txt.orig performance.txt

The grader may attempt to run the above commands with the environment variable
``PRELOAD_MKL`` set, if ``validate.py`` doesn't pass without it.  Your work is
considered pass in either case.

``make bench`` should reproduce the file ``performance.txt``.  The reproduced
performance numbers may differ.  Grader will check if they are reasonable.

``performance.txt`` may looks like:

.. code-block:: log

  Start multiply_naive (repeat=5), take min = 1.2351381540065631 seconds
  Start multiply_mkl (repeat=5), take min = 0.056176671001594514 seconds
  MKL speed-up over naive: 21.9867 x


.. vim: set ft=rst ff=unix fenc=utf8 et sw=2 ts=2 sts=2:
