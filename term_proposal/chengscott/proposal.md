# nsd PageRank

## Background

PageRank was an algorithm used by Google to rank the importance of resulting web pages from their search engine. PageRank simulates Web surfers, starting at a random page and randomly browsing around, would tend to converge and aggregated, and this process is allowed to iterate many times.

The algorithm accepts a directed graph as the input, where vertices and edges represents web pages and links between, respectively; It outputs PageRank scores for each vertices, which represents the likelihood of a page being visited by random web surfing. The problem can be reformulated in matrix calculations:

Let $d_i$ be the out-degree of vertex $i$. Then, define a column stochastic matrix $M$ as
$$
M_{ji}=
\begin{cases}
\frac{1}{d_i},\ \text{if } i \to j \\
0,\ \text{otherwise.}
\end{cases}
$$
In addition, define Google Matrix as $A = \beta M + (1 - \beta) \left[\frac{1}{N}\right]_{N\times N}$, where $N$ is the number of vertices, and $\beta\in(0,1)$ is a probability to teleport away.

The PageRank equation (matrix version) is:
$$
r = A \cdot r
$$

The calculation requires several passes of iterations to approximate PageRank values.

## Project

In this project, I am going to implement PageRank algorithm with different kinds of BLAS. I expect to get familiar with several BLAS libraries, and compare performance among mine implementations. The code will all be written in modern C++, and be tested with small and large datasets (refer to reference).

The input is a file containing lines of a vertex pair (i.e., an direct edge), which represents a directed graph. The output is a file containing lines of (vertex, PageRank) pair, which is sorted by PageRank number from high to low.

Scenario: An user wants to know the most popular websites, so he describe the website links as an directed edge in a file. By the output file of this program, he will directly know the top-k popular websites by a simple glance.

## Reference

- [PageRank Algorithm](https://en.wikipedia.org/wiki/PageRank)
- [Gnutella Dataset](https://snap.stanford.edu/data/p2p-Gnutella04.html)


This chrome plugin may help to read in GitHub: [MathJax Plugin for Github](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima?hl=en).

