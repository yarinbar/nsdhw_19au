nsd Term Proposal
===
GitHub@stvreumi

## Motivation
- When I interned in NCHC last summer, I've been asked to use some tools like [hyperopt](https://github.com/hyperopt/hyperopt) to tune the hyperparameters of deep learning models. I know how to use the tools to run the experiments to produce possible better hyperparameters, but I don't know what's the math behind these tools.
- **Bayesian Optimization**, is the keyword showing up in almost every optimization tool that I've used or interested in. I want to understand what Bayesian Optimization is by research and implementation in the final project.

## Status
- Actually, I have no idea of Bayesian Optimization. And I've forgotten most things in the probability course.
- But I have found sources called "[Think in Bayes](http://www.greenteapress.com/thinkbayes/thinkbayes.pdf)". The goal of this book is trying to teach readers the concepts of Bayesian Optimization through Python code.
  (Although I only finish preface when I writing this proposal)

## Plan

- My goal of the term project is trying to finish the materials, and implement some core functions of Bayesian Optimization in C++.
- I want to bind the core functions with `pybind11`, and I will write a Python program to call it to generate data and do some analysis on it.
- I will try to package my code to a python package called `simpleBayes`. There is a function in it called `maximize`, which can take search space , objective function and iterate numbers as parameters. The `maxmize` function will try to maximize the value of the objective function. The code may looks like:
  ```python
  from simpleBayes import maximize
  
  # objective function
  def f(x,y):
      return -((x-2)*(x-2) + (y+3)*(y+3))
  
  
  # search space is a `dict` in the form: 
  # parameter name of objective function: (lower bound, upper bound)
  # the type of lower bound and upper bound is `double`
  
  search_space = {
      'x': (-50.0, 50.0),
      'y': (-50.0, 50.0)
  }
  
  """
  try to maximize the value of objective function
  parameters:
      search_space: search space definition. 
      if the search space doesn't meet the spec, it will throw an error.
      f: the objective function we want to optimize
      num: the iterate number of of the experiment
  
  return:
      best: the best value of each parameter(`dict`)
      history: the history of the experiment (suggested parameters, objective function)(`list`)
  """
  best, history = maximize(search_space=search_space, f=f, num=50)
  ```
- I will use some simple funcions like `f` in the above code to test if it really can find values close to maxmimal. 
- After this course:
  - I may try to research on how to find better results in fewer iterations.
  - I could research on the related GitHub repos in the future. If possible, I may try to contribute some PRs to their repos.

## Resources
- [Think in Bayes](http://www.greenteapress.com/thinkbayes/thinkbayes.pdf)
  - [python2 repo](https://github.com/AllenDowney/ThinkBayes)
  - [python3 repo](https://github.com/AllenDowney/ThinkBayes2)
- [Automated Machine Learning Hyperparameter Tuning in Python](https://towardsdatascience.com/automated-machine-learning-hyperparameter-tuning-in-python-dfda59b72f8a)
- [Bayesian optimization, Martin Krasser](http://krasserm.github.io/2018/03/21/bayesian-optimization/)