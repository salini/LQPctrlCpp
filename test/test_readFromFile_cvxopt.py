#!/usr/bin/python
#author=joseph salini
#date=30 april 2012

from numpy import array, dot

from cvxopt import solvers, matrix
from cvxopt.solvers import qp as qpsolver

f = open("/home/joe/myProblem.txt", "r")
problemtext = f.read()
f.close()

print problemtext

Problist = problemtext.split()

print Problist

nProblem = int(Problist[Problist.index('nProblem:')+1])
nEquality = int(Problist[Problist.index('nEquality:')+1])
nInequality = int(Problist[Problist.index('nInequality:')+1])

print nProblem, nEquality, nInequality

P = array([float(v) for v in Problist[Problist.index('P:')+1 : Problist.index('P:')+ nProblem*nProblem +1]]).reshape(nProblem,nProblem)
q = array([float(v) for v in Problist[Problist.index('q:')+1 : Problist.index('q:')+ nProblem +1]]).reshape(nProblem)

A = array([float(v) for v in Problist[Problist.index('A:')+1 : Problist.index('A:')+ nEquality*nProblem +1]]).reshape(nEquality,nProblem)
b = array([float(v) for v in Problist[Problist.index('b:')+1 : Problist.index('b:')+ nEquality +1]]).reshape(nEquality)

G = array([float(v) for v in Problist[Problist.index('G:')+1 : Problist.index('G:')+ nInequality*nProblem +1]]).reshape(nInequality,nProblem)
h = array([float(v) for v in Problist[Problist.index('h:')+1 : Problist.index('h:')+ nInequality +1]]).reshape(nInequality)

print "P:\n", P
print "q:\n", q
print "A:\n", A
print "b:\n", b
print "G:\n", G
print "h:\n", h

Pp = matrix(P)
qp = matrix(q)
Gp = matrix(G)
hp = matrix(h)
Ap = matrix(A)
bp = matrix(b)

solvers.options.update({'abstol': 1e-12, 'reltol': 1e-12, 'feastol': 1e-12,})

results = qpsolver(Pp, qp, Gp, hp, Ap, bp)

X_sol = array(results['x']).flatten()

print "cost:", .5*dot(X_sol, dot(P, X_sol)) + dot(q, X_sol)
print X_sol
