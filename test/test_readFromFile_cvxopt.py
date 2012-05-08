#!/usr/bin/python
#author=joseph salini
#date=30 april 2012

from numpy import array, dot, zeros

from cvxopt import solvers, matrix
from cvxopt.solvers import qp as qpsolver

import sys

_fileName = sys.argv[1]
f = open(_fileName, "r")
problemtext = f.read()
f.close()

print "==============> "+_fileName+" <=============="

Problist = problemtext.split()

nProblem = int(Problist[Problist.index('nProblem:')+1])
nEquality = int(Problist[Problist.index('nEquality:')+1])
nInequality = int(Problist[Problist.index('nInequality:')+1])

P = zeros((nProblem, nProblem))
q = zeros(nProblem)
if "P:" in Problist:
    P = array([float(v) for v in Problist[Problist.index('P:')+1 : Problist.index('P:')+ nProblem*nProblem +1]]).reshape(nProblem,nProblem)
if "q:" in Problist:
    q = array([float(v) for v in Problist[Problist.index('q:')+1 : Problist.index('q:')+ nProblem +1]]).reshape(nProblem)

A = zeros((nEquality, nProblem))
b = zeros(nEquality)
if "A:" in Problist:
    A = array([float(v) for v in Problist[Problist.index('A:')+1 : Problist.index('A:')+ nEquality*nProblem +1]]).reshape(nEquality,nProblem)
if "b:" in Problist:
    b = array([float(v) for v in Problist[Problist.index('b:')+1 : Problist.index('b:')+ nEquality +1]]).reshape(nEquality)

G = zeros((nInequality, nProblem))
h = zeros(nInequality)
if "G:" in Problist:
    G = array([float(v) for v in Problist[Problist.index('G:')+1 : Problist.index('G:')+ nInequality*nProblem +1]]).reshape(nInequality,nProblem)
if "h:" in Problist:
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


