# [Semidefinite Programming](@id semidefinite-programming)

Semidefinite Programming (SDP) is a powerful optimization technique that can be
used to solve a wide variety of problems in science and engineering, including
many problems in quantum physics [vandenberghe1996semidefinite](@cite). It is a
generalization of linear programming, where instead of optimizing a linear
function over a set of linear inequalities, we optimize a linear function over a
set of "semidefinite" constraints.

A semidefinite constraint is a condition on a matrix to be "positive
semidefinite". A matrix is positive semidefinite if it is symmetric and all its
eigenvalues are non-negative. This condition is a generalization of the concept
of a non-negative number to matrices.

SDP is particularly useful for problems involving sums of squares of
polynomials, like the one we discussed in the previous section. This is because
the condition that a polynomial is a sum of squares can be expressed as a
semidefinite constraint. This allows us to use powerful numerical algorithms to
solve these problems.

# [Moment Sum-of-Hermitian-Square Hierarchy](@id moment-sohs-hierarchy)

The "Moment Sum-of-Hermitian-Square Hierarchy" is a powerful technique for
approximating the solution to polynomial optimization problems
[lasserre2010Moments](@cite). The solution is approximated by a sequence of
easier problems. Each problem in the sequence is a semidefinite program, which
can be solved efficiently.

The basic idea is to relax the original problem by considering "moments" of the
variables. The moments are related to the average values of the powers of the
variables. In quantum physics, these moments are related to the expectation
values of the powers of the operators. The constraints on the moments can be
expressed as semidefinite constraints. By adding more and more moments, we get a
tighter and tighter relaxation of the original problem, and the solution to the
relaxed problem gets closer and closer to the solution of the original problem.

This technique has been applied to many problems in quantum information theory,
such as calculating the ground state energy of a quantum system or determining
whether a quantum state is entangled.

## [Moment Problem](@id moment-problem)

In physics, the "moment problem" is a classic mathematical puzzle that has
surprising connections to quantum mechanics. Imagine you have a physical system,
like a particle, and you can measure certain properties of that system, like its
position or momentum. The results of these measurements will have a certain
statistical distribution. The "moments" of this distribution are the average
values of the powers of the measured quantity. For example, the first moment is
the average value, the second moment is related to the variance (how spread out
the values are), and so on.

The moment problem asks the following question: if you know all the moments of a
distribution, can you uniquely determine the distribution itself? In other
words, if you know all the average values of the powers of a physical quantity,
do you know everything there is to know about the probability of measuring a
certain value?

It turns out that the answer is not always yes. There are cases where different
distributions can have the same set of moments. This has important implications
in quantum mechanics, where the state of a system is described by a wave
function, which is related to a probability distribution. The moment problem
helps us understand what information is needed to fully characterize a quantum
state.

## [Sum-of-Hermitian-Square Problem](@id sohs-problem)

The "Sum-of-Hermitian-Square Problem" is a mathematical problem that arises in
quantum information theory and quantum computing. In this context, we are often
interested in optimizing certain quantities that are described by polynomials in
non-commuting variables. These variables represent physical observables that do
not commute with each other, like position and momentum in quantum mechanics.

A key question is to determine whether a given polynomial is always
non-negative. One way to prove that a polynomial is non-negative is to show that
it can be written as a sum of squares of other polynomials. In the case of
non-commuting variables, we consider sums of "Hermitian squares". A Hermitian
square is a product of a polynomial and its "Hermitian conjugate", which is a
generalization of the complex conjugate for matrices.

The paper [Sums of hermitian squares and the BMV
conjecture](https://arxiv.org/abs/0710.1074) by Igor Klep and Markus
Schweighofer, explores the connection between this problem and the
Bessis-Moussa-Villani (BMV) conjecture, a famous open problem in mathematical
physics. The authors use techniques from semidefinite programming to make
progress on this conjecture.
