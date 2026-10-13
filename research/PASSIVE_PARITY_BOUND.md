# Passive parity measurements: access and sampling bounds

Exactly three spin-½ bosons occupy three distinct input modes. Their joint
spin state is arbitrary, and the detector reports only total occupation
parity at each output. A passive device applies the same single-particle
unitary to both spin components. This account proves a three-mode
identifiability obstruction and a four-mode sampling lower bound. The fourth
input mode is vacuum.

The four-mode bound permits an arbitrary product of single-spin unitaries
**before** the passive device. Settings and these prerotations may depend on
the classical record from earlier independent copies. Additional particles
or spatial modes, spin-resolved detection, intermediate measurements within
a copy, interleaved site-dependent spin controls, nonunitary input filters,
and quantum memory between copies are outside this class. Any sequence of
spin-independent passive layers without intervening measurements is one
allowed spatial unitary. Every unknown input preparation counts.

Use the [repository convention](../THEORY.md)
$\chi=\mathbf S_1\cdot(\mathbf S_2\times\mathbf S_3)$, with
$\mathbf S=\boldsymbol\sigma/2$, and set

$$
h=\frac{\sqrt3}{4},\qquad
\Pi_{\rm m}=\Pi_++\Pi_-,\qquad
\operatorname{Tr}\chi^2=\frac34.
$$

The principal four-mode conclusion is

$$
\boxed{\sup_\rho n\operatorname{Var}_\rho(\widehat\chi)
\ge\frac{13}{32}}
$$

for every fixed-$n$, universally unbiased estimator in the declared class.
This is a lower bound, not the exact passive optimum. Section 6 gives a
finite-sample mean-square-error bound that also permits biased estimators.
The [finite-repulsion certificate](CERTIFICATE.md) supplies strictly smaller
all-input upper bounds, including a specified noisy channel with exact
response calibration. These are input-copy comparisons, not apparatus-time
or calibration-cost comparisons.

## 1. Compressing a spatial parity observable

Let $U$ be the spatial unitary, $S$ a subset of output modes, and

$$
W_S=U^\dagger\operatorname{diag}\bigl((-1)^{1_{j\in S}}\bigr)U.
$$

The product of the output parities in $S$ is the second quantization of
$W_S$. Compress it to one particle in each of the three occupied input
modes, and let $K_S$ be the corresponding $3\times3$ principal submatrix.
The resulting spin observable is

$$
B_S=\sum_{\pi\in S_3}
\left(\prod_{i=1}^3(K_S)_{i,\pi(i)}\right)P_\pi,
\qquad
P_\pi|s_1s_2s_3\rangle
=|s_{\pi(1)}s_{\pi(2)}s_{\pi(3)}\rangle.
\tag{1}
$$

To prove (1), expand the three transformed creation operators. A matrix
element between singly occupied input states must match the three outgoing
input positions to the three incoming positions bijectively. Each matching
is a permutation, contributes the displayed product of matrix entries, and
permutes the spin labels accordingly. Bosonic contractions have positive
sign. There is no occupation factorial at these distinct input positions.

For Hermitian $K$, the identity coefficient and the three transposition
coefficients are real. The two cycle coefficients are $z$ and $\bar z$,
where $z=K_{12}K_{23}K_{31}$. Since
$C-C^\dagger=4i\chi$, their chirality coefficient is

$$
c(B_S)=-4\operatorname{Im}(K_{12}K_{23}K_{31}).
\tag{2}
$$

Reversing the cycle convention reverses both signs and leaves all bounds
unchanged. For three qubits,
$C+C^\dagger=P_{12}+P_{23}+P_{13}-I$, so every $B_S$ belongs to
$\operatorname{span}_{\mathbb R}\{I,P_{12},P_{23},P_{13},\chi\}$.

## 2. Three modes cannot identify chirality

For a single output, $W=I-2vv^\dagger$. Its off-diagonal triangle product is
$-8|v_1v_2v_3|^2$, which is real. Equation (2) therefore gives zero
chirality coefficient. With three output modes and exactly three particles,
the product of all three parities is $-1$. A two-output parity is the
negative of the complementary single-output parity. All parity effects
therefore lie in

$$
\operatorname{span}_{\mathbb R}\{I,P_{12},P_{23},P_{13}\}.
\tag{3}
$$

Product-unitary prerotations do not remove this obstruction. An identity or
swap has Pauli weight at most two, whereas $\chi$ has Pauli weight three;
product unitaries preserve each Pauli support. Thus every rotated effect
in (3) remains Hilbert–Schmidt orthogonal to $\chi$.

The states $\rho_\pm=\Pi_\pm/2$ have chirality $\pm h$ and identical
outcome laws for every such setting. Their complete adaptive data laws
are also identical, by conditioning successively on the preceding record.
For any estimator $T$, including biased and nonlinear estimators,

$$
\max\{\mathbb E_+(T-h)^2,\mathbb E_-(T+h)^2\}
\ge\mathbb ET^2+h^2\ge\frac3{16}.
$$

Reporting zero has error at most $h^2$ on every input. Hence the exact
three-mode minimax mean-square error is $3/16$, independently of the number
of copies. This statement concerns parity detection; it does not apply to
particle-number-resolved Fourier readout discussed in
[the attribution account](../REFERENCES.md).

For example, the balanced three-mode Fourier device has

$$
E_{111}=\frac13\Pi_{\rm s},\qquad
E_{100}=E_{010}=E_{001}=\frac13I-\frac19\Pi_{\rm s}.
$$

Both chiral sectors consequently give probabilities $(0,1/3,1/3,1/3)$.

## 3. Four-mode chirality coefficients lie in a bounded cube

For two selected outputs, $W_S$ is a Hermitian unitary with two eigenvalues
$+1$ and two eigenvalues $-1$. The three-dimensional occupied-input
subspace intersects each of these two-dimensional eigenspaces. Its
compression $K_S$ therefore has eigenvalues $1,t,-1$, with $|t|\le1$.

For a Hermitian matrix $K=V\operatorname{diag}(\lambda_1,\lambda_2,
\lambda_3)V^\dagger$,

$$
\left|\operatorname{Im}(K_{12}K_{23}K_{31})\right|
=|J(V)|\prod_{i<j}|\lambda_i-\lambda_j|,
\quad
J(V)=\operatorname{Im}(V_{11}V_{22}\bar V_{12}\bar V_{21}).
\tag{4}
$$

One direct proof expands the left side as a cubic in the eigenvalues. It
vanishes whenever two eigenvalues coincide, because the remaining matrix
is a scalar identity plus a rank-one projector. It is therefore a constant
times the degree-three Vandermonde product. Expanding its
$\lambda_1^2\lambda_2$ coefficient identifies that constant as $J(V)$,
up to the ordering sign removed in (4).

The bound $|J(V)|\le1/(6\sqrt3)$ is elementary. Column phases make the
first row $(a,b,c)$ real and nonnegative. Write the second row as
$\alpha x+\beta y$, where $x,y$ are a real orthonormal basis of its
orthogonal plane and $|\alpha|^2+|\beta|^2=1$. Then

$$
|J(V)|=ab\,|\operatorname{Im}(\beta\bar\alpha)|
\,|x_1y_2-x_2y_1|
\le\frac{abc}{2}\le\frac1{6\sqrt3}.
$$

Here $|x_1y_2-x_2y_1|=c$ follows from $x\times y=\pm(a,b,c)$;
the last step is arithmetic–geometric mean with $a^2+b^2+c^2=1$.
Equations (2) and (4) now imply

$$
|c(B_S)|\le8(1-t^2)|J(V)|
\le B:=\frac4{3\sqrt3}.
\tag{5}
$$

Write an observed parity as signs $z_j=(-1)^{p_j}$, so
$\prod_{j=0}^3z_j=-1$. Fourier inversion on these eight outcomes uses
the identity, four single-output parities, and three pair parities
$B_{01},B_{02},B_{03}$. Let their chirality coefficients be $b_1,b_2,b_3$.
Single-output coefficients vanish, and therefore

$$
c_z=\frac18(r_1b_1+r_2b_2+r_3b_3),\qquad
r_j=z_0z_j,\qquad r_1r_2r_3=-1,
\quad b\in[-B,B]^3.
\tag{6}
$$

Each of the four possible $r$ occurs twice. Convexity in $b$ reduces the
maximum of $\sum_z|c_z|$ to cube vertices. There its eight terms have
magnitudes $3B/8$ twice and $B/8$ six times. Consequently

$$
h\sum_z|c_z|\le\frac12.
\tag{7}
$$

For a product-unitary prerotation $L=L_1\otimes L_2\otimes L_3$,
the rotated coefficients are $c'_z=g(L)c_z$, where

$$
g(L)=\frac43\operatorname{Tr}(\chi L^\dagger\chi L),
\qquad |g(L)|\le1.
$$

The identity and swap terms contribute zero by their Pauli support;
Hilbert–Schmidt Cauchy–Schwarz proves the bound on $g$. Thus the rotated
coefficients still have (6), with $b$ replaced by $g(L)b$. The rotated
effects need not preserve $\Pi_{\rm m}$; the following positivity argument
does not assume that they do.

## 4. A normalized Fisher-information bound

Use the valid submodel

$$
\rho_\theta=\frac{\Pi_{\rm m}}4+\frac43\theta\chi,
\qquad |\theta|<h,\qquad
\operatorname{Tr}(\rho_\theta\chi)=\theta.
$$

For either unrotated or prerotated effects write
$p_z(\theta)=q_z+\theta c_z$. Positivity at the endpoints gives
$q_z\ge h|c_z|$; completeness gives $\sum_zq_z=1$ and
$\sum_zc_z=0$. An outcome with $q_z=0$ also has $c_z=0$ and is removed.

At zero, the Fisher information is $F(0)=\sum_zc_z^2/q_z$.
For fixed $c$, this is convex in $q$. Its maximum under the relaxed
constraints $q_z\ge h|c_z|$ and $\sum q_z=1$ therefore assigns every
probability its minimum, except one outcome $j$ receiving the entire
remainder. The resulting value is

$$
G_j(b)=\frac{A_j(b)}h+
\frac{c_j(b)^2}{1-hA_j(b)},\qquad
A_j(b)=\sum_{i\ne j}|c_i(b)|.
$$

The function $H(A,x)=A/h+x^2/(1-hA)$ is jointly convex, by the
quadratic-over-linear perspective, and increasing in $A$.
Thus $G_j$ is convex in $b$, because $A_j$ is convex and $c_j$ is affine.
Its denominator is positive by (7). The maximum over the cube is attained
at a cube vertex.

At a vertex the minimum probabilities are $1/8$ twice and $1/24$ six
times, leaving remainder $1/2$. Assigning it to a small-coefficient outcome
gives

$$
F(0)\le\frac83-\left(\frac29-\frac2{117}\right)
=\boxed{\frac{32}{13}}.
\tag{8}
$$

Assigning it to a large-coefficient outcome gives only $32/15$.
The cube and positivity constraints are relaxations of physical devices;
(8) does not assert that a passive unitary attains equality.

The same argument gives a useful pointwise bound. At nonzero $\theta$,
the extreme-probability value is

$$
H_\theta(A,x)=\frac{hA+\theta x}{h^2-\theta^2}
+\frac{x^2}{1-hA+\theta x}.
$$

It remains jointly convex and increasing in $A$, with positive denominator.
At cube vertices, the largest value places the extra probability on a
small coefficient with the same sign as $\theta$. With $t=|\theta|/h$,

$$
\boxed{F(\theta)\le I(t):=
\frac{16(6+t)}{3(1-t^2)(13+t)}}.
\tag{9}
$$

Indeed the minimum-probability sum is $8/[3(1-t^2)]$, and assigning the
remainder subtracts $8/[3(1+t)(13+t)]$. This also proves (9) at zero.
$I(t)$ increases on $0\le t<1$: both $(6+t)/(13+t)$ and
$1/(1-t^2)$ are nonnegative and increasing.

## 5. Adaptation and unbiased estimators

Condition on each preceding classical record. The next spatial unitary and
product prerotation define an allowed measurement, so (8)–(9) still apply.
Conditional likelihood scores have conditional mean zero. Their cross
terms vanish in the total Fisher information, giving

$$
F_n(0)\le\frac{32n}{13},\qquad
F_n(\theta)\le nI(|\theta|/h).
$$

The likelihood argument remains regular for arbitrary adaptive setting
kernels. At each fixed preceding record, setting selection is independent
of the unknown parameter, and a measurement factor is $q_z+\theta c_z$.
If $q_z=0$, then $c_z=0$ and that outcome vanishes throughout the submodel.
Otherwise, on any compact interval $|\theta|\le a<h$,

$$
\left|\partial_\theta\log(q_z+\theta c_z)\right|
=\frac{|c_z|}{q_z+\theta c_z}
\le\frac1{h-a}.
$$

Thus the full fixed-$n$ score is bounded by $n/(h-a)$, and likelihood
ratios between two parameters in this interval are bounded by
$\exp(n|\theta-\theta_0|/(h-a))$. These bounds justify differentiating
estimator expectations whenever their second moment is finite; an infinite
variance already satisfies the lower bound. Randomized settings may have
continuous distributions: their parameter-independent conditional kernels
add no score, and the same bounds apply to the joint observed record.

Any universally unbiased estimator $T$ satisfies
$\mathbb E_\theta T=\theta$ on this submodel. At zero,
$\mathbb E_0[T\,\partial_\theta\log p_\theta]=1$.
Cauchy–Schwarz yields

$$
\operatorname{Var}_{\rho_0}(T)\ge\frac1{F_n(0)}
\ge\frac{13}{32n}.
$$

This covers arbitrary unbiased processing, not merely a sample mean or a
fixed score. The lower-bound state is $\rho_0=\Pi_{\rm m}/4$, not a
promise known to the estimator. State-independent reference preparations
carry no Fisher information about $\theta$; they may select later settings,
which the adaptive bound already permits. Discarding outcomes cannot
increase information per original preparation. The statement assumes a
fixed number of unknown copies; it does not assert the same formula for
arbitrary stopping rules.

## 6. Finite mean-square error, including biased estimators

For $0<a<h$, introduce the normalized auxiliary prior

$$
\pi_a(\theta)=\frac1a\cos^2\!\left(\frac{\pi\theta}{2a}\right),
\qquad -a\le\theta\le a.
$$

It vanishes at the endpoints and has Fisher information $\pi^2/a^2$.
For the joint density of parameter and data, integration by parts gives
$\mathbb E[(T-\theta)\partial_\theta\log(p_\theta\pi_a)]=1$.
The data score has conditional mean zero, so its cross term with the prior
score vanishes. Cauchy–Schwarz therefore gives, without unbiasedness,

$$
\boxed{\sup_\rho\mathbb E_\rho(T-\langle\chi\rangle)^2
\ge\frac1{nI(a/h)+\pi^2/a^2}}.
\tag{10}
$$

This is the scalar van Trees argument, reproduced here; see [GL95] in
[References](../REFERENCES.md). It requires uniform risk over arbitrary
inputs, not risk at one known state. The fixed-sample likelihood bounds
above and the finite prior information justify the differentiation and
integration by parts for uniformly finite-risk estimators.

In particular, let $R_n$ be the infimum of this worst-case MSE over all
allowed fixed-$n$ protocols and estimators, including biased estimators and
adaptive settings. For $n>1$, choose $a_n=h n^{-1/4}$ in (10). Then

$$
nR_n\ge
\left[I(n^{-1/4})+\frac{\pi^2}{h^2\sqrt n}\right]^{-1},
\qquad
\boxed{\liminf_{n\to\infty}nR_n\ge\frac{13}{32}}.
$$

Here the prior shrinks slowly enough that its information is $o(n)$,
while $I(n^{-1/4})\to I(0)=32/13$. This establishes the asymptotic
minimax lower bound without an unbiasedness or attainability assumption.

For any $0<t_0<1$, a uniform MSE requirement $\varepsilon^2$ implies

$$
n\ge\frac{\varepsilon^{-2}-16\pi^2/(3t_0^2)}{I(t_0)}
\tag{11}
$$

when the numerator is positive, followed by an integer ceiling. At
$\varepsilon=0.001$ and $t_0=1/10$,
$I(t_0)=97600/38907$ and (11) requires **396539 unknown copies**.
This is a necessary count, not an attained passive optimum.

For a comparison valid at a common fixed sample size, choose $t_0=1/5$.
Then $I(t_0)=775/297$ and the prior information is $400\pi^2/3$.
An exactly calibrated native estimator with coefficient below $0.367$
has smaller worst-case MSE than the lower bound (10) whenever

$$
n>\frac{0.367(4000/3)}{1-0.367(775/297)}
=\frac{5813280}{503}<12000,
$$

where $\pi^2<10$ was used conservatively. Thus **every integer
$n\ge12000$** gives a strict comparison even against biased passive
estimators. The [noisy certificate](CERTIFICATE.md) supplies the required
native coefficient under its exact-channel assumptions. Finite empirical
calibration must additionally satisfy the
[calibration and cost conditions](FINITE_CALIBRATION.md); it is not removed
by this target-copy calculation.

## 7. An explicit four-mode construction

The extra vacuum mode makes chirality identifiable without interactions.
For $0<x<1$, put $v=\sqrt x$ and $r=\sqrt{(1-x)/2}$, and use columns
$(0,1,2,3)$ with column 0 vacuum:

$$
U(x)=\begin{pmatrix}
v/\sqrt2&r/\sqrt2&1/\sqrt2&r/\sqrt2\\
r&(-v+i)/2&0&(-v-i)/2\\
r&(-v-i)/2&0&(-v+i)/2\\
v/\sqrt2&r/\sqrt2&-1/\sqrt2&r/\sqrt2
\end{pmatrix}.
\tag{12}
$$

Its rows are orthonormal by direct multiplication. Equation (1) gives four
nonzero chirality coefficients of magnitude $\sqrt x(1-x)/2$, in two
complementary pairs. The two singleton outcomes have
$q_{\rm m}=(1-x)(x+3)/16$; the two three-odd-site outcomes have
$q_{\rm m}=(1-x)(3x+1)/16$, where
$q_{\rm m}=\operatorname{Tr}(\Pi_{\rm m}E)/4$. Hence

$$
F(x)=\sum_p\frac{c_p^2}{q_{{\rm m},p}}
=\frac{32x(1-x^2)}{(x+3)(3x+1)}.
$$

Choose $a_p=c_p/[q_{{\rm m},p}F(x)]$. Uniformly permute the three occupied
input modes before the device, and multiply the reported score by the
permutation sign. The signed permutation average of its mean operator
retains only $\chi$, with coefficient $\sum_pa_pc_p=1$. The ordinary
permutation average of the squared-score operator is

$$
T=t_{\rm s}\Pi_{\rm s}+t_{\rm m}\Pi_{\rm m},\quad
t_{\rm m}=\frac{(x+3)(3x+1)}{32x(1-x)(x+1)},\quad
t_{\rm s}=\frac{(3x+1)^2}{16(1-x)(x+1)^2}.
$$

These formulas follow by substituting the four coefficients into (1);
the singleton symmetric-sector probabilities are $x(1-x)/2$, and the two
informative three-odd-site probabilities vanish in that sector. Moreover
$t_{\rm s}/t_{\rm m}=2x(3x+1)/[(x+1)(x+3)]\le1$.
The state $\Pi_{\rm m}/4$ has zero chirality and attains second moment
$t_{\rm m}$, so this estimator's exact worst-case variance coefficient is
$1/F(x)$.

For the simple algebraic choice $x=2/5$, the coefficient is $935/1344$.
The optimum within (12) occurs at the unique root in $(0,1)$ of
$3x^4+20x^3+12x^2-3=0$; uniqueness follows from its strictly positive
derivative for $x>0$. Its coefficient is approximately
$0.69513538981$. Thus a proved, non-tight bracket for the full passive
unbiased minimax coefficient is

$$
\frac{13}{32}\le V_{\rm passive}
\le\frac{935}{1344}.
$$

The improved decimal comes from the explicitly stated one-parameter family;
it is not a global optimality or interval-certification claim.

## 8. The same uniform ring without interactions

For $u=0$, let $x=2\tau\cos(\phi/4)$ and
$y=-2\tau\sin(\phi/4)$. Fourier diagonalization gives the circulant
single-particle propagator with first row

$$
\left(a,\frac{iX+Y}{2},d,\frac{iX-Y}{2}\right),\quad
a=\frac{\cos x+\cos y}{2},\quad d=\frac{\cos x-\cos y}{2},
\quad X=\sin x,\quad Y=\sin y.
$$

Substitution into (1) yields

$$
c_p=(-1)^{p_0+p_2}\kappa,\qquad
\kappa=\frac{\sin(2x)\sin(2y)[\cos(2x)-\cos(2y)]}{8}.
\tag{13}
$$

Equivalently, the only pair-parity chirality responses are $8\kappa$ for
$\{0,2\}$ and $-8\kappa$ for its complement. To bound (13), put
$m=[\cos(2x)+\cos(2y)]/2$ and
$d_0=|\cos(2x)-\cos(2y)|/2$. Then

$$
8|\kappa|=2d_0
\sqrt{[1-(m+d_0)^2][1-(m-d_0)^2]}
\le2d_0(1-d_0^2)\le\frac4{3\sqrt3}.
$$

The first inequality follows by expanding the squared radicand as
$(1-d_0^2)^2-2m^2(1+d_0^2)+m^4$ and using
$|m|\le1-d_0$. Thus $|\kappa|\le1/(6\sqrt3)$.
Product prerotations multiply all eight coefficients by the same factor
of magnitude at most one. Normalization and $q_p\ge h|c_p|$ then imply

$$
F(0)\le\frac{7|\kappa|}{h}
+\frac{\kappa^2}{1-7h|\kappa|}
\le\frac{80}{51}.
$$

The right side increases with $|\kappa|$ in the displayed range. Arbitrary
adaptive choices of flux and duration therefore still obey

$$
\boxed{\sup_\rho n\operatorname{Var}_\rho(\widehat\chi)
\ge\frac{51}{80}}
$$

for universally unbiased estimation in this uniform noninteracting-ring
class. This is a lower bound, not a claimed attained ring optimum. At the
finite-repulsion protocol's nominal settings $\tau=\pi$,
$\phi=\pm2\pi/3$, one has $y=\mp\pi$, so $\kappa=0$ exactly: that
noninteracting setting is chirality-blind, including with the allowed
prerotations.

## 9. Extra vacuum modes change the comparison

The four-mode restriction is essential. A passive construction with more
vacuum inputs can recover part of the number-resolved Fourier statistic
using only terminal total-site parity. The underlying Fourier spectral
measurement is established in [N26](../REFERENCES.md), including arbitrary
correlated internal states. The splitter construction below gives an
explicit precision reference outside the four-mode class.

Put $\omega=e^{2\pi i/3}$ and use the three-mode Fourier matrix
$F_{jk}=\omega^{-jk}/\sqrt3$, with $j,k\in\{0,1,2\}$. For occupation
$\mathbf n$ after this interferometer, let

$$
f(\mathbf n)=\sum_{j=0}^2 jn_j\pmod3,\qquad
\lambda_f=\frac12\sin\!\left(\frac{2\pi f}{3}\right)
\in\{0,h,-h\}.
$$

The count statistic is a spectral measurement of $C$ with eigenvalue
$\omega^f$, and hence of $\chi$ with eigenvalue $\lambda_f$, in the
repository's cycle convention. This follows directly by applying a cyclic
input shift before $F$: it becomes the output phase $\omega^f$.
Equivalently, if $E_{\mathbf n}$ is the spin effect for resolved counts,
then $E_{\mathbf n}$ is supported on the corresponding cyclic eigenspace.
This is an operator statement; it does not assume product spin inputs.

Assign the three Fourier outputs uniformly at random to three balanced
splitter groups with $L_0,L_1,L_2\ge1$ arms. Three cyclic assignments
suffice. Each splitter has one populated input and vacuum in its other
inputs, so the entire apparatus is a passive unitary on
$M=L_0+L_1+L_2$ spatial modes with $M-3$ vacuum inputs. Record the assignment.

Known particle number three makes **three odd output parities** an exact
certificate of one particle in each of three distinct arms. On such an
event, summing the odd bits within each group reconstructs the original
Fourier occupation and its $f$. Nonzero $f$ occurs only for permutations
of $(2,1,0)$. If the doubly occupied Fourier output enters an $L$-arm
splitter, its particles reach distinct arms with probability $1-1/L$.
This probability is independent of their internal spin state: both
particles occupy the same input spatial orbital. Distinct group counts
have disjoint output occupations, so there are no cross terms between
different Fourier count sectors. Averaging the recorded random assignment
therefore retains each informative count effect with the same factor

$$
s=1-\frac13\left(\frac1{L_0}+\frac1{L_1}+\frac1{L_2}\right).
$$

For $s>0$, report $\lambda_f/s$ on three-odd events and zero on every
other event. Every trial counts; there is no conditioning of the sample
size on a detection event. The averaged first and second score operators
are exactly

$$
A=\chi,\qquad
T=\frac{\chi^2}{s}=\frac{h^2}{s}\Pi_{\rm m}.
\tag{14}
$$

These identities hold for arbitrary correlated spin inputs and remain
valid after tensoring with an external reference. Since
$\rho_0=\Pi_{\rm m}/4$ has zero chirality, it attains the upper bound on
the single-copy variance. Thus the exact worst-case variance coefficient
of this estimator is $h^2/s$.

| Total modes | Splitter arms $(L_0,L_1,L_2)$ | Informative survival $s$ | Attained worst-case coefficient |
|---:|---|---:|---:|
| 6 | $(2,2,2)$ | $1/2$ | $3/8=0.375$ |
| 7 | $(3,2,2)$ | $5/9$ | $27/80=0.3375$ |
| 8 | $(3,3,2)$ | $11/18$ | $27/88\approx0.306818$ |
| 9 | $(3,3,3)$ | $2/3$ | $9/32=0.28125$ |

The eight-mode construction has a smaller coefficient than the present
four-mode interacting score, approximately $0.328698$. It uses the same
kind of terminal parity detector but twice as many spatial modes. This
makes the resource boundary explicit: the interacting advantage established
here is at fixed four-mode and detector resources. The table neither
identifies the passive optimum at any listed mode count nor proves that
eight modes are necessary to outperform the interacting score.

## Verification

[check_passive_parity.py](check_passive_parity.py) independently compares
the permutation formula with bosonic creation-operator amplitudes, checks
three-mode blindness, the convex-bound arithmetic and local-rotation
identity, the explicit four-mode score, the finite-risk constants, and the
uniform-ring formula against the full bosonic model. A fifth
group checks (14) for all four splitter constructions using independent
creation-operator amplitudes in the complete spinful three-particle Fock
space, together with exact rational coefficients. Its finite examples
cross-check the analytical proofs; they do not prove a global optimum by
numerical search. Reports are written only to new paths outside the
repository.
