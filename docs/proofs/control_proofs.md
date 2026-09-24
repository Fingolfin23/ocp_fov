# From Visibility Regions to Certifiable Trajectory and Braking Constraints

This note formalizes the geometric output, Frenet planning bands, segment-based braking distance, and historical blind-region mechanism in Chapter 8 of the report. It provides independent mathematical derivations; it does not claim that the current implementation or experiments have verified every assumption. The original report and source code are not modified by this note.

The input to this chapter is a **correct visibility region or a certified inner approximation** supplied by the geometric module. The preceding chapter establishes how that set is obtained from events on a common ray and the common road stations on the two boundaries. The results below determine when a correct geometric input also yields valid planning constraints.

## C.1 Notation, Model, and Scope of Certification

Let $\Omega\subset\mathbb R^2$ denote the road region and $q$ the observation point. Define the purely geometric visibility region by

$$
V(q)=\{x\in\Omega:[q,x]\subset\Omega\}.
\tag{C.1}
$$

This is a two-dimensional model in which the road boundary occludes visibility and the road interior contains no additional occluders. Independent occluders, when present, must be removed from the region through which a line of sight may pass. Sensor orientation, field-of-view angle, and range restrictions must also be imposed by intersecting the set on the right-hand side with the corresponding sensor coverage. The convention on whether a line of sight may touch an occluding boundary must agree with the preceding chapter. For robust guarantees, positive-margin erosion should replace critical tangency.

**Geometric visibility does not imply that a region has been confirmed obstacle-free.** To make safety claims, introduce a set $K_k\subseteq V(q_k)$ that observations at time $t_k$ certify as obstacle-free. If the environment contains only the known road boundaries and no other obstacles, the model permits $K_k=V(q_k)$. In an environment with unknown obstacles, Eq. (C.1) alone does not certify that all of $V(q_k)$ is free.

For a finite forward query window, distinguish $K_k^{\rm front}\subseteq V(q_k)$ from the known-free set used for full-vehicle constraints,
$K_k^{\rm known}=K_k^{\rm front}\cup K_k^{\rm near}$.
Here, $K_k^{\rm near}$ is a separately certified free neighborhood around and behind the current vehicle whose certificate remains valid; it may extend beyond the current forward window. In the full-vehicle erosions below, $K_k$ is to be understood as $K_k^{\rm known}$. For a point vehicle considered only within the window, one may take $K_k=K_k^{\rm front}$. An artificial entrance cap must not be treated as a collision wall, and free space must not be enlarged without certification.

Let $p$ be the vehicle reference point and $B(\psi)$ its occupied shape relative to that point. The vehicle occupies $p+B(\psi)$. If a circular envelope $B_\rho=\{b:\|b\|\le\rho\}$ satisfies $B(\psi)\subseteq B_\rho$, define a certified reference-point set by

$$
C_k\subseteq K_k\ominus B_\rho
=\{p:p+B_\rho\subseteq K_k\}.
\tag{C.2}
$$

If $q$ lies at the entrance of a finite forward window and the vehicle has nonzero size, $K_k$ must include a separately certified free neighborhood around the vehicle body and behind it before applying the erosion in Eq. (C.2). Eroding a window that contains only the region ahead would exclude the vehicle at its entrance. Artificial window caps must not be treated as physical obstacles.

Thus, $p\in C_k$ is sufficient for the entire vehicle to lie in known free space. It is generally not necessary, because the circular envelope and the inner approximation $C_k$ may be conservative. Set $\rho=0$ for a point vehicle. If $K_k$ is nonconvex, one may first select a convex subset and then erode that subset.

Assume that the reference curve $\gamma:I\to\mathbb R^2$ is a regular $C^2$ curve parameterized by arc length, with unit tangent $T(s)$, left unit normal $N(s)$, and curvature $\kappa(s)$. The Frenet map is

$$
\Phi(s,n)=\gamma(s)+nN(s),\qquad
-w_R(s)\le n\le w_L(s).
\tag{C.3}
$$

Throughout this chapter, the map is assumed to be **globally injective** on the road strip in use, with $1-\kappa(s)n>0$. The latter condition ensures only local nondegeneracy; it does not by itself imply global injectivity. Hairpin bends are allowed, but distinct road cross-sections must not overlap within the coordinate domain and receive ambiguous stations.

## C.2 Blind-Region Removal, Convex Decomposition, and Frenet Intervals

### Theorem C.1: From a Certified Partition to Planning Constraints

Suppose the geometric algorithm returns a candidate region $D_k\subseteq\Omega$ and blind regions $U_{k,j}$. If

$$
F_k=D_k\setminus\bigcup_j U_{k,j}\subseteq V(q_k),
\tag{C.4}
$$

has been proved, and the part of $F_k$ used for planning has been certified free, then for any $C_k\subseteq F_k\ominus B_\rho$, the constraint $p\in C_k$ ensures that the entire vehicle lies in the known-free region certified by that observation.

**Proof.** By the definition of erosion, $p\in C_k$ implies $p+B_\rho\subseteq F_k$. Since $B(\psi)\subseteq B_\rho$, it follows that $p+B(\psi)\subseteq F_k$. Combining this inclusion with Eq. (C.4) and the free-space certificate proves the claim. $\square$

The certification condition in Eq. (C.4) is essential. For example, choosing $D_k$ as the intersection of the road with one side of an ITP–FPV line requires the preceding chapter to establish that the half-space and the blind-region list account for every invisible part. A single line does not automatically have this property on an arbitrary nonconvex road.

**Convexity.** Even if $D_k$ and every $U_{k,j}$ are convex, $F_k$ is generally nonconvex. Removing blind regions does not, by itself, establish convexity. The following two constructions are valid:

1. For a polygonal $F_k$, compute a verified convex decomposition $F_k=\bigcup_{m=1}^{M_k}C_{k,m}$, where each $C_{k,m}$ is a convex polygon. Select one component, or several components with verified connections, for planning.
2. Without requiring two-dimensional convexity, compute the actual feasible set on each Frenet cross-section and retain one connected interval of that set.

A finite polygonal region can be decomposed into convex cells using a boundary arrangement or triangulation. Whether boundary points are feasible must follow the convention used in Eq. (C.1); positive-margin erosion avoids critical boundary cases. Selecting one convex component is conservative. Only the union of all components can potentially retain the full feasible set.

### Theorem C.2: Exact Cross-Sectional Representation

For any certified reference-point set $C$, define at a fixed station $s$

$$
E_C(s)=\{n\in[-w_R(s),w_L(s)]:\Phi(s,n)\in C\}.
\tag{C.5}
$$

The reference point $\Phi(s,n)$ lies in $C$ if and only if $n\in E_C(s)$. If $C$ is a finite polygonal set, $E_C(s)$ is a finite union of intervals and possibly isolated points. If $C$ is convex, $E_C(s)$ is a single interval, a singleton, or the empty set.

**Proof.** The equivalence follows directly from the definition of the preimage. At fixed $s$, the map $n\mapsto\Phi(s,n)$ is affine. The preimage of a convex set under an affine map is convex, and convex subsets of $\mathbb R$ are intervals. The boundary of a finite polygonal set has only finitely many types of intersection events with the cross-section line. Membership along the line can change only at those events, yielding a finite union of intervals. $\square$

Equation (C.5) is a **preimage or intersection at a fixed cross-section**. Orthogonally projecting an entire two-dimensional set onto a normal axis may mix points from different stations and cannot replace Eq. (C.5).

If $E_C(s)=[a,b]\cup[c,d]$ with $b<c$, imposing $a\le n\le d$ would incorrectly fill the infeasible gap. A safe single-band representation must select one component, such as $[a,b]$, or use a discrete component-selection variable to represent the interval union.

### Theorem C.3: Closed-Form Planning Bands for Convex Polygons

Let a certified convex set be defined by half-spaces:

$$
C=\bigcap_{i=1}^{m}\{x:a_i^\top x\le b_i\}.
\tag{C.6}
$$

At a given station $s$, define

$$
c_i(s)=a_i^\top N(s),\qquad
d_i(s)=b_i-a_i^\top\gamma(s).
$$

If every constraint with $c_i(s)=0$ satisfies $d_i(s)\ge0$, let

$$
L(s)=\max\left\{-w_R(s),\max_{c_i(s)<0}\frac{d_i(s)}{c_i(s)}\right\},
\tag{C.7}
$$

$$
U(s)=\min\left\{w_L(s),\min_{c_i(s)>0}\frac{d_i(s)}{c_i(s)}\right\}.
\tag{C.8}
$$

An empty inner maximum is $-\infty$, and an empty inner minimum is $+\infty$. The cross-sectional feasible set is $[L(s),U(s)]$ if and only if $L(s)\le U(s)$. If any constraint has $c_i=0$ and $d_i<0$, or if $L>U$, the cross-section is infeasible.

**Proof.** Substituting Eq. (C.3) into each half-space gives $c_i n\le d_i$. According to whether $c_i$ is positive, negative, or zero, this yields an upper bound, a lower bound, or a feasibility condition. Intersecting these constraints with the road-width interval gives the stated formulas. $\square$

This is an exact realization of Eqs. (8.2) and (8.8) of the report within the selected convex component. **An empty set cannot be repaired by averaging the lower and upper bounds.** Setting both bounds to their mean does not establish that the resulting point satisfies the original constraints.

### Corollary C.3a: Direct Incorporation of Vehicle Shape into Half-Spaces

For a circular vehicle envelope, $p+B_\rho\subseteq C$ if and only if

$$
a_i^\top p\le b_i-\rho\|a_i\|\quad(1\le i\le m).
\tag{C.9}
$$

Indeed, $\max_{b\in B_\rho}a_i^\top b=\rho\|a_i\|$. For the exact vehicle shape, replace $\rho\|a_i\|$ by the support function $\max_{b\in B(\psi)}a_i^\top b$. The latter generally makes the constraints depend on heading, so they need not remain constant upper and lower bounds.

## C.3 Position-Dependent Forward Distance

### Definition C.4: First-Exit Distance in Certified Free Space

Let $p\in C$, and let $h$ be a unit forward direction. Define

$$
d_C(p,h)=\sup\{d\ge0:p+\lambda h\in C\ \text{for all}\ \lambda\in[0,d]\}.
\tag{C.10}
$$

This measures the forward distance over which feasibility is maintained **continuously from the current position**. A ray may leave and later reenter a nonconvex set. A visible fragment beyond the first exit, or a farther intersection, does not extend the distance in Eq. (C.10).

### Theorem C.4: First-Exit Distance from a Convex Safe Set

If $C$ is given by Eq. (C.6) and $p\in C$, then

$$
d_C(p,h)=\min_{i:a_i^\top h>0}
\frac{b_i-a_i^\top p}{a_i^\top h},
\tag{C.11}
$$

with value $+\infty$ if there is no positive denominator. A finite sensor range or road length, when applicable, must be imposed as an additional constraint.

**Proof.** Along the ray, the $i$th constraint becomes $a_i^\top p+\lambda a_i^\top h\le b_i$. If $a_i^\top h\le0$, feasibility at the starting point prevents this constraint from being the first one violated in the positive direction. If the denominator is positive, the corresponding ratio is the largest admissible $\lambda$ for that constraint. The endpoint of the largest interval on which all constraints hold is their minimum. $\square$

For **fixed $s$, fixed $h$, and a fixed convex set $C$**, let $p(n)=\gamma(s)+nN(s)$. Then

$$
d_C(n)=\min_{i:a_i^\top h>0}
\left(\frac{b_i-a_i^\top\gamma(s)}{a_i^\top h}
-n\frac{a_i^\top N(s)}{a_i^\top h}\right).
\tag{C.12}
$$

Thus, $d_C(n)$ is the lower envelope of affine functions and is piecewise affine and concave. This provides an explicit, computable, certified realization of the report's $d_{\mathrm{front}}(n)$.

### Theorem C.5: Ray Intersection with a Single ITP–FPV Segment

Let $A$ and $B$ be the endpoints of the frontier segment, set $e=B-A$, and write the two-dimensional cross product as $[u,v]=u_xv_y-u_yv_x$. If $[h,e]\ne0$, the equation

$$
p+t h=A+u e
$$

has the unique solution

$$
t=\frac{[A-p,e]}{[h,e]},\qquad
u=\frac{[A-p,h]}{[h,e]}.
\tag{C.13}
$$

The forward ray intersects the segment if and only if $t\ge0$ and $0\le u\le1$. Because $h$ is a unit vector, $t$ is the Euclidean travel distance.

**Proof.** Taking the cross product of $t h-u e=A-p$ with $e$ and with $h$ gives Eq. (C.13). The nonzero denominator ensures uniqueness. The two parameter conditions enforce the forward ray and the finite segment, respectively. $\square$

With $A,B,h,s$ fixed, $t(n)$ is affine in $n$. However, it is the exact $d_{\mathrm{front}}(n)$ only if all of the following hold:

- The ray intersects the finite segment, rather than merely its supporting line.
- Every point $p(n)+\lambda h$, for $0\le\lambda\le t(n)$, lies in the certified-free reference-point set.
- The frontier is the first exit boundary along the ray.

If a road side or an Unknown Zone is encountered earlier, that earlier event must be used. If only $[0,t]\subseteq C$ is established, without immediate exit thereafter, $t$ is a lower bound on the certified safe distance and still suffices for conservative braking.

If $[h,e]=0$, the ray and segment are parallel. They have no intersection if they are not collinear; if they are collinear, their one-dimensional overlap must be treated directly instead of using Eq. (C.13). If the direction $h=h(n,\xi)$ varies, or if moving the viewpoint requires recomputing $A,B$, Eq. (C.13) remains applicable state by state, but neither affine dependence nor persistence of the same event branch follows automatically.

**Viewpoint relocation.** When the viewpoint $q=q(n)$ depends on an optimization variable, keeping the frontier fixed does not automatically establish that it still represents the instantaneous visibility region from the new viewpoint. Two rigorous interpretations are available:

1. Recompute and certify $C(q(n))$ for every candidate viewpoint; or
2. Interpret the fixed $C$ as previously observed free space whose validity persists, and establish its validity period without claiming that it equals the current visibility region from the new viewpoint.

## C.4 Braking Distance: Necessity and Sufficiency

### Theorem C.6: Exact Braking Condition for an Ideal Longitudinal Model

Suppose the vehicle travels along a certified-free path with available path length $D\ge0$. Consider the longitudinal model

$$
\dot\ell=v,\qquad \dot v=-a,\qquad
0\le a\le a_*\ (a_*>0),\qquad v\ge0.
\tag{C.14}
$$

Ignore actuation delay and allow the vehicle to maintain its target speed, or remain at rest, after reaching it. Given $v_0\ge0$ and a target upper bound $v_{\mathrm{safe}}\ge0$, reducing speed to at most $v_{\mathrm{safe}}$ within distance $D$ is possible if and only if

$$
v_0^2\le v_{\mathrm{safe}}^2+2a_*D.
\tag{C.15}
$$

**Proof.** If $v_0\le v_{\mathrm{safe}}$, the condition is already met at the initial point. Otherwise, $v>0$ before the first time the target speed is reached, so

$$
\frac{d(v^2)}{d\ell}=2\dot v=-2a\ge-2a_*.
$$

If the target is reached after distance $d\le D$, integration gives
$v_0^2-v_{\mathrm{safe}}^2=2\int_0^d a(\ell)\,d\ell\le2a_*D$,
proving necessity. Conversely, if Eq. (C.15) holds, constant maximum deceleration $a=a_*$ reaches the target after distance $(v_0^2-v_{\mathrm{safe}}^2)/(2a_*)\le D$, proving sufficiency. $\square$

Equation (8.1) of the report is precisely Eq. (C.15) with $D=d_{\mathrm{front}}(n)$. Two practical claims must nevertheless be distinguished:

- If the physical system is **guaranteed to achieve** deceleration $a_*$ along the entire braking trajectory, Eq. (C.15) is sufficient for that system. The system may be capable of stronger braking, so the condition need not be necessary.
- If $a_*$ is merely an upper bound on achievable deceleration, Eq. (C.15) is only necessary: it does not establish that the upper bound can always be attained. Cornering friction constraints, road grade, and actuator state must be included in the realizability analysis.

A positive $v_{\mathrm{safe}}$ guarantees only that the speed is reduced to that value before the boundary. **It does not, by itself, prevent collision with an unknown obstacle located immediately beyond the boundary.** To guarantee stopping within known free space against arbitrary unknown obstacle locations, set $v_{\mathrm{safe}}=0$ and use vehicle-shape erosion, error margins, and coverage of the actual braking path to ensure that the entire stopping occupancy remains free. The curvature-based expression $v_{\mathrm{safe}}=\sqrt{a_{y,\max}/|\kappa|}$ is a lateral-acceleration speed limit under the relevant curvature and steady-state assumptions; it is not equivalent to safety against unknown obstacles.

**Path consistency.** When using the ray distance $D=d_C(p,h)$, the theorem requires the braking occupancy to follow that fixed direction, or a separate proof that the actual braking trajectory's swept region lies in certified free space and has available arc length at least $D$. The initial tangent-ray length alone does not certify a path that turns during braking. For a prescribed curved braking path $\chi(\ell)$, the appropriate distance is

$$
D_\chi=\sup\{d\ge0:\chi(\ell)+B(\psi(\ell))\subseteq K
\ \forall\ell\in[0,d]\}.
\tag{C.16}
$$

The proof of Eq. (C.15) is unchanged when $D=D_\chi$.

### Corollary C.6a: A Conservative Condition with Reaction Delay

Suppose the delay is at most $\tau$, longitudinal acceleration during the delay is at most $\alpha\ge0$, and deceleration $a_*$ is guaranteed thereafter. Define

$$
d_\tau=v_0\tau+\tfrac12\alpha\tau^2,
\qquad v_\tau=v_0+\alpha\tau.
$$

The following condition suffices to reach the target speed within the certified distance:

$$
d_\tau+\frac{(v_\tau^2-v_{\mathrm{safe}}^2)_+}{2a_*}\le D.
\tag{C.17}
$$

The proof adds the maximum distance traveled during the delay to the subsequent braking distance from Theorem C.6. The swept region over both the delay and braking phases must still be certified as in Eq. (C.16).

### Corollary C.6b: Convex Speed Constraints under Fixed Geometry

Let $v_{\mathrm{safe}}$ and $a_*$ be fixed, and let $d_C(n)$ be given by Eq. (C.12), with $C$, $s$, and $h$ fixed as required there. Within a convex lateral feasible interval, the set

$$
\{(n,v):v\ge0,\quad v^2-v_{\mathrm{safe}}^2-2a_*d_C(n)\le0\}
\tag{C.18}
$$

is convex. Indeed, $v^2$ and $-d_C(n)$ are convex, so their sum is convex and its zero sublevel set is convex. This conclusion concerns only this constraint set under fixed geometry; the complete vehicle-dynamics OCP generally remains nonconvex.

Writing Eq. (C.12) as $d_C(n)=\min_i(A_i+B_i n)$, Eq. (C.18) is equivalently imposed by the following constraint for every relevant face:

$$
v^2-v_{\mathrm{safe}}^2\le2a_*(A_i+B_i n)\quad\text{for all }i.
\tag{C.18a}
$$

This expresses exactly the same feasible set without directly evaluating the nonsmooth minimum. Moving all terms of each inequality to the left yields a convex quadratic inequality. Frontier changes appear as changes in the active constraint, so the distance minimum need not be artificially smoothed.

### Proposition C.7: Correct Soft Constraints and the Scope of Their Guarantees

Let $g=v^2-v_{\mathrm{safe}}^2-2a_*D$. The correct constraint allowing a nonnegative violation $r\ge0$ is

$$
g-r\le0,\qquad r\ge0,
\quad J=J_0+\lambda\sum_k r_k^2,\quad\lambda>0.
\tag{C.19}
$$

For fixed values of all other variables, the optimal slack is $r=\max(0,g)$: feasibility requires $r\ge\max(0,g)$, and the squared penalty is increasing on the nonnegative half-line. The $+r$ in Eqs. (8.9) and (8.11) of the report should therefore be $-r$; $g+r\le0$ only tightens the original condition.

A finite quadratic penalty weight does not guarantee $r=0$. Consequently, the softened problem does not inherit the hard safety guarantee of Theorem C.6. A hard safety claim requires a separate nonrelaxable stopping constraint, verification that the final solution has $r=0$ with numerical margin, or an independent safety controller that handles violating states.

## C.5 Historical Visibility, Reobservation, and Blind-Region Memory

### Theorem C.8: Historical Intersection for a Fixed Set of Valid Constraints

At station $s_p$, suppose every active historical viewpoint $k\in\mathcal A_p$ supplies a certified interval
$I_k(s_p)=[L_k(s_p),U_k(s_p)]$,
and those historical certificates remain valid at the current time. The exact interval satisfying all historical constraints is

$$
\overline L(s_p)=\max_{k\in\mathcal A_p}L_k(s_p),
\qquad
\overline U(s_p)=\min_{k\in\mathcal A_p}U_k(s_p),
\tag{C.20}
$$

provided that $\overline L\le\overline U$; otherwise, the intersection is empty.

**Proof.** A real number belongs to every interval if and only if it is no smaller than every lower bound and no larger than every upper bound. $\square$

Equations (8.14) and (8.15) of the report use $\mathcal A_p=\{0,\ldots,p\}$ and correctly express the conservative policy of satisfying every historical planning band simultaneously. That policy, however, does not remove old restrictions upon reobservation. If any historical interval excluded a point, the permanent intersection excludes it permanently. The formula therefore cannot directly be interpreted as prohibiting a blind region only until it is observed again.

### Theorem C.9: Reobservation Updates in a Static Environment

Suppose road obstacles are stationary, every $K_k$ is genuinely free, and spatial registration is correct. Define the accumulated known-free set and unknown set by

$$
K_p^{\mathrm{known}}=\bigcup_{k=0}^{p}K_k,
\qquad U_p=\Omega\setminus K_p^{\mathrm{known}}.
\tag{C.21}
$$

Then $K_p^{\mathrm{known}}$ remains free at time $t_p$, and

$$
U_p=U_{p-1}\setminus K_p.
\tag{C.22}
$$

**Proof.** In a static environment, every previously free point remains free, so their union remains free. Equation (C.22) follows directly from De Morgan's law for complements. $\square$

Equation (C.22) precisely implements removal of an unknown region once it has been observed to be free. It does not require previously visible regions to remain visible in every frame. If the design explicitly requires continued direct observation, a historical intersection or the current visibility region should be used deliberately, rather than conflating these sets with unknown-region memory.

To convert Eq. (C.21) into a single planning band, apply Theorem C.2 anew to the accumulated free set. A union may produce multiple lateral components; taking only the smallest lower bound and the largest upper bound is not generally valid. An active-history rule expressed through $\mathcal A_p$ is also possible, but each discarded historical constraint must be shown to have been replaced by new free-space evidence. A reobservation flag alone does not justify unconditional deletion of an entire geometric constraint.

### Theorem C.10: Historical Certificates under Bounded Obstacle Motion

Suppose each obstacle reference-point trajectory $o(t)$ satisfies
$\|o(t)-o(t')\|\le w|t-t'|$.
Assume that no unmodeled sudden appearance or creation of obstacles occurs within the domain, and that obstacle size has been handled by the appropriate spatial inflation. If $K_k$ contains no such obstacle point at time $t_k$, then

$$
\widetilde K_{k\to p}=K_k\ominus B_{w(t_p-t_k)}
\tag{C.23}
$$

still contains no obstacle point at $t_p$. Consequently, $\bigcup_{k\le p}\widetilde K_{k\to p}$ is also currently free.

**Proof.** Suppose an obstacle reaches $y\in\widetilde K_{k\to p}$ at $t_p$. Its position $z$ at $t_k$ satisfies $\|z-y\|\le w(t_p-t_k)$. By the definition of erosion, $y+B_{w(t_p-t_k)}\subseteq K_k$, hence $z\in K_k$, contradicting the absence of obstacle points in $K_k$ at that time. $\square$

To certify that the region remains free throughout a future braking interval, replace $t_p$ in Eq. (C.23) by the latest future time being certified, or certify the braking trajectory in space and time pointwise. **Current freedom from obstacles does not automatically imply freedom from obstacles throughout future braking.** This theorem supplies a conservative sufficient certificate, not a necessary characterization of traversable space among moving obstacles.

## C.6 Precise Meaning of the Unknown-Zone Soft Penalty

Let each prohibited region $U_j$ be closed, and define the continuous signed-distance function

$$
\sigma_j(x)=\operatorname{dist}(x,U_j^c)
-\operatorname{dist}(x,U_j).
\tag{C.24}
$$

It is positive in the interior of $U_j$, negative outside, and zero on the boundary. Equation (8.16) of the report can be defined precisely as

$$
J_{\mathrm{UZ}}=\lambda_{\mathrm{UZ}}
\sum_{k,j}\left[\max\{0,\sigma_j(\Phi(s_k,n_k))\}\right]^2.
\tag{C.25}
$$

The penalty is nonnegative and vanishes at all checked nodes if and only if no node enters the interior of any $U_j$. For a convex polygon $U_j=\{x:a_{ji}^\top x\le b_{ji}\ \forall i\}$ with unit normal vectors, one may also use

$$
\widehat\sigma_j(x)=\min_i(b_{ji}-a_{ji}^\top x)
\tag{C.26}
$$

as a function with the same inside/outside sign. Outside the polygon, it is generally not the negative Euclidean distance to the polygon, but it suffices for interior detection. Summing positive penalties for individual faces would penalize some points outside the polygon and cannot be claimed equivalent to Eq. (C.25).

A finite weight again does not guarantee that the final penalty vanishes. The report's use of this term as soft guidance is therefore appropriate. If a smooth approximation is used, its zero set and sign must be checked separately before assigning them the same exact meaning.

## C.7 From Discrete Nodes to Continuous Trajectories

Imposing planning-band or speed constraints only at OCP grid nodes does not establish safety between the nodes. The following are directly applicable sufficient certificates.

### Theorem C.11: Straight Connections within a Common Convex Cell

If adjacent reference points $p_k,p_{k+1}$ belong to the same certified convex set $C$, and their actual connection is the Cartesian linear interpolation
$p(\theta)=(1-\theta)p_k+\theta p_{k+1}$,
$\theta\in[0,1]$, then the entire connection lies in $C$.

**Proof.** This is the definition of convexity. $\square$

The conclusion does not automatically apply to linear interpolation in Frenet coordinates along a curved reference line, because $\Phi(s,n(s))$ is generally not a Cartesian straight line.

### Theorem C.12: Linear Interpolation within a Single Frenet Band

Suppose the valid planning band on a grid interval $[s_k,s_{k+1}]$ is $L(s)\le n\le U(s)$, where $L$ is convex and $U$ is concave. Let $n(s)$ be the linear interpolation of the endpoint values. If both endpoints satisfy the band constraints, the entire interval satisfies them.

**Proof.** For $s=(1-\theta)s_k+\theta s_{k+1}$, convexity and endpoint feasibility give

$$
L(s)\le(1-\theta)L(s_k)+\theta L(s_{k+1})
\le(1-\theta)n_k+\theta n_{k+1}=n(s).
$$

The upper bound follows similarly from concavity of $U$. $\square$

The condition holds naturally if $L,U,n$ are affine on the same pieces; every breakpoint must be included in the verification grid. General OCP interpolation polynomials or integrated trajectories require independent continuous verification.

### Theorem C.13: A Node-Margin Certificate from Derivative Bounds

For an actual continuous candidate trajectory, suppose every constraint residual $g_j(s)$ is absolutely continuous on a grid interval and satisfies $|g'_j(s)|\le M_j$ almost everywhere. If the interval length is $\Delta s$ and both endpoints satisfy

$$
g_j(s_k),g_j(s_{k+1})\le-\frac12M_j\Delta s,
\tag{C.27}
$$

then $g_j(s)\le0$ throughout the interval.

**Proof.** Every interior point is at most $\Delta s/2$ from its nearest endpoint. By the derivative bound, the residual can increase by at most $M_j\Delta s/2$ relative to that endpoint. This increase is offset by the endpoint margin, leaving the residual nonpositive. $\square$

Here, $M_j$ must be a verifiable upper bound; an estimate from finitely many samples does not constitute a proof. Geometric branch changes may be verified piecewise, or a valid Lipschitz bound may replace a classical derivative bound.

## C.8 A Combined Result for the Manuscript

**Theorem C.14 (A Combined Certificate from Geometric Frontiers to Safe Braking).** Assume that:

1. The geometric classification module returns a correct visibility region or inner approximation, with all relevant blind regions excluded.
2. The free-space certificates remain valid during vehicle passage and braking; vehicle shape and errors are covered by erosion or swept-volume verification.
3. Frenet coordinates are single-valued on the operating domain. The selected planning bands are constructed according to Theorems C.2–C.3 and hold along the entire continuous trajectory.
4. $D$ is a certified lower bound on the continuously free length available along an actually realizable braking path.
5. Longitudinal braking has the guaranteed capability in Theorem C.6. Delay is included through Eq. (C.17) when necessary, and the hard stopping constraint uses $v_{\mathrm{safe}}=0$ with no positive slack.
6. The braking policy is realizable under the vehicle's full dynamics, steering limits, and tire constraints.

Then every state satisfying these constraints admits a realizable trajectory that stops within certified free space.

**Proof.** Assumptions 1–3 validate the input free set and the reference-point or full-vehicle constraints. Assumption 4 makes the distance certificate cover the entire braking path. Assumption 5 and Theorem C.6 or Corollary C.6a ensure that speed reaches zero within that certified distance. Assumption 6 lifts the ideal longitudinal policy to a realizable policy for the full system. Assumption 2 ensures that every occupied vehicle region remains free at its corresponding time, so the stopping trajectory is collision-free. $\square$

The theorem establishes the **existence of a safe stopping policy**. Proving that a specified closed-loop controller always chooses safe actions, that the next optimization problem remains feasible, or that safety holds over an infinite horizon additionally requires a controller-selection or backup policy, recursive feasibility, and assumptions on environmental changes. These properties do not follow merely from the existence of an OCP solution at each instant.

## C.9 Exact Correspondence with Chapter 8 of the Report

| Report component | Rigorous formulation in this note | Required qualifications |
|---|---|---|
| Section 8.2, Eq. (8.1): $d_{\mathrm{front}}(n)$ and the speed limit | Theorems C.4–C.6, Eqs. (C.12)–(C.18) | First exit; finite segment; actual braking path; guaranteed deceleration |
| Section 8.3, Eq. (8.2), and Section 8.4, Eq. (8.8): convex planning bands | Theorems C.1–C.3 | Set subtraction does not ensure convexity; fixed-cross-section preimages; empty-set handling |
| Figure 17: Unknown Zones and the planning region | Eq. (C.4), convex components, and interval unions | Completeness of blind-region coverage follows from the geometric theorems in the preceding chapter |
| Eqs. (8.9), (8.11)–(8.13): slack | Proposition C.7 | Use $-r$; a finite soft penalty does not guarantee hard safety |
| Section 8.6.1, Eqs. (8.14)–(8.15): historical envelopes | Theorems C.8–C.10 | Permanent intersection, release upon reobservation, and dynamic validity are distinct mechanisms |
| Section 8.6.2, Eq. (8.16): blind-region penalty | Eqs. (C.24)–(C.26) | Inside/outside sign for the whole polygon; soft guidance only |
| Section 8.7: discrete OCP results | Theorems C.11–C.13 | Node feasibility requires additional continuous-trajectory certification before claiming continuous safety |

These corrections do not invalidate the geometric method of aligning common stations on both road boundaries along a single ray, then classifying and combining local events. They give separate, explicit assumptions and conclusions for geometric correctness, exact road visibility, conservative planning bands, and full-vehicle safety.

## Appendix C.A: Exact Visibility Has a Single Cross-Sectional Interval in a Hole-Free Road

This appendix proves a stronger result than the general statement in Theorem C.2. It directly supports representing the **exact point-vehicle visibility region of the original geometric model** by a single Frenet planning band, without first proving convexity of the two-dimensional visibility region or selecting a convex component to obtain a single interval.

### Theorem C.15: The Visible-Chord Property of a Jordan Road Domain

Let $\Omega\subset\mathbb R^2$ be a closed Jordan domain: a simple closed curve together with its bounded interior. Thus, the road has no holes and $\mathbb R^2\setminus\Omega$ is connected. Let $P\in\Omega$, and define visibility with boundary contact allowed:

$$
V(P)=\{q\in\Omega:[P,q]\subseteq\Omega\}.
\tag{C.A1}
$$

If $q_1,q_2\in V(P)$ and $[q_1,q_2]\subseteq\Omega$, then

$$
[q_1,q_2]\subseteq V(P).
\tag{C.A2}
$$

**Proof.** First suppose $P,q_1,q_2$ are not collinear, and let $\Delta$ be their closed triangle. Visibility of the endpoints gives $[P,q_1],[P,q_2]\subseteq\Omega$, and the hypothesis gives the third side $[q_1,q_2]\subseteq\Omega$. Hence $\partial\Delta\subseteq\Omega$.

Suppose there exists $z\in\operatorname{int}\Delta\setminus\Omega$. The set $\mathbb R^2\setminus\Omega$ is connected and, because $\Omega$ is bounded, contains points outside $\Delta$. It is disjoint from $\partial\Delta$. It would therefore meet both the interior and exterior components of $\mathbb R^2\setminus\partial\Delta$, contradicting connectedness. Consequently, $\Delta\subseteq\Omega$.

For every $q\in[q_1,q_2]$, the segment $[P,q]\subseteq\Delta\subseteq\Omega$, so $q\in V(P)$. If the three points are collinear, every such $q$ satisfies $[P,q]\subseteq[P,q_1]\cup[P,q_2]$ and is again visible. $\square$

### Corollary C.15a: An Exact Single Frenet Band

In the globally embedded road-strip model of the main note, every full cross-section

$$
F_s=\{\Phi(s,n):-w_R(s)\le n\le w_L(s)\}
$$

is a convex line segment contained in $\Omega$. By Theorem C.15,

$$
V(P)\cap F_s
\quad\text{is a closed line segment, a singleton, or the empty set.}
\tag{C.A3}
$$

Therefore, the exact lateral visibility set

$$
E_V(s)=\{n\in[-w_R(s),w_L(s)]:\Phi(s,n)\in V(P)\}
\tag{C.A4}
$$

has the following lossless representation whenever it is nonempty:

$$
E_V(s)=[L_V(s),U_V(s)],\qquad
L_V(s)=\min E_V(s),\quad U_V(s)=\max E_V(s).
\tag{C.A5}
$$

**Proof.** The entire chord between any two visible points on the cross-section lies in $\Omega$ and is therefore visible by Theorem C.15, proving convexity. Closedness must also be established. If $q_m\in V(P)$ and $q_m\to q$, then for every $\lambda\in[0,1]$,
$(1-\lambda)P+\lambda q_m\in\Omega$.
Because $\Omega$ is closed, the limit $(1-\lambda)P+\lambda q$ also lies in $\Omega$, hence $q\in V(P)$. Thus $V(P)$ is closed, and its intersection with the compact cross-section is a compact convex subset. In the affine coordinate $n$, it is exactly the closed interval in Eq. (C.A5). $\square$

For any additional fixed convex set $G$, such as the report's FOV half-space $\mathcal H$,

$$
(V(P)\cap G)\cap F_s=(V(P)\cap F_s)\cap(G\cap F_s)
\tag{C.A6}
$$

is still a single interval, a singleton, or the empty set, because the right-hand side is the intersection of two one-dimensional convex sets. This does not require $V(P)\cap G$ to be convex in two dimensions.

### Direct Implications for the Report

If the geometric module proves that the union of all Unknown Zones is exactly the invisible part of the road, removing those regions yields precisely $V(P)$. Intersecting this set with a fixed half-space still admits an **exact single-band representation** through Eq. (C.A5). Under the original road assumptions, nonconvexity of the two-dimensional visibility region therefore does not by itself require multiple Frenet bands. The report's phrase "the now convex set" can be replaced by the more precise statement:

> In a hole-free, globally embedded road strip, the intersection of the exact visibility region with each full Frenet cross-section is a convex interval. Lateral feasibility can therefore be represented exactly by one lower and one upper bound, without requiring the two-dimensional visibility region itself to be convex.

This does not imply that the full feasible set in $(s,n)$ is convex. As $s$ varies, $L_V(s)$ and $U_V(s)$ need not be convex and concave, respectively. Nor does it imply convexity of the OCP with vehicle dynamics.

### Scope and Relation to the General Results Above

- The theorem requires the region through which sight lines may pass to have no holes. Independent obstacles inside the road must be removed from that region. Its complement is then generally disconnected, and the triangle-filling argument no longer applies.
- The theorem concerns exact $V(P)$ or its intersection with a convex set. An arbitrary conservative subset $C\subseteq V(P)$ need not retain the property: artificially removing a middle portion of a cross-section, for example, produces a union of intervals. Theorem C.2 is therefore still needed for general $C$.
- Vehicle-shape erosion $V(P)\ominus B_\rho$, contraction of historical certificates, and other additional constraints cannot be assumed to preserve cross-sectional convexity merely from Eq. (C.A3). The transformed set requires a new proof or an explicit check. If additional constraints are known to give intervals on the same cross-sections, their intersection does preserve the interval property.
- Changing the visibility convention for boundary contact, introducing a nonconvex sensor field of view, or using an approximate overestimate of blind regions requires rechecking the assumptions. An explicit inner approximation always supports conservative planning, but its single-band representation cannot be called exact without further conditions.

This corollary is a structural result about road geometry. Once a complete geometric partition has been proved to yield exact $V(P)$, it can be passed directly to an exact point-vehicle Frenet band, strengthening the interface between the geometric proof and Chapter 8 of the report.
