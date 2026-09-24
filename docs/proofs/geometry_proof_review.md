# Proof Review of Two-Stage Event Comparison in a Frenet Road Strip

Date: 2026-09-22. Independent mathematical review retaining the complete multiple-event treatment in Section 6.3 of the report. Same-side events generate shadow pockets, opposite-side events generate cross-road fronts, and all opposite-side events are then compared. The original report and source code have not been modified.

## 1. Model and Three Distinct Levels of the Argument

Let the centerline $\gamma:[a,b]\to\mathbb R^2$ be parameterized by arc length, with $\tau=\gamma'$ and left normal $N$. Both road widths are continuous and strictly positive. Define

$$
K=\{(s,n):a\le s\le b,\ -w_R(s)\le n\le w_L(s)\},
\qquad X(s,n)=\gamma(s)+nN(s).
$$

Assume that $X$ is a homeomorphism from $K$ onto the closed road strip $\overline\Omega$ and satisfies $1-\kappa(s)n>0$ wherever smooth. Thus, the road does not intersect itself, and $\Omega$ is a Jordan domain. More importantly, each station fiber

$$
F_s=X(\{s\}\times[-w_R(s),w_L(s)])
$$

is a convex straight line segment. This is stronger than assuming an arbitrary topological road polygon; counterexamples involving arbitrary polygons cannot replace an analysis under this assumption.

The viewpoint $P$ has station $s_0$, and the events under consideration satisfy $s>s_0$. Throughout the following argument, set $a=s_0$, truncate the road at $s_0$, place $P$ in the relative interior of the starting cap, and consider only directions entering the road from that cap. Choose $r$ sufficiently small that $\overline\Omega\cap B(P,r)$ is a closed half-disk meeting neither another boundary portion nor any event cut. Choose $P_\varepsilon$ in its interior. Every reference to the “observer side” means the interior connected component containing $P_\varepsilon$; the actual viewpoint remains $P$, not $P_\varepsilon$. Finiteness of the events and their disjointness from $P$ allow a single such $r$ to be chosen.

Define $V(P)=\{z\in\overline\Omega:[P,z]\subset\overline\Omega\}$ and $V_{\rm reg}(P)=\operatorname{closure}_{\overline\Omega}(\operatorname{int}V(P))$. Since $V(P)$ is closed, $V_{\rm reg}\subset V$. The basic general-position version below excludes isolated, zero-width grazing sightlines, in which case the two sets coincide. If degeneracies are retained, the region-closure statement concerns $V_{\rm reg}$, and the event rays themselves must be classified separately.

The original track may be closed, but the theorem concerns a finite embedded window after unwrapping. If the initial and terminal caps of an entire lap are identified to form an annulus, the Jordan-disk conclusions cannot be reused without further justification.

Let a candidate tangency point be $T=P+\rho_T d$ and its hit point be $U=P+\rho_U d$, with $0<\rho_T<\rho_U$ and unit ray direction $d$. A **proper radial crosscut** is $Q=[T,U]$, where $T$ and $U$ lie on the road boundary and the relative interior of $Q$ lies entirely in $\Omega$. To avoid zero-clearance degeneracies under closed-domain visibility, the nondegenerate conclusions below additionally require $U$ to be a transverse exit: a short ray interval immediately beyond $U$ lies outside $\overline\Omega$.

Three different claims must be distinguished:

1. A local angular extremum generates a candidate event.
2. A same-side/opposite-side station comparison returns a proper crosscut.
3. Combining the shadows of all crosscuts exactly recovers the visible region.

There are strong positive theorems for the second and third claims, but each has its own assumptions. The existence of a local extremum alone does not establish a complete safety certificate.

## 2. Station Monotonicity along a Straight Segment

**Lemma 1 (convex-fiber lemma).** If $L$ is any straight line segment contained in $\overline\Omega$, then the station function $s$ along $L$ is monotone, allowing constant intervals.

**Proof.** The function $s$ is continuous. Every level set $\{z\in L:s(z)=c\}=L\cap F_c$ is the intersection of two convex line segments and is therefore empty or a connected interval. If a continuous real-valued function is not monotone, there are three ordered parameters $t_1<t_2<t_3$ for which the middle value is strictly greater than both endpoint values or strictly less than both. Choose a level between the middle value and the two endpoint values. The intermediate value theorem gives a preimage of that level on each side of the middle parameter, with a point outside the level set between them. The level set is therefore disconnected, a contradiction. ∎

This proof includes the degenerate case in which $L$ is collinear with a normal fiber; that case need not be excluded.

**Corollary.** If $Q=[T,U]$ is a proper crosscut with $s_U>s_T$, then station is nondecreasing on $Q$, with $\max_Q s=s_U$ and $\min_Q s=s_T$.

## 3. First Stage: When the First Station Hit Is the Geometric First Hit

**Theorem 2 (sufficient conditions for a station-first certificate).** Starting from $T$ in direction $d$, assume that:

- A short open ray interval immediately beyond $T$ lies entirely in $\Omega$.
- Station strictly increases on this short interval.
- $H$ is the first subsequent boundary contact, and the ray exits the domain immediately beyond $H$.
- The search enumerates all exact boundary intersections on the ray satisfying $\rho>\rho_T$ and $s>s_T$.

Then $H$ is precisely the intersection with the smallest station among these intersections. Consequently, finding the smallest station on each of the same and opposite boundaries and comparing them returns the geometric first hit; $[T,H]$ is a proper crosscut.

**Proof.** Since $[T,H]\subset\overline\Omega$, Lemma 1 and the local increase imply that station is nondecreasing there and $s_H>s_T$. Suppose there were a further ray intersection $W$ beyond $H$ with $s_T<s_W\le s_H$. The intermediate value theorem gives $B\in[T,H]$ such that $s_B=s_W$. Points $B$ and $W$ belong to the same convex normal fiber $F_{s_W}$, so $[B,W]\subset F_{s_W}\subset\overline\Omega$. However, $[B,W]$ contains the interval immediately beyond $H$ where the ray leaves the domain, a contradiction. Every further intersection with $\text{station}>s_T$ therefore has station greater than $s_H$. Thus $H$ uniquely minimizes station. ∎

If the first exit occurs at the terminal cap, it is treated as a search-window endpoint event. The same argument shows that it must not be replaced by an in-window sidewall intersection selected later along the ray beyond that cap.

**Checkable form.** At a regular Frenet point, the ray satisfies

$$
\frac{ds}{d\rho}=\frac{d\cdot\tau(s)}{1-\kappa(s)n}.
$$

Since the denominator is positive, $d\cdot\tau(s_T)>0$ guarantees local forward station motion. Combined with a local test that the ray enters $\Omega$ after the tangency point, this allows the theorem to be applied. Prior visibility of $T$ is not required; hidden tangency points may still participate in the complete event partition.

If $[P,T]$ is visible, $s_T>s_0$, and $T$ is a tangency that actually generates a shadow window, these local conditions hold automatically: station monotonicity along the visible ray provides the forward ordering, and the shadow window itself provides the inward continuation.

The theorem clarifies the distinction from general polygons. In a general polygon, station ordering and radial distance ordering may conflict. In an embedded Frenet strip, the stated local conditions and convex-fiber structure exclude the conflicts that would invalidate the first-stage comparison.

## 4. The Jordan-Separation Meaning of Same-Side and Opposite-Side Events

**Theorem 3 (shadow property of a radial crosscut).** Let $Q=[T,U]$ be a proper radial crosscut, with $P\notin Q$ and $P,T,U$ collinear, and impose the transverse-exit condition at $U$ stated above. By the Jordan crosscut theorem, $\Omega\setminus Q$ has exactly two components. Write $A_Q$ for the component containing $P_\varepsilon$ and $S_Q$ for the other component. No point in $S_Q$ is visible from $P$.

**Proof.** Suppose $z$ were visible. Since $z$ is an interior road point, the initial portion of its sightline from $P$ enters the cut-free local half-disk defined above and belongs to the same interior component as $P_\varepsilon$. To reach $S_Q$, the sightline would therefore have to cross $Q$. If $z$ is not on the supporting line $\ell=\operatorname{aff}(P,T,U)$, then $[P,z]$ meets $\ell$ only at $P$. Since $P\notin Q$, the sightline cannot cross $Q$, a contradiction. If $z$ lies on $\ell$ but before $Q$, then $[P,z]$ still does not meet $Q$, and the same separation argument applies. If $z$ lies beyond $U$, then $[P,z]$ contains the exit interval immediately beyond $U$, so $z$ is not visible. The cut $Q$ itself is a component boundary and does not belong to $S_Q$. ∎

**Degenerate cases.** If the transverse-exit condition at $U$ is removed while boundary grazing is still counted as visible, noncollinear points remain invisible, but collinear zero-clearance rays require separate treatment. A robust-visibility definition or a general-position convention may be used; this qualification cannot simply be omitted.

**Same-side case.** If $T$ and $U$ are on the same sidewall, the crosscut and the sidewall arc between its endpoints that passes through neither cap enclose a local shadow pocket. If $s_T,s_U>s_0$, this pocket does not contain $P$. This is the topological meaning of the report's same-side→unknown-zone classification. The boundary of the actual region must include the road-boundary arc; it is not a triangle formed by three collinear points.

**Opposite-side case.** If $T$ and $U$ lie on the left and right sidewalls, respectively, $Q$ separates the starting cap from the terminal cap and is therefore a cross-road front. Its far-side component is invisible. This conclusion does not depend on proving $T$ visible in advance; the complete method may retain hidden events as redundant subdivisions.

## 5. Second Stage: Why $\min(s_U)$ Selects the Earliest Cross-Road Front

**Lemma 4 (consistent endpoint order of disjoint cross-road cuts).** For two disjoint proper left-to-right crosscuts $Q_1,Q_2$, the station order of their endpoints on the left sidewall is the same as the station order of their endpoints on the right sidewall.

**Proof.** If the orders were reversed, the four endpoints of the two cuts would alternate along the Jordan boundary. The two boundary arcs separated by one cut would then contain the two endpoints of the other cut separately. The latter would have to cross the former, contradicting disjointness. ∎

Distinct rays from a common $P$ have no intersection other than $P$, and none of the cuts contains $P$. The lemma therefore applies to all proper cross-road cuts on distinct rays. Degenerate repeated events on the same ray can first be merged or handled using the first transverse exit.

**Theorem 5 (second-stage station comparison).** For each opposite-side cut, denote its left and right endpoint stations by $(l_j,r_j)$. The original method requires the hit to satisfy $s_U>s_T$, so

$$
s_U=\max(l_j,r_j).
$$

The cuts form a chain in their topological front-to-back order. Lemma 4 shows that both $l_j$ and $r_j$ increase along this chain, so $\max(l_j,r_j)$ also increases along it. Thus $\min_j s_U$ selects exactly the complete cross-road cut nearest the starting cap.

This is not a numerical heuristic, and it does not require interpreting sidewall station as a global Euclidean ray distance. It is the precise topological reason why the second comparison in the report is valid.

## 6. Cross-Road Fronts and the Maximum Visible Station

**Theorem 6 (upper bound and equality from one valid cross-road cut).** Let $Q=[T,U]$ be a proper left-to-right radial crosscut with $s_U>s_T$, and suppose its far-side component is invisible. Then every visible point $z$ satisfies $s(z)\le s_U$. If $U$ is also visible, the maximum visible station is exactly $s_U$.

**Proof.** Lemma 1 gives the station range $[s_T,s_U]$ on $Q$. The entire substrip $\{s>s_U\}$ is disjoint from $Q$, and every point in that substrip can be connected to the terminal cap within the same substrip. Hence the entire substrip lies in the far-side component and is invisible, establishing the upper bound. If $U$ is visible, it attains that bound. ∎

In particular, $[P,T]\subset\overline\Omega$ together with properness of $Q$ guarantees $[P,U]\subset\overline\Omega$, so $U$ is visible. There is no need to search for extrema of the station gradient in the two-dimensional road interior: the convex straight-fiber structure already excludes a point in the near-side component whose station exceeds $s_U$.

A hidden cross-road cut alone provides an upper bound. Completeness of the windows or a visible-prefix certificate is required to strengthen that bound to equality.

## 7. Closing the Argument for the Complete Multiple-Event Partition

The following are sufficient conditions for a main theorem. Suppose the sidewalls are piecewise smooth and the events are finite. Include all relevant smooth tangencies and occluding polygonal corners; exclude degeneracies such as simultaneous multiple tangencies or handle them explicitly. Generate a proper radial cut satisfying Theorems 2 and 3 for every event, and remove its shadow component away from $P$. The essential completeness condition is:

Every point in $\operatorname{int}(\Omega)\cap\partial V$ lies on at least one generated cut.

**The completeness argument can be made precise as follows.** For an interior visibility-boundary point $z$, closed-domain visibility gives $[P,z]\subset\overline\Omega$. If $[P,z]$, apart from $P$, stays clear of the road boundary, a small neighborhood of $z$ remains visible, contradicting $z\in\partial V$. Hence the segment has a boundary contact. Take $T$ to be its last contact before $z$. This contact must be a smooth tangency or a suitable nonconvex corner; a transverse contact would send part of the ray outside the road, contradicting visibility. The segment from $T$ to its first subsequent boundary hit $U$ is a proper radial window containing $z$. A complete list of tangency and corner events must enumerate $T$.

All cuts are nonintersecting radial segments from $P$. Their far-side shadow components form a laminar family: any two are either disjoint or one contains the other. To handle invisible arcs of the original sidewall correctly, define the output explicitly:

$$
\begin{aligned}
C_P&=\text{the interior connected component of }
\\
&\qquad\Omega\setminus\Bigl(\bigcup_j Q_j\Bigr)\text{ containing }P_\varepsilon,\\
K_{\rm vis}&=\operatorname{closure}_{\overline\Omega}(C_P).
\end{aligned}
$$

It is insufficient to remove open shadow components from the closed road and then retain all original boundary points indiscriminately; that would incorrectly retain invisible sidewall arcs. By definition, $C_P$ is connected, and it contains no interior visibility boundary covered by the complete event list. It contains visible interior points near $P$. Therefore it cannot also contain invisible interior points, since an interior path connecting the two kinds of points would have to cross the visibility boundary. Conversely, Theorem 3 proves that every separated component is invisible, so no visible interior point is removed. Consequently $K_{\rm vis}=V_{\rm reg}$, and in the basic nondegenerate model $K_{\rm vis}$ also equals $V$.

A hidden tangency point need not be removed in advance. If it generates a certified cut, it only subdivides an already invisible region and does not invalidate the conclusion.

**Short proof that the earliest opposite-side cut is automatically visible.** Let $Q^*$ be the earliest opposite-side cut. For any other cut $R$, disjointness implies that the relative interior of $Q^*$ lies entirely on one side of $R$. If $R$ is a same-side cut, the closure of its far-side pocket meets only one sidewall and cannot contain $Q^*$, which joins both sidewalls. If $R$ is a sidewall-to-terminal-cap window, its far side meets only one sidewall and that cap and likewise cannot contain $Q^*$. If $R$ is another opposite-side cut, Theorem 5 places $Q^*$ on its near side. Hence no other shadow excludes $Q^*$. The cuts are finite and pairwise disjoint, so every interior point of $Q^*$ has a small half-neighborhood on its near side contained in $C_P$. Therefore $Q^*\subset\operatorname{closure}(C_P)=V_{\rm reg}$. Closedness places its endpoints in $V_{\rm reg}\subset V$ as well, so $T^*$ and $U^*$ are both visible. Together with Theorems 5 and 6, this proves that $\min(s_U)$ from the complete event algorithm equals the maximum visible station, without an additional assumption that the earliest tangency point is visible.

The argument requires general-position conditions such as the absence of shared event endpoints. Coincident events must be handled under a consistent degeneracy convention. If there is no opposite-side cut, the event list is complete, and terminal-cap windows are handled correctly, then the visible core reaches the search terminal cap, and the maximum station is $b$.

## 8. What Can Be Concluded about the Report and Current Code

The following positive statements can be proved:

- The same-side/opposite-side classification has a precise shadow-pocket/cross-road-front interpretation.
- Under the Frenet structure and checkable local conditions, station-first selection agrees with the geometric first hit.
- Taking $\min(s_U)$ over all opposite-side events selects the earliest cross-road front exactly.
- Correct and complete partitioning by all events recovers the visible region and identifies that minimum with the maximum visible station.

The following conditions still need to be stated explicitly in the algorithm or checked:

- Use exact intersections or intersections with certified error bounds; an arbitrary direction-threshold sample hit cannot directly be treated as a proper cut.
- Check inward continuation and forward station motion after the tangency point, or certify properness of the entire cut directly.
- Use the actual boundary arc for a same-side shadow pocket.
- Include all true visibility windows in the candidate list. Heuristic filters based on $\kappa$, $\theta$, or `min_arc` do not inherit completeness without proof.
- Establish consistent conventions for grazing, multiple tangencies, overlap, the horizon cap, and station unwrapping on closed tracks.

This review provides no counterexample in which both comparison stages fail while all the stated sufficient conditions hold. On the contrary, both stages can be proved rigorously under those conditions. It also does not claim that a code snapshot whose compliance with the conditions has not yet been verified is already a certified implementation of the theorem.
