# FOV Tangency Events, Visibility Partitioning, and a Certified Algorithm: Conditional Proofs

This note formulates the chain “local bearing extrema—tangent rays—same-side/opposite-side intersections—visibility partitioning” as verifiable mathematical statements. The proofs concern a static, planar geometric model with opaque road boundaries. Geometric visibility is not identified with reliable sensor detection or vehicle dynamic safety. The original report and code are not modified.

## 1. Model, Regularity, and Required Conventions

Let the centerline $\gamma:[a,b]\to\mathbb R^2$ be parameterized by arc length, with unit tangent and normal $t(s),n(s)$ and curvature $\kappa(s)$. Let the left and right widths satisfy $w_L,w_R>0$. The closed road is

$$
\overline\Omega=\{\Phi(s,z)=\gamma(s)+z n(s):a\le s\le b,\ -w_R(s)\le z\le w_L(s)\}.
$$

Assume that $\Phi$ is an injective embedding of this strip and $1-\kappa(s)z>0$. The road is therefore homeomorphic to a rectangle, with two side boundaries and two caps, and every point has a unique longitudinal coordinate $S(q)$. Write the left and right boundaries as

$$
B_\eta(s)=\gamma(s)+h_\eta(s)n(s),\quad h_L=w_L,\quad h_R=-w_R.
$$

Using second derivatives requires the boundaries to be piecewise $C^2$. For example, $\gamma\in C^3$ and $C^2$ widths are sufficient. A centerline that is only $C^2$ does not by itself guarantee $C^2$ normal-offset boundaries; piecewise $C^2$ boundaries may instead be assumed directly. The formulas below hold on each smooth patch, with junctions treated as separate events.

Fix a viewpoint $P\in\Omega$, with $S(P)=s_0$ and $a<s_0<b$. To simplify boundary cases, place the viewpoint in the open road interior and subsequently restrict the output to $S(q)\ge s_0$. A finite search cap truncates the computational domain; it does not represent a physical obstacle. For a closed track, first select and explicitly specify an embedded forward window that does not wrap around. This note makes no claim to have analyzed the region beyond that window.

The default convention **allows ideal grazing**:

$$
V(P)=\{q\in\overline\Omega:[P,q]\subseteq\overline\Omega\}.
$$

Thus a ray that touches a tangency point may remain visible beyond it if it stays inside the road. A tangency point is not necessarily an occluding endpoint. If strictly positive clearance is desired, visibility can instead be computed in an eroded road $\Omega_\varepsilon$; this is clearer than conflating the two definitions.

A finite-event algorithm requires a further condition: all tangency roots and boundary intersections of each ray under consideration must admit finite enumeration. $C^2$ regularity alone does not guarantee this; smooth functions can have infinitely many accumulating zeros. Any of the following additional sets of conditions may be used:

- Use polynomial, rational, or suitable algebraic representations on finitely many curve patches, with certified root solving, and exclude identically zero relevant equations.
- Explicitly assume that the relevant roots are finite and have all been certified.
- Require every tangency root on the compact parameter interval to be nondegenerate, with endpoint degeneracies handled separately. Together with interval exclusion, this establishes finiteness.

The general-position version excludes entire boundary segments that coincide with a ray. Such a segment must not be silently skipped as if it were an ordinary isolated root; it requires separate one-dimensional interval geometry.

## 2. Bearing Derivatives and Tangency Events

For either boundary $B(s)$, write $r(s)=B(s)-P$. Use a continuously unwrapped bearing angle

$$
\theta(s)=\operatorname{unwrap}\operatorname{atan2}(r_y(s),r_x(s))-\psi_P.
$$

A fixed reference heading $\psi_P$ adds or subtracts only a constant and does not change the locations of extrema. A branch jump of the raw $\operatorname{atan2}$ at $\pm\pi$ is not a geometric extremum.

**Proposition 1 (bearing derivative).** On a smooth patch with $r(s)\ne0$,

$$
\theta'(s)=\frac{g(s)}{\|r(s)\|^2},\qquad g(s)=r(s)\times B'(s).
$$

Here $x\times y=x_1y_2-x_2y_1$. Consequently, $\theta'(s_T)=0$ if and only if the line from the vehicle to the boundary point is parallel to the boundary tangent.

**Proof.** Differentiating $\operatorname{atan2}(y,x)$ gives $(xy'-yx')/(x^2+y^2)$; substituting $r$ gives the formula. When $B'\ne0$ and $r\ne0$, a zero cross product is equivalent to parallel vectors. ∎

The Frenet boundary derivative is

$$
B_\eta'(s)=(1-\kappa(s)h_\eta(s))t(s)+h_\eta'(s)n(s).
$$

This also shows why a centerline-curvature threshold alone cannot determine boundary tangency events. Width variation changes the boundary tangent, and small local centerline curvature does not prove the absence of occlusion events.

**Proposition 2 (nondegenerate extrema).** If $g(s_T)=0$ and

$$
g'(s_T)=r(s_T)\times B''(s_T)\ne0,
$$

then $s_T$ is a strict local extremum, and

$$
\theta''(s_T)=\frac{r(s_T)\times B''(s_T)}{\|r(s_T)\|^2}.
$$

**Proof.** We have $g'=B'\times B'+r\times B''=r\times B''$. When differentiating $g/\|r\|^2$ at a root, the term involving the derivative of the denominator vanishes because it is multiplied by $g=0$. A nonzero second derivative gives a strict local extremum and also makes the root isolated. ∎

If $r(s_T)=\lambda B'(s_T)$, then

$$
g'(s_T)=\lambda\,[B'(s_T)\times B''(s_T)]
=\lambda\,\kappa_B(s_T)\|B'(s_T)\|^3.
$$

A tangency event at zero boundary curvature may therefore be degenerate. The condition $g=0$ guarantees only a stationary bearing or tangency, not a strict extremum; a higher-order stationary point may leave visibility unchanged. Both types of simple extrema should enter the candidate set and then undergo visibility and direction classification. For a general complex road, a rule that keeps only maxima on the left and minima on the right cannot replace a completeness proof.

## 3. A Precise Version of “Every Local Extremum Implies a Blind Region”

The following concepts must be distinguished:

1. **Tangency candidate:** $g(s_T)=0$. This is a potential visibility-change event.
2. **Forward candidate:** additionally satisfies $s_T>s_0$ and the selected sensor angular range, if such a range is imposed. Without a sensor-angle restriction, do not add the latter condition.
3. **Visible tangency point:** $[P,T]\subseteq\overline\Omega$. This permits earlier grazing at other collinear tangency points.
4. **Exposed tangency point:** $[P,T)\subset\Omega$, the root is nondegenerate, and there is no other degenerate event at the same angle. This excludes events dominated by a nearer occluder.
5. **Shadow-front event:** the first-exit boundary branches actually differ on the two sides of the tangent ray. This can be certified by the complete event classification in the next section.

**Proposition 3 (an exposed nondegenerate tangency creates a local shadow).** Let $T$ be an exposed nondegenerate tangency point, and consider only targets in the road interior. Then some neighborhood beyond $T$ contains a nonempty open set of invisible road points, and the tangent ray forms part of the boundary of this local shadow.

**Proof sketch.** Take $T$ as the origin and the ray from the vehicle through the tangency point as the $x$-axis, with $P=(-\rho,0)$. A nondegenerate tangency gives a local boundary graph $y=f(x)$, where $f(0)=f'(0)=0$ and $f''(0)\ne0$. Because $[P,T)$ approaches $T$ from inside the road, the tangent line locally lies inside the road on both sides of the tangency, rather than outside. Rotating the ray slightly toward the occluded side makes it cross the boundary into the exterior before reaching targets beyond the tangency. Rotating slightly toward the other side allows it to continue into a small neighborhood beyond the tangency. Compact subsegments of $[P,T)$ away from $T$ have positive clearance, so sufficiently small directional perturbations are not first blocked by another obstacle. Hence there are nonempty local regions that are invisible on one side and visible on the other. ∎

This cannot be strengthened to “each extremum corresponds to an independent blind region.” Hidden tangencies may be entirely redundant; collinear tangencies can merge into one event; local shadows can be nested or overlap; degenerate stationary points may cause no state change; and events outside the sensor range do not affect the current field of view. The precise statement is: **all relevant tangency events form a candidate set for visibility changes, and certified exposed nondegenerate events create local occlusion fronts.**

## 4. First Intersection and First Exit

For a unit direction $u_\alpha$, let $R_\alpha(\rho)=P+\rho u_\alpha$, $\rho\ge0$. Define

$$
\rho_{\rm exit}(\alpha)=\sup\{r\ge0:R_\alpha([0,r])\subseteq\overline\Omega\}.
$$

This defines the **initial visible prefix** from the viewpoint. Points on a ray that subsequently leaves and reenters the road remain invisible.

- In a generic direction with no tangency or vertex event, the smallest positive boundary intersection equals $\rho_{\rm exit}$.
- In a tangency-event direction, the nearest boundary intersection may only be a grazing contact. Check whether the ray subsequently enters the exterior. If it is inside the road on both sides of the contact, continue to the actual exit.
- Once the ray has passed through the exterior, a farther tangency or intersection cannot restore visibility.

The ray-intersection equation is

$$
F_\alpha(s)=(B(s)-P)\times u_\alpha=0,
\qquad (B(s)-P)\cdot u_\alpha>0.
$$

An angular-alignment test with a tolerance does not establish an exact intersection. Finding one local root with an optimizer does not establish that it is the first intersection.

## 5. Complete Event Partitioning and an Exact Certified Algorithm

### 5.1 Event Set

Construct a finite set of directions $E$ containing:

- Ray directions of all tangency roots on every smooth boundary patch.
- Directions to patch endpoints, road-cap corners, and nonsmooth junctions.
- The two limiting directions of the sensor angular range, if a range restriction is imposed.

Merge coincident directions while recording every event on each direction. Sort directions around the viewpoint; consecutive event directions define open angular sectors. Only finitely many event cells are sampled here, rather than using arbitrary dense sampling of the entire field of view.

### 5.2 An Implementable Certified Version

```text
CertifiedVisibility(P, boundary patches):
  1. Certifiably enumerate every root of g(s)=0; add patch-endpoint
     and corner directions to obtain E.
  2. Sort E cyclically and merge events with the same angle.
  3. For each open sector I:
       Choose a direction alpha* certified to lie inside I.
       Enumerate and sort all positive intersections of this ray
       with all boundary patches.
       Identify the continuous intersection branch containing the
       nearest boundary intersection to P.
       Represent or certifiably continue this branch throughout I
       as rho_first(alpha).
       Output the sector's visible set 0 <= rho <= rho_first(alpha).
  4. For each event ray:
       Enumerate all positive intersections and sort them by rho.
       Certify road membership between consecutive intersections.
       Starting at P, retain the entire prefix up to the first
       exterior interval; allow grazing contacts that do not enter
       the exterior to be passed.
  5. Combine the sector and event-ray outputs, then clip by forward
     station and the declared sensor range.
```

In Step 3, the numerical radius from a single direction cannot be used for the entire sector. Store the actual continuous intersection branch, or a parameter interval and equation sufficient to evaluate that branch accurately in any direction within the sector.

**Theorem 4 (completeness of the event algorithm).** Under the embedding, finite-root, and certified root-enumeration assumptions of Section 1, the algorithm above returns exactly the defined $V(P)$, with the declared range clipping.

**Proof.**

1. Within an open sector, every ray intersection with a smooth boundary is transverse. If $\partial_sF_\alpha=B'(s)\times u_\alpha=0$, the point also satisfies the tangency equation and its direction should already be an event direction, a contradiction.
2. By the implicit function theorem, every intersection continues continuously through the sector. An intersection cannot appear or disappear without a tangency, cannot pass a patch endpoint without an endpoint event, and cannot pass through the viewpoint because $P\in\Omega$. Compactness of the road excludes changes involving intersections at infinity within a finite sector.
3. Distinct intersection branches cannot exchange their radial order. An exchange would require the same direction and distance, hence the same spatial point. Boundary embedding excludes self-intersections of nonadjacent patches, and patch junctions have already been added to the event set.
4. The identity of the first-intersection branch is therefore fixed in each open sector. Since the viewpoint is in the interior, this first intersection is the first exit. Points at smaller radii are exactly the visible points, while points beyond it are invisible because the sightline passes outside the road.
5. Between consecutive intersections on an event ray, no boundary is present, so road membership is constant and can be certified once. The longest initial prefix containing no exterior interval is exactly the visible portion under the grazing-allowed definition.
6. Every direction belongs to an open sector or an event ray. These two classes omit no direction, proving the result. ∎

Event rays must be handled separately rather than recovered solely as the closure of open-sector outputs. Multiple collinear grazing contacts can create limiting visible sets of zero angular width and zero clearance. Such sets may deliberately be excluded in applications requiring strict clearance, but the model must state this explicitly.

### 5.3 Equivalent Boundary-Segment Formulation

Add **all intersections** between every event ray and the boundary as subdivision points. This partitions the boundary into finitely many open arcs. The angle on each arc is monotone and lies in one open sector, and its intersection branch has a fixed rank in the radial ordering. Its visible/invisible status is consequently constant. Choosing one interior representative point per arc and applying a first-exit test certifies the entire arc. Event endpoints are handled separately.

This gives a conservative complete partition that is straightforward to prove. Retaining only necessary prime occluders may reduce the partition, but the removal rule itself requires proof. For published background on the general “tangency—ray reintersection—segmentation” construction, see Section 3.1 of Elber et al. (2005). The proof here is independently based on transversality, the implicit function theorem, and invariance of root ordering.

## 6. Road Ordering: Why Same-Side/Opposite-Side Classification Can Replace General Verification

This section develops a structure more specialized than a general curve-visibility algorithm. It identifies the usable assumptions and explains why longitudinal road station cannot be applied unconditionally to arbitrary free curves.

The preceding sections assume $P\in\Omega$. The proposed forward window may instead begin at $s_0$, placing $P$ in the relative interior of the starting cap; the complete definition and proof adjustments appear in Section 6.3. In that version, the “interior cell containing the viewpoint” always means the cell containing a forward interior marker $P_\varepsilon$ sufficiently close to $P$. The actual origin of every sightline remains $P$, not $P_\varepsilon$.

**Lemma 5 (straight lines and Frenet fibers).** If a line $\ell$ does not coincide with the supporting line of any road-normal fiber, then $S$ is strictly monotone on every connected line segment $Q\subset\overline\Omega\cap\ell$.

**Proof.** If two distinct points of $Q$ had the same $S=s$, both would lie on $\ell$ and on the normal line $\gamma(s)+\mathbb R n(s)$. Two lines with two intersections must coincide, contradicting the assumption. Therefore $S|_Q$ is continuous and injective. A continuous injective function on an interval is strictly monotone. ∎

Directions coinciding with a normal fiber must be treated separately; the strict-monotonicity conclusion does not apply to them. The noncoincidence condition may be localized to the longitudinal interval actually used by the algorithm.

For an in-road ray $q(\rho)=P+\rho u$, differentiation of the inverse Frenet coordinates gives

$$
\frac{dS}{d\rho}=\frac{u\cdot t(S)}{1-\kappa(S)z}.
$$

With a positive denominator, $u\cdot t(s_T)>0$ is a direct local test for increasing station beyond the tangency. Together with fiber injectivity, it prevents reversal of order along the entire continuous interior ray segment. A candidate on a branch that does not satisfy this forward condition must not be processed directly by the forward minimum-station rule.

**Corollary 6 (sufficient conditions for finding the first exit by road station).** Suppose that immediately beyond a tangency point $T$, the selected ray remains inside the road and station increases along its initial interior path. Let $U$ be the subsequent first true exit. Assume there is no unhandled intermediate grazing contact and that Lemma 5's noncoincidence condition holds. Then $s_U>s_T$, and no farther boundary intersection on the ray can have $S\in(s_T,s_U)$. Consequently, selecting the smallest longitudinal station among all relevant intersections with $s>s_T$ is equivalent to selecting the first exit beyond the tangency.

**Proof.** On $[T,U]$, $S$ is continuous and strictly increasing, covering $[s_T,s_U]$. If a farther intersection had a longitudinal coordinate in the interior of this interval, its normal fiber would already have intersected the same ray once on $[T,U]$. This would give two intersections, contradicting Lemma 5. ∎

An intermediate grazing contact, if present, must be processed as a further event. Under the grazing-allowed model, it cannot automatically be treated as an exit. If $[P,T]\subseteq\overline\Omega$, $s_T>s_0$, and no degenerate collinear event intervenes, the same fiber lemma supplies the forward condition at the tangency. If $T$ is already hidden, inward continuation and forward motion immediately beyond the tangency must be checked independently.

### 6.1 A True Radial Proper Crosscut

A segment $Q=[T,U]$ is a proper crosscut of the road if and only if its endpoints lie on the boundary, its open segment $(T,U)\subset\Omega$, and it has no other boundary contact. Additionally require $P,T,U$ to be collinear, with $T,U$ on the same positive ray from $P$, and $P\notin Q$. A **radial window** used here to certify occlusion also requires $U$ to be an actual exit: a short ray segment beyond $U$ must lie outside the road. The first transverse boundary intersection satisfies this condition. The general-position assumptions exclude a ray tangent at both $T$ and $U$.

- If $T,U$ lie on opposite side boundaries, $Q$ splits the rectangular strip into two components containing the starting and terminal caps, respectively.
- If both endpoints lie on the same side, $Q$ and the arc between those endpoints on that side that passes through no cap form a Jordan curve. They cut off a lateral pocket containing neither cap. It must be checked that $P$ lies in the other component.

**Proposition 7 (invisibility beyond a radial window).** For either type of radial window above with an actual exit endpoint, every point in the interior of the component not containing $P$ is invisible.

**Proof.** If such a point $q$ were visible, $[P,q]$ would be a path in the closed road from one component to the other and would have to cross $Q$. If the line through $P,q$ differs from the line through $P,T,U$, the lines meet only at $P$, but $P\notin Q$, a contradiction. If they are the same line and the point is beyond $U$, the ray first enters the exterior immediately after the actual exit; subsequent reentry cannot restore visibility. Boundary and collinear degenerate cases are handled using the first-exit definition of Section 4. ∎

A general proper crosscut without a guarantee of actual exit beyond $U$ yields a Jordan-separation argument that directly excludes only noncollinear points. A collinear ray allowing zero-clearance grazing may continue through an endpoint. This is why the actual-exit condition is essential.

This explains the geometric content of the proposed classification: same-side events can form local shadow pockets, whereas opposite-side events create a front cutting across the entire road corridor. The construction can certify downstream invisibility without requiring $T$ itself to be visible. However, calling $U$ the actual farthest visible point additionally requires certification of the visible prefix $[P,T]$.

### 6.2 Two-Stage Station Comparison

Let $Q_k$ be a family of opposite-side proper crosscuts on distinct rays. Distinct rays intersect only at $P$, and $P\notin Q_k$, so these crosscuts are pairwise disjoint. For disjoint left-to-right crosscuts, Jordan separation gives the same endpoint order on both sidewalls: if $Q_i$ precedes $Q_j$, then

$$
s_L^{(i)}<s_L^{(j)},\qquad s_R^{(i)}<s_R^{(j)}.
$$

**Corollary 8 (sufficient conditions for the proposed min-hit rule).** If every candidate is generated by a forward hit satisfying $s_U^{(k)}>s_T^{(k)}$, then

$$
s_U^{(k)}=\max\{s_L^{(k)},s_R^{(k)}\}.
$$

Since both endpoint sequences have the same order, $\max$ also strictly preserves that order. Hence

$$
k_*\in\arg\min_k s_U^{(k)}
$$

selects the earliest crosscut, even if the hit endpoints of different candidates lie on different sidewalls.

This is a **conditional but substantive road-specific theorem**: the first stage uses station order to determine and classify the exit among same-side and opposite-side candidates for each tangency; the second stage minimizes hit station over opposite-side candidates. It does not apply directly to uncertified approximately aligned points. Forward hits, proper crosscuts, a shared longitudinal parameter, and candidate completeness cannot be omitted.

If the complete prefix $P\to T\to U$ of this earliest cut is visible and attention is restricted to its upstream road component, fiber monotonicity further gives a maximum visible station of $s_U$. Station on the cut does not exceed $s_U$; the substrip with $s>s_U$ connects to the terminal cap without crossing the cut and therefore lies in the invisible component; and $U$ itself is visible. Without a verified visible prefix, $s_U$ is only an upper bound, and equality cannot be asserted.

Completeness still requires proof that the candidate extractor retains every exposed event capable of creating the earliest cutoff or a lateral pocket. Angle and curvature thresholds, minimum tangency spacing, and local optimization that returns only one root can invalidate this step unless a certified fallback is provided.

### 6.3 A Direct Completeness Proof Using All Same-Side and Opposite-Side Windows

The proposed method need not be reformulated as a sector-by-sector solver. The following argument directly shows that, under general-position conditions, certified enumeration of all tangency events, generation of their same-side/opposite-side windows, and removal of the corresponding shadow regions produce the complete visible region. The method in Section 5 can serve as a degeneracy fallback, verifier, or baseline.

**Interior and cap viewpoints.** This subsection may retain the open road $\Omega$ and interior viewpoint $P$ from Section 1, enumerating events in all directions. To consider only the proposed forward window, instead use

$$
\overline\Omega_+=\{\Phi(s,z):s_0\le s\le b,\ -w_R(s)\le z\le w_L(s)\},
\qquad \Omega_+=\operatorname{int}\overline\Omega_+,
$$

with $P=\Phi(s_0,z_0)$ in the relative interior of the starting cap. After choosing this version, suppress the subscript $+$: $\Omega$ always denotes the open interior and $\overline\Omega$ the closed road including its caps. Visibility remains defined by $[P,q]\subseteq\overline\Omega$. Choose $P_\varepsilon=P+\varepsilon t(s_0)$ as an interior cell marker. Since $z_0$ is in the interior of the width interval and the event set is finite, a single sufficiently small $\varepsilon>0$ can be chosen so that the marker is visible and lies on the observer side of every partition cut. In the interior-viewpoint version, use $P$ itself as the marker.

In the cap version, an exposed tangency requires $(P,T)\subset\Omega$, allowing the true ray origin $P$ to remain on the cap. Every visible ray reaching the road interior initially points strictly inward, so the Frenet ordering lemma ensures increasing forward station for its subsequent relevant tangency and hit points. In the interior-viewpoint version, computing the complete omnidirectional visible region requires more than the events with $s>s_0$; the forward two-stage pseudocode below specifically concerns the cap version.

Define the regularized output

$$
V_{\rm reg}(P)=\overline{\operatorname{int}_{\mathbb R^2}V(P)}^{\,\overline\Omega}.
$$

This is the relative closure of the final open visible cell. It does not mean adding back every original road-boundary point or every candidate-window boundary.

Impose the following general-position assumptions: relevant tangency roots are nondegenerate; no positive ray has two boundary tangencies; corner directions are handled separately; the ray stays inside the road beyond an exposed tangency; and its next boundary contact is a transverse exit or an explicitly handled search cap. Let $\mathcal W$ contain the radial window from each exposed tangency $T$ to that first exit $U$.

**Theorem 9 (tangency windows completely describe the interior visibility boundary).** Under the conditions above and apart from separately handled corner events,

$$
\partial V(P)\cap\Omega
=\bigcup_{Q\in\mathcal W} \bigl(Q\cap\Omega\bigr).
$$

**Proof.**

1. The set $V(P)$ is closed. If $q_j\to q$ and every $[P,q_j]\subset\overline\Omega$, then the limit point at each fixed segment fraction remains in the closed set $\overline\Omega$, giving $[P,q]\subset\overline\Omega$.
2. Take $x\in\partial V(P)\cap\Omega$. Closedness of $V$ makes $x$ visible. If $P$ is interior and $[P,x]$ has no boundary contact, the compact segment has positive clearance, making a small neighborhood of $x$ visible, a contradiction. If $P$ is on the starting cap, a visible ray to an interior point $x$ is strictly inward at $P$: an outward direction exits immediately, while a direction along the cap's supporting line cannot reach a visible interior point from the cap's relative interior. Hence a short prefix has an inward property stable under small perturbations of the target. After removing that prefix, the remaining compact segment has uniform positive clearance if it has no boundary contact; this again makes a neighborhood of $x$ visible. Thus, in either version, $(P,x)$ must contain a boundary contact $T$. This does not assert that the entire open segment beyond the cap origin has a uniform positive distance from the boundary.
3. The contact cannot be transverse, since that would place part of the segment outside the road, contradicting visibility of $x$. It is therefore a tangency or a separately listed corner event. General position gives a unique nondegenerate tangency $T$ with $(P,T)\subset\Omega$; the interior version also has $P\in\Omega$. Point $x$ lies after $T$ and before the first actual exit $U$, so it belongs to the associated window. Forward ordering in the cap version follows from the fiber-order conditions stated above.
4. Conversely, take an interior point $x$ of an exposed window. The segment $[P,x]$ has no boundary contact except $T$ and possibly the cap origin. If the origin is on the cap, use a short strictly inward prefix for stability. The other compact subsegments away from $T$ have positive clearance. By Proposition 3's local nondegenerate tangency structure, small target perturbations toward one side are occluded and perturbations toward the other side are visible. Hence $x\in\partial V$. ∎

**Lemma 10 (shadow regions are nested or disjoint).** Windows on distinct rays are disjoint. For each window, denote by $H_Q$ the open component that it cuts off away from the observer marker, which is $P$ in the interior version and $P_\varepsilon$ in the cap version. Then any two regions $H_Q,H_R$ are disjoint or one contains the other.

**Proof.** Two radial lines meet only at $P$, and the windows do not contain $P$. In a topological disk, a proper crosscut divides the road into two components. The interior of a second disjoint crosscut must lie entirely in one of them. The region it further cuts off away from the same observer marker is consequently a subregion of the first region, a region containing it, or a disjoint region; partial crossing overlap is impossible. Equivalently, map the domain to a rectangle and apply Jordan separation. ∎

**Corollary 11 (the observer cell after removing all shadows is fully visible).** Define $C_P$ explicitly as the connected component of the open set $\Omega\setminus\bigcup_{Q\in\mathcal W}Q$ containing the observer marker. Then

$$
C_P=\operatorname{int}_{\mathbb R^2} V(P),\qquad
V_{\rm reg}(P)=\overline{C_P}^{\,\overline\Omega}.
$$

Every other open cell is invisible. The output is precisely this relative closure; all original road boundaries or hidden-candidate boundaries must not be added unconditionally. In general position without additional isolated grazing rays, the output equals $V(P)$. If degeneracies are allowed and pointwise exact visibility under the original grazing definition in Section 1 is required, use a first-exit test to check lower-dimensional sets on the event rays separately.

**Proof.** Proposition 7 guarantees that every $H_Q$ is invisible. Conversely, Theorem 9 guarantees that $C_P$ contains no point of $\partial V$. In the interior version, a neighborhood of $P$ is visible. In the cap version, $P_\varepsilon$ has a visible neighborhood ensured by a short strictly inward prefix at the origin and clearance along the remaining compact segment. If connected $C_P$ contained an invisible point, a path from this neighborhood to that point would cross $\partial V$, a contradiction. Thus all of $C_P$ is visible. Every other cell is separated from the observer side by at least one window and lies in that window's $H_Q$, hence is invisible. Taking the relative closure according to the regularized-output definition gives the second identity. ∎

This yields the complete certification procedure for the proposed forward-window method with a cap viewpoint:

```text
RoadOrderedVisibility:
  A. Certifiably enumerate all relevant tangency and corner events.
  B. For each event, compare forward stations on both boundaries
     to find and certify the first exit. Retain a local shadow
     pocket for a same-side event and a cross-road shadow for an
     opposite-side event.
  C. Certify every record by proper-window and actual-exit checks.
     Use an exact event fallback for failed, degenerate, or
     uncertifiably ordered events.
  D. Remove the windows from the open road and select the open
     cell C_P containing P_epsilon. Output
     closure_in_closed_road(C_P), without automatically restoring
     any other boundary points.
  E. Apply Corollary 8's min-hit rule to opposite-side windows
     to select the forward cutoff. Retain the same-side shadow
     records and event history preceding that cutoff.
```

Hidden candidates need not all be discarded beforehand. If a candidate truly generates a window satisfying Proposition 7, its shadow remains invisible, although it may already be contained in an earlier shadow and thus be redundant. Its ray boundary must not be added to the visible set without verification. Whenever a curvature or peak threshold removes a boundary patch, an additional certificate must ensure that at least every true front in $\mathcal W$ has been retained.

If the algorithm outputs only the earliest opposite-side front and the unknown zones before it, the output should explicitly be described as the currently certified forward corridor. Events beyond the opposite-side cutoff may be skipped because they are already invisible. In a window with no opposite-side cutoff, however, all relevant same-side events up to the cap must be retained. Failure to find an FPV must not automatically be interpreted as the absence of unknown zones.

## 7. Complexity: Which Claims Are Currently Justified

Let $N$ be the number of curve patches, $K$ the number of isolated tangency roots, $C$ the number of corner/junction directions, and $E\le K+C$ the number of distinct event directions. Let $I(N)$ denote the cost of certified intersection of one ray with all curve patches, $M$ the maximum number of intersections of such a ray, and $T_g(N)$ the cost of certifying all tangency roots.

The structural complexity of the naive event-partitioning method is

$$
O\big(T_g(N)+E\log E+E[I(N)+M\log M]+T_{\rm continuation}\big).
$$

The bit complexity of exact algebraic-number comparison and root isolation cannot be treated as constant time. This expression counts geometric operations. Outputting continuous curve boundaries also requires accounting for branch representation or certified continuation.

If the curve patches are polynomials of fixed degree $d$, the ray-intersection polynomial has degree at most $d$. The tangency polynomial $(B-P)\times B'$ has degree at most $2d-2$: its highest-degree term cancels because a vector's cross product with itself is zero. For a rational curve $B=(X/W,Y/W)$, $W\ne0$, let $Q=(X-P_xW,Y-P_yW)$. The tangency numerator is $Q\times Q'$, with the same degree bound. Thus the numbers of isolated roots satisfy $K=O(Nd)$ and $M=O(Nd)$. For fixed degree, temporarily disregarding coefficient bit complexity, naively intersecting every patch in every sector is roughly quadratic; including sorting gives $O(N^2\log N)$.

For an exact normal offset of a centerline spline, the boundary is generally not a polynomial of the same degree. Degree bounds for a cubic boundary must therefore not be applied directly; the count must be derived from the actual boundary representation.

If the two-stage road method scans both arrays of length $M$ from scratch for every tangency, its worst-case cost remains $O(KM)$; it is not inherently linear. An improvement such as $O(K+M)$ requires proof that pointers or root brackets can be reused monotonically across events with a bounded total number of advances. Early termination, few events in practice, and smaller constants are also valuable, but must be distinguished from asymptotic complexity.

## 8. Certifiable Conditions for Numerical Implementation and Discretization Error

### 8.1 Tangency Localization

For an isolated root, if $|g'|\ge\gamma>0$ on a certified bracket and the computation/model error satisfies $|\widehat g-g|\le\varepsilon_g$, then the parameter error on the correct root branch is bounded by

$$
|\widehat s_T-s_T|\le\varepsilon_g/\gamma
$$

plus the width of the remaining numerical root bracket. This requires the root branch to have been identified and all other roots to have been enumerated. As a root approaches degeneracy, $\gamma\to0$, and the bound naturally deteriorates.

Noiseless uniform sampling on a known unimodal interval can localize the sample maximum to within approximately one step of the true peak. Without that prior knowledge, any fixed step can skip a narrow pair of extrema. **Dense sampling is not a general completeness condition.** Completeness can instead be established by isolating all roots on every spline patch or proving $0\notin g(I)$ for each interval declared root-free.

### 8.2 Alignment and Intersection Localization

For the intersection equation $F(s)=(B(s)-P)\times u=0$, suppose transversality on a bracket satisfies

$$
|F'(s)|=|B'(s)\times u|\ge\beta>0,
$$

Then a normal residual $\varepsilon_F$ causes parameter error of at most approximately $\varepsilon_F/\beta$; an interval formulation can make this bound rigorous. If the candidate is within distance $R$ of the tangency point, an angular-alignment tolerance $\varepsilon_{\rm dir}$ gives only the normal-residual bound $R\sin\varepsilon_{\rm dir}$. Converting this to longitudinal error still requires division by a lower bound on transversality. Near a tangency, $\beta$ is small, so a fixed angular tolerance may correspond to a large index or station error.

The ordering of multiple intersections must be certified using disjoint root or position intervals. If each candidate station has error at most $\varepsilon_s$, a separation greater than $2\varepsilon_s$ is required to preserve ordering directly. Overlapping intervals require further refinement, merged degenerate events, or conservative output; a tiny floating-point difference must not decide the order.

### 8.3 Polyline Approximation

If $B\in C^2$, the parameter step is $h$, and $\sup\|B''\|\le M_2$, the positional error of linear interpolation satisfies

$$
\|B(s)-B_{\rm chord}(s)\|\le M_2h^2/8.
$$

This is only a geometric approximation error. It does not automatically imply a uniform bound on visible distance, because visibility can change discontinuously at narrow passages, near bitangencies, and limiting grazing configurations.

A provably conservative approach is to establish a free-space signed-distance error $|d_{\Omega}-d_{\widehat\Omega}|\le\delta$ and compute in the eroded approximate road $\{d_{\widehat\Omega}\ge\delta\}$. This set is contained in the true free space, so its visible set is conservative as well. A boundary Hausdorff bound alone, without certified topology and inside/outside orientation, does not justify this inclusion.

### 8.4 Filtering and Fallback in the Original Method

- Minimum-angle, centerline-curvature, and tangency-spacing thresholds may be used for heuristic acceleration, but every rejected interval requires an independent certificate that it contains no effective occlusion event.
- When merging unknown zones, use an outer approximation of the occluded set or verify directly that the resulting planning corridor lies in the visible region. Geometric smoothing must not erase a real blind region.
- Failure to find an FPV does not mean that the road ahead is fully visible. Return a certified visible prefix, use the computational-window cap, or invoke a complete-event fallback.
- One practical architecture is: a fast two-stage road algorithm generates candidates; root isolation and necessary ray checks certify them; frames that fail certification use the complete event algorithm. This preserves the opportunity for efficiency in the specialized method while providing an explicit correctness fallback.

## 9. Scope of Claims Suitable for a Paper

Once both the proof and implementation satisfy the conditions above, the following claim is justified:

> In an embedded Frenet road strip with finitely many nondegenerate tangency events and forward proper crosscuts, two-stage event comparison uses the shared longitudinal order of the left and right boundaries to identify the earliest cross-road visibility cutoff while retaining same-side local shadow pockets. Complete root enumeration or a certified fallback guarantees completeness of the output.

The following should not currently be stated without qualification: every extremum creates an independent blind region; arbitrary C² curves have finitely many events; approximate alignment identifies the first intersection; enlarging the window never changes the result; or the method is inherently linear. These are propositions that need additional assumptions or counterexample testing, not reasons to dismiss the method's value.

Published geometric background is cited only to position the method: Elber, Sayegh, Barequet, and Martin, *Two-Dimensional Visibility Charts for Continuous Curves*, 2005, Section 3.1, DOI 10.1109/SMI.2005.48. Author-provided full text: https://gershon.cs.technion.ac.il/papers/crv_vis.pdf . General event partitioning has prior precedent. The aspects that particularly merit separate proof and evaluation here are the order comparisons supported by Frenet fibers, the two occlusion classes, and a certified implementation.
