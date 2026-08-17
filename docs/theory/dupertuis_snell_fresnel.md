# Generalization of complex Snell–Descartes and Fresnel laws

> **Transcription note (for GOAD theory work).** This is a faithful markdown
> transcription of the scientific paper listed below, prepared to support the
> branch-selection / Fresnel-coefficient discussion in the GOAD far-field and
> refraction work. Equation numbers are the paper's own, so later work can cite
> "Dupertuis eq. (N)". The full reference list is omitted; figure captions are
> given briefly. Uncertain glyph readings are marked `[?]` inline and collected
> at the end.
>
> - Source file: `GeneralizationComplexSnellFresnel.pdf`
> - Citation: M. A. Dupertuis, M. Proctor, and B. Acklin, "Generalization of
>   complex Snell–Descartes and Fresnel laws," *J. Opt. Soc. Am. A* **11**(3),
>   1159–1166 (March 1994).
> - Transcribed: 2026-08-17

---

**Authors**

M. A. Dupertuis and M. Proctor
Institut de Micro- et Optoélectronique, Département de Physique, Ecole
Polytechnique Fédérale de Lausanne, PHB-Ecublens, CH-1015 Lausanne, Switzerland

B. Acklin
Institut de Microtechnique de l'Université, Rue A. L. Bréguet 2, CH-2000
Neuchâtel, Switzerland

Received March 19, 1993; revision received August 13, 1993; accepted August 27, 1993

## Abstract

Generalized Snell–Descartes and Fresnel laws are derived for harmonic
inhomogeneous plane waves that are incident upon a static interface between two
continuous absorbing dielectric media that are macroscopically characterized by
their electric and magnetic permittivities and their conductivities. A
coordinate-free formalism based on complex vector algebra is used to carry out
all discussions. Surprisingly, the usual complex Snell–Descartes laws for
reflection and refraction and Fresnel laws for polarization are recovered only in
the special case in which the vector characterizing the direction of
inhomogeneity is in the plane of incidence. In the more general case a new
deflection angle between planes of incidence and refraction has to be introduced.
An experiment is proposed to test this prediction. A generalized form of the TE
and TM modes (with respect to the interface), which are elliptically polarized
and which are called parallel electric and parallel magnetic modes, also emerges.

## Introduction

The study of light propagation has been central to the emergence of ancient and
modern physical sciences and has closely depended on and stimulated the
development of mathematics. The authors recall the invention of the optical
lunette and the optical microscope early in the seventeenth century, Descartes's
discovery of the basic laws of geometric optics in 1637, and Huygens and
Newton's theory of colors in 1704. Most solid foundations of modern classical
optics were set up by Fresnel and Maxwell (1820 and 1864, respectively). By the
end of the nineteenth century the study of light propagation in moving dielectric
media culminated in Einstein's special relativity (1905). In the twentieth
century many authors continued to study the classical nonrelativistic laws of
refraction, reflection, energy balance, and reciprocity for plane waves incident
upon an interface; a partial but already extensive list of papers is provided
(Refs. 2–23).

This paper is concerned with a generalization of the basic laws of refraction and
reflection at a static interface, that is, generalized Snell–Descartes and Fresnel
laws, which, surprisingly, the authors "have not been able to find in the
literature." The study deals with the reflection and refraction of an arbitrarily
polarized harmonic inhomogeneous plane wave (HIPW) at an interface between two
homogeneous isotropic absorbing (or amplifying) dielectric media.

Several review articles on electromagnetic propagation of plane waves through an
interface are cited (König in the *Handbuch der Physik* Vol. XX, 1928; papers by
Šantavý and Knittl on reversibility in thin-film optics; a review by von
Fragstein summarizing a recurrent controversy about the Poynting vector; work by
Pincherle and Mahan, a quantum treatment of evanescent waves, a derivation based
on molecular optics, modern analytical investigations of reciprocity and angular
momentum balance, and supercurrent investigations that include diffraction
effects due to finite apertures). Finally, in a very recent paper Chen and Nelson
derived Snell–Descartes and Fresnel laws without resorting to boundary
conditions, using a "wave-vector space method," an asymptotic method similar to
those used in quantum scattering. One of the strongest motivations for many
recent studies has been the optics of thin films, which plays an increasingly
important role in modern optoelectronic technology.

HIPW's are stationary solutions of Maxwell's equations for an unbounded medium and
take the following form:

$$\mathbf{E} = \tfrac{1}{2}\{\mathbf{E}_0\exp[-i(\mathbf{k}\cdot\mathbf{r}-\omega t)] + \text{c.c.}\}, \tag{1}$$

where $\omega$ is real and $\mathbf{k}$ is the complex wave vector. The wave
amplitudes may be exponentially decreasing or increasing in a direction
perpendicular to the direction of propagation (in a nonabsorbing medium). As a
consequence of the general formulation, a generalized form of the
Snell–Descartes and Fresnel laws is found (Ref. 24), which are derived and given
below. HIPW's are very natural and useful in a number of problems, such as the
description of evanescent and surface waves and the well-known exact Sommerfeld
solution to the diffraction of a plane wave by a conducting half-plane. HIPW's are
also often used in two-dimensional boundary-value problems. Their usefulness stems
from the fact that the path of integration of a usual Fourier transform of an
arbitrary bounded distribution can be modified in the complex plane by use of the
Cauchy theorem. The bounded distribution can then be viewed as a linear
superposition of HIPW's.

The observability of new effects issuing from the laws derived in this paper is
dependent on the possibility of generating experimentally bounded HIPW's that, in
a region of space, would be close to idealization, like ordinary homogeneous
plane waves. This point is discussed in the Conclusions. **Generalized energy
balance and reciprocity relations for HIPW will be investigated in a subsequent
paper.²⁹ In that paper arbitrary stacks of absorbing (or amplifying) dielectric
layers will be considered because of their importance for the optics of thin
films.** Because some of the general properties of HIPW's do not seem to be well
known, they are presented below in detail.

## Harmonic Inhomogeneous Plane Waves

Maxwell's equations for the complex components of a HIPW in a homogeneous
isotropic medium, characterized by its electric and magnetic permittivities
$\epsilon$ and $\mu$, respectively, and its conductivity $\sigma$, can be written
as

$$\mathbf{k}\cdot\mathbf{E}_0 = \mathbf{k}\cdot\mathbf{H}_0 = 0, \tag{2}$$

$$\mathbf{k}\wedge\mathbf{E}_0 = \omega\mu_0\mu\,\mathbf{H}_0, \tag{3}$$

$$\mathbf{k}\wedge\mathbf{H}_0 = -\omega\epsilon_0\bar\epsilon\,\mathbf{E}_0, \tag{4}$$

where $\bar\epsilon = \epsilon - i\sigma/(\omega\epsilon_0)$ is the complex
dielectric constant. Here the complex vector scalar product is to be interpreted
as it would be for a real vector. All the usual real and complex vector
identities that are extensively used throughout this formulation are presented in
Appendix A. Note the differences with respect to the frequent scalar product
$\mathbf{a}\circ\mathbf{b} = \mathbf{a}^*\cdot\mathbf{b}$.

Equations (2)–(4) include the complex transversality (2) and imply the complex
perpendicularity $\mathbf{E}_0\cdot\mathbf{H}_0 = 0$. They do not, however, imply
that HIPW's are purely transverse waves or that $\mathbf{E}(\mathbf{r},t)$ is
always perpendicular to $\mathbf{H}(\mathbf{r},t)$. Substituting Eq. (3) into
Eq. (4) and using Eqs. (A4) and (2), one can find the dispersion relation

$$\mathbf{k}^2 - \omega^2\mu_0\epsilon_0\,\mu\bar\epsilon = 0. \tag{5}$$

Relation (5) suggests a relation to disentangle the complex wave vector
$\mathbf{k}$:

$$\mathbf{k} = \tilde n\,k_0\,\mathbf{n}, \tag{6}$$

where $\tilde n = \sqrt{\mu\bar\epsilon}$ is the complex refractive index,
$k_0 = \omega/c$, and $\mathbf{n}$ is a complex unit vector satisfying
$\mathbf{n}\cdot\mathbf{n} = 1$. The vector $\mathbf{n}$ can be defined with the
inhomogeneity parameter $\beta$ and two real, perpendicular unit vectors:

$$\mathbf{n} = \hat n_\parallel\cosh\beta + i\,\hat n_\perp\sinh\beta. \tag{7}$$

The parameter $\beta$ is restricted to positive values to avoid redundancy (by
reversing the sign of both $\beta$ and $\hat n_\perp$). It is a parameter whose
nature is geometric and whose values may be imposed by boundary conditions at
infinity or at an interface (e.g., in the case of evanescent or surface waves).
Planes of constant phase are simply defined by $\mathrm{Re}(\mathbf{k})\cdot
\mathbf{r} = 0$, and planes of constant amplitude are defined by
$\mathrm{Im}(\mathbf{k})\cdot\mathbf{r} = 0$. It is well known that planes of
constant amplitude and constant phase coincide only when the wave is homogeneous
($\beta = 0$). It is easy to see that this case is equivalent to

$$\mathbf{k}\wedge\mathbf{k}^* = |\tilde n|^2 k_0^2(2i\sinh\beta\cosh\beta)(\hat n_\perp\wedge\hat n_\parallel) = 0, \tag{8}$$

which, with the help of Eq. (A3), implies that

$$\mathbf{k}^*\cdot\mathbf{E}_0 = -\frac{1}{\omega\epsilon_0\bar\epsilon}\mathbf{k}^*\cdot(\mathbf{k}\wedge\mathbf{H}_0) = -\frac{1}{\omega\epsilon_0\bar\epsilon}(\mathbf{k}^*\wedge\mathbf{k})\cdot\mathbf{H}_0 = 0. \tag{9}$$

However, the inverse is not true, because using Eqs. (A6) and (2) one can show
that

$$(\mathbf{k}\wedge\mathbf{k}^*)\cdot(\mathbf{E}_0\wedge\mathbf{E}_0^*) = -|\mathbf{k}\cdot\mathbf{E}_0^*|^2. \tag{10}$$

Similar considerations can be phrased in terms of the magnetic field
$\mathbf{H}_0$ because of the symmetry of Maxwell's equations (2)–(4), but unless
they are necessary for clarity such symmetric considerations are omitted. The
magnetic-field component can always be obtained from $\mathbf{E}_0$ by the use of
Eq. (3).

## Polarization Properties of HIPW

HIPW's have the following unusual polarization properties:

**1.** HIPW's are not purely transverse waves: a longitudinal component
$E_{0\parallel}\hat n_\parallel$ is implied by the complex transversality
condition (2). If the electric field is kept along a right-handed coordinate
system $\{\hat e_1, \hat e_2, \hat n_\parallel\}$,

$$\mathbf{E}_0 = E_{01}\hat e_1 + E_{02}\hat e_2 + E_{0\parallel}\hat n_\parallel, \tag{11}$$

the longitudinal component is equal to

$$E_{0\parallel} = (-i\tanh\beta)(E_{01}\cos\psi + E_{02}\sin\psi), \tag{12}$$

where the direction of the inhomogeneity of the wave is specified by the angle
$\psi$ such that $\hat n_\perp = \hat e_1\cos\psi + \hat e_2\sin\psi$. However the
field vector still oscillates in a plane and along an ellipse. Rewriting Eq. (1)
as

$$\mathbf{E} = \{\mathrm{Re}(\mathbf{E}_0)\cos[\mathrm{Re}(\mathbf{k})\cdot\mathbf{r}-\omega t] + \mathrm{Im}(\mathbf{E}_0)\sin[\mathrm{Re}(\mathbf{k})\cdot\mathbf{r}-\omega t]\}\exp[-\mathrm{Im}(\mathbf{k})\cdot\mathbf{r}], \tag{13}$$

one sees that the normal to the oscillation plane of the electric-field vector is
given by

$$\hat n_E \sim \mathrm{Re}(\mathbf{E}_0)\wedge\mathrm{Im}(\mathbf{E}_0) = -\frac{1}{4i}\mathbf{E}_0\wedge\mathbf{E}_0^*. \tag{14}$$

Of course the normal remains undefined for $E$-linearly polarized waves precisely
characterized by

$$\mathbf{E}_0\wedge\mathbf{E}_0^* = 0. \tag{15}$$

This relation stems from the fact that in an $E$-linearly polarized wave the
electric-field derivative is always parallel to the electric-field vector itself
[$\mathbf{E}\wedge\partial\mathbf{E}/\partial t = 0$ is equivalent to Eq. (15)].

**2.** Interestingly, HIPW's cannot exhibit simultaneous linear polarization for
$\mathbf{E}_0$ and $\mathbf{H}_0$. To demonstrate this property, examine the
vectorial product $\mathbf{k}\wedge\mathbf{k}^*$ and substitute the following
relationship [obtained from Eqs. (A4) and (2)]:

$$\mathbf{E}_0^2\,\mathbf{k} = -\mathbf{E}_0\wedge(\mathbf{E}_0\wedge\mathbf{k}). \tag{16}$$

With the help of Eqs. (3) and (A7) one finds that

$$(\mathbf{E}_0^2\mathbf{E}_0^{*2})(\mathbf{k}\wedge\mathbf{k}^*) = (\omega\mu_0\mu)^2(\mathbf{E}_0\wedge\mathbf{H}_0)\wedge(\mathbf{E}_0^*\wedge\mathbf{H}_0^*) = (\omega\mu_0\mu)^2\{[(\mathbf{H}_0\wedge\mathbf{H}_0^*)\cdot\mathbf{E}_0]\mathbf{E}_0^* + [(\mathbf{E}_0\wedge\mathbf{E}_0^*)\cdot\mathbf{H}_0]\mathbf{H}_0^*\}. \tag{17}$$

Thus, if $E$ and $H$ are simultaneously linearly polarized, the right-hand side
vanishes and the wave is necessarily homogeneous [cf. Eq. (8)].

**3.** On the other hand a circularly polarized HIPW is always simultaneously
circular for $E$ and $H$. A circular polarization is characterized by a constant
$|\mathbf{E}|$, which implies (using Eq. (1)) that

$$\mathbf{E}_0^2 = 0. \tag{18}$$

Using Eqs. (4) and (2) it is easy to see that this is also equivalent to

$$\mathbf{H}_0^2 = 0 \iff \mathbf{E}_0\wedge\mathbf{H}_0 = 0, \tag{19}$$

where the last equivalence follows directly from Eq. (16). To prove the statement
concerning circular polarization, it suffices to notice that
$\mathbf{H}_0^2 = 0 \iff |\mathbf{H}| = \text{constant}$.

**4.** Because of the longitudinal component, the plane of oscillation of
$\mathbf{E}$ and $\mathbf{H}$ is not perpendicular to the direction of propagation
$\hat n_\parallel$ for a HIPW. Moreover the magnetic field $\mathbf{H}$ is not in
general perpendicular to the electric field unless one polarization is linear:

$$\mathbf{E}\cdot\mathbf{H} = \tfrac{1}{2}\,\mathrm{Re}(\mathbf{E}_0\cdot\mathbf{H}_0^*) \sim (\tilde n\,\mathbf{n} - \tilde n^*\mathbf{n}^*)\cdot(\mathbf{E}_0\wedge\mathbf{E}_0^*). \tag{20}$$

Equation (20) also shows that the fields are always perpendicular in the usual
degenerate case, in which the wave is homogeneous and the medium is nonabsorbing.

**5.** Finally, the polarization of a HIPW is characterized with Jones's
parameterization:

$$E_{01} = A\exp(i\delta)(\cos\varphi\cos\varepsilon - i\sin\varphi\sin\varepsilon), \tag{21}$$

$$E_{02} = A\exp(i\delta)(\sin\varphi\cos\varepsilon + i\cos\varphi\sin\varepsilon). \tag{22}$$

Here $\varphi$ is the azimuth and $\varepsilon$ the ellipticity angle of the field
with respect to the direct reference frame $\hat e_1, \hat e_2, \hat n_\parallel$.
$A$ and $\delta$ are the amplitude and the absolute phase of the field,
respectively.

**The complex Poynting vector for HIPW's, the Poynting theorem, and the mixed
Poynting vector will be discussed in a subsequent paper.²⁹**

## Reflection and Refraction of a HIPW at an Interface

To study the transmission of an arbitrary HIPW that is incident upon an interface
between two homogeneous isotropic media, the parameterization is carefully
completed. The interface is supposed to include the origin, its direction being
specified by the real unit vector $\hat s$ normal to the surface. The following
notations are adopted for each wave: the incident wave is denoted by the subscript
$i$, the reflected wave by the subscript $r$, and the transmitted wave by the
subscript $t$. A polar parameterization is used for the complex refractive indices
of the two media, $\tilde n_j = |\tilde n_j|\exp(i\eta_j)$, where $j = i, r, t$. To
characterize the respective orientations of the waves and the interface the angles
of incidence, reflection, and refraction are introduced by

$$\cos\theta_j = \hat s\cdot\hat n_{\parallel j}, \qquad \sin\theta_j = \hat e_{1j}\cdot(\hat s\wedge\hat n_{\parallel j}), \qquad j = i, r, t. \tag{23}$$

This definition finally fixes the direction of the real unit vector $\hat e_1$ [as
well as the respective right-handed coordinate system of the incident wave
$\{\hat e_{1i}, \hat e_{2i}, \hat n_{\parallel i}\}$; see Eq. (11)]. $\hat e_{1i}$
is the direction perpendicular to the plane of incidence, and it is chosen
$\hat e_{1i} = \hat e_{1r}$. This choice is important for the description of
inhomogeneous waves and the polarization. To fix the respective orientation of
$\hat e_{1t}$ [and the respective right-handed coordinate system of the refracted
wave $\{\hat e_{1t}, \hat e_{2t}, \hat n_{\parallel t}\}$; see Eq. (11)], a
preliminary study has shown that a new angle $\Delta$ must be introduced, defined
by

$$\cos\Delta = \hat e_{1i}\cdot\hat e_{1t}, \qquad \sin\Delta = \hat s\cdot(\hat e_{1i}\wedge\hat e_{1t}). \tag{24}$$

The appearance of this angle is new in the context of this problem, and it is
called the **deflection angle** because it is the angle between the planes of
incidence and refraction. In Fig. 1 a pictorial representation of the angles
involved is given.

[Figure 1: Coordinate systems, vectors, and angles introduced in the text to
describe the incident, reflected, and transmitted HIPW's at the interface. A
figurative exponential tail is drawn along $\hat n_{\parallel}$ to represent the
HIPW. $\hat n_{\perp j}$ is always perpendicular to $\hat n_{\parallel j}$ and
makes an angle $\psi_j$ with $\hat e_{1j}$ ($j = i, r, t$).]

## Kinematics: Generalized Snell–Descartes Laws

The phase-matching condition for the waves at the interface requires

$$\mathbf{k}_i\cdot\mathbf{r}_\parallel = \mathbf{k}_r\cdot\mathbf{r}_\parallel = \mathbf{k}_t\cdot\mathbf{r}_\parallel, \tag{25}$$

where $\mathbf{r}_\parallel$ is such that $\mathbf{r}_\parallel = -\hat s\wedge(\hat s\wedge\mathbf{r}_\parallel)$. Substituting in Eq. (25) the complex decomposition of $\mathbf{k}$ in perpendicular and parallel components $\mathbf{k} = (\hat s\cdot\mathbf{k})\hat s - \hat s\wedge(\hat s\wedge\mathbf{k})$ and using Eq. (A3), one obtains

$$(\hat s\wedge\mathbf{k}_i)\cdot(\hat s\wedge\mathbf{r}_\parallel) = (\hat s\wedge\mathbf{k}_r)\cdot(\hat s\wedge\mathbf{r}_\parallel) = (\hat s\wedge\mathbf{k}_t)\cdot(\hat s\wedge\mathbf{r}_\parallel). \tag{26}$$

Because Eq. (26) is valid for any $\mathbf{r}_\parallel$ one sees that the complex
in-plane components of $\mathbf{k}$ must be conserved:

$$\hat s\wedge\mathbf{k}_i = \hat s\wedge\mathbf{k}_r = \hat s\wedge\mathbf{k}_t. \tag{27}$$

Together with the dispersion relation (5), Eq. (27) also implies for the normal
component of the reflected wave vector,

$$\hat s\cdot\mathbf{k}_r = -\hat s\cdot\mathbf{k}_i, \tag{28}$$

as expected.

The apparent simplicity of Eq. (27) hides a generalization of the complex
Snell–Descartes laws. To show this generalization one projects these equations on
the in-plane basis $\{\hat e_{1i}, \hat s\wedge\hat e_{1i}\}$ and extensively uses
the vector identities of Appendix A and the parameterization of the HIPW's
detailed above. After lengthy calculations one finds for reflection

$$\cosh\beta_r\sin\theta_r = \cosh\beta_i\sin\theta_i, \tag{29}$$

$$\sinh\beta_r\sin\psi_r\cos\theta_r = \sinh\beta_i\sin\psi_i\cos\theta_i, \tag{30}$$

$$\sinh\beta_r\cos\psi_r = \sinh\beta_i\cos\psi_i. \tag{31}$$

For refraction one finds

$$|\tilde n_t|\cosh\beta_t\sin\theta_t\cos\Delta = |\tilde n_i|[\cos(\eta_i-\eta_t)\cosh\beta_i\sin\theta_i + \sin(\eta_i-\eta_t)\sinh\beta_i\sin\psi_i\cos\theta_i], \tag{32}$$

$$|\tilde n_t|\cosh\beta_t\sin\theta_t\sin\Delta = -|\tilde n_i|\sin(\eta_i-\eta_t)\sinh\beta_i\sin\psi_i\cos\theta_i, \tag{33}$$

$$|\tilde n_t|(\sinh\beta_t)[\cos\psi_t\sin\theta_t\sin\Delta + \sin\psi_t\cos\theta_t\cos\Delta] = |\tilde n_i|[-\sin(\eta_i-\eta_t)\cosh\beta_i\sin\theta_i + \cos(\eta_i-\eta_t)\sinh\beta_i\sin\psi_i\cos\theta_i], \tag{34}$$

$$|\tilde n_t|(\sinh\beta_t)(\cos\psi_t\cos\Delta - \sin\psi_t\cos\theta_t\sin\Delta) = |\tilde n_i|\cos(\eta_i-\eta_t)\sinh\beta_i\cos\psi_i. \tag{35}$$

The generalized Snell–Descartes laws shown in Eqs. (29)–(31) for reflection and
in Eqs. (32)–(35) for refraction differ considerably from the well-known complex
Snell–Descartes laws:

$$\sin\tilde\theta_r = \sin\tilde\theta_i, \tag{36}$$

$$\tilde n_t\sin\tilde\theta_t = \tilde n_i\sin\tilde\theta_i, \tag{37}$$

$$\tilde\theta_j = \theta_j - i\beta_j, \qquad j = i, r, t. \tag{38}$$

These complex laws are recovered only in the special case in which

$$\psi_i = \frac{\pi}{2} \implies \psi_t = -\psi_r = \frac{\pi}{2}. \tag{39}$$

In vectorial form this condition may also be expressed as

$$\hat s\cdot(\mathbf{n}\wedge\mathbf{n}^*) = 0 \quad\text{or}\quad \hat s\cdot(\hat n_\parallel\wedge\hat n_\perp) = 0. \tag{40}$$

The simplified form of the Snell–Descartes laws in this case is due to the high
symmetry resulting from the invariance under lateral spatial translations. Another
special case would be given by $\hat s\cdot\hat n_\perp = 0$, but apparently it
does not lead to great simplifications.

A number of other consequences and special cases of the general laws (29)–(31)
and (32)–(35) are worth discussing:

**(a)** For the first time the deflection angle $\Delta$ appears, which is the
angle between the planes of incidence and refraction on each side of the
interface.

**(b)** The reflection laws (29)–(31) determine uniquely the alternative
reflection laws

$$\theta_r = \pi - \theta_i, \tag{41}$$

$$\psi_r = -\psi_i, \tag{42}$$

$$\beta_r = \beta_i. \tag{43}$$

The proof is as follows: One can substitute Eqs. (29) and (31) into Eq. (30) by
the use of trigonometric and hyperbolic sum rules. Solving the resulting
second-order equation by division (with the trivial solution $\sinh^2\beta_r =
\sinh^2\beta_i$) yields

$$\sinh^2\beta_r = \sinh^2\beta_i\sin^2\theta_i\cos^2\psi_i - \cos^2\psi_i. \tag{44}$$

Because all quantities are real this equation implies that $\beta_r \le \beta_i$.
But, because Eqs. (29)–(31) are symmetric with respect to the interchange
$i \leftrightarrow r$, a symmetric relation $\beta_i \le \beta_r$ also holds.
Therefore, $\beta_r$ is equal to $\beta_i$, which implies Eqs. (41)–(43). A last
remark is that the unusual law $\theta_r = \pi - \theta_i$ (instead of $\theta_r =
-\theta_i$) comes from the definition (23).

**(c)** The conditions for the occurrence of the unusual case $\Delta \ne 0$ are
investigated. From Eq. (33) one deduces that this unusual case occurs only when
the right-hand side of Eq. (33) is nonzero. Therefore, for purely inhomogeneous
waves, when Eq. (39) fails, and when the arguments of the complex refractive
indices of the two media are different, the deflection angle does not vanish.

**(d)** Finally, Eqs. (34) and (35) give a hint about the possible ways to
generate purely inhomogeneous waves from a homogeneous incident wave. If
$\beta_i = 0$, one can find with a little algebra that

$$|\tilde n_t|\sinh\beta_t\sin\psi_t\cos\theta_t = -|\tilde n_i|\sin(\eta_i-\eta_t)(\sin\theta_i\cos\Delta). \tag{45}$$

This means that generally $\beta_t \ne 0$ only if $(\eta_i - \eta_t) \ne 0$. The
special cases are when the transmitted wave is evanescent ($\theta_t = \pi/2$) or
when $\psi_t = 0$ or $\Delta = \pi/2$.

**Numerical solution of the transmitted wave (branch/selection statement).**
Numerically the direct solution of the generalized laws (32)–(35) in terms of the
parameters of the transmitted wave involves finding the unique solution of four
simultaneous nonlinear equations. However, the multivaluedness of the inverse
trigonometric functions and the many different conditions for signs and bounds of
the unknown variables makes direct solution complicated. The authors preferred to
obtain the transmitted parameters by disentangling the transmitted wave vector,
which can be simply calculated using the dispersion relation (5) and the
conservation of the tangential components (27). The deflection angle can then be
obtained by forming the scalar products involved in its definition (24). The
generalized laws (32)–(35) were checked to be satisfied by this numerical
solution.

## Dynamics: Generalized Polarizations and Generalized Fresnel Laws

At the interface the integral form of Maxwell's equations allows one to derive the
continuity or discontinuity equations for the fields between the two media:

$$\hat s\wedge\mathbf{E}_i + \hat s\wedge\mathbf{E}_r - \hat s\wedge\mathbf{E}_t = 0, \tag{46}$$

$$\mu_i\,\hat s\cdot(\mathbf{H}_i + \mathbf{H}_r) - \mu_t\,\hat s\cdot(\mathbf{H}_t) = 0, \tag{47}$$

$$\hat s\wedge\mathbf{H}_i + \hat s\wedge\mathbf{H}_r - \hat s\wedge\mathbf{H}_t = 0, \tag{48}$$

$$\epsilon_i\,\hat s\cdot(\mathbf{E}_i + \mathbf{E}_r) - \epsilon_t\,\hat s\cdot(\mathbf{E}_t) = \frac{\tilde\Sigma}{2}. \tag{49}$$

In the derivation of Eqs. (46)–(49) the fact is used that the surface current
density vanishes if the quantities $\epsilon$, $\mu$, and $\sigma$ remain constant
infinitely close to the interface. As a consequence, the parallel component of
$\mathbf{H}$ remains continuous. On the other hand the surface charge density
$\tilde\Sigma$, which was retained in Eq. (49), does not vanish whenever

$$\epsilon_i\sigma_t - \epsilon_t\sigma_i \ne 0. \tag{50}$$

This is easily checked by examining the equation for the conservation of charge and
Ohm's law on each side of the interface (see, for example, Ref. 32, p. 555).

Generalized TE and TM polarizations with respect to the interface can be defined
only by use of complex orthogonal unit vectors,

$$\mathbf{q}_1 = \frac{\hat s\wedge\mathbf{n}}{[(\hat s\wedge\mathbf{n})^2]^{1/2}}, \qquad \mathbf{q}_2 = \mathbf{n}\wedge\mathbf{q}_1, \tag{51}$$

satisfying $\mathbf{q}_1^2 = \mathbf{q}_2^2 = 1$, $\mathbf{q}_1\cdot\mathbf{q}_2 = 0$
(but $\mathbf{q}_1\cdot\mathbf{q}_2^* \ne 0$), and $\mathbf{q}_1 = (\mathbf{q}_2
\wedge\hat s)/(\hat s\cdot\mathbf{n})$. The projections $(E_{PM}, E_{PE})$ and
$(H_{PM}, H_{PE})$ of the electric and magnetic fields $\mathbf{E}_0$ and
$\mathbf{H}_0$ along $(\mathbf{q}_1, \mathbf{q}_2)$, respectively, are defined by

$$\mathbf{E}_{PE} = E_{PE}\mathbf{q}_1, \qquad \mathbf{H}_{PE} = H_{PE}\mathbf{q}_2, \tag{52}$$

$$\mathbf{E}_{PM} = E_{PM}\mathbf{q}_2, \qquad \mathbf{H}_{PM} = H_{PM}\mathbf{q}_1. \tag{53}$$

The two sets of modes satisfy independently Maxwell's equations (2)–(4) and are
not coupled by the interface. One can also write the complex projections as

$$E_{PE} = \mathbf{q}_1\cdot\mathbf{E}_0, \qquad H_{PE} = \mathbf{q}_2\cdot\mathbf{H}_0, \tag{54}$$

$$E_{PM} = \mathbf{q}_2\cdot\mathbf{E}_0, \qquad H_{PM} = \mathbf{q}_1\cdot\mathbf{H}_0. \tag{55}$$

Because these generalized modes possess a longitudinal component one cannot
meaningfully use the designation generalized TE modes. They are called instead
parallel electric (PE) modes because the polarization ellipse of the electric
field $\mathbf{E}_{PE} = E_{PE}\mathbf{q}_1$ always remains parallel to the
interface plane ($\hat s\cdot\mathbf{E}_{PE} = 0$). Similar considerations hold of
course for the parallel magnetic (PM) modes. Using Jones's parameterization one can
easily see that these modes are in general elliptically polarized, e.g.,

$$\mathbf{E}_{PE} = A_{PE}\exp(i\delta_{PE})[\mathbf{e}_{PE} + (-i\tanh\beta)(\hat n_\perp\cdot\mathbf{e}_{PE})\hat n_\parallel], \tag{56}$$

where

$$\mathbf{e}_{PE} = \frac{\hat e_1 + \chi_{PE}\hat e_2}{(1 + |\chi_{PE}|^2)^{1/2}}, \tag{57}$$

$$\chi_{PE} = \frac{i\tanh\beta\cos\theta\cos\psi}{\sin\theta - i\tanh\beta\cos\theta\sin\psi}. \tag{58}$$

Solving the boundary conditions (46)–(49) one derives the generalized Fresnel laws
for PE-polarized HIPW's:

$$E_{PEr} = r_{PE}E_{PEi} \quad\text{with}\quad r_{PE} = \left(\frac{\hat s\cdot\mathbf{k}_i/\mu_i - \hat s\cdot\mathbf{k}_t/\mu_t}{\hat s\cdot\mathbf{k}_i/\mu_i + \hat s\cdot\mathbf{k}_t/\mu_t}\right), \tag{59}$$

$$E_{PEt} = t_{PE}E_{PEi} \quad\text{with}\quad t_{PE} = \left(\frac{2\hat s\cdot\mathbf{k}_i/\mu_i}{\hat s\cdot\mathbf{k}_i/\mu_i + \hat s\cdot\mathbf{k}_t/\mu_t}\right). \tag{60}$$

For PM-polarized HIPW's,

$$H_{PMr} = r_{PM}H_{PMi} \quad\text{with}\quad r_{PM} = \left(\frac{\hat s\cdot\mathbf{k}_i/\bar\epsilon_i - \hat s\cdot\mathbf{k}_t/\bar\epsilon_t}{\hat s\cdot\mathbf{k}_i/\bar\epsilon_i + \hat s\cdot\mathbf{k}_t/\bar\epsilon_t}\right), \tag{61}$$

$$H_{PMt} = t_{PM}H_{PMi} \quad\text{with}\quad t_{PM} = \left(\frac{2\hat s\cdot\mathbf{k}_i/\bar\epsilon_i}{\hat s\cdot\mathbf{k}_i/\bar\epsilon_i + \hat s\cdot\mathbf{k}_t/\bar\epsilon_t}\right). \tag{62}$$

These laws are kept in condensed form because developing the complex scalar
products in terms of the parameterization is cumbersome. Because most often one
deals with electric fields, it is noted here that Eq. (61) holds also for $E_{PM}$
because $E_{PM} = -[\tilde n/(\omega\epsilon_0\mu)]\,H_{PM}$ [?], and that Eq. (62)
can be replaced by the following nonsymmetrical expression:

$$E_{PMt} = \left[\frac{2\tilde n_t(\hat s\cdot\mathbf{k}_i)}{\tilde n_i(\bar\epsilon_i/\bar\epsilon_t)(\hat s\cdot\mathbf{k}_i) + \tilde n_t(\hat s\cdot\mathbf{k}_t)}\right]E_{PMi}. \tag{63}$$

From Eq. (49) the surface charge density at the interface can be evaluated. It is
noted that PE waves do not generate a surface charge density. It has been checked
that the expressions of Fresnel coefficients for PM waves $r_{PM}$ and $t_{PM}$
were compatible with the continuity equation for charges at the interface cited
above [relation (50)].

To conclude this section it is mentioned first that the transmission of an
arbitrarily polarized wave can be calculated using the principle of linear
superposition in PE and PM modes with respect to this interface. Second, **the
Fresnel coefficients $r_{PE}$, $t_{PE}$, $r_{PM}$, $t_{PM}$ given by Eqs.
(59)–(62), respectively, also play an important role for the formulation of energy
balance and reciprocity relations, which are extensively discussed in a subsequent
paper.²⁹**

## Proposal for an Experiment

An experimental verification of the generalized laws can be made under various
circumstances. In Fig. 2 a case is proposed that could be realized in practice to
generate a HIPW and demonstrate a deflection. Equation (33) indicates that such a
deflection occurs only under restricted conditions, that is, when the right-hand
side of Eq. (33) is nonzero [e.g., at an interface between an absorbing and a
transparent medium, when Eq. (39) fails]. Therefore one might expect to observe a
deflection at the second interface of an absorbing wedge for particular
orientations of an incoming ordinary plane wave. The implicit assumption here is
that the lateral spatial extent of the incoming beam is so small that a
superposition of the multiply reflected waves in the wedge at the second interface
is avoided. In Fig. 2 the different media are air, rutile (TiO₂), indium tin oxide
(ITO:SnO₂/In₂O₃), and air, which have the following refractive indices at a
wavelength of 633 nm: 1, 2.6, 2, and 1, respectively. Only the ITO wedge is
slightly absorbing; its absorption coefficient $\alpha = -2k_0\,\mathrm{Im}(\tilde
n)$ is equal to 0.214 μm⁻¹. The role of the rutile–ITO interface is to create the
HIPW, and the role of the ITO–air interface (which is tilted with respect to the
normal out of the plane of the page) is to show the peculiarities of the
generalized laws. The calculated angles (angles of incidence, reflection,
refraction, deflection, etc.) are displayed in Table 1. Near normal incidence the
deflection angle $\Delta$ has an appreciable magnitude of 27°, and the wave has a
slightly elliptical polarization.

[Figure 2: Proposed experiment to test the generalized Snell–Descartes and Fresnel
laws in an ITO wedge on a rutile slab. A linearly polarized PE homogeneous plane
wave is incident upon the rutile slab, coupled onto the rutile–ITO interface where
the HIPW is generated. The orientation of the ITO–air final interface is such that
$\psi \ne \pi/2$; therefore only the generalized laws are valid for that interface.
The different angles and parameters are indicated in Table 1 and are not to scale
in the figure.]

**Table 1. Angles and Parameters of the HIPW and of Each of the Layers in the
Proposed Experiment** (angles in degrees; material parameters at 633 nm).

| Layer | $n$ | $\alpha$ (μm⁻¹) | $\theta_i$ | $\psi$ | $\beta$ | $\varphi$ | $\varepsilon$ | $\Delta$ |
|---|---|---|---|---|---|---|---|---|
| Rutile | 2.6 | 0 | 47 | 90 | 0 | 0 | 0 | 0 |
| ITO (w.r.t. first interface) | 2.0 | 0.214 | 71.92 | −90 | $1.57\times10^{-2}$ | 0 | 0 | 0 |
| ITO (w.r.t. second interface) | 2.0 | 0.214 | 0.01 | 22 | $1.57\times10^{-2}$ | 0 | 0 | 0 |
| Air | 1.0 | 0 | $1.86\times10^{-2}$ | −5.16 | $3.13\times10^{-2}$ | 27.16 | $2.4\times10^{-3}$ | 27.16 |

## Conclusions

A generalization of complex Snell–Descartes and Fresnel laws has been derived for
the most general kind of plane wave, that is, the HIPW. The polarization
properties of HIPW's and the form of generalized PE and PM modes have been
determined. The theory does not include anisotropic materials (birefringence) or
partial polarization states. It includes arbitrary values of all the wave
quantities specified by the parameterization angles $\theta_i$, $\theta_r$,
$\theta_t$, $\psi_i$, $\psi_r$, $\psi_t$, $\Delta$; the plane-wave inhomogeneities
$\beta_i$, $\beta_r$, $\beta_t$; and homogeneous isotropic dissipative media
characterized by the complex refractive indices $\tilde n_i$, $\tilde n_t$. **The
theory can describe numerous known phenomena, such as total internal reflection,
polarization by reflection (Brewster), metallic-type reflection, refraction of
metals, the skin effect, and surface waves at the interface.** New phenomena are
predicted by the present theory, which cannot be described by standard reflection
and refraction laws, when the incident wave includes an inhomogeneity in the
transverse direction [when Eq. (39) fails]. In particular a deflection of the HIPW
by the interface [characterized by the deflection angle $\Delta$, Eqs. (24)] is
predicted when the condition (50) is satisfied. These effects should be observable
in the limited region of space where edge diffraction effects that are due to the
finite lateral spatial extension of the incoming beam are negligible. To the best
of the authors' knowledge, these new phenomena do not seem to have been noticed
before (Refs. 27, 28, 33), because in all usual problems starting from ordinary
homogeneous plane waves that are incident upon an arbitrary set of parallel
interfaces Eq. (39) remains valid. An experiment has been proposed to demonstrate
the generalized laws and the deflection. Finally, it should also be mentioned that
once a HIPW is generated the standard Snell–Descartes and Fresnel laws fail also
at an interface between pure dielectrics even if the wave is not deflected. It is
not yet clear whether such effects could be usefully exploited in optical devices.
Thin-film technology may one day face layers with nonparallel interfaces, where the
phenomena discussed here might become important. Further generalizations for static
or moving interfaces between birefringent or optically active media should be
possible in the classical and the relativistic cases (Refs. 30, 1).

## Appendix A. Complex Vectorial Algebra

Because complex vector identities are not usually well known and play an important
role in this scheme, the main results for the two kinds of scalar product are
presented. Let $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$, $\mathbf{d}$ be four
arbitrary complex vectors. In this paper the scalar product and the vectorial
product are defined as

$$\mathbf{a}\cdot\mathbf{b} = \sum_{i=1}^{3} a_i b_i = \mathbf{b}\cdot\mathbf{a}, \tag{A1}$$

$$\mathbf{a}\wedge\mathbf{b} = \sum \epsilon_{ijk} a_j b_k = -\mathbf{b}\wedge\mathbf{a}, \tag{A2}$$

respectively, where $\epsilon_{ijk}$ is the fully antisymmetric tensor of rank two
[sic; rank three]. In this case all the usual vector identities that hold for real
vectors are valid:

$$\mathbf{a}\cdot(\mathbf{b}\wedge\mathbf{c}) = \mathbf{c}\cdot(\mathbf{a}\wedge\mathbf{b}) = \mathbf{b}\cdot(\mathbf{c}\wedge\mathbf{a}), \tag{A3}$$

$$\mathbf{a}\wedge(\mathbf{b}\wedge\mathbf{c}) = (\mathbf{a}\cdot\mathbf{c})\mathbf{b} - (\mathbf{a}\cdot\mathbf{b})\mathbf{c}, \tag{A4}$$

$$\mathbf{a}\wedge(\mathbf{b}\wedge\mathbf{c}) + \mathbf{b}\wedge(\mathbf{c}\wedge\mathbf{a}) + \mathbf{c}\wedge(\mathbf{a}\wedge\mathbf{b}) = 0, \tag{A5}$$

$$(\mathbf{a}\wedge\mathbf{b})\cdot(\mathbf{c}\wedge\mathbf{d}) = (\mathbf{c}\cdot\mathbf{a})(\mathbf{d}\cdot\mathbf{b}) - (\mathbf{d}\cdot\mathbf{a})(\mathbf{c}\cdot\mathbf{b}), \tag{A6}$$

$$(\mathbf{a}\wedge\mathbf{b})\wedge(\mathbf{c}\wedge\mathbf{d}) = [(\mathbf{c}\wedge\mathbf{d})\cdot\mathbf{a}]\mathbf{b} - [(\mathbf{c}\wedge\mathbf{d})\cdot\mathbf{b}]\mathbf{a} = [(\mathbf{a}\wedge\mathbf{b})\cdot\mathbf{d}]\mathbf{c} - [(\mathbf{a}\wedge\mathbf{b})\cdot\mathbf{c}]\mathbf{d}. \tag{A7}$$

A last very useful identity allows one to perform the decomposition of an arbitrary
real or complex vector $\mathbf{x}$ along three noncollinear vectors:

$$[(\mathbf{b}\wedge\mathbf{c})\cdot\mathbf{a}]\mathbf{x} = [(\mathbf{b}\wedge\mathbf{c})\cdot\mathbf{x}]\mathbf{a} + [(\mathbf{c}\wedge\mathbf{a})\cdot\mathbf{x}]\mathbf{b} + [(\mathbf{a}\wedge\mathbf{b})\cdot\mathbf{x}]\mathbf{c}. \tag{A8}$$

On the other hand the scalar product for complex vectors is sometimes defined as

$$\mathbf{a}\circ\mathbf{b} = \sum_{i=1}^{3} a_i^* b_i = \mathbf{a}^*\cdot\mathbf{b}, \tag{A9}$$

because the real norm of a complex vector is then properly given by the scalar
product of the vector with itself. This definition is somewhat awkward because most
identities (A3)–(A8) do not remain valid. A simple way that was found to solve this
problem is to redefine the vectorial product:

$$\mathbf{a}\overset{\wedge}{\circ}\mathbf{b} = \sum_{i=1}^{3}\epsilon_{ijk}a_j^* b_k^* = \mathbf{a}^*\wedge\mathbf{b}^* = -\mathbf{b}\overset{\wedge}{\circ}\mathbf{a}. \tag{A10}$$

The identities (A3)–(A8) can now be used. It must, however, be emphasized that the
order of the arguments in Eqs. (A3)–(A8) has been carefully chosen and must now be
absolutely respected because the scalar product is not commutative anymore:
$\mathbf{a}\circ\mathbf{b} \ne \mathbf{b}\circ\mathbf{a}$.

## Acknowledgments

The authors thank R. M. A. Azzam, F. K. Reinhart, and E. Wolf for useful
discussions; C. K. Carniglia, who provided an extensive bibliography on the subject
of refraction between two absorbing media and the mixed Poynting vector (Ref. 34);
and Ron Miles for critical readings of the manuscript.

## References

Full reference list omitted. Key entries relevant to GOAD work:

- Ref. 24: M. A. Dupertuis and M. Proctor, "Generalization of complex
  Snell–Descartes and Fresnel laws," *Opt. Photon. News* **1**(9), p. A-126
  (1990). (Prior conference note by the same authors.)
- Ref. 29 (the "subsequent paper" / companion paper on energy): M. A. Dupertuis,
  B. Acklin, and M. Proctor, "Generalized energy balance and reciprocity relations
  for thin-film optics," *J. Opt. Soc. Am. A* **11**, 1167–1174 (1994). (Immediately
  following article, same issue.)
- Ref. 34: C. K. Carniglia, "Reflection and transmission at a boundary between two
  absorbing media: a partial bibliography," presented at the Thin Films Technical
  Group meeting, Optical Society of America Annual Meeting, November 8, 1990,
  Boston, Mass.

---

## Collected uncertain glyph readings `[?]`

- **Eq. (44)** RHS `$\sinh^2\beta_i\sin^2\theta_i\cos^2\psi_i - \cos^2\psi_i$`:
  the exact grouping of factors and the trailing `$-\cos^2\psi_i$` term is at the
  limit of legibility at the rendered resolution. The physical argument in the
  text (that the equation forces $\beta_r \le \beta_i$, hence by symmetry
  $\beta_r = \beta_i$) is clear, but the precise RHS should be re-checked against a
  higher-resolution scan before being quoted numerically.
- **Eq. (56)**: the `+` joining `$\mathbf{e}_{PE}$` and the longitudinal term
  `$(-i\tanh\beta)(\hat n_\perp\cdot\mathbf{e}_{PE})\hat n_\parallel$` was inferred
  from the structure of Eq. (12); the printed bracket separator was hard to read.
- **Prefactor $E_{PM} = -[\tilde n/(\omega\epsilon_0\mu)]H_{PM}$** (unnumbered,
  page 5 right column): the exact denominator (whether it is
  $\omega\epsilon_0\mu$, $\omega\epsilon_0\bar\epsilon$, or similar) is uncertain
  and looks dimensionally odd as read; verify against a cleaner scan.
- **Eq. (A2)**: the text describes $\epsilon_{ijk}$ as "the fully antisymmetric
  tensor of rank two" — this is a printed error in the original (it is rank three);
  transcribed verbatim with a `[sic]`.
- Media subscripts on `$|\tilde n|$`, `$\beta$`, and `$\theta$` in Eqs. (32)–(35)
  (i vs t) were read carefully but a couple sit near the resolution limit; the
  physically-required assignment (LHS = transmitted quantities, RHS = incident)
  was used where a raw glyph was ambiguous.
