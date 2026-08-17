# Dupertuis, Acklin & Proctor — Generalized energy balance and reciprocity relations for thin-film optics (transcription)

> **Transcription for GOAD theory work.** This is a faithful markdown transcription of Dupertuis, Acklin & Proctor, "Generalized energy balance and reciprocity relations for thin-film optics," *J. Opt. Soc. Am. A* **11**(3), 1167–1174 (March 1994). It is the companion "paper II" to their generalized Snell–Fresnel paper (transcribed separately as `dupertuis_snell_fresnel.md`, referenced here as "paper I").
>
> Source file: `/Users/ixguard/Downloads/GeneralizedEnergyOSAAV11-394b.pdf`
> Transcription date: 2026-08-17
>
> The paper's own equation numbers are preserved so later work can cite "Dupertuis II eq. (N)". The reference list is omitted (numeric superscripts are kept as citation markers). Figure captions: this article contains no figures. Uncertain glyph readings are marked `[?]` inline and collected at the end.

---

# Generalized energy balance and reciprocity relations for thin-film optics

**M. A. Dupertuis**
*Institut de Micro- et Optoélectronique, Département de Physique, Ecole Polytechnique Fédérale de Lausanne, PHB-Ecublens, CH-1015 Lausanne, Switzerland*

**B. Acklin**
*AT&T Bell Labs, Room 4G-520, 101 Crawfords Corner Road, Holmdel, New Jersey 07733-3030*

**M. Proctor**
*Institut de Micro- et Optoélectronique, Département de Physique, Ecole Polytechnique Fédérale de Lausanne, PHB-Ecublens, CH-1015 Lausanne, Switzerland*

Received March 19, 1993; revised manuscript received August 13, 1993; accepted August 27, 1993

*Vol. 11, No. 3/March 1994/J. Opt. Soc. Am. A, pp. 1167–1174*

## Abstract

Energy balance and reciprocity relations are studied for harmonic inhomogeneous plane waves that are incident upon a stack of continuous absorbing dielectric media that are macroscopically characterized by their electric and magnetic permittivities and their conductivities. New cross terms between parallel electric and parallel magnetic modes are identified in the fully generalized Poynting vector. The symmetry and the relations between the general Fresnel coefficients are investigated in the context of energy balance at the interface. The contributions of the so-called mixed Poynting vector are discussed in detail. In particular a new transfer matrix is introduced for energy fluxes in thin-film optics based on the Poynting and mixed Poynting vectors. Finally, the study of reciprocity relations leads to a generalization of a theorem of reversibility for conducting and dielectric media.

---

## 1. INTRODUCTION

Generalized Snell–Descartes and Fresnel laws for arbitrary harmonic inhomogeneous plane waves (HIPW's) were derived in a preceding paper,¹ which we call paper I for simplicity. Previously published papers in this field during the twentieth century were briefly reviewed in paper I, and a tentative explanation of the surprising fact that such a discovery was so late was put forward in the Conclusion. The main purpose of the present paper is to illuminate the consequences of the general formulation of energy balance and reversibility for thin-film optics.

Energy balance is one of the oldest and most basic laws of physics. Its formulation in the frame of electromagnetic theory is from Poynting (1884). In the most simple case of a plane wave at normal incidence upon an interface between two purely dielectric media, Fresnel coefficients derived from Maxwell equations are real ($r_\pm, t_\pm \in \Re$) and satisfy the Stokes reversibility principle (1849):

$$r_+ r_+ + t_+ t_- = 1, \tag{1}$$

$$r_+ t_+ + t_+ r_- = 0, \tag{2}$$

$$r_- r_- + t_- t_+ = 1, \tag{3}$$

$$r_- t_- + t_- r_+ = 0. \tag{4}$$

Here and in the remainder of this paper the indices $+$ and $-$ denote quantities related to forward and backward waves with respect to the normal to the interface, respectively.

Stokes proposed a generalization of Eqs. (1)–(4):

$$r_+ r_+^* + t_+ t_-^* = 1, \tag{5}$$

$$r_+ t_+^* + t_+ r_-^* = 0, \tag{6}$$

$$r_- r_-^* + t_- t_+^* = 1, \tag{7}$$

$$r_- t_-^* + t_- r_+^* = 0, \tag{8}$$

but unfortunately complex Fresnel coefficients do not in general satisfy these relations. These relations look like an energy balance relation, but in 1962 Knittl outlined most beautifully how this connection was in fact loose.² The historical misconception of these relationships as energy balance relations remains in most textbooks; it is even argued often that energy balance can be deduced from the Stokes principle of reversibility of light beams. With the help of three examples, Knittl showed that

1. Reversibility does not unambiguously imply energy balance and, vice versa,
2. Fresnel coefficients cannot be deduced from reversibility even when completed with an additional condition of energy conservation.

Two of Knittl's conclusions are especially relevant: "The classical boundary conditions of the electromagnetic theory, with its double reality of **E** and **H** vectors, cannot be circumvented in deriving Fresnel coefficients," and "The energetic interpretation, if possible, is a matter of the terms inserted into the first reversibility equation, rather than of the structure of the equation itself." Chen and Nelson have, however, recently shown that a new wave-vector space method,³,⁴ reminiscent of quantum scattering theory, allows one to derive Fresnel coefficients without using boundary conditions. Nevertheless we largely agree with Knittl's statement that the electromagnetic theory cannot be circumvented in deriving Fresnel coefficients. Therefore the only path to a fully justified energetic interpretation is to start from Maxwell's equation and use the Poynting theorem.

We must, however, recognize that, in purely dielectric media, reversibility and energy balance often seem to be the same thing. For example, Ou and Mandel⁵ recently and conversely showed that, for the special case of a beam splitter, reciprocity relations could be derived from energy balance. Care should still be taken, even in the restricted purely dielectric case. For example, in total reflection the Fresnel coefficients satisfy both the energy balance relation²

$$r_+ r_+^* = 1 \tag{9}$$

and the reciprocity relations (1)–(4) with $r_\pm, t_\pm$ complex.⁶ Obviously Eq. (9) strongly differs from Eqs. (1)–(4) with complex coefficients and from Eqs. (5)–(8).

Similar nontrivial considerations apply to absorbing media: it is known that at normal incidence upon a simple dielectric–metal boundary relations (1)–(4) are satisfied instead of Eqs. (5)–(8). Another important conclusion from Knittl's research is that reciprocity relations could not be found for metal–metal interfaces.⁷ Soon thereafter, Šantavý showed directly from the symmetry of Maxwell's equations completed with a phenomenological conductivity term that some kind of reciprocity relations still hold in a generalized sense when the operation ħ [?], representing the change of a sign of the conductivities in both media, is introduced.⁶ Šantavý was able to unify all the special cases of reciprocity relations derived earlier in a common framework. Although we show below that Šantavý's reversibility theorem is not fully general, his research is fundamental as it takes a new and correct approach to the problem. This was recognized by Knittl, who also underlined some of its consequences and showed how to use these new reciprocity relations.²

During this century another problem linked with energy balance at an interface between an absorbing and a dielectric medium has been recurring. When a wave is incident from the absorbing side, a new nonintuitive term, called the mixed Poynting vector, appears in the energy balance. This term is due to absorption of the standing-wave interference between the incident and reflected waves. The presence of this term is one of the best examples of the failure of the energetic interpretation of the first reversibility equation (1). Born and Ladenburg⁸ discovered this effect in 1911, but it was only in 1950 that the effect was fully interpreted.⁷,⁹ Since then, a number of researchers have still disregarded this term or treated it incorrectly. It is important only in the so-called strongly absorbing media [see inequality (70) below]. A number of related papers are reviewed by Fragstein.¹⁰

We also discuss situations in which the mixed Poynting vector becomes important, and we present this vector in a unified general framework by the use of the coordinate-free vector formalism of paper I. Thus we hope to fulfill Fragstein's wish of a modern final presentation¹⁰ of the mixed Poynting vector. In addition we show the relevance of this concept by outlining a complex power transfer matrix formalism for thin-film optics that might be useful for researchers who desire to work in terms of power fluxes.

The above considerations on energy balance and reversibility justify the approach taken in this paper to study each issue from first principles, that is, Maxwell's equations. The connection between them will then become clearer and will be outlined when appropriate. As a result of our general approach new features will also appear, however.

In Section 2 we study the Poynting vector of an arbitrary HIPW. The problem of energy balance at an interface is settled in Section 3 and leads in Section 4 to the discovery that the Poynting vector is not necessarily separable into pure parallel electric (PE) and parallel magnetic (PM) contributions. In Sections 5 and 6 we study power reflection and transmission, as well as time-averaged intensity, and derive identities for generalized Fresnel coefficients that can be interpreted energetically. Complex power transmission coefficients are naturally suggested in this approach, and in Section 7 they are further developed into a complex power transfer matrix technique for studying energy fluxes in thin-film optics. The mixed Poynting vector appears naturally in this context. In Sections 8 and 9 we study reversibility and give the first fully correct formulation of the general theorem of reversibility for thin-film optics in Section 9. In particular Fresnel coefficients satisfy this theorem and illustrate a number of special cases.

## 2. POYNTING VECTOR FOR HIPW's

Let us briefly recall Poynting's theorem for harmonic fields in a complex form. The spatially dependent harmonic electric and magnetic fields are expanded as follows:

$$\mathbf{E}(\mathbf{r}, t) = \tfrac{1}{2}[\mathbf{E}(\mathbf{r})\exp(i\omega t) + \text{c.c.}], \tag{10}$$

$$\mathbf{H}(\mathbf{r}, t) = \tfrac{1}{2}[\mathbf{H}(\mathbf{r})\exp(i\omega t) + \text{c.c.}], \tag{11}$$

where $\omega$ is real. The time average over one period of product vectors that are quadratic in the fields such as the Poynting vector can be written as

$$\langle \mathbf{S}(\mathbf{r}, t)\rangle = \langle \mathbf{E}(\mathbf{r}, t) \wedge \mathbf{H}(\mathbf{r}, t)\rangle = \tfrac{1}{2}\text{Re}[\mathbf{E}(\mathbf{r}) \wedge \mathbf{H}^*(\mathbf{r})] = \text{Re}[\mathbf{S}(\mathbf{r})]. \tag{12}$$

$\mathbf{S}(\mathbf{r})$ is called the **complex Poynting vector**.

The time-averaged power that is dissipated is proportional to $\int_V \mathbf{J}(\mathbf{r}, t)\cdot \mathbf{E}(\mathbf{r}, t)\,\mathrm{d}^3 x$, where $\mathbf{J}(\mathbf{r}, t)$ is the current density. The complex version of the Poynting theorem can be derived in local form from Maxwell's equations and the Gauss theorem:

$$\nabla \cdot \mathbf{S}(\mathbf{r}) + \tfrac{1}{2}\mathbf{J}^*(\mathbf{r})\cdot\mathbf{E}(\mathbf{r}) - 2i\omega[w_e(\mathbf{r}) - w_m(\mathbf{r})] = 0, \tag{13}$$

where the electric and magnetic energy densities are expressed as

$$w_e(\mathbf{r}) = \tfrac{1}{4}\mathbf{E}(\mathbf{r})\cdot\mathbf{D}^*(\mathbf{r}), \tag{14}$$

$$w_m(\mathbf{r}) = \tfrac{1}{4}\mathbf{B}(\mathbf{r})\cdot\mathbf{H}^*(\mathbf{r}), \tag{15}$$

respectively. The real part of Eq. (13) expresses the conservation of energy for time-averaged quantities, and the imaginary part relates to the reactive stored energy and its alternating flow. The second term represents the power dissipated per unit volume.

As in paper I we concentrate on the case of HIPW's and homogeneous media characterized by constant electric and magnetic permittivities ($\epsilon, \mu$) and constant conductivity $\sigma$. We can therefore write

$$\mathbf{E}(\mathbf{r}) = \mathbf{E}\exp(-i\mathbf{k}\cdot\mathbf{r}), \tag{16}$$

$$\mathbf{D}(\mathbf{r}) = \epsilon_0\bar\epsilon\,\mathbf{E}(\mathbf{r}), \tag{17}$$

$$\mathbf{H}(\mathbf{r}) = \mathbf{H}\exp(-i\mathbf{k}\cdot\mathbf{r}), \tag{18}$$

$$\mathbf{B}(\mathbf{r}) = \mu_0\mu\,\mathbf{H}(\mathbf{r}), \tag{19}$$

where $\mathbf{k}$ is the complex wave vector

$$\mathbf{k} = \bar{n}k_0\mathbf{n}, \tag{20}$$

which satisfies the standard dispersion relation for HIPW, which was derived in Eq. (5) of paper I. In Eq. (20) $\bar{n} = \sqrt{\mu\bar\epsilon}$ is the complex refractive index, $\bar\epsilon = \epsilon - i\sigma/(\omega\epsilon_0)$ is the complex dielectric constant, $k_0 = \omega/c$ is the vacuum wave vector, and $\mathbf{n}$ is a complex unit vector satisfying $\mathbf{n}\cdot\mathbf{n} = 1$. The complex Poynting vector for HIPW now reads

$$\mathbf{S}(\mathbf{r}) = \mathbf{S}\exp[2\,\text{Im}(\mathbf{k})\cdot\mathbf{r}] = \frac{1}{2\mu_0\mu}[|\mathbf{E}|^2\mathbf{k}^* - (\mathbf{E}\cdot\mathbf{k}^*)\mathbf{E}^*]\exp[2\,\text{Im}(\mathbf{k})\cdot\mathbf{r}], \tag{21}$$

where $\mathbf{S} = \tfrac{1}{2}\mathbf{E}\wedge\mathbf{H}^*$. The second term is clearly due to the inhomogeneity [see Eq. (9) of paper I]. The instantaneous and time-averaged electric and magnetic energy densities, Eqs. (14) and (15), can be reexpressed as

$$w_e = \tfrac{1}{4}\epsilon_0\bar\epsilon|\mathbf{E}|^2\exp[2\,\text{Im}(\mathbf{k})\cdot\mathbf{r}], \tag{22}$$

$$w_m = \tfrac{1}{4}\mu_0\mu|\mathbf{H}|^2\exp[2\,\text{Im}(\mathbf{k})\cdot\mathbf{r}]. \tag{23}$$

They are in general neither independent of spatial coordinates nor equal¹¹ because for a general HIPW we can write

$$|\mathbf{H}|^2 = \frac{2\mathbf{k}\cdot\mathbf{S}}{\omega\mu_0\mu} = \frac{1}{\omega\mu_0\mu}(|\mathbf{k}|^2|\mathbf{E}|^2 - |\mathbf{k}\cdot\mathbf{E}^*|^2). \tag{24}$$

A symmetric expression for $|\mathbf{E}|^2$ can be found by the use of the relation $|\mathbf{E}|^2 = 2\mathbf{k}^*\cdot\mathbf{S}/(\omega\epsilon_0\bar\epsilon^*)$. Let us now express again the Poynting theorem for HIPW's in local form:

$$\nabla\cdot\mathbf{S}(\mathbf{r}) = \frac{i\omega}{2}(\epsilon_0\bar\epsilon^*|\mathbf{E}|^2 - \mu_0\mu|\mathbf{H}|^2)\exp[2\,\text{Im}(\mathbf{k})\cdot\mathbf{r}], \tag{25}$$

where we see clearly that the real part depends only on the dissipation (the imaginary part of $\bar\epsilon$).

To complete the review of HIPW's carried out in paper I we use Eq. (21) to derive two additional properties linked with the complex Poynting vector:

1. The vector

$$\mathbf{k}\wedge\mathbf{S} = \frac{1}{2\mu_0\mu}[|\mathbf{E}|^2(\mathbf{k}\wedge\mathbf{k}^*) - (\mathbf{E}\cdot\mathbf{k}^*)(\mathbf{k}\wedge\mathbf{E}^*)], \tag{26}$$

by virtue of two typical properties of HIPW's [Eqs. (8) and (9), paper I], vanishes for a homogeneous wave. Conversely, because $\mathbf{k}\wedge\mathbf{k}^*$ and $\mathbf{k}\wedge\mathbf{E}^*$ are never parallel for purely inhomogeneous waves we deduce that $\mathbf{k}\wedge\mathbf{S}$ vanishes if, and only if, the wave is homogeneous. The same can be shown for $\mathbf{k}^*\wedge\mathbf{S} = 0$. Thus the planes of constant phase or constant amplitude characterized by $(\mathbf{k}\pm\mathbf{k}^*)\cdot\mathbf{r} = 0$ are perpendicular to the complex Poynting vector **S** if and only if the wave is purely inhomogeneous.

2. The complex Poynting vector is complex perpendicular to **E** and $\mathbf{H}^*$ by construction. It is easy to show that it is complex perpendicular to $\mathbf{E}^*$ and **H** only if the polarization is circular.

3. This last property is not incompatible with $\mathbf{S}(\mathbf{r}, t)$ always being perpendicular to $\mathbf{E}(\mathbf{r}, t)$ and $\mathbf{H}(\mathbf{r}, t)$ (by construction).

Having these concepts well in hand, let us turn to the detailed study of the energy flow at an interface between two homogeneous media that are characterized by their permittivities and conductivities. As in paper I, the incident wave is denoted by the subscript $i$, the reflected wave is denoted by the subscript $r$, and the transmitted wave is denoted by the subscript $t$; the normal to the interface is denoted by $\hat{s}$; the quantities related to the first and the second media are indexed by $(i, t)$, except in Eq. (27), where we use $(L, R)$ because of possible confusion.

## 3. ENERGY BALANCE AT AN INTERFACE

The proper energy balance relations at an interface are derived from the integral form of the Poynting theorem (13), using a small Gaussian pillbox surrounding the interface. The perpendicular component of the Poynting vector is conserved on each side of the interface:

$$S_L^\perp(\mathbf{r}) = \hat{s}\cdot\mathbf{S}_L(\mathbf{r}) = \hat{s}\cdot\mathbf{S}_R(\mathbf{r}) = S_R^\perp(\mathbf{r}). \tag{27}$$

Here we indexed each side by $(L, R)$ to avoid confusion with the Poynting vectors for the incident and transmitted waves, which we need to introduce below. The apparent simplicity of the statement (27) hides far-reaching consequences, which are discussed below. First, the interface allows one to separate arbitrary fields into PE and PM polarization components (see paper I). For the Poynting vector, however, a new cross term has to be introduced. Second, the mixed Poynting vector contributes at the interface and needs to be reviewed carefully. The general formulation of the energy balance can be viewed as a condition on the Fresnel coefficients, which is close to (but distinct from) the condition proposed by Stokes,² which was based on reversibility. A detailed comparison is possible after the discussion of the consequences of reversibility carried out in Section 9.

## 4. SEPARABILITY OF THE POYNTING VECTOR INTO PE AND PM COMPONENTS

At a planar interface an arbitrary HIPW can be decomposed into PE and PM components (see paper I). The uniquely defined components satisfy

$$\mathbf{E}_\text{PE} = E_\text{PE}\mathbf{q}_1, \qquad \mathbf{H}_\text{PE} = H_\text{PE}\mathbf{q}_2, \tag{28}$$

$$\mathbf{E}_\text{PM} = E_\text{PM}\mathbf{q}_2, \qquad \mathbf{H}_\text{PM} = H_\text{PM}\mathbf{q}_1, \tag{29}$$

where the complex unit orthogonal orientation vectors $\mathbf{q}_1, \mathbf{q}_2$ were discussed in paper I:

$$\mathbf{q}_1 = \frac{\hat{s}\wedge\mathbf{n}}{[(\hat{s}\wedge\mathbf{n})^2]^{1/2}}, \qquad \mathbf{q}_2 = \mathbf{n}\wedge\mathbf{q}_1. \tag{30}$$

The electric-field and magnetic-field amplitudes can be obtained simply by projection:

$$E_\text{PE} = \mathbf{q}_1\cdot\mathbf{E}, \qquad H_\text{PE} = \mathbf{q}_2\cdot\mathbf{H}, \tag{31}$$

$$E_\text{PM} = \mathbf{q}_2\cdot\mathbf{E}, \qquad H_\text{PM} = \mathbf{q}_1\cdot\mathbf{H}, \tag{32}$$

and Maxwell's equations imply the following relation between the amplitudes,

$$E_\text{PM} = -\frac{\bar{n}}{\bar\epsilon}\sqrt{\frac{\mu_0}{\epsilon_0}}H_\text{PM}, \tag{33}$$

$$H_\text{PE} = \frac{\bar{n}}{\mu}\sqrt{\frac{\epsilon_0}{\mu_0}}E_\text{PE}. \tag{34}$$

The decomposition of the Poynting vector is then

$$\mathbf{S} = \tfrac{1}{2}(\mathbf{E}\wedge\mathbf{H}^*) = \mathbf{S}_\text{PE} + \mathbf{S}_\text{PM} + \mathbf{S}_X, \tag{35}$$

where

$$\mathbf{S}_\text{PE} = \tfrac{1}{2}E_\text{PE}H_\text{PE}^*(\mathbf{q}_1\wedge\mathbf{q}_2^*), \tag{36}$$

$$\mathbf{S}_\text{PM} = \tfrac{1}{2}E_\text{PM}H_\text{PM}^*(\mathbf{q}_2\wedge\mathbf{q}_1^*), \tag{37}$$

$$\mathbf{S}_X = \tfrac{1}{2}[E_\text{PE}H_\text{PM}^*(\mathbf{q}_1\wedge\mathbf{q}_1^*) + E_\text{PM}H_\text{PE}^*(\mathbf{q}_2\wedge\mathbf{q}_2^*)]. \tag{38}$$

The first two terms represent the standard PE and PM energetic contributions. The last term $\mathbf{S}_X$ is a new cross energy flux term between the PE and PM parts of the HIPW with respect to the given interface. The spatial dependence of $\mathbf{S}_X(\mathbf{r})$ is the same as for the other terms. Let us investigate when this new term vanishes. Using complex vectorial identities (see Appendix A of paper I), we can show that

$$\mathbf{q}_1\wedge\mathbf{q}_1^* = \frac{(\mathbf{n}\wedge\mathbf{n}^*)\cdot\hat{s}}{|(\hat{s}\wedge\mathbf{n})^2|}\hat{s}, \tag{39}$$

$$\mathbf{q}_2\wedge\mathbf{q}_2^* = [(\mathbf{n}\wedge\mathbf{n}^*)\cdot\mathbf{q}_1^*]\mathbf{q}_1 + [(\mathbf{q}_1\wedge\mathbf{q}_1^*)\cdot\mathbf{n}^*]\mathbf{n}. \tag{40}$$

We have shown in paper I that the standard complex Snell–Descartes and Fresnel laws are recovered if the arbitrarily polarized HIPW satisfies the condition

$$\hat{s}\cdot(\mathbf{n}\wedge\mathbf{n}^*) = 0. \tag{41}$$

Therefore, in terms of $\mathbf{q}_1$ and $\mathbf{q}_2$, we obtain in this case

$$\mathbf{q}_1\wedge\mathbf{q}_1^* = 0, \qquad \mathbf{q}_2\wedge\mathbf{q}_2^* = 0. \tag{42}$$

Therefore the condition (41) is a sufficient condition for a vanishing cross Poynting vector $\mathbf{S}_X$. A closer look reveals that conditions (41) or (42) are necessary and sufficient conditions. Condition (41) also implies the interesting identity $\mathbf{q}_1\cdot\mathbf{n}^* = 0$.

## 5. POWER REFLECTION AND TRANSMISSION AND THE MIXED POYNTING VECTOR

The energy balance equation (27) may now be disentangled in terms of incident, reflected, and transmitted HIPW waves:

$$S_i^\perp + S_r^\perp + S_m^\perp = S_t^\perp, \tag{43}$$

where we have factored out the spatial dependence along the interface by the use of the generalized Snell–Descartes laws discussed in paper I. In addition to the standard energy flux terms,

$$\mathbf{S}_j = \tfrac{1}{2}(\mathbf{E}_j\wedge\mathbf{H}_j^*), \qquad j = i, r, t, \tag{44}$$

Equation (43) contains the mixed Poynting vector

$$\mathbf{S}_m^\perp = \tfrac{1}{2}\hat{s}\cdot(\mathbf{E}_i\wedge\mathbf{H}_r^* + \mathbf{E}_r\wedge\mathbf{H}_i^*). \tag{45}$$

This term originates from the standing-wave pattern between the incident and the reflected waves in a dissipative medium. The mixed Poynting vector vanishes under specific conditions that we examine in detail once we have expressed the complex Poynting vector in terms of the Fresnel coefficients and the incident energy flux.

## 6. INTENSITY REFLECTION AND TRANSMISSION COEFFICIENTS AND FRESNEL COEFFICIENTS

We may be tempted to define the intensity reflection and transmission coefficients as

$$S_r^\perp = \bar{\mathcal{R}}S_i^\perp, \tag{46}$$

$$S_t^\perp = \bar{\mathcal{T}}S_i^\perp. \tag{47}$$

However, because $S^\perp$ is in general complex, we show that $\bar{\mathcal{R}}$ and $\bar{\mathcal{T}}$ must also be defined as complex (although in many cases $\bar{\mathcal{R}}$ remains real), and we call $\bar{\mathcal{R}}$ and $\bar{\mathcal{T}}$ the complex reflection and transmission coefficients for the energy flow, respectively. The complex energy balance equation (27) is therefore not fully appropriate to define the intensity reflection and transmission coefficients. Instead we must adopt the usual definition for the reference incident intensity, as the time-averaged intensity,

$$\phi_i = \text{Re}(S_i^\perp) = \frac{1}{2\omega\mu_0\mu_i}\text{Re}[|\mathbf{E}_i|^2\mathbf{k}_i^\perp - (\mathbf{E}_i\cdot\mathbf{k}_i^\perp)E_i^{*\perp}[?]], \tag{48}$$

and define the real intensity reflection and transmission coefficients, $\mathcal{R}$ and $\mathcal{T}$, respectively, by

$$\phi_r = \text{Re}(S_r^\perp) = -\mathcal{R}\phi_i, \tag{49}$$

$$\phi_t = \text{Re}(S_t^\perp) = \text{Re}(\bar{\mathcal{T}}S_i^\perp) = \mathcal{T}\phi_i. \tag{50}$$

This is consistent with the fact that only time-averaged quantities are observable. The mixed Poynting vector also obeys an equation of the type

$$S_m^\perp = \bar{\mathcal{M}}S_i^\perp, \tag{51}$$

where $\bar{\mathcal{M}}$ is complex and gives rise to an interference flux,

$$\phi_m = \text{Re}(S_m^\perp) = \text{Re}(\bar{\mathcal{M}}S_i^\perp) = \mathcal{M}\phi_i, \tag{52}$$

which can be positive or negative and affects the flux balance:

$$\phi_i + \phi_r + \phi_m = \phi_t. \tag{53}$$

The physical reality of the interference flux $\phi_m$ has been a subject of controversy that has been discussed by many authors (Ref. 10 and references therein) because this term may be positive and looks like a source of energy at the interface, seemingly contradicting our intuition of energy conservation in a dissipative medium. In agreement with Fragstein, we advocate the inclusion of the interference part of the flux in the incident flux. This forms a new interpretation of the left-hand side of Eq. (53) that is compatible with our intuition of causality. As expected, the total Poynting vector never leads to a positive divergence of the total energy flow,¹² and thus no breakdown of orthodox concepts occurs. Despite this a few authors committed mistakes of interpretation and sometimes even provided solutions that were not compatible with Maxwell's equations at the boundary (in particular Vasicek¹³). Clearly the interference flux is related to an interference between the forward and backward waves in the first medium, and we analyze this in more detail below.

Using equations (59)–(62) in paper I as well as the definitions (46) and (51), we have derived the different components of the Poynting vector of an arbitrarily polarized HIPW in terms of the Fresnel coefficients. For a PE wave we obtain

$$\bar{\mathcal{R}}_\text{PE} = \mathcal{R}_\text{PE} = |r_\text{PE}|^2, \tag{54}$$

$$\bar{\mathcal{T}}_\text{PE} = |t_\text{PE}|^2\frac{\hat{s}\cdot\mathbf{k}_t^*/\mu_t}{\hat{s}\cdot\mathbf{k}_i^*/\mu_i}, \tag{55}$$

as well as the purely imaginary mixed Poynting vector contribution

$$\bar{\mathcal{M}}_\text{PE} = r_\text{PE} - r_\text{PE}^* = 2i\,\text{Im}(r_\text{PE}). \tag{56}$$

The corresponding relationships for PM waves are obtained by replacing $\mu$ by $\bar\epsilon^*$. The PE–PM cross-term coefficients give rise to

$$\bar{\mathcal{R}}_X = -r_\text{PE}^*r_\text{PM}, \tag{57}$$

$$\bar{\mathcal{T}}_X = t_\text{PE}^*t_\text{PM}\frac{|\hat{s}\cdot\mathbf{k}_t|^2}{|\hat{s}\cdot\mathbf{k}_i|^2}\frac{\mu_i\bar\epsilon_i}{\mu_t\bar\epsilon_t}, \tag{58}$$

$$\bar{\mathcal{M}}_X = -(r_\text{PE}^* + r_\text{PM}). \tag{59}$$

The complex energy balance condition at the interface can then be written separately for the PE-, PM-, or $X$-term components of the Poynting vector as

$$\bar{\mathcal{R}}_J + \bar{\mathcal{T}}_J = 1 + \bar{\mathcal{M}}_J, \qquad J = \text{PE}, \text{PM}, X. \tag{60}$$

These relationships can be understood as identities for generalized Fresnel coefficients or for complex reflection and transmission coefficients.

Some comments are necessary. We first discuss the mixed Poynting vector contributions in detail and then the new cross terms.

Equation (60) differs from the Stokes relation (5) by the contribution of the mixed Poynting vector $\bar{\mathcal{M}}J$. Let us analyze the conditions under which this term vanishes. The conditions differ for the cases of PE, PM, or arbitrarily polarized modes. For PE waves we obtain easily

$$(\hat{s}\cdot\mathbf{k}_i)(\hat{s}\cdot\mathbf{k}_r^*) = (\hat{s}\cdot\mathbf{k}_i^*)(\hat{s}\cdot\mathbf{k}_r), \tag{61}$$

which is a generalization of known identities.¹⁴ For PM waves Eq. (61) is replaced by

$$\frac{\hat{s}\cdot\mathbf{k}_i}{\bar\epsilon_i}\frac{\hat{s}\cdot\mathbf{k}_r^*}{\bar\epsilon_r^*} = \frac{\hat{s}\cdot\mathbf{k}_i^*}{\bar\epsilon_i^*}\frac{\hat{s}\cdot\mathbf{k}_t}{\bar\epsilon_t} \tag{62}$$

and reduces to the previous case for normal incidence and homogeneous waves. If standard Fresnel–Descartes laws are valid [cf. condition (41)] we can write $\hat{s}\cdot\mathbf{k}_i = \bar{n}_i\cos(\bar\theta_i)$, where $\bar\theta_i$ is the imaginary angle defined by equation (39) in paper I. To illustrate the physical interpretation of the mixed Poynting vector, and the so-called interference flux, let us consider the standard case of normal incidence and homogeneous waves. In this case $\cos(\bar\theta_i) = \cos(\bar\theta_t) = 1$ and we get $\bar{n}_i\bar{n}_i^* = \bar{n}_t\bar{n}_t^*$. If the arguments of the complex refractive indices ($\eta_i$ and $\eta_t$ in paper I) are equal, this is equivalent to

$$n_i\kappa_i = n_t\kappa_t, \tag{63}$$

when $n_i$, $\kappa_i$ and $n_t$, $\kappa_t$ denote the real and imaginary parts of the complex refractive indices, respectively. Then the reflection and transmission coefficients for fluxes (assuming $\mu_i = \mu_t$ for simplicity) become

$$\mathcal{R} = \left|\frac{\bar{n}_i - \bar{n}_t}{\bar{n}_i + \bar{n}_t}\right|^2 = \frac{(n_i - n_t)^2 + (\kappa_i - \kappa_t)^2}{(n_i + n_t)^2 + (\kappa_i + \kappa_t)^2}, \tag{64}$$

$$\mathcal{T} = \frac{1}{n_i}\text{Re}(\bar{\mathcal{T}}\bar{n}_i^*) = \frac{1}{n_i}\text{Re}\left(\frac{4\bar{n}_i\bar{n}_t^*}{|\bar{n}_i + \bar{n}_t|^2}\bar{n}_i^*\right) = \frac{4n_in_t\left(1 + \dfrac{\kappa_i^2}{n_i^2}\right)}{(n_i + n_t)^2 + (\kappa_i + \kappa_t)^2}, \tag{65}$$

$$\mathcal{M} = \frac{4\kappa_i(n_i\kappa_t - n_t\kappa_i)}{n_i[(n_i + n_t)^2 + (\kappa_i + \kappa_t)^2]}. \tag{66}$$

We see that $\mathcal{M}$ vanishes not only when Eq. (63) is satisfied but also when $\kappa_i = 0$. Therefore in this important case the interference flux vanishes when the mixed Poynting vector does not.

The physical origin of the interference flux is, however, best understood when one examines the form of the mixed Poynting vector inside the medium of incidence. From the general form

$$\begin{aligned}
S_m^\perp(\mathbf{r}) &= \tfrac{1}{2}(\mathbf{E}_i\wedge\mathbf{H}_r^*)\exp[-i(\mathbf{k}_i - \mathbf{k}_r^*)\cdot\mathbf{r}] + \tfrac{1}{2}(\mathbf{E}_r\wedge\mathbf{H}_i^*)\exp[-i(\mathbf{k}_r - \mathbf{k}_i^*)\cdot\mathbf{r}] \\
&= \tfrac{1}{2}(\mathbf{E}_i\wedge\mathbf{H}_r^*)\exp\{-2i\,\text{Re}(\hat{s}\cdot\mathbf{k}_i)(\hat{s}\cdot\mathbf{r}) - 2\,\text{Im}[\hat{s}\wedge(\hat{s}\wedge\mathbf{k}_i)]\cdot\mathbf{r}\} \\
&\quad + \tfrac{1}{2}(\mathbf{E}_r\wedge\mathbf{H}_i^*)\exp\{+2i\,\text{Re}(\hat{s}\cdot\mathbf{k}_i)(\hat{s}\cdot\mathbf{r}) - 2\,\text{Im}[\hat{s}\wedge(\hat{s}\wedge\mathbf{k}_i)]\cdot\mathbf{r}\},
\end{aligned} \tag{67, 68}$$

we find that

$$S_m^\perp(\mathbf{r}) = \tfrac{1}{2}\{r_\text{PE}\exp[-2in_ik_0(\hat{s}\cdot\mathbf{r})] - \text{c.c.}\} \tag{69}$$

under our simplifying hypotheses. Obviously Eq. (69) describes a dephased standing-wave-like pattern with a vanishing spatial average. $S_m^\perp(\mathbf{r})$ is purely imaginary and describes the redistribution of energy flux resulting from the spatially varying energy dissipation in the standing-wave pattern. Depending on the phase of $r_\text{PE}$, this term gives a negative or positive contribution at the interface that is precisely the mixed Poynting vector. This term has a purely local character, and its contribution is usually small because it occurs at a scale small compared with the wavelength in the material. This term becomes important in the case of a strongly absorbing medium¹⁴,¹⁵ characterized by ($\mu_i = \mu_t$ for simplicity):

$$\frac{\kappa_i}{n_i} > \frac{\kappa_t}{n_t} > 0.1. \tag{70}$$

In Section 6 we have presented general coefficients for the transmission of energy flow. We have shown that in the general case there exists a new cross PE–PM term that obeys a separate transmission law. We have emphasized the distinction between two types of reflection, transmission, and mixed Poynting vector coefficients: (a) complex coefficients for the energy flow (Poynting vector) and (b) real coefficients for time-averaged intensities (the only physically observable quantity).

## 7. ENERGY BALANCE IN A STACK OF PARALLEL THIN ABSORBING FILMS

Because the optical properties of thin films become increasingly important in modern optical technology, we illustrate briefly how the complex energy flux transfer coefficients and the mixed Poynting vector concepts can be used to study energy balance in a stack of parallel thin films.

We assume that the field in each layer is given by

$$\mathbf{E}^{(n)}(\mathbf{r}) = [\mathbf{E}_{\text{PE}+}^{(n)} + \mathbf{E}_{\text{PM}+}^{(n)}]\exp[-i\mathbf{k}_+^{(n)}\cdot(\mathbf{r} - \mathbf{r}_{n-1})] + [\mathbf{E}_{\text{PE}-}^{(n)} + \mathbf{E}_{\text{PM}-}^{(n)}]\exp[-i\mathbf{k}_-^{(n)}\cdot(\mathbf{r} - \mathbf{r}_{n-1})], \tag{71}$$

$$\mathbf{H}^{(n)}(\mathbf{r}) = [\mathbf{H}_{\text{PE}+}^{(n)} + \mathbf{H}_{\text{PM}+}^{(n)}]\exp[-i\mathbf{k}_+^{(n)}\cdot(\mathbf{r} - \mathbf{r}_{n-1})] + [\mathbf{H}_{\text{PE}-}^{(n)} + \mathbf{H}_{\text{PM}-}^{(n)}]\exp[-i\mathbf{k}_-^{(n)}\cdot(\mathbf{r} - \mathbf{r}_{n-1})]. \tag{72}$$

The layer index $n$ runs from 0 to $N + 1$ for a set of $N$ layers sandwiched between arbitrary semi-infinite media. In other cases, such as an infinite Bragg dielectric mirror, the index might run from $-\infty$ to $+\infty$. In each layer the material constants $\bar\epsilon^{(n)}, \mu^{(n)}$ are known, but the wave vector and the field amplitudes must be calculated.

If the forward wave vector $\mathbf{k}^{(m)}$ is known in the $m$th layer, the forward and backward wave vectors in an arbitrary layer, $\mathbf{k}_+^{(n)}$ and $\mathbf{k}_-^{(n)}$, can be constructed from two constraints:

1. The tangential component of the complex wave vector **k** is conserved [Snell–Descartes law in vectorial form, equation (27) in paper I]:

$$\hat{s}\wedge\mathbf{k}_+^{(n)} = \hat{s}\wedge\mathbf{k}_-^{(n)} = \hat{s}\wedge\mathbf{k}^{(m)}. \tag{73}$$

2. The complex square of the wave vector $[\mathbf{k}_\pm^{(n)}]^2 = (\bar{n}^{(n)})^2 k_0^2$ depends on the complex refractive index $\bar{n}^{(n)}$ and the vacuum wave vector $k_0$ only. Thus the perpendicular component is

$$\hat{s}\cdot\mathbf{k}_+^{(n)} = -\hat{s}\cdot\mathbf{k}_-^{(n)} = \{[\mathbf{k}^{(n)}]^2 - [\hat{s}\wedge(\hat{s}\wedge\mathbf{k}^{(m)})]^2\}^{1/2}. \tag{74}$$

The standard approach to calculate the field amplitude in each layer is to use a transfer matrix formalism.¹⁶ Maxwell's boundary conditions (46)–(49) from paper I are written at each interface for the total fields:

$$\mathbf{E}_{\text{PE}+}^{(n)}\exp[-i\hat{s}\cdot\mathbf{k}_+^{(n)}d_n] + \mathbf{E}_{\text{PE}-}^{(n)}\exp[-i\hat{s}\cdot\mathbf{k}_-^{(n)}d_n] = \mathbf{E}_{\text{PE}+}^{(n+1)} + \mathbf{E}_{\text{PE}-}^{(n+1)}, \tag{75}$$

$$\frac{1}{\mu^{(n)}}\{[\hat{s}\cdot\mathbf{k}_+^{(n)}]\mathbf{E}_{\text{PE}+}^{(n)}\exp[-i\hat{s}\cdot\mathbf{k}_+^{(n)}d_n] + [\hat{s}\cdot\mathbf{k}_-^{(n)}]\mathbf{E}_{\text{PE}-}^{(n)}\exp[-i\hat{s}\cdot\mathbf{k}_-^{(n)}d_n]\} = \frac{1}{\mu^{(n+1)}}\{[\hat{s}\cdot\mathbf{k}_+^{(n+1)}]\mathbf{E}_{\text{PE}+}^{(n+1)} + [\hat{s}\cdot\mathbf{k}_-^{(n+1)}]\mathbf{E}_{\text{PE}-}^{(n+1)}\}. \tag{76}$$

Solving for $\mathbf{E}_{\text{PE}+}^{(n)}$ and $\mathbf{E}_{\text{PE}-}^{(n)}$, we obtain the following transfer matrix:

$$\begin{pmatrix}\mathbf{E}_{\text{PE}+}^{(n)} \\ \mathbf{E}_{\text{PE}-}^{(n)}\end{pmatrix} = \begin{bmatrix}\mathcal{A}_\text{PE}^{(n)} & \mathcal{B}_\text{PE}^{(n)} \\ \mathcal{C}_\text{PE}^{(n)} & \mathcal{D}_\text{PE}^{(n)}\end{bmatrix}\begin{pmatrix}\mathbf{E}_{\text{PE}+}^{(n+1)} \\ \mathbf{E}_{\text{PE}-}^{(n+1)}\end{pmatrix}, \tag{77}$$

where

$$\mathcal{A}_\text{PE}^{(n)} = \tfrac{1}{2}\exp(i\hat{s}\cdot\mathbf{k}_+^{(n)}d_n)\left[1 + \frac{\mu^{(n)}}{\mu^{(n+1)}}\frac{\hat{s}\cdot\mathbf{k}_+^{(n+1)}}{\hat{s}\cdot\mathbf{k}_+^{(n)}}\right], \tag{78}$$

$$\mathcal{B}_\text{PE}^{(n)} = \tfrac{1}{2}\exp(i\hat{s}\cdot\mathbf{k}_+^{(n)}d_n)\left[1 - \frac{\mu^{(n)}}{\mu^{(n+1)}}\frac{\hat{s}\cdot\mathbf{k}_-^{(n+1)}}{\hat{s}\cdot\mathbf{k}_+^{(n)}}\right], \tag{79}$$

$$\mathcal{C}_\text{PE}^{(n)} = \tfrac{1}{2}\exp(-i\hat{s}\cdot\mathbf{k}_-^{(n)}d_n)\left[1 - \frac{\mu^{(n)}}{\mu^{(n+1)}}\frac{\hat{s}\cdot\mathbf{k}_+^{(n+1)}}{\hat{s}\cdot\mathbf{k}_-^{(n)}}\right], \tag{80}$$

$$\mathcal{D}_\text{PE}^{(n)} = \tfrac{1}{2}\exp(-i\hat{s}\cdot\mathbf{k}_-^{(n)}d_n)\left[1 + \frac{\mu^{(n)}}{\mu^{(n+1)}}\frac{\hat{s}\cdot\mathbf{k}_-^{(n+1)}}{\hat{s}\cdot\mathbf{k}_+^{(n)}}\right], \tag{81}$$

where $d_n$ is the distance between each interface. Similar equations can be obtained for PM modes with the simple replacement of **E** by **H** and the replacement of the permittivity $\mu$ by $\bar\epsilon$. Using such transfer matrix equations for stacks with arbitrary layers, one can easily calculate the complex reflection and transmission coefficients for field amplitudes.¹⁶

We now outline how to formulate the problem directly in terms of the Poynting vector. The mixed Poynting vector appears naturally in this context, and so do the complex power reflection and transmission coefficients. The basic idea is to express, at each interface, the output perpendicular components of the Poynting vector in terms of the input components:

$$S_{\text{PE}+}^{\perp(n)} = \bar{\mathcal{T}}_{\text{PE}+}^{(n+1)}S_{\text{PE}+}^{\perp(n)} - \bar{\mathcal{R}}_{\text{PE}-}^{(n+1)}S_{\text{PE}-}^{\perp(n+1)} + \bar{\mathcal{R}}_{\text{PEm}+}^{(n)}S_{\text{PEm}+}^{\perp(n)} - \bar{\mathcal{R}}_{\text{PEm}-}^{(n)}S_{\text{PEm}-}^{\perp(n+1)}, \tag{82}$$

$$S_{\text{PE}-}^{\perp(n)} = -\bar{\mathcal{T}}_{\text{PE}-}^{(n)}S_{\text{PE}-}^{\perp(n+1)} - \bar{\mathcal{R}}_{\text{PE}+}^{(n)}S_{\text{PE}+}^{\perp(n)} - \bar{\mathcal{R}}_{\text{PEm}+}^{(n)}S_{\text{PEm}+}^{\perp(n)} + \bar{\mathcal{R}}_{\text{PEm}-}^{(n)}S_{\text{PEm}-}^{\perp(n+1)}, \tag{83}$$

$$S_{\text{PEm}+}^{\perp(n)} = -\bar{\mathcal{M}}_{\text{PE}+}^{(n)}S_{\text{PE}+}^{\perp(n)} + \bar{\mathcal{M}}_{\text{PE}-}^{(n)}S_{\text{PE}-}^{\perp(n+1)}\bar{\mathcal{M}}_{\text{PEm}+}^{(n)}S_{\text{PEm}+}^{\perp(n)}, \tag{84}$$

$$S_{\text{PEm}-}^{\perp(n)} = (S_{\text{PEm}+}^{\perp(n)})^*, \tag{85}$$

where generalized complex energy flux transfer coefficients $\bar{\mathcal{T}}_{\text{PE}+}^{(n+1)}, \ldots, \bar{\mathcal{M}}_{\text{PEm}+}^{(n)}$ have been implicitly introduced. We notice that it is necessary to include backward-type coefficients because $S_{\text{PE}-}^{\perp(n+1)}$ must be included among the input components. Identical relations hold for PM modes, and similar expressions can be found for the crossed PE–PM components.

We may now solve Eqs. (82)–(85) in terms of the components at each interface and express the resulting equations in terms of transfer matrices for the perpendicular components of the Poynting vector in each layer:

$$\begin{pmatrix}S_{\text{PE}+}^{\perp(n)} \\ S_{\text{PEm}+}^{\perp(n)} \\ S_{\text{PEm}-}^{\perp(n)} \\ S_{\text{PE}-}^{\perp(n)}\end{pmatrix} = \begin{bmatrix}\mathcal{E}_\text{PE}^{(n)} & \mathcal{F}_\text{PE}^{(n)} & \mathcal{G}_\text{PE}^{(n)} & \mathcal{H}_\text{PE}^{(n)} \\ \mathcal{I}_\text{PE}^{(n)} & \mathcal{J}_\text{PE}^{(n)} & \mathcal{K}_\text{PE}^{(n)} & \mathcal{L}_\text{PE}^{(n)} \\ \mathcal{N}_\text{PE}^{(n)} & \mathcal{O}_\text{PE}^{(n)} & \mathcal{P}_\text{PE}^{(n)} & \mathcal{Q}_\text{PE}^{(n)} \\ \mathcal{T}_\text{PE}^{(n)} & \mathcal{U}_\text{PE}^{(n)} & \mathcal{V}_\text{PE}^{(n)} & \mathcal{W}_\text{PE}^{(n)}\end{bmatrix}\begin{pmatrix}S_{\text{PE}+}^{\perp(n+1)} \\ S_{\text{PEm}+}^{\perp(n+1)} \\ S_{\text{PEm}-}^{\perp(n+1)} \\ S_{\text{PE}-}^{\perp(n+1)}\end{pmatrix}, \tag{86}$$

where the matrix elements $\mathcal{E}_\text{PE}^{(n)}, \ldots, \mathcal{W}_\text{PE}^{(n)}$ can be evaluated in terms of the complex energy flux transfer coefficients $\bar{\mathcal{T}}_\text{PE}^{(n+1)}, \ldots, \bar{\mathcal{M}}_\text{PEm}^{(n)}$. We omit such lengthy expressions here, as well as expressions of the transfer coefficients in terms of Fresnel coefficients. Nowadays computer software such as MATHEMATICA (Wolfram Inc.) can reliably provide and evaluate such expressions. We propose to use Eq. (86) as a working tool for all situations in which one is directly interested in a detailed understanding of energy fluxes inside a stack of thin films. In this picture the mixed Poynting vector finds its most natural and general use.

## 8. SYMMETRY OF MAXWELL'S EQUATIONS AND GENERALIZED RECIPROCITY RELATIONS

As outlined in Section 1, we discuss reciprocity relations from first principles and independently of energy balance. Maxwell's equations for isotropic conducting media are obeyed with the fields considered in Section 7:

$$\nabla\wedge\mathbf{E} = -\mu_0\mu\partial_t\mathbf{H}, \qquad \nabla\cdot(\epsilon_0\epsilon\mathbf{E}) = \rho,$$
$$\nabla\wedge\mathbf{H} = \sigma\mathbf{E} + \epsilon_0\epsilon\partial_t\mathbf{E}, \qquad \nabla\cdot(\mu\mathbf{H}) = 0. \tag{87}$$

It is obvious that these equations are invariant with respect to time reversal $T: t \to t' = -t$ provided that the following transformations are applied simultaneously:

$$\mathbf{E}' = \mathbf{E}, \qquad \mathbf{H}' = -\mathbf{H}, \qquad \sigma' = -\sigma. \tag{88}$$

Negative conductivities phenomenologically describe media that provide optical gain.

We recall also that Maxwell equations are separately invariant under charge conjugation $C: \rho \to \rho'' = -\rho$ and spatial inversion $P: \mathbf{r} \to \mathbf{r}''' = -\mathbf{r}$. Although charge conjugation can be applied without consequences in the present framework, we note that the invariance under spatial inversion is broken by both the interface and the implicit spatial asymptotic behavior present in Eq. (89).

Šantavý⁶ showed that time-reversal invariance was at the root of generalized reversibility relations that could be used successfully in many circumstances. We reformulate this theorem in our framework and show that Šantavý's formulation does not describe the most general case. This fact can be easily demonstrated by looking at the transformation of the complex wave vector [Eq. (6) of paper I] for a HIPW. If we denote the mathematical operation of inversing the conductivity sign by ħ [?], and the operation of inversing the sign of the inhomogeneity of the wave by ♭ [?], we find the identity

$$\mathbf{k}^{♭ħ} = \bar{n}^ħ k_0\mathbf{n}^♭ = \bar{n}^*k_0\mathbf{n}^* = \mathbf{k}^*. \tag{89}$$

This type of identity also applies to the Fresnel coefficients, which depend only on the complex wave vector and (for the PM case) on the complex dielectric permittivity. Dealing essentially with the case of homogeneous waves and normal incidence, Šantavý did not need to introduce the ♭ [?] operation.

Let us suppose that two PE waves $\mathbf{E}_{\text{PE}+}^{(0)}$ and $\mathbf{E}_{\text{PE}+}^{(N+1)}$ are incident upon each side of the stack of Section 7 (a completely parallel discussion also holds for PM waves of course). In Section 8 the $N + 1$th layer width is assumed to be large, so that the incident amplitudes represent the electric field infinitely closely to the interfaces on both sides. It is possible to express the output fields as a linear superposition of the input fields:

$$\begin{pmatrix}\mathbf{E}_{\text{PE}-}^{(0)} \\ \mathbf{E}_{\text{PE}+}^{(N+1)}\end{pmatrix} = \begin{bmatrix}\mathfrak{R}_{\text{PE}+} & \Theta_{\text{PE}-} \\ \Theta_{\text{PE}+} & \mathfrak{R}_{\text{PE}-}\end{bmatrix}\begin{pmatrix}\mathbf{E}_{\text{PE}+}^{(0)} \\ \mathbf{E}_{\text{PE}-}^{(N+1)}\end{pmatrix} = \mathfrak{M}_\text{PE}\begin{pmatrix}\mathbf{E}_{\text{PE}+}^{(0)} \\ \mathbf{E}_{\text{PE}-}^{(N+1)}\end{pmatrix}. \tag{90}$$

We now apply transformation rules (88) to find time-reversed solutions of Maxwell's equations:

$$\begin{pmatrix}\mathbf{E}_{\text{PE}-}^{(0)ħ*} \\ \mathbf{E}_{\text{PE}+}^{(N+1)ħ*}\end{pmatrix} = \mathfrak{M}_\text{PE}^{♭ħ*}\begin{pmatrix}\mathbf{E}_{\text{PE}+}^{(0)ħ*} \\ \mathbf{E}_{\text{PE}-}^{(N+1)ħ*}\end{pmatrix}. \tag{91}$$

We note that we recover the usual antilinearity of the time-reversal transformation on complex amplitudes. Because of the symmetry of the complex wave vector (89), the waves described by the complex conjugate amplitudes on the right-hand side of this equation are outgoing waves because

$$\{\exp[-i(\mathbf{k}\cdot\mathbf{r} - \omega t)]\}^{♭ħ*} = \exp[+i(\mathbf{k}\cdot\mathbf{r} + \omega t)]. \tag{92}$$

Because the relationship between output and input wave amplitudes are linear and unique Eqs. (90) and (91) are equivalent, and the following matrix identity holds:

$$\mathfrak{M}_\text{PE}\mathfrak{M}_\text{PE}^{♭ħ*} = 1. \tag{93}$$

This relationship can be rewritten in terms of the scattering amplitudes $\mathfrak{R}_\text{PE\pm}$ and $\Theta_\text{PE\pm}$, for which we find the generalized form of the reciprocity relations:

$$\mathfrak{R}_{\text{PE}+}\Theta_{\text{PE}-}^{♭ħ*} + \Theta_{\text{PE}-}\mathfrak{R}_{\text{PE}+}^{♭ħ*} = 1, \tag{94}$$

$$\mathfrak{R}_{\text{PE}+}\Theta_{\text{PE}+}^{♭ħ*} + \Theta_{\text{PE}+}\mathfrak{R}_{\text{PE}+}^{♭ħ*} = 0, \tag{95}$$

$$\mathfrak{R}_{\text{PE}-}\Theta_{\text{PE}+}^{♭ħ*} + \Theta_{\text{PE}+}\mathfrak{R}_{\text{PE}-}^{♭ħ*} = 1, \tag{96}$$

$$\mathfrak{R}_{\text{PE}-}\Theta_{\text{PE}-}^{♭ħ*} + \Theta_{\text{PE}-}\mathfrak{R}_{\text{PE}+}^{♭ħ*} = 0. \tag{97}$$

Other forms of these relationships can be obtained by redistributing $*$, ħ, or ♭ [?], using the fact that these operations are commutative and idempotent. Šantavý noticed that the following interesting identities,

$$\Theta_{\text{PE}+}\Theta_{\text{PE}-}^{♭ħ*} = \Theta_{\text{PE}-}\Theta_{\text{PE}+}^{♭ħ*}, \tag{98}$$

$$\mathfrak{R}_{\text{PE}+}\Theta_{\text{PE}-}^{♭ħ*} = \mathfrak{R}_{\text{PE}-}\Theta_{\text{PE}+}^{♭ħ*}, \tag{99}$$

are a consequence of the generalized reversibility equations.

## 9. CONSEQUENCES OF THE GENERALIZED RECIPROCITY RELATIONS

The generalized reciprocity relations (94)–(97) can be used to recover the known reciprocity rules for the special cases detailed in the Introduction. Instead let us first state some simple general consequences:

(a) For a single interface the coefficients $\mathfrak{R}_\text{PE\pm}$ and $\Theta_\text{PE\pm}$ reduce to the Fresnel coefficients $r_\text{PE\pm}$ and $t_\text{PE\pm}$. Because of the symmetry of the wave vector and of the Fresnel coefficients (89) the nonstandard Stokes relations (1)–(4) with complex coefficients are always obeyed.

(b) If more than one interface is considered, the general scattering amplitudes $\mathfrak{R}_\text{PE\pm}$ and $\Theta_\text{PE\pm}$ involve complex propagators of the form $\exp(-i\hat{s}\cdot\mathbf{k}d_n)$ that are not invariant [see Eq. (89)] and transform according to Eq. (92). Therefore the generalized reciprocity relations do not in general reduce to simple special cases like Eqs. (5)–(8) or Eqs. (1)–(4).

(c) If more than one interface is considered and in purely dielectric media $\sigma = 0$, the ħ operation is unoperative ($r^ħ = r$). For homogeneous incident waves above the angle for total internal reflection (in all layers), the ♭ [?] operation is also inoperative. The scattering amplitudes, however, transform into complex conjugates because of the complex propagators [see Eq. (92)], so that the standard generalized complex Stokes relations (5)–(8) are recovered.

These considerations generalize all previous restricting conditions linked with earlier forms of the real or complex Stokes relations [Eqs. (1)–(4) or Eqs. (5)–(8)]. As a final remark we mention that for the case of total reflection in purely dielectric media (subcase of case a above) Šantavý had to resort to *ad hoc* rules (see Ref. 6) because he did not include the possibility of inhomogeneous waves in his theorem. In the present formulation this case is covered by the very general result (a), which is a consequence of the generalized reversibility theorem, Eqs. (94)–(97).

Another important application of the generalized reversibility equations concerning the use of recurrent calculations was pointed out by Knittl.² We do not reproduce here the calculations of section 3 of Knittl's paper. It is sufficient to point out that our generalized formulation is also suited for such methods, especially in the case of nonperpendicular incidence, in which the wave inhomogeneity becomes crucial.

## 10. CONCLUSIONS

We have discussed in detail energy balance and reciprocity relations for HIPW's that are important for applications in thin-film optics. In summarizing previous papers on this subject we have emphasized how loose the connection is between these seemingly related issues. The consideration of arbitrary inhomogeneous plane waves allowed us to give a (hopefully) final description of energy balance and reciprocity. Among the new features are the presence of a cross PE–PM term in the Poynting vector and a more general formulation of the reversibility theorem for both dielectric and conducting layered media.

The study of energy balance led to detailed considerations on the role of the mixed Poynting vector and the transfer of the different kinds of energy fluxes through the interface. A new transfer matrix formalism for energy fluxes was presented as a possible generalized use of the mixed Poynting vector concept.

The general theorem on reversibility allows us to find quickly and test the symmetry relations obeyed by the scattering amplitudes in arbitrary stacks of conducting dielectric media. It also brings new insight in the most simple cases, such as total reflection at the boundary between pure dielectrics.

## ACKNOWLEDGMENT

We thank D. F. Nelson for useful discussions on this subject and for sending us a preprint of Refs. 3 and 4.

## REFERENCES

*[Reference list omitted per transcription scope. Key references cited above: Ref. 1 = Dupertuis, Proctor & Acklin, "Generalization of complex Snell–Descartes and Fresnel laws," J. Opt. Soc. Am. A 11, 1159–1166 (1994) (= "paper I"); Ref. 2 = Z. Knittl, "The principle of reversibility and thin film optics," Opt. Acta 9 (1962); Refs. 3, 4 = B. Chen and D. F. Nelson, wave-vector space method; Ref. 6 = I. Šantavý, "On the reversibility of light beams in conducting media," Opt. Acta 8, 301–307 (1961); Ref. 8 = Born and Ladenburg (1911); Refs. 9, 10 = C. von Fragstein (1950, 1987) on the mixed Poynting vector; Ref. 16 = Yariv and Yeh, Optical Waves in Crystals (Wiley, 1984).]*

---

## Collected uncertain glyph readings [?]

- **ħ operation** (Sections 1, 8, 9; Eqs. 89–99): the symbol for "inversing the conductivity sign" — rendered here as `ħ`. In the original it is a script/German-Fraktur-like letter attached as a superscript operation. Its role (flip $\sigma \to -\sigma$) is unambiguous from context; only the exact glyph is uncertain.
- **♭ operation** (Sections 8, 9; Eqs. 89–99): the symbol for "inversing the sign of the inhomogeneity of the wave" — rendered here as `♭` (flat/b-like). Role unambiguous ($\mathbf{n} \to \mathbf{n}^*$, i.e. flip the sign of the imaginary/inhomogeneous part), exact glyph uncertain. Together $♭ħ*$ is the combined operation appearing in the reciprocity relations.
- **Eq. (48)**: the second term in the bracket, transcribed as $(\mathbf{E}_i\cdot\mathbf{k}_i^\perp)E_i^{*\perp}$, has an ambiguous final factor in the scan (whether it is $E_i^{*\perp}$ or $\mathbf{E}_i^{*\perp}$, scalar vs vector). Structure follows Eq. (21).
- **Eq. (84)**: the third term $\bar{\mathcal{M}}_{\text{PEm}+}^{(n)}S_{\text{PEm}+}^{\perp(n)}$ appears to be missing a preceding operator sign ($+$ or $-$) in the scan; transcribed as printed.
- **Eqs. (67)/(68)**: printed in the original as two separately numbered lines for one multiline expression; combined here under label "(67, 68)".
