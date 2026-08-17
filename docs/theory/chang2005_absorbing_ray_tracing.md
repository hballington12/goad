# Chang et al. (2005) — Ray tracing in absorbing media

> **Transcription for GOAD theory work.**
> Source file: `Chang 2005 Ray tracing in absorbing media.pdf`
> Transcribed: 2026-08-16
>
> This is a faithful markdown/LaTeX transcription of the mathematical content and
> section structure of the paper. Equation numbers match the original so that
> later work can cite e.g. "Chang eq. (12)". Purely bibliographic content (full
> reference list) is omitted and figure captions are given briefly. Uncertain
> readings are marked `[?]` and collected at the bottom of the file.

---

**Ray tracing in absorbing media**

Peter C.Y. Chang, J.G. Walker, K.I. Hopcraft

- School of Electrical and Electronic Engineering, University of Nottingham, University Park, Nottingham, NG7 2RD, UK
- Theoretical Mechanics Division, School of Mathematical Sciences, University of Nottingham, University Park, Nottingham, NG7 2RD, UK

*Journal of Quantitative Spectroscopy & Radiative Transfer* **96** (2005) 327–341.
Received 8 October 2004; accepted 23 December 2004. doi:10.1016/j.jqsrt.2005.01.001

## Abstract

In ray tracing, it is standard practice to consider solely the real part of any medium's refractive index for determining the changes in a ray's direction at the interface between two media. Any absorption is accounted for by reducing the intensity or weight of the ray in accordance with its propagation distance. This practice is adequate at optical wavelengths where absorption is generally very low but in the infrared and other regimes absorption can be significant and refraction angles and reflectances can deviate substantially from lossless cases. We calculate the quantitative consequence of a more rigorous approach to ray tracing for a lossy wedge and explore the potential differences in results from a Monte Carlo simulation of reflection from a rough, lossy surface. Comparisons between the traditional and complex ray tracing approaches are made and show that there can be significant discrepancies.

*Keywords:* Ray tracing; Lossy media; Inhomogeneous plane waves; Rough surface scattering

---

## 1. Introduction

Ray tracing is used commonly as a tool for calculating the response of a scene to light. The scenarios can include computer-generated imagery [1], single scattering and multiple scattering situations. In a single scattering scenario, the representative object can be a composite of smaller objects and/or have a complicated geometry [2–4]. Multiple scattering effects are important for propagation through turbid media, from structured and rough surfaces [5–7], and for thermal systems [8,9].

The usual assumptions underlying ray tracing are that geometrical optics applies, all intervening media are uniform and dielectric and all interfaces are planar. The first assumption implies that the wavelength of light is short compared with other scales in the scene. The second restricts the media to be homogeneous and non-absorbing. Note that classically dielectrics are defined as non-conductors and have no internal loss mechanism (their electric susceptibility is real) at frequencies far from atomic and molecular resonances. These first two assumptions allow us to use rays to represent the propagation of light. The third implies that curved surfaces are approximated to be locally flat at points where rays intersect them. Together these restrictions permit the use of Snell's law and Fresnel coefficients for calculating a ray's properties on its interaction with any surface in the scene. If there are inhomogeneities in the media or the surfaces encountered have large curvatures, higher order corrections to the ray have to be applied using the eikonal equation [10].

Strong absorption can occur in electrically conductors and semiconductors [11], biological systems such as human tissue and industrial systems such as engines and furnaces where radiative transfer could be the dominant means of energy transport [9].

In absorbing media, there are three separate problems that need to be addressed. The familiar Snell's law and Fresnel formulae are no longer valid and a ray is attenuated along its path. It is common for the first two problems to be neglected by using the real part of the refractive indices that describe the media alone. The path loss is accounted for by a simple exponential law. This will be shown to be inadequate because they rely on ad hoc simplifications. We shall refer to this traditional approach as the 'standard' method: see, for example, Refs. [4,7].

The following sections will detail the physics for a more rigorous approach, examine the low absorption case, illustrate the differences with a quantitative study of light passing through a lossy wedge and finally the scattering from a rough lossy surface is studied with a Monte Carlo ray tracing model.

---

## 2. Physical basis of ray tracing in absorbing media

Geometrical optics is based on the solution to Maxwell's equations for electromagnetic radiation of short wavelength. For a homogeneous medium, it approximates light propagation with a ray whose direction coincides with the normal to planes of constant phase of a travelling plane wave.

The problem of solving Maxwell's equations for plane waves travelling through an interface is canonical. However, when there is absorption or loss in one or both media separated by the interface, it is only recently that the problem has been completely settled [3,12–22]. The definitive paper is by Dupertuis et al. [21] which clarified some of the issues surrounding propagation between two absorbing media and pointed out flaws in previous work. We use this paper's results to derive the more rigorous ray tracing method and discuss the predicted discrepancies between it and the traditional approach for two readily realisable instances. There is also a "high-frequency" solution derived by Yang and Liou [2]. Owing to the complex nature of the refractive index, the wave vector is also complex and hence the plane waves supported by lossy media are inhomogeneous — that is, the normals to the planes of constant phase and constant amplitude do not coincide.

The inhomogeneity leads to longitudinal components of the electric and magnetic field vectors and this makes it more difficult to find formulae equivalent to the Fresnel formulae for dielectric media. Many textbooks and papers gloss over the problem of two absorbing media and even present a complex version of Snell's law which is correct only under certain conditions.

The rest of this section reviews the derivation of the generalized Snell's law and Fresnel formulae. Results are presented that form the basis of a more rigorous ray tracing methodology.

### 2.1. Inhomogeneous plane wave solution to Maxwell's equations

Maxwell's equations for a homogeneous, isotropic, linear, stationary and charge-free system are

$$\nabla \wedge \mathbf{E}(\mathbf{r},t) = -\mu\,\partial_t \mathbf{H}(\mathbf{r},t), \qquad \nabla \cdot \mathbf{E}(\mathbf{r},t) = 0,$$

$$\nabla \wedge \mathbf{H}(\mathbf{r},t) = \varepsilon\,\partial_t \mathbf{E}(\mathbf{r},t) + \sigma \mathbf{E}(\mathbf{r},t), \qquad \nabla \cdot \mathbf{H}(\mathbf{r},t) = 0. \tag{1}$$

An absorbing or lossy medium can be modelled with a complex refractive index $\tilde{n} = n + \mathrm{i}k$.[?1] For plane wave solutions, we impose $\mathbf{E}(\mathbf{r},t) = \mathbf{E}\exp(\mathrm{i}\mathbf{k}\cdot\mathbf{r} - \mathrm{i}\omega t)$ and $\mathbf{H}(\mathbf{r},t) = \mathbf{H}\exp(\mathrm{i}\mathbf{k}\cdot\mathbf{r} - \mathrm{i}\omega t)$, and obtain

$$\mathbf{k} \wedge \mathbf{E} = \mu\omega\mathbf{H}, \qquad \mathbf{k}\cdot\mathbf{E} = 0,$$

$$\mathbf{k} \wedge \mathbf{H} = -\tilde{\varepsilon}\omega\mathbf{E}, \qquad \mathbf{k}\cdot\mathbf{H} = 0, \tag{2}$$

where $\tilde{\varepsilon} = \varepsilon + \mathrm{i}\sigma/\omega = \varepsilon_0\tilde{n}^2$.

Assuming the medium is non-magnetic ($\mu = \mu_0$), the wave vector $\mathbf{k} = k_0(N\mathbf{e} + \mathrm{i}K\mathbf{f})$ is complex for $k_0 = \omega/c$ and $c^2 = 1/\varepsilon_0\mu_0$. From Eqs. (2), it can be shown that the wave vector's components obey

$$N^2 - K^2 = n^2 - \kappa^2 \quad\text{and}\quad NK\mathbf{e}\cdot\mathbf{f} = n\kappa, \tag{3}$$

where $\mathbf{e}$ and $\mathbf{f}$ are unit normals to the planes of constant phase and amplitude, respectively. Note that $\mathbf{E}$ and $\mathbf{H}$ are also complex vectors. Complex vectors are known as bivectors [23] and obey the usual identities involving scalar and vector products.

As a consequence of Eqs. (3), plane waves can be inhomogeneous when the medium is lossy and the normals to the planes of constant phase and the planes of constant amplitude no longer coincide. The associated ray is complex with its real part perpendicular to surfaces of constant phase and its imaginary part perpendicular to surfaces of constant amplitude.

$N$ and $K$ can be interpreted as effective or *apparent* refractive indices. Defining $\cos\alpha = \mathbf{e}\cdot\mathbf{f}$,

$$N^2 = \tfrac{1}{2}\!\left[\,n^2 - \kappa^2 + \sqrt{(n^2 - \kappa^2)^2 + 4\!\left(\frac{n\kappa}{\cos\alpha}\right)^{\!2}}\,\right] \tag{4}$$

and $K$ can be obtained from (3). Eq. (4) determines $N$ when given the refractive index $\tilde{n}$ and the inhomogeneity angle $\alpha$.

An alternative formulation is to define $\mathbf{k} = k_0\tilde{n}\mathbf{g}$, where $\mathbf{g}\cdot\mathbf{g} = 1$. The complex unit vector $\mathbf{g}$ can be resolved into real and imaginary components

$$\mathbf{g} = \cosh\beta\,\mathbf{g}_\parallel + \mathrm{i}\sinh\beta\,\mathbf{g}_\perp \tag{5}$$

with $\mathbf{g}_\parallel\cdot\mathbf{g}_\parallel = \mathbf{g}_\perp\cdot\mathbf{g}_\perp = 1$ and $\mathbf{g}_\parallel\cdot\mathbf{g}_\perp = 0$. The inhomogeneity parameter $\beta$ acts in a similar manner to the angle $\alpha$ in Eq. (4). It can be seen from

$$NK\mathbf{e}\wedge\mathbf{f} = \tfrac{1}{2}(n^2 + \kappa^2)\sinh 2\beta\,\mathbf{g}_\parallel\wedge\mathbf{g}_\perp \tag{6}$$

that a wave is homogeneous when $\beta = 0$. In fact, we can derive a direct relationship between $\alpha$ and $\beta$:

$$\sinh 2\beta = \frac{2n\kappa}{n^2 + \kappa^2}\tan\alpha \tag{7}$$

by using (6) and (3) together with the fact that $\alpha \geqslant 0$ and the convention that $\beta \geqslant 0$. Note that in a lossless medium, $\kappa = 0$, there are two possibilities: $N = n$ and $K = 0$ or $\alpha = \pi/2$ (where $\beta$ is not determined by (7)). The first is the homogeneous case whilst the second case represents $\mathbf{f}$ being perpendicular to $\mathbf{e}$.

### 2.2. Planar interface

Consider two half spaces separated by a plane with normal $\mathbf{n}$ and let the two media have refractive indices $n_1,\kappa_1,n_2,\kappa_2$ with $\mathbf{k}_\mathrm{i},\mathbf{k}_\mathrm{r}$ and $\mathbf{k}_\mathrm{t}$ denoting the incident, reflected and transmitted wave vectors. The geometry is shown in Fig. 1: $\theta_\mathrm{i}$ and $\psi_\mathrm{i}$ are the incident angles for the real and imaginary parts of the wave vector. Without loss of generality, we can assume $\mathbf{f}_\mathrm{i}$ has an azimuthal angle $\phi_\mathrm{i}$.

*[Figure 1: geometry of the interface between two lossy media; the real part (with unit vector $\mathbf{e}$, angle $\theta$) and imaginary part (with unit vector $\mathbf{f}$, angles $\psi$ and azimuth $\phi$) of the incident and transmitted wave vectors are shown in separate diagrams (a) and (b).]*

#### 2.2.1. Snell's law

Matching the phases of plane waves at the interface gives the boundary conditions

$$\mathbf{n}\wedge\mathbf{k}_\mathrm{i} = \mathbf{n}\wedge\mathbf{k}_\mathrm{r} = \mathbf{n}\wedge\mathbf{k}_\mathrm{t} \tag{8}$$

from which a generalized Snell's law follows:

$$N_1\sin\theta_\mathrm{i} = N_2\sin\theta_\mathrm{t}, \qquad K_1\sin\psi_\mathrm{i} = K_2\sin\psi_\mathrm{t}, \tag{9}$$

where the new apparent refractive indices can be calculated from the quartic equation

$$N_2^4 - \left[N_s^2 + K_s^2 + (n_2^2 - \kappa_2^2)\right]N_2^2 + N_s^2 K_s^2 + (n_2^2 - \kappa_2^2)N_s^2 - (n_2\kappa_2 - N_s K_s\cos\phi_\mathrm{i})^2 = 0 \tag{10}$$

with $N_s \equiv N_1\sin\theta_\mathrm{i}$, $K_s \equiv K_1\sin\psi_\mathrm{i}$. The normal components of the incident and reflected wave vectors satisfy

$$\mathbf{n}\cdot\mathbf{k}_\mathrm{r} = -\mathbf{n}\cdot\mathbf{k}_\mathrm{i} \tag{11}$$

and hence $\theta_\mathrm{r} = \pi - \theta_\mathrm{i}$, $\psi_\mathrm{r} = \pi - \psi_\mathrm{i}$ and $\phi_\mathrm{r} = \phi_\mathrm{i}$. Note that $\phi_\mathrm{t} = \phi_\mathrm{i}$ throughout.

Therefore, the real and imaginary components of the wave vector have reflected and transmitted directions given by

$$\mathbf{e}_\mathrm{r} = \mathbf{e}_\mathrm{i} - 2\cos\theta_\mathrm{i}\,\mathbf{n}, \qquad \mathbf{e}_\mathrm{t} = \frac{N_1\mathbf{e}_\mathrm{i}}{N_2} + \left(\cos\theta_\mathrm{t} - \frac{N_1\cos\theta_\mathrm{i}}{N_2}\right)\mathbf{n},$$

$$\mathbf{f}_\mathrm{r} = \mathbf{f}_\mathrm{i} - 2\cos\psi_\mathrm{i}\,\mathbf{n}, \qquad \mathbf{f}_\mathrm{t} = \frac{K_1\mathbf{f}_\mathrm{i}}{K_2} + \left(\cos\psi_\mathrm{t} - \frac{K_1\cos\psi_\mathrm{i}}{K_2}\right)\mathbf{n}, \tag{12}$$

where $\cos\theta_\mathrm{t}$ and $\cos\psi_\mathrm{t}$ have the same signs as $\cos\theta_\mathrm{i} = \mathbf{n}\cdot\mathbf{e}_\mathrm{i}$ and $\cos\psi_\mathrm{i} = \mathbf{n}\cdot\mathbf{f}_\mathrm{i}$, respectively. These equations allow us to trace the progress of a ray through any number of interfaces between absorbing media.

#### 2.2.2. Fresnel coefficients

The integral form of Maxwell's equations give rise to continuity conditions at the interface. These are that the tangential components of $\mathbf{E}$ and $\mathbf{H}$ are continuous and that the normal components of $\tilde{\varepsilon}\mathbf{E}$ and $\mu\mathbf{H}$ are discontinuous due to surface charges and currents, respectively. To apply these conditions, we need to decompose the field vectors into components tangential and normal to the interface. Following Dupertuis et al. [21], we define a set of complex orthonormal vectors:

$$\mathbf{h}_1 = \frac{\mathbf{n}\wedge\mathbf{k}}{\sqrt{(\mathbf{n}\wedge\mathbf{k})^2}}, \qquad \mathbf{h}_2 = \frac{\mathbf{k}\wedge\mathbf{h}_1}{k_0\tilde{n}}. \tag{13}$$

The projections of the field vectors on $\mathbf{h}_1$ and $\mathbf{h}_2$ correspond to the TE (perpendicular or s) and TM (parallel or p) polarizations, respectively, since they match definitions with the homogeneous case, where $\mathbf{k}$ is real or $\mathbf{f}$ is parallel to $\mathbf{e}$.

Solving the boundary conditions gives the generalized Fresnel formulae for reflection and transmission coefficients

$$r_\mathrm{s} = \frac{k_\mathrm{i} - k_\mathrm{t}}{k_\mathrm{i} + k_\mathrm{t}}, \qquad t_\mathrm{s} = \frac{2k_\mathrm{i}}{k_\mathrm{i} + k_\mathrm{t}},$$

$$r_\mathrm{p} = \frac{\tilde{n}_2^2 k_\mathrm{i} - \tilde{n}_1^2 k_\mathrm{t}}{\tilde{n}_2^2 k_\mathrm{i} + \tilde{n}_1^2 k_\mathrm{t}}, \qquad t_\mathrm{p} = \frac{2\tilde{n}_1\tilde{n}_2 k_\mathrm{i}}{\tilde{n}_2^2 k_\mathrm{i} + \tilde{n}_1^2 k_\mathrm{t}}, \tag{14}$$

where $k \equiv \mathbf{k}\cdot\mathbf{n} = k_0(N\cos\theta + \mathrm{i}K\cos\psi)$.[?2]

It should be emphasized that in the absorbing medium, $E_\mathrm{s} = \mathbf{E}\cdot\mathbf{h}_1$ is not the transverse electric field component and similarly $E_\mathrm{p}$ is not the TM component. However, in air or vacuo $E_\mathrm{s}$ and $E_\mathrm{p}$ become the usual quantities on exiting lossy media by transformation according to the formulae in Eq. (14).

#### 2.2.3. Path loss

A ray attenuates along a path in an absorbing medium. With a complex wave vector, the normal to planes of constant amplitude is at an angle $\alpha$ to the normal to planes of constant phase. So along the ray's path, $\mathbf{r} = l\mathbf{e}$,

$$\mathbf{E}(l\mathbf{e}) = \mathbf{E}\exp(\mathrm{i}\,\mathbf{k}\cdot\mathbf{r})$$
$$= \mathbf{E}\exp(-k_0 K l\cos\alpha)\exp(\mathrm{i}k_0 N l), \tag{15}$$

where $l$ is the path length.[?3]

Comparing the above with the usual loss factor $\exp(-k_0\kappa l)$, we see that the absorption coefficient in a lossy medium is dependent on the angle of entry into that medium.

---

## 3. Low absorption case

The formulae in the previous section allow us to explore the behaviour of rays at low absorption.

For simplicity, consider a ray travelling across a planar interface from air or vacuum to a medium of $\tilde{n} = n + \mathrm{i}\kappa$, where $n > 1$ and $x = (\kappa/n) \ll 1$. By using the binomial expansion in $x$, we find that, to second order,

$$N_2 \approx n\left[1 + \frac{\sin^2\theta_\mathrm{i}}{2(n^2 - \sin^2\theta_\mathrm{i})}x^2\right],$$

$$K_2 \approx \frac{n\kappa}{\sqrt{n^2 - \sin^2\theta_\mathrm{i}}}\left[1 - \frac{n^2\sin^2\theta_\mathrm{i}}{2(n^2 - \sin^2\theta_\mathrm{i})^2}x^2\right],$$

$$\sin\theta_\mathrm{t} \approx \frac{\sin\theta_\mathrm{i}}{n}\left[1 - \frac{\sin^2\theta_\mathrm{i}}{2(n^2 - \sin^2\theta_\mathrm{i})}x^2\right]. \tag{16}$$

These demonstrate that the lossless case ($N_2 = n$, $K_2 = 0$ and $\sin\theta_\mathrm{t} = \sin\theta_\mathrm{i}/n$) is reached in the limit of $\kappa = 0$.

Also, the real part of the apparent refractive index and the sine of the transmitted angle deviate from the lossless case by a squared dependency on the imaginary part of the refractive index. This second-order correction together with the maximum value of the coefficients ($1/(2(n^2 - 1))$) implies that $x$ has to be of reasonable size, say $> 0.15$ for $n = 1.5$, for the more rigorous method to show differences from the standard method greater than 1%.

The loss factor is determined to be

$$K_2\cos\theta_\mathrm{t} \approx \kappa\left[1 - \frac{\sin^2\theta_\mathrm{i}}{2(n^2 - \sin^2\theta_\mathrm{i})}x^2\right]. \tag{17}$$

Thus path loss has a second-order correction too. However, the reflection and transmission coefficients have first-order corrections. For example,

$$r_\mathrm{s} \approx \frac{\cos\theta_\mathrm{i} - \sqrt{n^2 - \sin^2\theta_\mathrm{i}}}{\cos\theta_\mathrm{i} + \sqrt{n^2 - \sin^2\theta_\mathrm{i}}}\left[1 + \frac{2\mathrm{i}}{n^2 - 1}\frac{n^2\cos\theta_\mathrm{i}}{\sqrt{n^2 - \sin^2\theta_\mathrm{i}}}x\right]. \tag{18}$$

Hence the standard Fresnel coefficients are only a good approximation for very low absorptions. However, the reflectivities and transmittivities have second-order corrections because they are defined using the absolute squares of the coefficients:

$$|r_\mathrm{s}|^2 \approx \left(\frac{\cos\theta_\mathrm{i} - \sqrt{n^2 - \sin^2\theta_\mathrm{i}}}{\cos\theta_\mathrm{i} + \sqrt{n^2 - \sin^2\theta_\mathrm{i}}}\right)^{\!2} \times \left[1 + \frac{2\left(2n^4 - (n^2 + 1)\sin^2\theta_\mathrm{i}\right)n^2\cos\theta_\mathrm{i}}{(n^2 - 1)^2(n^2 - \sin^2\theta_\mathrm{i})^{3/2}}x^2\right]. \tag{19}$$

The maximum value of the coefficient of $x^2$ in the reflectivity is $4n^3/(n^2 - 1)^2$ and is much greater than the corresponding coefficients in the approximations for the $\sin\theta_\mathrm{t}$ and the path loss. So the amount of absorption allowed for a given limit on the deviation from the standard reflectivity is much less than previously calculated. For example, for differences less than 1%, $x < 0.034$ when $n = 1.5$.

Yet another difference is that the reflectivity of light polarized parallel to the plane of incidence does not vanish at any angle for lossy media. That is, at Brewster's angle, as defined by $\tan(\theta_\mathrm{B}) = n$ [22], $r_\mathrm{p}$ is finite. In fact,

$$|r_\mathrm{p}(\theta_\mathrm{B})|^2 \approx \frac{(n^2 - 1)^2}{4n^4}x^2, \tag{20}$$

which shows the zero-order term, the standard reflectivity in lossless media, vanishes but that there is a positive second-order term.

Therefore, we can conclude that rays propagating through low absorption media can be approximated accurately by the standard approach. The only exception to this is that the polarization state is altered significantly as shown by Eq. (18) in a manner not accounted for using the standard approach. This is because the imaginary part of $r_\mathrm{p}r_\mathrm{s}^*$ is directly proportional $x$ to first order. So any absorption in the medium will lead to linearly polarized states being reflected in an elliptically polarized state.

---

## 4. Refraction by a lossy wedge

We now illustrate the differences between the standard and more rigorous methods by considering a wedge of absorbing material of angle $\gamma$ as shown in Fig. 2. Various quantities can be calculated analytically for this simple geometry using the results of the previous section.

Snell's law (9) with Eqs. (3) and (4) determines all the incident and transmitted angles and hence the output angle. The surrounding medium is assumed to be air with $N_1 = 1$, $K_1 = 0$ and $\mathbf{e}_\mathrm{i} = \mathbf{f}_\mathrm{i}$. The dependence of the output angle upon $\theta_\mathrm{i}$ and $\tilde{n}$ is shown in Fig. 3, which also demonstrates that total internal reflection (TIR) occurs only in the lossless case [22] when

$$\sin\theta_\mathrm{TIR} = \sqrt{n^2 - 1}\,\sin\gamma - \cos\gamma.$$

For the case shown, $\gamma = 30^\circ$ and $n = 1.5$ so there is no transmission when $\theta_\mathrm{i} < \theta_\mathrm{TIR} = -17.88^\circ$.

*[Figure 2: an illustration of refraction by a lossy wedge of internal angle $\gamma$; the complex refractive index is assumed to be $\tilde{n} = 1.5 + \mathrm{i}\kappa$. The incident ray enters the vertical face at angle $\theta_\mathrm{i}$, travels a horizontal distance $L$ from the apex, and exits the sloping face at output angle $\theta_\mathrm{o}$.]*

*[Figure 3: plots of output angle $\theta_\mathrm{o}$ vs input angle $\theta_\mathrm{i}$ for $\kappa = 0.0, 0.1, 0.2, 0.5, 1.0, 1.5$ on a wedge of internal angle $\gamma = 30^\circ$.]*

Examining the output apparent refractive index in Fig. 4 shows that the transmitted rays are homogeneous only when $\kappa = 0$. This is not surprising since the rays in the lossy wedge retain their inhomogeneity on refraction by virtue of Snell's law (9). The inhomogeneity that occurs for $\kappa = 0$ when $\theta_\mathrm{i} < \theta_\mathrm{TIR}$ is shown by $K > 0$ even though there is no transmission. This is due to the evanescent surface wave on the upper surface of the wedge, i.e. $\theta_\mathrm{o} = 90^\circ$.

*[Figure 4: plots of apparent refractive index $(N, K)$ vs input angle $\theta_\mathrm{i}$ for $\kappa = 0.0, 0.1, 0.2, 0.5, 1.0, 1.5$.]*

Applying (14), allows us to calculate the transmittances (initially neglecting path loss) for rays that pass through the wedge using $T = 1 - |r|^2$. The plots for s and p polarizations, shown in the left and right graphs in Fig. 5, respectively, exhibit the same TIR for zero absorption. In the conventional TIR zone, the two transmittances rise to maxima with increasing $\kappa$ before falling to zero. These maxima occur below $\kappa = 0.6$. This behaviour occurs because of the interplay between decreasing $\theta_\mathrm{o}$ leading to higher transmittance and the increasing $K$ which lowers transmittance.

*[Figure 5: plots of s and p polarized transmittance ($T_\mathrm{s}$, $T_\mathrm{p}$), not including path loss, vs input angle for $\kappa = 0.0, 0.1, 0.2, 0.5, 1.0, 1.5$.]*

To study the attenuation due propagation losses, let the ray's entry point be a distance $L$ from the apex of the wedge. If $D$ is the attenuation factor then we define the dimensionless attenuation coefficient $d$ through

$$D = \exp(-2k_0 L d). \tag{21}$$

In Fig. 6, the dotted curves are the coefficients for standard treatment and the solid curves represent the more rigorous approach. The dotted curves always lie above the corresponding solid ones. This implies that the losses are overestimated by the standard approach. This is particularly significant in the TIR zone $\theta_\mathrm{i} < -17.88^\circ$, when the ray has a longer path through lossy medium.

*[Figure 6: plots of the dimensionless attenuation coefficient $d$ vs input angle for all input angles and $\kappa = 0.1, 0.2, 0.5, 1.0, 1.5$; dotted lines = standard method, solid lines = more rigorous method.]*

By specifying the distance from the ray's entry point to the wedge's apex as $L = 2\lambda/\pi$, the total transmittance can be calculated. Both polarization states show measurable differences in the plots of Fig. 7 between the two ray tracing approaches. The standard method never has any transmission in the TIR zone and always has higher p transmittance when $\theta_\mathrm{i} > \theta_\mathrm{TIR}$. The more rigorous method shows that the p polarized transmittance behaves differently to the s polarized transmittance in the TIR zone — it possesses a point of inflection.

*[Figure 7: plots of the total s and p transmittance coefficient vs input angle for all input angles and $\kappa$ values, for an entry point of $L = 2\lambda/\pi$; dotted lines = standard method, solid lines = more rigorous method.]*

---

## 5. Monte Carlo simulation

A simulation of reflections from a slab that has an upper Gaussian rough surface and a lower plane surface over a half-space of air was modified to incorporate the more rigorous approach. The Monte Carlo simulation is based on micro-facets [5–7] and allows rays to penetrate the surface. The lower plane is situated a distance $H$ below the mean upper surface and any rays that impinge on it are reflected with proper weighting and polarization state. Rays of different polarization are launched at a given incident angle $\theta_\mathrm{i}$ and the polarimetric bidirectional reflectance distribution function (BRDF) is calculated. For each of the four entry polarization states, 160 million rays were launched and over 2440 different surface realizations were generated for an ensemble average of the results. The rays that return from the slab are binned in divisions of $2^\circ$ and $10^\circ$ in polar and azimuthal angles, respectively.

A refractive index of $\tilde{n} = 2.694 + \mathrm{i}0.509$ was chosen together with an incident angle of $40^\circ$. This value of the refractive index was chosen to match that of glass at $\lambda = 10\,\mu\mathrm{m}$ [11, p. 762]. The surface properties were an RMS slope angle of $18.4^\circ$ and RMS height $\lambda/2$ thus giving a correlation length of $3\lambda$. The lower surface was offset by $H = 2\lambda$.

Fig. 8 shows the unpolarized responses as a function of the output polar angle in the plane of incidence. The more rigorous approach gives a response greater by about 6.6% than the standard method. This is because the reflectivity is higher for an absorbing medium. For a smaller value of absorption, $\kappa = 0.2$, Fig. 9 demonstrates that the more rigorous approach yields results that are around 1.1% greater than standard method.

Fig. 10 shows plots that summarize the angular response for all Mueller matrix elements. Overall, the standard method gives Mueller matrix elements that are smaller in magnitude than the more rigorous approach. The main qualitative differences shown by the more rigorous approach is illustrated in Fig. 11. The $M_{24}$, $M_{34}$, $M_{42}$ and $M_{43}$ elements were found to be zero or close to zero in the standard method. The top two plots in Fig. 11 demonstrates that right circularly polarized light reflected by the surface possesses non-zero portion of light linearly polarized at $0^\circ$ and $135^\circ$. The bottom two plots shows that any linearly polarized light has, on reflection, an elliptical polarization state. The ratio $M_{43}/M_{11}$ was found to vary up to nearly 0.2 for $\tilde{n} = 2.694 + \mathrm{i}0.509$; for $\kappa = 0.2$, it had a maximum of almost 0.08.

*[Figure 8: plot of the Mueller matrix element $M_{11}$ in the plane of incidence vs output polar angle for the standard and more rigorous methods; light incident at $\theta_\mathrm{i} = 40^\circ$ on a rough surface of RMS slope angle $18.435^\circ$, RMS height $\lambda/2$ and refractive index $\tilde{n} = 2.694 + \mathrm{i}0.509$.]*

*[Figure 9: similar plot of $M_{11}$ as Fig. 8 but with refractive index $\tilde{n} = 2.694 + \mathrm{i}0.2$.]*

*[Figure 10: plots of all Mueller matrix elements for the standard approach for a selection of azimuthal angles $\varphi = 0, 30, 60, 90$.]*

*[Figure 11: plots of some Mueller matrix elements ($M_{24}$, $M_{34}$, $M_{42}$, $M_{43}$) for the more rigorous approach for a selection of azimuthal angles $\varphi = 0, 30, 60, 90$.]*

---

## 6. Summary and conclusion

A review of the propagation of plane waves between two lossy media was undertaken. The generalized forms of Snell's law and Fresnel formulae were given. A more rigorous method of ray tracing was detailed and then applied to two very different situations.

The low absorption case was also examined and, generally, the more rigorous method shows a deviation from the standard approach by a second order term in $\kappa/n$. This deviation was more pronounced for the reflectivities and transmittivities than the angle of refraction or the path loss coefficient. However, the polarization changes a ray suffers (as encoded in $M_{34}$ and $M_{43}$) from traversing across an interface between lossless and lossy media are directly proportional to $\kappa/n$ in the low absorption limit.

A lossy wedge was used as an example to study the effects of the more rigorous approach and to point out features of including absorption in a correct manner. Total internal reflection was found to occur only for lossless media. The path loss was clearly shown to be greater for the standard method.

A model of a lossy slab with a rough upper surface and a plane lower surface was examined with a polarized Monte Carlo ray tracing computer program. The standard approach and the more rigorous method were implemented separately. The results showed that the new method has greater unpolarized reflectance. In the case of glass with an entry angle of $40^\circ$, the more rigorous response was found to have an increase of about 6% when compared with the standard approach. Also the more rigorous method has non-zero $M_{24}$, $M_{34}$, $M_{42}$ and $M_{43}$ elements of the Mueller matrix because the single reflections has non-zero $M_{42}$ and $M_{43}$ elements. These elements are zero in the standard approach as the Fresnel coefficients are purely real.

We conclude with recommendations for ray tracing involving absorbing media. If polarization information is needed, it is essential to use the more rigorous method for, say, $\kappa/n > 0.01$. For unpolarized ray tracing, the standard approach works well when $\kappa/n < 0.07$. These guidelines are rough estimates for obtaining results within around 1% of the more rigorous method.

This paper has concentrated on reflection and transmission within lossy media but further work that may be interesting is the case of multiple scattering in a lossy turbid media. The Mie solution for light scattering from spheres has been extended to spheres in absorbing media. This solution could be used with the more rigorous ray tracing method detailed herein to enable new results to be found.

## Acknowledgements

This work is supported by the UK Engineering and Physical Science Research Council.

---

## Uncertain readings

- **[?1]** In Eq. (1)-context the complex refractive index is written $\tilde{n} = n + \mathrm{i}k$ in the scan. Here the imaginary index is denoted $k$, but everywhere else in the paper the imaginary part of the refractive index is $\kappa$ (kappa). Almost certainly $\tilde{n} = n + \mathrm{i}\kappa$; the lowercase $k$ appears to be an OCR/typography clash with the wave-number symbol $k$.
- **[?2]** Eq. (14): the s-transmission coefficient is transcribed as $t_\mathrm{s} = 2k_\mathrm{i}/(k_\mathrm{i}+k_\mathrm{t})$ by analogy with the standard Fresnel form and the definition $k \equiv \mathbf{k}\cdot\mathbf{n}$. The scan glyph could be read as an upper-case $K_\mathrm{i}$ in the numerator; the surrounding symbols ($r_\mathrm{s}$, $k_\mathrm{i}$, $k_\mathrm{t}$) all use lower-case $k$, so $2k_\mathrm{i}$ is the consistent reading.
- **[?3]** Eq. (15): the first line is transcribed as $\mathbf{E}(l\mathbf{e}) = \mathbf{E}\exp(\mathrm{i}\,\mathbf{k}\cdot\mathbf{r})$ with $\mathbf{r} = l\mathbf{e}$. The scan appears to abbreviate the exponent argument (possibly printed as $\mathbf{k}\cdot\mathbf{e}$ or $\mathbf{k}\cdot\mathbf{r}$); the second line $\exp(-k_0 K l\cos\alpha)\exp(\mathrm{i}k_0 N l)$ is unambiguous and follows from $\mathbf{k}\cdot\mathbf{r} = l\,k_0(N + \mathrm{i}K\cos\alpha)$ using $\mathbf{e}\cdot\mathbf{f} = \cos\alpha$.
