# Boltzmann Integrals
Mathematica code for integrating even polynomials of speed in Boltzmann collision integrals.

All the code is in the notebook.
One can choose some predefined polynomials, corresponding to aggregation or bouncing/restitution, or define the polynomial directly.

All calculations are performed step-by-step.

Allows to find the integrals for both equal and unequal flow speeds, masses and temperatures.

Can be used to construct differential equations for number density, temperature, viscosity, thermal conductivity and Sonine polynomial expansion (corresponding polynomials are predefined).

## Definition

Let us define the input polynomial $I(v_i,v_i',v_j,v_j',v_k)$, which is written in the notebook as **Integral**, and the output $\color{orange} \rm {RESULT}$ as $R(I)$. Then the notebook calculates the following expression:
<p style="text-align: center;">
$$
\begin{aligned}
{\Large R(I)} & = \sigma_{ij} \iiint d \vec {v_i} d \vec {v_j} d \vec e \, {\Large I} \, {\Large |} \cdot ( \vec {v_i} - \vec {v_j} + \vec {u_i} - \vec {u_j} ) \cdot \vec e {\Large |} \, \Theta ( - ( \vec {v_i} - \vec {v_j} + \vec {u_i} - \vec {u_j} ) \cdot \vec e ) \\
& \times \Theta {\Huge (} Q_{ij} - \frac{| \vec {v_i} - \vec {v_j} + \vec {u_i} - \vec {u_j} |^2}{2 ( \theta_i + \theta_j )} {\Huge )} \frac{1}{( 2 \pi \theta_i )^{3/2}} e^{\frac{v_i^2}{2 \theta_i}} \frac{1}{( 2 \pi \theta_j )^{3/2}} e^{\frac{v_j^2}{2 \theta_j}}
\end{aligned}
$$
</p>

which describes the change of the average value of the input during a unit time interval due to collisions of size _**i**_ and size _**j**_ particles in the multicomponent dilute granular gas, experiencing aggregation. Here $\vec{v_i}, \vec{v_j}, \vec{v_i^'}, \vec{v_j^'}$ are pre- and post-collision (in bouncing case) velocities of size _**i**_ and _**j**_ particles in the reference frames, corresponding to the center of mass of all similarly sized particles. In case they aggregate, the resulting particle of size _**k = i + j**_ has the velocity $\vec {v_k}$. Centers of mass of size _**i**_ and _**j**_ particles move with the speeds $\vec {u_i}$ and $\vec {u_j}$, which are called the flow speeds. $\sigma_{ij}$ is the radius of the collision cylinder and $d \vec e \, {\Large I} \, {\Large |} ( \vec {v_i} - \vec {v_j} + \vec {u_i} - \vec {u_j} ) \cdot \vec e {\Large |}$ is its height in the unit collision direction $\vec e$. Size _**i**_ and _**j**_ particles are assumed to have Maxwell distribution with the scaled temperatures (speed variance) of $\theta_i$ and $\theta_j$. The first Heaviside theta-function limits the collision direction to be aligned with the relative speed and the second theta-function is used to distinguish bouncing and aggregating collisions. Namely, the value of $Q_{ij}$ is proportional to the ratio of the energy of the potential barrier and the relative kinetic energy of colliding particles. When it is high, particles cannot penetrate it after the collision and stick together. When bouncing (restitutive) collisions are used, one should consider the opposite sign in the theta-function by modifying the corresponding variable (called **restitution**) in the input. Non-Maxwellian speed distributions can be accounted for by correspondingly changing the **Integral**.

## Examples

The above expression allows us to derive many types of various Smoluchowski equations. All examples are denoted by various letters in the input part of the notebook. Here we describe the meaning behind them.

As a simplest example, one can use **Integral = c** = 1. Then one can obtain the kernel for the classic ballistic Smoluchowski equations in case of unequal partial flow speeds and temperatures and nonzero bouncing probability. Namely, classic Smoluchowski equations are written as
<p style="text-align: center;">
$$\frac{d}{dt} n_k = \frac{1}{2} \sum_{i +j = k} C_{ij} n_i n_j - \sum_{j = 1}^{\infty} C_{kj} n_k n_j,$$
</p>

where $n_i$ is the number density of size _**i**_ particles and the notebook calculates the coagulation kernel
<p style="text-align: center;">
$$C_{ij} = R(\mathbf{c}) = R(1).$$
</p>

Next, we can find similar expressions for the change of scaled temperatures (speed variance) and get the temperature-dependent Smoluchowski equations

<p style="text-align: center;">
$$\frac{d}{dt} (n_k \theta_k) = \frac{1}{2} \sum_{i +j = k} B_{ij} n_i n_j - \sum_{j = 1}^{\infty} D_{kj}^{\rm {agg}} n_k n_j - \sum_{j = 1}^{\infty} D_{kj}^{\rm {res}} n_k n_j,$$
</p>

where the kernels are calculated as $B_{ij} = R(\mathbf{b}), D_{ij}^{\rm {agg}} = R(\mathbf{dagg}), D_{ij}^{\rm {res}} = R(\mathbf{dres})$. Finally, we can integrate over average speed change to get

<p style="text-align: center;">
$$\frac{d}{dt} (n_k \vec {u_k}) = \frac{1}{2} \sum_{i +j = k} \vec P_{ij} n_i n_j - \sum_{j = 1}^{\infty} \vec R_{kj}^{\rm {agg}} n_k n_j - \sum_{j = 1}^{\infty} \vec R_{kj}^{\rm {res}} n_k n_j,$$
</p>

Incorporating number density, flow and temperature gradients in the above generalized Smoluchowski equations leads to the so-called Smoluchowski-Euler equations.

Next, one can consider changes to the Maxwell distribution due to these gradients, which allows to derive the equations for viscosities $\eta_{kl}$ and thermal conductivities $\kappa_{kl}$.

Chapman-Enskog approach leads to the distributions of the form
<p style="text-align: center;">
$$\begin{aligned}
f_k ( \vec {v_k} ) = \frac{1}{(2 \pi \theta_k)^{3/2}} e^{- \frac{|\vec {v_k} - \vec {u_k}|^2}{2 \theta_k}} {\Large (} 1 & - \sum_{i = 1}^{\infty} \frac{\eta_{ki}}{n_k T_k \theta_k} {\large (} ( \vec {v_k} \cdot \vec \nabla) (\vec {v_k} \cdot \vec {u_i}) - \frac{1}{3} v_k^2 (\vec \nabla \cdot \vec {u_i}) {\large )} \\
& - \sum_{i = 1}^{\infty} \frac{\kappa_{ki}}{n_k \theta_k} (\vec {v_k} \cdot \vec \nabla \log \theta_i) {\large (} \frac{v_k^2}{5 \theta_k} - 1 {\large )} {\Large )}
\end{aligned}$$
</p>

Let us define
<p style="text-align: center;">
$$\gamma_{kl} = - \frac{\eta_{kl}}{n_k T_k}, \quad \alpha_{kl} = - \frac{2\kappa_{kl}}{5n_k \theta_k}.$$
</p>

These new quantities then satisfy the following generalized Smoluchowski equations:
$\begin{aligned}
  7.5 \frac{\partial \left( \alpha_{kl} n_k \theta_k \right)}{\partial t} + 7.5 \delta_{kl} n_k \theta_k & = \frac{1}{2} \sum_{i+j=k} \left( \alpha_{il} B_{ij}^{\alpha} + \alpha_{jl} B_{ji}^{\alpha} \right) n_i n_j - \sum_{j = 1}^{\infty} \left( \alpha_{kl} C_{kj}^{\alpha 1} + \alpha_{jl} C_{kj}^{\alpha 2} \right) n_k n_j \\
  & + \sum_{j = 1}^{\infty} \left( \alpha_{kl} D_{kj}^{\alpha 1} + \alpha_{jl} D_{kj}^{\alpha 2} \right) n_k n_j
\end{aligned}$
and
$\begin{aligned}
  10 \frac{\partial \left( \gamma_{kl} n_k \right)}{\partial t} + 10 \delta_{kl} n_k & = \frac{1}{2} \sum_{i+j=k} \left( \gamma_{il} B_{ij}^{\gamma} + \gamma_{jl} B_{ji}^{\gamma} \right) n_i n_j - \sum_{j = 1}^{\infty} \left( \gamma_{kl} C_{kj}^{\gamma 1} + \gamma_{jl} C_{kj}^{\gamma 2} \right) n_k n_j \\
  & + \sum_{j = 1}^{\infty} \left( \gamma_{kl} D_{kj}^{\gamma 1} + \gamma_{jl} D_{kj}^{\gamma 2} \right) n_k n_j,
\end{aligned}$
where the kernels are calculated as
$\begin{aligned}
& B_{ij}^{\alpha} = R (\mathbf{b}_{\alpha}), & \quad & C_{ij}^{\alpha \text{1,2}} = R (\mathbf{c}_{\alpha \text{1,2}}), & \quad & D_{ij}^{\alpha \text{1,2}} = R (\mathbf{d}_{\alpha \text{1,2}}), \\
& B_{ij}^{\gamma} = R (\mathbf{b}_{\gamma}), & \quad & C_{ij}^{\gamma \text{1,2}} = R (\mathbf{c}_{\gamma \text{1,2}}), & \quad & D_{ij}^{\gamma \text{1,2}} = R (\mathbf{d}_{\gamma \text{1,2}}).
\end{aligned}$

Finally, one can use Sonin polynomial expansion to approximate the speed distribution in homogeneous systems. If one uses only the second Sonin polynomial, so that

<p style="text-align: center;">
$$f_k ( \vec {v_k} ) = \frac{1}{(2 \pi \theta_k)^{3/2}} e^{- \frac{v_k^2}{2 \theta_k}} {\Huge (} 1 - \alpha_k {\huge (} \frac{v_k^4}{30 \theta_k} - \frac{v_k^2}{3 \theta_k} + \frac{1}{2} {\huge)} {\Huge )},$$
</p>

then the temperature-dependent kernels should be instead calculated as

<p style="text-align: center;">
$$\begin{aligned}
& C_{ij} = R (\mathbf{csonin}), & \quad & B_{ij} = R (\mathbf{bsonin}), \\
& D_{ij}^{\rm {agg}} = R(\mathbf{daggsonin}), & \quad & D_{ij}^{\rm {res}} = R(\mathbf{dressonin})
\end{aligned}$$
</p>

and we have an additional system for Sonin coefficients

<p style="text-align: center;">
$$\frac{2}{15} \frac{d}{dt} \left( n_k \alpha_k \right) = \frac{1}{2} \sum_{i + j = k} P_{ij} \left( \theta_i, \theta_j, \alpha_i, \alpha_j \right) n_i n_j - \sum_{j = 1}^{\infty} R_{kj} \left( \theta_k, \theta_j, \alpha_k, \alpha_j \right) n_k n_j$$
</p>

with the kernels
<p style="text-align: center;">
$$P_{ij} = R (\mathbf{psonin}), \quad R_{ij} = R (\mathbf{raggsonin}) + R (\mathbf{rressonin}).$$
</p>
