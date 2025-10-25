# Project VI: Growth of Budworms in a forest

Ludwig and coworkers proposed a simple model to quantify the population of spruce budworms in balsam fir forest. The state of the forest is characterized by the variable $S(t)$ - the average height of the trees in the forest and the general energy reserve in the forest $E(t)$ - this is an effective measurement of how healthy the forest actually is. In the presence of constant Budworm population $B$, the following equations govern $S$ and $E$:

$$
\begin{aligned}
\dot{S} &= r_S S \left(1 - \frac{K_E S}{K_S E}\right) \\
\dot{E} &= r_E E \left(1 - \frac{E}{K_E}\right) - P \frac{B}{S}
\end{aligned}
$$

where $r_S$, $r_E$, $K_S$, $K_E$, and $P$ are all positive constants.

---

## Project Objectives

(a) Nondimensionalize the model. How many dimensionless groups are needed?

(b) Find out all the fixed points.

(c) Solve the above equations numerically using a 4th order Runge-Kutta method (or any other appropriate numerical method). You are required to use **MATLAB** or **Python** for solving the problem. You are **not** allowed to use any in-built solvers (e.g., `ode45` in MATLAB) to solve the ODEs.

Plot $E(t)$ vs $S(t)$ (using numerical solutions) - this is called as the **Phase portrait** of the system and it helps us understand the long term behaviour of the system. Assume various ranges of the parameter values. Focus particularly on large Budworm population.
