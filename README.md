Incremental MPM
======

**Accompanying article**: [notes and examples for the material point method](https://nialltl.neocities.org/articles/mpm_guide.html)

<img src="https://nialltl.neocities.org/articles/img/mpm_guide/mpm_neohookean.gif" width="50%" /><img src="https://nialltl.neocities.org/articles/img/mpm_guide/mpm_fluid_constitutive_model.gif" width="50%" />

**Unity version:** 2018.3.10f1 — Untested in other versions.

overview
=======

this is an MIT-licensed commented implementation of the MLS-MPM algorithm for simulating soft bodies and fluids. it's intended to accompany [this article](https://nialltl.neocities.org/articles/mpm_guide.html). 

this project is made in **Unity**, and makes use of the **Burst compiler** and **High-Performance C#**. it's partially parallelised to run at interactive rates, but this is intended more as a learning reference than a final implementation.

all examples are single-file and self-contained. this package contains 3 scenes:

* a stripped-back example demonstrating how particle-grid transfers work in MLS-MPM

* non-linear elasticity using a Neo-Hookean model

* liquid simulation using a constitutive model for isotropic newtonian fluids

if you have any questions, they might be answered by [the article](https://nialltl.neocities.org/articles/mpm_guide.html), but if not feel free to contact me [@nialltl](https://twitter.com/nialltl) :o)
