#incremental MPM

![Example 1, no material model](https://nialltl.neocities.org/articles/img/mpm_guide/mpm_first_pass.gif)

![Example 2, Neo-Hookean elasticity](https://nialltl.neocities.org/articles/img/mpm_guide/mpm_neohookean.gif)

![Example 3, fluids](https://nialltl.neocities.org/articles/img/mpm_guide/mpm_fluid_constitutive_model.gif)

#overview

this is an MIT-licensed commented implementation of the MLS-MPM algorithm for simulating soft bodies and fluids. it's intended to be accompanied by this article: [notes and examples for the material point method](https://nialltl.neocities.org/articles/mpm_guide.html).

all examples are single self-contained files with comments. this package contains 3 scenes:

* a stripped-back example demonstrating how particle-grid transfers work in MLS-MPM

* non-linear elasticity using a Neo-Hookean model

* liquid simulation using a constitutive model for isotropic newtonian fluids

if you have any questions, they might be answered by [the article]((https://nialltl.neocities.org/articles/mpm_guide.html), but if not feel free to contact me [@nialltl](https://twitter.com/nialltl) :o)

**Unity version:** 2018.3.10f1 — Untested in other versions.

#license

MIT License

Copyright (c) 2019 niall tl

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.