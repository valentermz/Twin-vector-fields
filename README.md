# Twin vector fields

This is a `SageMath` code used to compute twin vector fields. These computations follow the ideas presented in the following paper:

* Ramirez, V. Twin vector fields and independence of spectra for quadratic vector fields. *J. Dynam. Control Syst.* **23** (2017), no. 3, 623–633.

DOI: [10.1007/s10883-016-9344-5](http://doi.org/10.1007/s10883-016-9344-5), arXiv: [1508.02413](https://arxiv.org/abs/1508.02413)

## The theorem in question

We work with quadratic vector fields 

```
P(x,y)\frac{\partial}{\partial x} + Q(x,y)\frac{\partial}{\partial y}
```

on `\mathbb{C}^2` having isolated singularities. 

**Definition:** We say that two vector fields *v* and *v'* are *twin vector fields* if they have the same singular set, and for each singular point they the linearization matrices have the same eigenvalues.

**Theorem:** A generic quadratic vector field has exactly one twin.
