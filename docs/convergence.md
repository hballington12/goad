# Convergence

When computing orientation averaged scattering, it's not usually known beforehand exactly how many orientations are required to converge on the desired result. GOAD's solution to this is called a `Convergence`, which uses Westford's algorithm to track the mean and variance of one or more prescribed convergence variables. The simulation runs until the convergence criteria are met, or some maximum number of orientations is reached. A simple example:

{{code_block('examples/convergence', 'basic')}}

which produces the following output:

```bash
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GOAD: [Convergence]  ▰▱▱▱▱▱▱

[Orientations: 158 (max 10000)] [0.010 sec/orientation] [Minimum orientations check: ✓]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Asymmetry       0.7726 ± 0.0154 ━━━━━━━━━━━━━━━━━━━━━━━━━ 100% [SEM: 2.00% / 2.00%]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```
