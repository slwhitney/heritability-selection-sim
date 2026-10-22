# Heritability Is Not Enough: Simulating (No) Response to Selection

A population-genetics simulation in R showing that a trait can be **heritable** *and* **correlated with fitness**, yet still **fail to evolve**.

The classic breeder's equation (R = h²S) predicts that a heritable trait under selection will shift over generations. This project builds a scenario designed to *break* that prediction, which shows why heritability is a necessary but not sufficient condition for evolutionary change.

📄 **Full write-up with code and figures:** [Final.pdf](Final.pdf)

---

## The model

Each individual's phenotype is the sum of three components:

```
P = X − N + E
```

| Component | Meaning | Variance | Inherited? |
|---|---|---|---|
| **X** | Additive genetic value | 0.40 | ✅ Yes |
| **N** | Nutrition | 0.40 | ❌ No (fresh each generation) |
| **E** | Residual environmental noise | 0.20 | ❌ No |

This gives a heritability of **h² = 0.4**. The catch is that **fitness depends only on nutrition (N)**, which is not passed on to offspring. A population of 2,000 individuals is followed for 20 generations of fitness-weighted reproduction under two scenarios:

1. **Linear selection:** fitness `w = 1 + 0.3·N`
2. **Stabilizing selection:** fitness is a Gaussian centered on an optimal nutrition level, `w = exp(−(N − 0.5)² / 2ω²)` with ω = 0.5

## Results

### Scenario 1: Linear selection

<p align="center"><img src="figures/linear_selection.png" width="85%" alt="Linear selection: heritability, phenotype-fitness covariance, and mean phenotype over 20 generations"></p>

- **(A)** Parent–offspring regression recovers the expected heritability (slope ≈ **0.39**).
- **(B)** Phenotype and fitness are clearly correlated: **Cov(P, w) = −0.115**.
- **(C)** The breeder's equation would predict the mean phenotype to fall by about 0.05 per generation (about −0.9 after 20 generations). Instead, the mean **just drifts around zero**.

**Why:** the phenotype–fitness covariance comes entirely from nutrition, which is not inherited. Selection picks well-fed parents, not parents with different genes, so there is no covariance between *genetic value* and fitness and therefore no evolutionary response.

### Scenario 2: Stabilizing selection

<p align="center"><img src="figures/stabilizing_heritability_covariance.png" width="85%" alt="Stabilizing selection: heritability and phenotype-fitness relationship"></p>
<p align="center"><img src="figures/stabilizing_no_response.png" width="70%" alt="Mean phenotype under stabilizing selection over 20 generations"></p>

- **(D)** Heritability is unchanged (slope ≈ **0.39**).
- **(E)** Fitness peaks at an intermediate phenotype, with **Cov(P, w) = −0.154**.
- **(F)** Again, **no directional response**. The mean phenotype fluctuates near zero across all 20 generations.

## Takeaway

Heritability tells you how much offspring resemble their parents. It does not tell you whether a population will evolve. A response to selection requires fitness to covary with the **heritable** part of the trait. When selection acts through environmental factors, even a strongly heritable trait stays put.

---

## Repository contents

```
├── selection_simulation.R   # Full simulation and plotting code
├── Final.pdf                # Knitted R Markdown report
└── figures/                 # Figures from the report
```

## How to run

Requires R with **tidyverse** and **patchwork**:

```r
install.packages(c("tidyverse", "patchwork"))
source("selection_simulation.R")
```

The simulation is stochastic, so exact values change slightly between runs. Add `set.seed()` at the top for reproducible output.
