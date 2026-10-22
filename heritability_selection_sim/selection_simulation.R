# Heritability Is Not Enough: simulating (no) response to selection
# Sarah Whitney, May 2025
#
# Phenotype P = X - N + E
#   X = additive genetic value (heritable)
#   N = nutrition (environmental, NOT heritable)
#   E = residual environmental noise
# Fitness depends only on nutrition N, so P is heritable AND correlated with
# fitness, yet selection acts on a non-heritable component.

library(tidyverse)
library(patchwork)

# ---- Parameters -------------------------------------------------------------
N      <- 2000   # population size
var_X  <- 0.40   # additive genetic variance
var_N  <- 0.40   # nutrition variance
var_E  <- 0.20   # residual environmental variance
s      <- 0.30   # strength of the fitness effect of nutrition
n_gen  <- 20     # generations to follow

h2 <- var_X / (var_X + var_N + var_E)   # realised heritability
h2

# ---- Helper functions -------------------------------------------------------
## one generation, returns a data.frame of individuals
draw_generation <- function(N, var_X, var_N, var_E) {
  tibble(
    X = rnorm(N, 0, sqrt(var_X)),
    N = rnorm(N, 0, sqrt(var_N)),
    E = rnorm(N, 0, sqrt(var_E))
  ) %>%
    mutate(
      P = X - N + E,
      w = 1 + s * N       # linear, always positive if s < ~3/SD(N)
    )
}

## mid-parent regression slope approx h2
midparent_offspring <- function(df) {
  # random pairing
  moms <- sample(nrow(df))
  dads <- sample(nrow(df))
  midP <- (df$P[moms] + df$P[dads]) / 2
  midX <- (df$X[moms] + df$X[dads]) / 2

  # offspring inherit only X; nutrition & E are fresh draws
  X_off <- rnorm(length(midX), mean = midX, sd = sqrt(var_X / 2))
  N_off <- rnorm(length(midX), 0, sqrt(var_N))
  E_off <- rnorm(length(midX), 0, sqrt(var_E))
  P_off <- X_off - N_off + E_off

  tibble(midP, P_off)
}

# ---- Scenario 1: linear selection on nutrition ------------------------------
gen0  <- draw_generation(N, var_X, var_N, var_E)
regdf <- midparent_offspring(gen0)
slp   <- coef(lm(P_off ~ midP, data = regdf))[2]

pA <- ggplot(regdf, aes(midP, P_off)) +
  geom_point(alpha = 0.3, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "darkblue") +
  labs(x = "Mid-parent phenotype",
       y = "Offspring phenotype",
       title = "(A) Heritability") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6,
           label = sprintf("slope approx h2 = %.2f", slp)) +
  theme_minimal(base_size = 11)

cov_pw <- cov(gen0$P, gen0$w)

pB <- ggplot(gen0, aes(P, w)) +
  geom_point(alpha = 0.3, colour = "tomato") +
  geom_smooth(method = "lm", se = FALSE, colour = "firebrick") +
  labs(x = "Phenotype",
       y = "Fitness",
       title = "(B) Covariance between P and w") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6,
           label = sprintf("Cov(P , w) = %.3f", cov_pw)) +
  theme_minimal(base_size = 11)

mean_P <- numeric(n_gen + 1)
gen    <- gen0
mean_P[1] <- mean(gen$P)

for (g in 1:n_gen) {
  ## selection proportional to w
  sel_idx <- sample(nrow(gen), N, replace = TRUE, prob = gen$w)
  parents <- gen[sel_idx, ]

  ## offspring genetics come only from X (nutrition/env. are new each time)
  X_off <- rnorm(N, parents$X, sd = sqrt(var_X / 2))
  N_off <- rnorm(N, 0, sqrt(var_N))
  E_off <- rnorm(N, 0, sqrt(var_E))

  gen <- tibble(X = X_off, N = N_off, E = E_off,
                P = X_off - N_off + E_off,
                w = 1 + s * N_off)
  mean_P[g + 1] <- mean(gen$P)
}

pC <- tibble(Generation = 0:n_gen, MeanP = mean_P) %>%
  ggplot(aes(Generation, MeanP)) +
  geom_line(colour = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "(C) Mean phenotype across 20 generations",
       x = "Generation", y = expression(bar(P))) +
  theme_minimal(base_size = 11)

(pA | pB) / pC

# ---- Scenario 2: stabilizing selection on nutrition -------------------------
omega <- 0.5   # strength of stabilizing selection
N_opt <- 0.5   # optimal nutrition value

stab_w <- function(N) {
  exp(-((N - N_opt)^2) / (2 * omega^2))
}

draw_generation_stab <- function(N, var_X, var_N, var_E) {
  tibble(
    X = rnorm(N, 0, sqrt(var_X)),
    N = rnorm(N, 0, sqrt(var_N)),
    E = rnorm(N, 0, sqrt(var_E))
  ) %>%
    mutate(
      P = X - N + E,
      w = stab_w(N)
    )
}

gen0_stab  <- draw_generation_stab(N, var_X, var_N, var_E)
regdf_stab <- midparent_offspring(gen0_stab)

slp_stab    <- coef(lm(P_off ~ midP, data = regdf_stab))[2]
cov_pw_stab <- cov(gen0_stab$P, gen0_stab$w)

pA_stab <- ggplot(regdf_stab, aes(midP, P_off)) +
  geom_point(alpha = 0.3, colour = "darkcyan") +
  geom_smooth(method = "lm", se = FALSE, colour = "darkblue") +
  labs(x = "Mid-parent Phenotype", y = "Offspring Phenotype",
       title = "(D) Heritability (Stabilizing)") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("slope approx h2 = %.2f", slp_stab)) +
  theme_minimal()

pB_stab <- ggplot(gen0_stab, aes(P, w)) +
  geom_point(alpha = 0.3, colour = "purple") +
  geom_smooth(method = "loess", se = FALSE, colour = "pink") +
  labs(x = "Phenotype", y = "Fitness",
       title = "(E) Phenotype-Fitness Covariance (Stabilizing)") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("Cov(P,w) = %.3f", cov_pw_stab)) +
  theme_minimal()

pA_stab | pB_stab

mean_P_stab <- numeric(n_gen + 1)
gen_stab    <- gen0_stab
mean_P_stab[1] <- mean(gen_stab$P)

for (g in 1:n_gen) {
  sel_idx <- sample(nrow(gen_stab), N, replace = TRUE, prob = gen_stab$w)
  parents <- gen_stab[sel_idx, ]

  X_off <- rnorm(N, parents$X, sd = sqrt(var_X / 2))
  N_off <- rnorm(N, 0, sqrt(var_N))
  E_off <- rnorm(N, 0, sqrt(var_E))

  gen_stab <- tibble(X = X_off, N = N_off, E = E_off,
                     P = X_off - N_off + E_off,
                     w = stab_w(N_off))

  mean_P_stab[g + 1] <- mean(gen_stab$P)
}

pC_stab <- tibble(Generation = 0:n_gen, MeanP = mean_P_stab) %>%
  ggplot(aes(Generation, MeanP)) +
  geom_line(colour = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "(F) No Response Under Stabilizing Selection",
       x = "Generation", y = "Mean Phenotype") +
  theme_minimal()

pC_stab

# ---- Combined figure --------------------------------------------------------
(pA | pB) / pC /
  (pA_stab | pB_stab) / pC_stab
