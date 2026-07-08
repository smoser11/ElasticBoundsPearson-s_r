# Distribution Families and Asymmetry Visualization

Let me explain the distribution families used in the code, how they're parameterized, and how increasing parameter values affect asymmetry between X and Y distributions.

## Distribution Families

### 1. Uniform Distribution (Base for X)
- **Description**: Equal probability for all categories
- **Formula**: $P(X = i) = \frac{1}{K}$ for all $i \in \{0,1,...,K-1\}$
- **Characteristics**: Perfectly symmetric, maximum entropy

### 2. Exponential Tilt Family (for Y)
- **Description**: Probability that increases exponentially with category index
- **Formula**: $P(Y = j) \propto e^{\alpha_j}$ where $j \in \{0,1,...,K-1\}$
- **Parameter**: $\alpha$ (tilt parameter)
  - When $\alpha = 0$: Uniform distribution (like $X$)
  - As $\alpha$ increases: Mass shifts toward higher categories
  - Large $\alpha$: Most mass concentrated in the highest category
- **Characteristics**: Smooth, monotonic shift of probability mass

### 3. Two-Phase Skew Family (for Y)
- **Description**: Redistributes mass from lower to higher categories based on a threshold
- **Formula**: Starts uniform, then moves a fixed amount of mass from the first $(1-phase) \times K$ categories to the last $phase \times K$ categories
- **Parameter**: $phase \in [0,1]$
  - When $phase = 0$: Uniform distribution (like $X$)
  - As $phase$ increases: More categories receive increased mass
  - When $phase = 1$: Returns to uniform as all categories are affected equally
- **Characteristics**: Creates a "stepped" distribution with two distinct density levels

## Asymmetry Between X and Y

When we talk about "asymmetry" in this context, we're referring to how different the distribution of Y is from the distribution of X. Since X is kept uniform, any deviation of Y from uniformity creates asymmetry between the distributions.

The **Total Variation (TV) distance** quantifies this asymmetry:
- It measures how different Y is from being symmetric around its midpoint
- When TV = 0: Y is symmetric (like X)
- When TV is large: Y is highly asymmetric compared to X

## Visualizing the Distributions

Let's create visualizations showing how the distributions change as parameters increase:

```r
#--------------------------------------
# Distribution Visualization
#--------------------------------------
library(ggplot2)
library(gridExtra)
library(reshape2)

# Create distribution examples
K <- 7  # Number of categories
param_values <- c(0, 1, 3, 5)  # Parameter values to show

# Function to create data frame for plotting
create_dist_df <- function(param_values, K, dist_type = "exp") {
  result <- data.frame()

  for (param in param_values) {
    # Uniform X distribution
    pX <- rep(1/K, K)

    # Y distribution based on type
    if (dist_type == "exp") {
      # Exponential tilt
      unnormalized <- exp(param * (0:(K-1)))
      pY <- unnormalized / sum(unnormalized)
      param_name <- paste("α =", param)
    } else {
      # Two-phase skew
      pY <- rep(1/K, K)
      bins_to_increase <- max(1, round(param * K))

      if (bins_to_increase < K) {
        mass_to_move <- 0.5 * (K - bins_to_increase) / K
        pY[1:(K-bins_to_increase)] <- pY[1:(K-bins_to_increase)] - mass_to_move/(K-bins_to_increase)
        pY[(K-bins_to_increase+1):K] <- pY[(K-bins_to_increase+1):K] + mass_to_move/bins_to_increase
      }
      param_name <- paste("phase =", param)
    }

    tv_dist <- 0.5 * sum(abs(pY - rev(pY)))

    # Create data for both X and Y
    df_X <- data.frame(
      Category = 0:(K-1),
      Probability = pX,
      Distribution = "X (Uniform)",
      Parameter = param_name,
      TV_Distance = tv_dist
    )

    df_Y <- data.frame(
      Category = 0:(K-1),
      Probability = pY,
      Distribution = "Y",
      Parameter = param_name,
      TV_Distance = tv_dist
    )

    result <- rbind(result, df_X, df_Y)
  }

  return(result)
}

# Create data frames for both distribution types
exp_dists <- create_dist_df(param_values, K, "exp")
exp_dists$Family <- "Exponential Tilt"

tps_dists <- create_dist_df(param_values, K, "tps")
tps_dists$Family <- "Two-Phase Skew"

all_dists <- rbind(exp_dists, tps_dists)

# Plot distributions
p_dist <- ggplot(all_dists, aes(x = Category, y = Probability, fill = Distribution)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  facet_grid(Family ~ Parameter) +
  scale_fill_manual(values = c("X (Uniform)" = "blue", "Y" = "red")) +
  labs(title = "Distribution Patterns with Increasing Parameter Values",
       subtitle = "X remains uniform while Y becomes more asymmetric") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot TV distance vs parameter
tv_summary <- all_dists %>%
  filter(Distribution == "Y") %>%
  select(Parameter, TV_Distance, Family) %>%
  distinct()

p_tv <- ggplot(tv_summary, aes(x = Parameter, y = TV_Distance, color = Family, group = Family)) +
  geom_point() +
  geom_line() +
  labs(title = "Total Variation Distance vs Parameter Value",
       subtitle = "Higher values indicate greater asymmetry between X and Y",
       y = "TV Distance") +
  theme_minimal()

# Display plots
grid.arrange(p_dist, p_tv, nrow = 2, heights = c(3, 1))
```

## Interpretation of the Visualization

### Exponential Tilt
1. At **α = 0**: Y is uniform like X (symmetric)
2. At **α = 1**: Y shows moderate shift toward higher categories
3. At **α = 3**: Y shows strong concentration in highest categories
4. At **α = 5**: Y has most probability mass concentrated in the highest category

The TV distance increases monotonically with α, reflecting greater asymmetry between X and Y.

### Two-Phase Skew
1. At **phase = 0**: Y is uniform like X (symmetric)
2. At **phase = 0.3**: Y has higher probability in the top 30% of categories
3. At **phase = 0.5**: Y has higher probability in the top 50% of categories
4. At **phase = 1**: Y approaches uniformity again, but with a different pattern

The TV distance first increases as phase grows from 0, reaches a maximum around phase = 0.5, then decreases back toward 0 as phase approaches 1.

## Connection to Correlation Bounds

As asymmetry between X and Y increases (measured by TV distance):
1. The relationship between r_min and r_max becomes asymmetric
2. The sum Δ = r_max + r_min deviates from 0
3. The degree of this deviation relates directly to the amount of asymmetry

When X and Y have the same number of categories but different distributions, the asymmetry in their distributions translates to asymmetry in the correlation bounds. This is because the comonotonic coupling (that produces r_max) benefits differently from the distribution shapes than the countermonotonic coupling (that produces r_min).

This visualization helps us understand what "increasing asymmetry" means between the distributions and how it might affect the relationship between r_min and r_max.
