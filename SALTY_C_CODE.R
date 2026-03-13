##########################################
# Organic Matter AFDM vs Salinity GLMM
# Periphyton is "algae" in-script
##########################################

# -------------------------------
# 1. Load libraries
# -------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 2. Read data
# -------------------------------
setwd("C:/Users/lisat/OneDrive/Desktop/Salty_Carbon_Data_and_Scripts")

df <- read.csv("SALTY_C_AFDM_YEAR_1_ALL.csv", check.names=FALSE) %>%
  mutate(
    Site  = factor(Site),
    OM    = factor(OM),  # e.g., CBOM, FBOM, Periphyton
    month = factor(Season, levels=c("February","May","August","October")),
    SC    = as.numeric(SC),
    AFDM  = as.numeric(AFDM)
  )

# -------------------------------
# 3. Plot aesthetics
# -------------------------------
OM_colors <- c("CBOM"="#E69F00", "FBOM"="#0072B2", "Periphyton"="#228B22")
OM_shapes <- c("CBOM"=16, "FBOM"=17, "Periphyton"=15)

month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2023)"
)

SC_seq <- seq(min(df$SC, na.rm=TRUE), max(df$SC, na.rm=TRUE), length.out=100)

# -------------------------------
# 4. Fit full model
# -------------------------------
mod_AFDM <- lmer(AFDM ~ SC * month * OM + (1 | Site), data=df)

# -------------------------------
# 5. Extract all SC slopes (table)
# -------------------------------
AFDM_SC_slopes_all <- emtrends(mod_AFDM, ~ month | OM, var="SC") %>%
  summary(infer=TRUE) %>%
  as.data.frame() %>%
  mutate(
    t_value = SC.trend / SE,
    p_value = 2 * pt(-abs(t_value), df),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    SC_slope  = round(SC.trend, 2),
    Std_Error = round(SE, 2),
    t_value   = round(t_value, 2),
    p_value   = signif(p_value, 3),
    lower.CL  = round(lower.CL, 2),
    upper.CL  = round(upper.CL, 2)
  ) %>%
  select(OM, month, SC_slope, Std_Error, df, t_value, p_value, lower.CL, upper.CL, Significance)

# -------------------------------
# Print full slopes table
# -------------------------------
print("=== AFDM ~ SC Slopes by OM and Month ===")
print(AFDM_SC_slopes_all)

# -------------------------------
# Save table for manuscript
# -------------------------------
write.csv(AFDM_SC_slopes_all, "AFDM_SC_slopes_by_OM_and_Month.csv", row.names=FALSE)

# -------------------------------
# 6. Extract only significant slopes for plotting
# -------------------------------
AFDM_SC_slopes_sig <- AFDM_SC_slopes_all %>% filter(Significance != "ns")

pred_sig <- emmeans(mod_AFDM, ~ SC | month * OM, at=list(SC=SC_seq), cov.reduce=FALSE) %>%
  as.data.frame() %>%
  semi_join(AFDM_SC_slopes_sig, by=c("month","OM"))

# -------------------------------
# 7. Plot significant slopes
# -------------------------------
p_AFDM_sig <- ggplot(df, aes(SC, AFDM, color=OM, shape=OM)) +
  geom_jitter(width=0.1, size=3, alpha=0.7) +
  geom_ribbon(data=pred_sig,
              aes(x=SC, ymin=lower.CL, ymax=upper.CL, group=interaction(month,OM)),
              inherit.aes=FALSE, fill="grey70", alpha=0.3) +
  geom_line(data=pred_sig,
            aes(x=SC, y=emmean, group=interaction(month,OM)),
            inherit.aes=FALSE, color="black", size=1.2) +
  facet_grid(OM ~ month, scales="free_y", labeller=labeller(month=month_labels)) +
  scale_color_manual(values=OM_colors) +
  scale_shape_manual(values=OM_shapes) +
  theme_bw(base_size=16) +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=14, face="bold")) +
  labs(x="Specific Conductivity (µS/m)", y="AFDM (g/m²)")

# -------------------------------
# Print figure
# -------------------------------
print(p_AFDM_sig)

##########################################
# Periphyton (algae in-script) CHLA vs Salinity GLMM 
##########################################

# -------------------------------
# 1. Load libraries
# -------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 2. Read and clean data
# -------------------------------
setwd("C:/Users/lisat/OneDrive/Desktop/Salty_Carbon_Data_and_Scripts")

df <- read.csv("SALTY_C_CHLA_ROCKS.csv", check.names=FALSE) %>%
  mutate(
    Site  = factor(Site),
    month = factor(Season, levels=c("February","May","August","October")),
    SC    = as.numeric(SC),
    CHLA  = as.numeric(CHLA)
  )

# -------------------------------
# 3. Prediction sequence
# -------------------------------
SC_seq <- seq(min(df$SC, na.rm=TRUE), max(df$SC, na.rm=TRUE), length.out=100)

# -------------------------------
# 4. Fit GLMM
# -------------------------------
mod_CHLA <- lmer(CHLA ~ SC * month + (1 | Site), data=df)

# -------------------------------
# 5. Extract SC slopes from model
# -------------------------------
CHLA_SC_slopes_all <- emtrends(mod_CHLA, ~ month, var="SC") %>%
  summary(infer=TRUE) %>%
  as.data.frame() %>%
  mutate(
    SC_slope  = signif(SC.trend, 3),
    Std_Error = signif(SE, 3),
    t_value   = signif(t.ratio, 3),
    p_value   = signif(p.value, 3),
    lower.CL  = signif(lower.CL, 3),
    upper.CL  = signif(upper.CL, 3),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(month, SC_slope, Std_Error, df, t_value, p_value, lower.CL, upper.CL, Significance)

# -------------------------------
# Print full slopes table
# -------------------------------
print("=== CHLA ~ SC Slopes by Month ===")
print(CHLA_SC_slopes_all)

# -------------------------------
# Save table for manuscript
# -------------------------------
write.csv(CHLA_SC_slopes_all, "CHLA_SC_slopes_by_month.csv", row.names=FALSE)

# -------------------------------
# 6. Extract only significant slopes for plotting
# -------------------------------
CHLA_SC_slopes_sig <- CHLA_SC_slopes_all %>% filter(Significance != "ns") %>% select(month)

pred_sig <- emmeans(mod_CHLA, ~ SC | month, at=list(SC=SC_seq), cov.reduce=FALSE) %>%
  as.data.frame() %>%
  semi_join(CHLA_SC_slopes_sig, by="month")

# -------------------------------
# 7. Plot significant slopes
# -------------------------------
algae_color <- "#228B22"  # Periphyton

month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2024)"
)

p_CHLA_sig <- ggplot(df, aes(x=SC, y=CHLA)) +
  geom_jitter(color=algae_color, size=3, width=0.1, alpha=0.7) +
  geom_ribbon(data=pred_sig,
              aes(x=SC, ymin=lower.CL, ymax=upper.CL),
              fill=algae_color, alpha=0.3, inherit.aes=FALSE) +
  geom_line(data=pred_sig, aes(x=SC, y=emmean),
            color=algae_color, size=1.2, inherit.aes=FALSE) +
  facet_wrap(~ month, scales="fixed", labeller=labeller(month=month_labels)) +
  theme_bw(base_size=16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=14, face="bold")
  ) +
  labs(
    x = "Specific Conductivity (µS/m)",
    y = expression("Periphyton chlorophyll " * italic(a) * " concentration (g/"*m^2*")")
  )

# -------------------------------
# Print figure
# -------------------------------
print(p_CHLA_sig)

##########################################
# Organic Matter Isotopes GLMM
##########################################

# -------------------------------
# 1. Load libraries
# -------------------------------
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 2. Set working directory
# -------------------------------
setwd("C:/Users/lisat/OneDrive/Desktop/Salty_Carbon_Data_and_Scripts")

# -------------------------------
# 3. Read data
# -------------------------------
df <- read.csv("SALTY_C_SI_YEAR_1_ALL.csv", check.names = FALSE)

# -------------------------------
# 4. Data preparation
# -------------------------------
df <- df %>%
  mutate(
    Season = factor(Season, levels = c("February", "May", "August", "October")),
    Site   = factor(Site),
    OM     = factor(OM),
    SC     = as.numeric(SC),
    d13C   = as.numeric(d13C),
    d15N   = as.numeric(d15N)
  )

# -------------------------------
# 5. Fit full interaction models
# -------------------------------
mod_d13C <- lmer(d13C ~ SC * Season * OM + (1 | Site), data = df)
mod_d15N <- lmer(d15N ~ SC * Season * OM + (1 | Site), data = df)

# -------------------------------
# 6. Function to extract SC slopes for each Season × OM
# -------------------------------
extract_slopes_all <- function(model, isotope_name) {
  emtrends(model, ~ Season | OM, var = "SC") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(
      t_value = SC.trend / SE,
      p_value = 2 * pt(-abs(t_value), df),
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      Isotope   = isotope_name,
      SC_slope  = signif(SC.trend, 3),
      Std_Error = signif(SE, 3),
      t_value   = signif(t_value, 3),
      p_value   = signif(p_value, 3)
    ) %>%
    select(Isotope, OM, Season, SC_slope, Std_Error, df, t_value, p_value, Significance)
}

# -------------------------------
# 7. Extract slopes for δ13C and δ15N
# -------------------------------
d13C_SC_slopes <- extract_slopes_all(mod_d13C, "d13C")
d15N_SC_slopes <- extract_slopes_all(mod_d15N, "d15N")

# -------------------------------
# 8. Print tables separately
# -------------------------------
cat("\n=== δ13C SC Slopes (All) ===\n")
print(d13C_SC_slopes)

cat("\n=== δ15N SC Slopes (All) ===\n")
print(d15N_SC_slopes)

# -------------------------------
# 9. Salinity prediction sequence
# -------------------------------
SC_seq <- seq(min(df$SC, na.rm = TRUE), max(df$SC, na.rm = TRUE), length.out = 100)

# -------------------------------
# 10. Predicted values by Season × OM
# -------------------------------
emm_d13C <- emmeans(mod_d13C, ~ SC | Season + OM, at = list(SC = SC_seq))
emm_d15N <- emmeans(mod_d15N, ~ SC | Season + OM, at = list(SC = SC_seq))

pred_fixed <- bind_rows(
  as.data.frame(emm_d13C) %>% mutate(Isotope = "d13C"),
  as.data.frame(emm_d15N) %>% mutate(Isotope = "d15N")
)

# Keep only significant Season × OM slopes for plotting
pred_fixed_sig <- pred_fixed %>%
  semi_join(bind_rows(
    d13C_SC_slopes %>% filter(Significance != "ns"),
    d15N_SC_slopes %>% filter(Significance != "ns")
  ), by = c("Season", "OM", "Isotope"))

# -------------------------------
# 11. Pivot raw data to long format
# -------------------------------
df_long <- df %>%
  pivot_longer(cols = c(d13C, d15N), names_to = "Isotope", values_to = "Value")

# -------------------------------
# 12. Plot aesthetics
# -------------------------------
OM_shapes <- c("CBOM" = 16, "FBOM" = 17, "Periphyton" = 15)
OM_colors <- c("CBOM" = "#E69F00", "FBOM" = "#0072B2", "Periphyton" = "#228B22")
month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2023)"
)

# -------------------------------
# 13. δ13C Plot
# -------------------------------
p_d13C <- ggplot(df_long %>% filter(Isotope == "d13C"),
                 aes(x = SC, y = Value, color = OM, shape = OM)) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  geom_ribbon(data = pred_fixed_sig %>% filter(Isotope == "d13C"),
              aes(x = SC, ymin = lower.CL, ymax = upper.CL,
                  fill = OM, group = interaction(Season, OM)),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_fixed_sig %>% filter(Isotope == "d13C"),
            aes(x = SC, y = emmean, color = OM, group = interaction(Season, OM)),
            size = 1.2, inherit.aes = FALSE) +
  scale_color_manual(values = OM_colors, guide = "legend") +
  scale_fill_manual(values = OM_colors, guide = "none") +
  scale_shape_manual(values = OM_shapes) +
  facet_wrap(~ Season, scales = "free_y", labeller = labeller(Season = month_labels)) +
  theme_bw(base_size = 16) +
  labs(x = "Specific Conductivity (µS/m)",
       y = expression(delta^{13} * "C"),
) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

print(p_d13C)

# -------------------------------
# 14. δ15N Plot
# -------------------------------
p_d15N <- ggplot(df_long %>% filter(Isotope == "d15N"),
                 aes(x = SC, y = Value, color = OM, shape = OM)) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  geom_ribbon(data = pred_fixed_sig %>% filter(Isotope == "d15N"),
              aes(x = SC, ymin = lower.CL, ymax = upper.CL,
                  fill = OM, group = interaction(Season, OM)),
              alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = pred_fixed_sig %>% filter(Isotope == "d15N"),
            aes(x = SC, y = emmean, color = OM, group = interaction(Season, OM)),
            size = 1.2, inherit.aes = FALSE) +
  scale_color_manual(values = OM_colors, guide = "legend") +
  scale_fill_manual(values = OM_colors, guide = "none") +
  scale_shape_manual(values = OM_shapes) +
  facet_wrap(~ Season, scales = "free_y", labeller = labeller(Season = month_labels)) +
  theme_bw(base_size = 16) +
  labs(x = "Specific Conductivity (µS/m)",
       y = expression(delta^{15} * "N"),
) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
print(p_d15N)

#########################################
# Organic Matter C:N Ratio GLMM
##########################################

# -------------------------------
# 1. Load libraries
# -------------------------------
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 2. Set working directory
# -------------------------------
setwd("C:/Users/lisat/OneDrive/Desktop/Salty_Carbon_Data_and_Scripts")

# -------------------------------
# 3. Read data
# -------------------------------
df <- read.csv("SALTY_C_CN_YEAR_1_ALL.csv", check.names = FALSE)

# -------------------------------
# 4. Data preparation
# -------------------------------
df <- df %>%
  mutate(
    Site = factor(Site),
    month = factor(Season, levels = c("February", "May", "August", "October")),
    OM = factor(OM),
    SC = as.numeric(SC),
    C_N = as.numeric(C_N)
  )

# -------------------------------
# 5. Fit Full GLMM
# -------------------------------
mod_CN <- lmer(
  C_N ~ SC * month * OM + (1 | Site),
  data = df
)

# -------------------------------
# 6. Extract SC slopes for each Month × OM
# -------------------------------
CN_SC_slopes <- emtrends(
  mod_CN,
  ~ month | OM,
  var = "SC"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    t_value = SC.trend / SE,
    p_value = 2 * pt(-abs(t_value), df),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    SC_slope = signif(SC.trend, 3),
    Std_Error = signif(SE, 3),
    t_value = signif(t_value, 3),
    p_value = signif(p_value, 3),
    lower.CL = signif(lower.CL, 3),
    upper.CL = signif(upper.CL, 3)
  ) %>%
  select(OM, month, SC_slope, Std_Error, df, t_value, p_value, lower.CL, upper.CL, Significance)

# -------------------------------
# 7. Print SC slopes table (all)
# -------------------------------
print("=== C:N ~ SC Slopes by Month × OM (All Slopes) ===")
print(CN_SC_slopes)
write.csv(CN_SC_slopes, "CN_SC_slopes_by_month_OM.csv", row.names = FALSE)

# -------------------------------
# 8. Salinity prediction sequence
# -------------------------------
SC_seq <- seq(
  min(df$SC, na.rm = TRUE),
  max(df$SC, na.rm = TRUE),
  length.out = 100
)

# -------------------------------
# 9. Predicted values for Month × OM
# -------------------------------
emm_CN <- emmeans(
  mod_CN,
  ~ SC | month + OM,
  at = list(SC = SC_seq)
)
pred_CN <- as.data.frame(emm_CN)

# -------------------------------
# 10. Keep only significant slopes for plotting
# -------------------------------
sig_CN_slopes <- CN_SC_slopes %>%
  filter(Significance != "ns") %>%
  select(OM, month)

pred_CN_sig <- pred_CN %>%
  semi_join(sig_CN_slopes, by = c("OM", "month"))

# -------------------------------
# 11. Plot settings
# -------------------------------
OM_shapes <- c("CBOM" = 16, "FBOM" = 17, "Periphyton" = 15)
OM_colors <- c("CBOM" = "#E69F00",
               "FBOM" = "#0072B2",
               "Periphyton" = "#228B22")

month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2024)"
)

# -------------------------------
# 12. Faceted Plot: OM × Month (only significant lines)
# -------------------------------
p_CN <- ggplot(df, aes(x = SC, y = C_N, color = OM, shape = OM)) +
  
  # Raw data points
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  
  # Only significant ribbons
  geom_ribbon(
    data = pred_CN_sig,
    aes(x = SC,
        ymin = lower.CL,
        ymax = upper.CL,
        fill = OM,
        group = interaction(OM, month)),
    alpha = 0.3,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  # Only significant lines
  geom_line(
    data = pred_CN_sig,
    aes(x = SC,
        y = emmean,
        color = OM,
        group = interaction(OM, month)),
    size = 1.2,
    inherit.aes = FALSE
  ) +
  
  # Manual scales
  scale_color_manual(values = OM_colors) +
  scale_fill_manual(values = OM_colors) +
  scale_shape_manual(values = OM_shapes) +
  
  # Facet grid: OM rows × Month columns
  facet_grid(
    OM ~ month,
    labeller = labeller(month = month_labels),
    scales = "free_y"
  ) +
  
  # Theme + labels
  theme_bw(base_size = 16) +
  labs(
    x = "Specific Conductivity (µS/m)",
    y = "C:N"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

# -------------------------------
# Print the faceted plot
# -------------------------------
print(p_CN)

##########################################
# Supplemental Organic Matter δ13C vs δ15N Biplot
##########################################

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 1. Load data
# -------------------------------
setwd("C:/Users/lisat/OneDrive/Desktop/Salty_Carbon_Data_and_Scripts")
df <- read.csv("SALTY_C_SI_YEAR_1_ALL.csv", check.names = FALSE)

# -------------------------------
# 2. Prepare data
# -------------------------------
df <- df %>%
  mutate(
    Season = factor(Season, levels = c("February", "May", "August", "October")),
    Site   = factor(Site),
    OM     = factor(OM),
    SC     = as.numeric(SC),
    d13C   = as.numeric(d13C),
    d15N   = as.numeric(d15N)
  )

# -------------------------------
# 3. Fit full interaction models
# -------------------------------
mod_d13C <- lmer(d13C ~ SC * Season * OM + (1 | Site), data = df)
mod_d15N <- lmer(d15N ~ SC * Season * OM + (1 | Site), data = df)

# -------------------------------
# 4. Extract SC slopes and significance
# -------------------------------
extract_slopes <- function(model, isotope) {
  emtrends(model, ~ Season | OM, var = "SC") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    mutate(
      t_value = SC.trend / SE,
      p_value = 2 * pt(-abs(t_value), df),
      Significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      Isotope = isotope,
      SC_slope = SC.trend
    ) %>%
    select(Isotope, OM, Season, SC_slope, SE, df, t_value, p_value, Significance)
}

slopes_d13C <- extract_slopes(mod_d13C, "d13C")
slopes_d15N <- extract_slopes(mod_d15N, "d15N")

# Combine slope tables
all_slopes <- bind_rows(slopes_d13C, slopes_d15N)
cat("\n=== SC Slopes for δ13C and δ15N ===\n")
print(all_slopes)

# -------------------------------
# 5. Salinity sequence for predictions
# -------------------------------
SC_seq <- seq(min(df$SC, na.rm = TRUE), max(df$SC, na.rm = TRUE), length.out = 100)

# -------------------------------
# 6. Predicted values (only significant slopes)
# -------------------------------
emm_d13C <- emmeans(mod_d13C, ~ SC | Season + OM, at = list(SC = SC_seq))
emm_d15N <- emmeans(mod_d15N, ~ SC | Season + OM, at = list(SC = SC_seq))

pred <- bind_rows(
  as.data.frame(emm_d13C) %>% mutate(Isotope = "d13C"),
  as.data.frame(emm_d15N) %>% mutate(Isotope = "d15N")
)

sig_slopes <- bind_rows(
  slopes_d13C %>% filter(Significance != "ns"),
  slopes_d15N %>% filter(Significance != "ns")
)

pred_sig <- pred %>% semi_join(sig_slopes, by = c("Season", "OM", "Isotope"))

# -------------------------------
# 7. Pivot predictions to wide for biplot
# -------------------------------
pred_wide <- pred_sig %>%
  pivot_wider(names_from = Isotope, values_from = emmean) 

# -------------------------------
# 8. Plot aesthetics
# -------------------------------
OM_shapes <- c("CBOM" = 16, "FBOM" = 17, "Periphyton" = 15)
OM_colors <- c("CBOM" = "#E69F00", "FBOM" = "#0072B2", "Periphyton" = "#228B22")
month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2024)"
)

# -------------------------------
# 9. Biplot δ13C vs δ15N with significant slopes
# -------------------------------
ggplot(df, aes(x = d13C, y = d15N, color = OM, shape = OM)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(data = pred_wide, aes(x = d13C, y = d15N, group = interaction(OM, Season)),
            size = 1.2, inherit.aes = FALSE) +
  geom_vline(xintercept= -28, color = "darkred", linetype= "dashed", size=1) + #can remove this if needed
  facet_grid(OM ~ Season, labeller = labeller(Season = month_labels)) +
  scale_color_manual(values = OM_colors) +
  scale_shape_manual(values = OM_shapes) +
  theme_bw(base_size = 16) +
  labs(
    x = expression(delta^{13}*"C"),
    y = expression(delta^{15}*"N"),
    color = "OM",
    shape = "OM"
  ) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

#########################################
# Supplemental Organic Matter Carbon (C) and Nitrogen (N) GLMM
#########################################

# -------------------------------
# 1. Load libraries
# -------------------------------
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 2. Set working directory
# -------------------------------
setwd("C:/Users/lisat/OneDrive/Desktop/Salty_Carbon_Data_and_Scripts")

# -------------------------------
# 3. Read data
# -------------------------------
df <- read.csv("SALTY_C_CN_YEAR_1_ALL.csv", check.names = FALSE)

# -------------------------------
# 4. Data preparation
# -------------------------------
df <- df %>%
  mutate(
    Site  = factor(Site),
    month = factor(Season, levels = c("February", "May", "August", "October")),
    OM    = factor(OM),
    SC    = as.numeric(SC),
    C     = as.numeric(C),
    N     = as.numeric(N)
  )

# -------------------------------
# 5. Salinity prediction sequence
# -------------------------------
SC_seq <- seq(
  min(df$SC, na.rm = TRUE),
  max(df$SC, na.rm = TRUE),
  length.out = 100
)

# -------------------------------
# 6. Plot aesthetics
# -------------------------------
OM_shapes <- c("CBOM" = 16, "FBOM" = 17, "Periphyton" = 15)
OM_colors <- c(
  "CBOM" = "#E69F00",
  "FBOM" = "#0072B2",
  "Periphyton" = "#228B22"
)

month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2024)"
)

# ============================================================
# ======================= CARBON (C) =========================
# ============================================================

# -------------------------------
# 7. Fit GLMM
# -------------------------------
mod_C <- lmer(
  C ~ SC * month * OM + (1 | Site),
  data = df
)

# -------------------------------
# 8. Extract SC slopes (table)
# -------------------------------
C_SC_slopes <- emtrends(
  mod_C,
  ~ month | OM,
  var = "SC"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    t_value = SC.trend / SE,
    p_value = 2 * pt(-abs(t_value), df),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(OM, month, SC.trend, SE, df,
         t_value, p_value, lower.CL, upper.CL, Significance)

print("=== Carbon (C) ~ SC Slopes by Month × OM ===")
print(C_SC_slopes)
write.csv(C_SC_slopes, "C_SC_slopes_by_month_OM.csv", row.names = FALSE)

# -------------------------------
# 9. Predicted values
# -------------------------------
emm_C <- emmeans(
  mod_C,
  ~ SC | month + OM,
  at = list(SC = SC_seq)
)

pred_C <- as.data.frame(emm_C)

# -------------------------------
# 10. Keep only significant slopes
# -------------------------------
sig_C <- C_SC_slopes %>%
  filter(Significance != "ns") %>%
  select(OM, month)

pred_C_sig <- pred_C %>%
  semi_join(sig_C, by = c("OM", "month"))

# -------------------------------
# 11. Plot Carbon
# -------------------------------
p_C <- ggplot(df, aes(x = SC, y = C, color = OM, shape = OM)) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  geom_ribbon(
    data = pred_C_sig,
    aes(x = SC, ymin = lower.CL, ymax = upper.CL,
        fill = OM, group = interaction(OM, month)),
    alpha = 0.3,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_line(
    data = pred_C_sig,
    aes(x = SC, y = emmean,
        color = OM, group = interaction(OM, month)),
    size = 1.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = OM_colors) +
  scale_fill_manual(values = OM_colors) +
  scale_shape_manual(values = OM_shapes) +
  facet_grid(
    OM ~ month,
    labeller = labeller(month = month_labels),
    scales = "free_y"
  ) +
  theme_bw(base_size = 16) +
  labs(
    x = "Specific Conductivity (µS/m)",
    y = "% Carbon"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

print(p_C)

# ============================================================
# ====================== NITROGEN (N) ========================
# ============================================================

# -------------------------------
# 12. Fit GLMM
# -------------------------------
mod_N <- lmer(
  N ~ SC * month * OM + (1 | Site),
  data = df
)

# -------------------------------
# 13. Extract SC slopes (table)
# -------------------------------
N_SC_slopes <- emtrends(
  mod_N,
  ~ month | OM,
  var = "SC"
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    t_value = SC.trend / SE,
    p_value = 2 * pt(-abs(t_value), df),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  select(OM, month, SC.trend, SE, df,
         t_value, p_value, lower.CL, upper.CL, Significance)

print("=== Nitrogen (N) ~ SC Slopes by Month × OM ===")
print(N_SC_slopes)
write.csv(N_SC_slopes, "N_SC_slopes_by_month_OM.csv", row.names = FALSE)

# -------------------------------
# 14. Predicted values
# -------------------------------
emm_N <- emmeans(
  mod_N,
  ~ SC | month + OM,
  at = list(SC = SC_seq)
)

pred_N <- as.data.frame(emm_N)

# -------------------------------
# 15. Keep only significant slopes
# -------------------------------
sig_N <- N_SC_slopes %>%
  filter(Significance != "ns") %>%
  select(OM, month)

pred_N_sig <- pred_N %>%
  semi_join(sig_N, by = c("OM", "month"))

# -------------------------------
# 16. Plot Nitrogen
# -------------------------------
p_N <- ggplot(df, aes(x = SC, y = N, color = OM, shape = OM)) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  geom_ribbon(
    data = pred_N_sig,
    aes(x = SC, ymin = lower.CL, ymax = upper.CL,
        fill = OM, group = interaction(OM, month)),
    alpha = 0.3,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_line(
    data = pred_N_sig,
    aes(x = SC, y = emmean,
        color = OM, group = interaction(OM, month)),
    size = 1.2,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = OM_colors) +
  scale_fill_manual(values = OM_colors) +
  scale_shape_manual(values = OM_shapes) +
  facet_grid(
    OM ~ month,
    labeller = labeller(month = month_labels),
    scales = "free_y"
  ) +
  theme_bw(base_size = 16) +
  labs(
    x = "Specific Conductivity (µS/m)",
    y = "% Nitrogen"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

print(p_N)