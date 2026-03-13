# OM Script
# Periphyton is "algae" in-script

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("")

df <- read.csv("SALTY_C_AFDM_YEAR_1_ALL.csv", check.names=FALSE) %>%
  mutate(
    Site  = factor(Site),
    OM    = factor(OM),  # e.g., CBOM, FBOM, Periphyton
    month = factor(Season, levels=c("February","May","August","October")),
    SC    = as.numeric(SC),
    AFDM  = as.numeric(AFDM)
  )

OM_colors <- c("CBOM"="#E69F00", "FBOM"="#0072B2", "Periphyton"="#228B22")
OM_shapes <- c("CBOM"=16, "FBOM"=17, "Periphyton"=15)

month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2023)"
)

SC_seq <- seq(min(df$SC, na.rm=TRUE), max(df$SC, na.rm=TRUE), length.out=100)

mod_AFDM <- lmer(AFDM ~ SC * month * OM + (1 | Site), data=df)

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

print("=== AFDM ~ SC Slopes by OM and Month ===")
print(AFDM_SC_slopes_all)

write.csv(AFDM_SC_slopes_all, "AFDM_SC_slopes_by_OM_and_Month.csv", row.names=FALSE)

AFDM_SC_slopes_sig <- AFDM_SC_slopes_all %>% filter(Significance != "ns")

pred_sig <- emmeans(mod_AFDM, ~ SC | month * OM, at=list(SC=SC_seq), cov.reduce=FALSE) %>%
  as.data.frame() %>%
  semi_join(AFDM_SC_slopes_sig, by=c("month","OM"))

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

print(p_AFDM_sig)

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("")

df <- read.csv("SALTY_C_CHLA_ROCKS.csv", check.names=FALSE) %>%
  mutate(
    Site  = factor(Site),
    month = factor(Season, levels=c("February","May","August","October")),
    SC    = as.numeric(SC),
    CHLA  = as.numeric(CHLA)
  )

SC_seq <- seq(min(df$SC, na.rm=TRUE), max(df$SC, na.rm=TRUE), length.out=100)

mod_CHLA <- lmer(CHLA ~ SC * month + (1 | Site), data=df)

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

print("=== CHLA ~ SC Slopes by Month ===")
print(CHLA_SC_slopes_all)

write.csv(CHLA_SC_slopes_all, "CHLA_SC_slopes_by_month.csv", row.names=FALSE)

CHLA_SC_slopes_sig <- CHLA_SC_slopes_all %>% filter(Significance != "ns") %>% select(month)

pred_sig <- emmeans(mod_CHLA, ~ SC | month, at=list(SC=SC_seq), cov.reduce=FALSE) %>%
  as.data.frame() %>%
  semi_join(CHLA_SC_slopes_sig, by="month")

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

print(p_CHLA_sig)

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("")

df <- read.csv("SALTY_C_SI_YEAR_1_ALL.csv", check.names = FALSE)

df <- df %>%
  mutate(
    Season = factor(Season, levels = c("February", "May", "August", "October")),
    Site   = factor(Site),
    OM     = factor(OM),
    SC     = as.numeric(SC),
    d13C   = as.numeric(d13C),
    d15N   = as.numeric(d15N)
  )

mod_d13C <- lmer(d13C ~ SC * Season * OM + (1 | Site), data = df)
mod_d15N <- lmer(d15N ~ SC * Season * OM + (1 | Site), data = df)

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

d13C_SC_slopes <- extract_slopes_all(mod_d13C, "d13C")
d15N_SC_slopes <- extract_slopes_all(mod_d15N, "d15N")

cat("\n=== δ13C SC Slopes (All) ===\n")
print(d13C_SC_slopes)

cat("\n=== δ15N SC Slopes (All) ===\n")
print(d15N_SC_slopes)

SC_seq <- seq(min(df$SC, na.rm = TRUE), max(df$SC, na.rm = TRUE), length.out = 100)

emm_d13C <- emmeans(mod_d13C, ~ SC | Season + OM, at = list(SC = SC_seq))
emm_d15N <- emmeans(mod_d15N, ~ SC | Season + OM, at = list(SC = SC_seq))

pred_fixed <- bind_rows(
  as.data.frame(emm_d13C) %>% mutate(Isotope = "d13C"),
  as.data.frame(emm_d15N) %>% mutate(Isotope = "d15N")
)

pred_fixed_sig <- pred_fixed %>%
  semi_join(bind_rows(
    d13C_SC_slopes %>% filter(Significance != "ns"),
    d15N_SC_slopes %>% filter(Significance != "ns")
  ), by = c("Season", "OM", "Isotope"))

df_long <- df %>%
  pivot_longer(cols = c(d13C, d15N), names_to = "Isotope", values_to = "Value")

OM_shapes <- c("CBOM" = 16, "FBOM" = 17, "Periphyton" = 15)
OM_colors <- c("CBOM" = "#E69F00", "FBOM" = "#0072B2", "Periphyton" = "#228B22")
month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2023)"
)

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

library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("")

df <- read.csv("SALTY_C_CN_YEAR_1_ALL.csv", check.names = FALSE)

df <- df %>%
  mutate(
    Site = factor(Site),
    month = factor(Season, levels = c("February", "May", "August", "October")),
    OM = factor(OM),
    SC = as.numeric(SC),
    C_N = as.numeric(C_N)
  )

mod_CN <- lmer(
  C_N ~ SC * month * OM + (1 | Site),
  data = df
)

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

print("=== C:N ~ SC Slopes by Month × OM (All Slopes) ===")
print(CN_SC_slopes)
write.csv(CN_SC_slopes, "CN_SC_slopes_by_month_OM.csv", row.names = FALSE)

SC_seq <- seq(
  min(df$SC, na.rm = TRUE),
  max(df$SC, na.rm = TRUE),
  length.out = 100
)

emm_CN <- emmeans(
  mod_CN,
  ~ SC | month + OM,
  at = list(SC = SC_seq)
)
pred_CN <- as.data.frame(emm_CN)

sig_CN_slopes <- CN_SC_slopes %>%
  filter(Significance != "ns") %>%
  select(OM, month)

pred_CN_sig <- pred_CN %>%
  semi_join(sig_CN_slopes, by = c("OM", "month"))

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

p_CN <- ggplot(df, aes(x = SC, y = C_N, color = OM, shape = OM)) +
  
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  
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
  
  geom_line(
    data = pred_CN_sig,
    aes(x = SC,
        y = emmean,
        color = OM,
        group = interaction(OM, month)),
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
    y = "C:N"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

print(p_CN)

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

# -------------------------------
# 1. Load data
# -------------------------------
setwd("")
df <- read.csv("SALTY_C_SI_YEAR_1_ALL.csv", check.names = FALSE)

df <- df %>%
  mutate(
    Season = factor(Season, levels = c("February", "May", "August", "October")),
    Site   = factor(Site),
    OM     = factor(OM),
    SC     = as.numeric(SC),
    d13C   = as.numeric(d13C),
    d15N   = as.numeric(d15N)
  )

mod_d13C <- lmer(d13C ~ SC * Season * OM + (1 | Site), data = df)
mod_d15N <- lmer(d15N ~ SC * Season * OM + (1 | Site), data = df)

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

all_slopes <- bind_rows(slopes_d13C, slopes_d15N)
cat("\n=== SC Slopes for δ13C and δ15N ===\n")
print(all_slopes)

SC_seq <- seq(min(df$SC, na.rm = TRUE), max(df$SC, na.rm = TRUE), length.out = 100)

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

pred_wide <- pred_sig %>%
  pivot_wider(names_from = Isotope, values_from = emmean) 

OM_shapes <- c("CBOM" = 16, "FBOM" = 17, "Periphyton" = 15)
OM_colors <- c("CBOM" = "#E69F00", "FBOM" = "#0072B2", "Periphyton" = "#228B22")
month_labels <- c(
  "February" = "February (Winter 2024)",
  "May"      = "May (Spring 2024)",
  "August"   = "August (Summer 2024)",
  "October"  = "October (Autumn 2024)"
)

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
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("")

df <- read.csv("SALTY_C_CN_YEAR_1_ALL.csv", check.names = FALSE)

df <- df %>%
  mutate(
    Site  = factor(Site),
    month = factor(Season, levels = c("February", "May", "August", "October")),
    OM    = factor(OM),
    SC    = as.numeric(SC),
    C     = as.numeric(C),
    N     = as.numeric(N)
  )

SC_seq <- seq(
  min(df$SC, na.rm = TRUE),
  max(df$SC, na.rm = TRUE),
  length.out = 100
)

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

mod_C <- lmer(
  C ~ SC * month * OM + (1 | Site),
  data = df
)

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

emm_C <- emmeans(
  mod_C,
  ~ SC | month + OM,
  at = list(SC = SC_seq)
)

pred_C <- as.data.frame(emm_C)

sig_C <- C_SC_slopes %>%
  filter(Significance != "ns") %>%
  select(OM, month)

pred_C_sig <- pred_C %>%
  semi_join(sig_C, by = c("OM", "month"))

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

mod_N <- lmer(
  N ~ SC * month * OM + (1 | Site),
  data = df
)

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

emm_N <- emmeans(
  mod_N,
  ~ SC | month + OM,
  at = list(SC = SC_seq)
)

pred_N <- as.data.frame(emm_N)

sig_N <- N_SC_slopes %>%
  filter(Significance != "ns") %>%
  select(OM, month)

pred_N_sig <- pred_N %>%
  semi_join(sig_N, by = c("OM", "month"))

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