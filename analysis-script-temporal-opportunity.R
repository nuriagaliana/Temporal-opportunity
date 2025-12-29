# ============================================================================
# COMPLETE TEMPORAL OPPORTUNITY ANALYSIS - SEASON 2007-2008
# ============================================================================

# LOAD LIBRARIES -------------------------------------------------------------
library(bipartite)
library(igraph)
library(AICcmodavg)
library(viridis)
library(ggplot2)
library(readxl)
library(lubridate)
library(dplyr)
library(tidyr)
library(vegan)
library(reshape2)
library(betapart)
library(gridExtra)
library(grid)
library(knitr)
library(extrafont)

# ============================================================================
# DATA IMPORT AND PREPARATION
# ============================================================================

# Metadata of the Sierras
sierras <- read.csv('./sierras2.csv', header = T, sep = ';')

# Treating the data for 2008
data_2008 <- read_excel("./Network 2008_temporal_opportunity.xlsx")[1:4]
data_2008$Date <- as.Date(data_2008$Date)
data_2008$Year <- year(data_2008$Date)
data_2008$Day <- yday(data_2008$Date)
data_2008$season <-'2007-2008'
data_2008$interaction <- paste0(data_2008$Plant_species,data_2008$Insect_species)

# Use only 2008 data
data_all <- na.omit(data_2008[,c('Sierra','Date','Plant_species','Insect_species','Year','Day','season','interaction')])

# Species richness and number of links per day
data_all_richness  <- data_all %>%
  dplyr::group_by(Sierra, Date, Year, Day) %>%
  dplyr::mutate(plant_richness = length(unique(Plant_species)),
                insect_richness = length(unique(Insect_species)),
                number_of_links = length(unique(interaction)))

# Remove 'No visit' and 'SV' from interactions
data_all_links <- data_all_richness[data_all_richness$Insect_species!='No visit' & 
                                      data_all_richness$Insect_species!='SV',]

# Overall species richness and links per sierra
total_links <- data_all_links %>%
  dplyr::group_by(Sierra, season) %>%
  dplyr::mutate(
    total_number_of_links = length(unique(interaction)),
    total_plant_richness = length(unique(Plant_species)),
    total_insect_richness = length(unique(Insect_species)),
    pollinator_plant_ratio = length(unique(Insect_species))/length(unique(Plant_species)),
    links.sp_insect = length(unique(interaction))/length(unique(Insect_species)),
    links.sp_plant = length(unique(interaction))/length(unique(Plant_species))
  )

# Meta-community metrics across all sierras
overall_data <- total_links %>%
  dplyr::group_by(season) %>%
  dplyr::mutate(
    meta_number_of_links = length(unique(interaction)),
    meta_plant_richness = length(unique(Plant_species)),
    meta_insect_richness = length(unique(Insect_species))
  )

# Merge with sierra metadata
overall_data_metadata <- merge(overall_data, sierras, by='Sierra')

# ============================================================================
# FIGURE 1 ANALYSES: Species Richness, Total Links, Links/Species vs Area
# ============================================================================

cat("\n=== FIGURE 1 ANALYSES ===\n")

# Prepare data
fig1_data <- overall_data_metadata %>%
  filter(season == '2007-2008') %>%
  select(Sierra, area, total_plant_richness, total_insect_richness, 
         total_number_of_links, links.sp_plant, links.sp_insect) %>%
  distinct() %>%
  mutate(
    is_difuntito = Sierra == 'Difuntito',
    total_species = total_plant_richness + total_insect_richness
  )

fig1_analysis <- fig1_data %>% filter(Sierra != 'Difuntito')

# Models for statistics (log-log for power-law)
lm_plant_rich <- lm(log10(total_plant_richness) ~ log10(area), data = fig1_analysis)
lm_poll_rich <- lm(log10(total_insect_richness) ~ log10(area), data = fig1_analysis)
lm_total_rich <- lm(log10(total_species) ~ log10(area), data = fig1_analysis)
lm_total_links <- lm(log10(total_number_of_links) ~ log10(area), data = fig1_analysis)
lm_links_plant <- lm(log10(links.sp_plant) ~ log10(area), data = fig1_analysis)
lm_links_poll <- lm(log10(links.sp_insect) ~ log10(area), data = fig1_analysis)

cat("\nPanel A: Species Richness vs Area\n")
cat("Plants: Slope =", round(coef(lm_plant_rich)[2], 3), 
    ", R² =", round(summary(lm_plant_rich)$r.squared, 3),
    ", p =", round(summary(lm_plant_rich)$coefficients[2,4], 4), "\n")
cat("Pollinators: Slope =", round(coef(lm_poll_rich)[2], 3),
    ", R² =", round(summary(lm_poll_rich)$r.squared, 3),
    ", p =", round(summary(lm_poll_rich)$coefficients[2,4], 4), "\n")
cat("Total: Slope =", round(coef(lm_total_rich)[2], 3),
    ", R² =", round(summary(lm_total_rich)$r.squared, 3),
    ", p =", round(summary(lm_total_rich)$coefficients[2,4], 4), "\n")

cat("\nPanel B: Total Links vs Area\n")
cat("Slope =", round(coef(lm_total_links)[2], 3),
    ", R² =", round(summary(lm_total_links)$r.squared, 3),
    ", p =", round(summary(lm_total_links)$coefficients[2,4], 4), "\n")

cat("\nPanel C: Links/Species vs Area\n")
cat("Plants: Slope =", round(coef(lm_links_plant)[2], 3),
    ", R² =", round(summary(lm_links_plant)$r.squared, 3),
    ", p =", round(summary(lm_links_plant)$coefficients[2,4], 4), "\n")
cat("Pollinators: Slope =", round(coef(lm_links_poll)[2], 3),
    ", R² =", round(summary(lm_links_poll)$r.squared, 3),
    ", p =", round(summary(lm_links_poll)$coefficients[2,4], 4), "\n")


# ============================================================================
# FIGURE 1: Area-Dependent Patterns in Network Structure
# ============================================================================

cat("\n=== GENERATING FIGURE 1 ===\n")

# Color scheme
color_pollinator <- "#E69F00"  
color_plant <- "#009E73"
color_total <- "grey30"

# ============================================================================
# Prepare data for Figure 1
# ============================================================================

fig1_data <- overall_data_metadata %>%
  filter(season == '2007-2008') %>%
  select(Sierra, area, total_plant_richness, total_insect_richness, 
         total_number_of_links, links.sp_plant, links.sp_insect) %>%
  distinct() %>%
  mutate(
    is_difuntito = Sierra == 'Difuntito',
    total_species = total_plant_richness + total_insect_richness
  )

# Data for analysis (excluding Difuntito)
fig1_analysis <- fig1_data %>% filter(Sierra != 'Difuntito')

# ============================================================================
# Statistical Models
# ============================================================================

# Models for statistics (log-log for power-law)
lm_plant_rich <- lm(log10(total_plant_richness) ~ log10(area), data = fig1_analysis)
lm_poll_rich <- lm(log10(total_insect_richness) ~ log10(area), data = fig1_analysis)
lm_total_rich <- lm(log10(total_species) ~ log10(area), data = fig1_analysis)
lm_total_links <- lm(log10(total_number_of_links) ~ log10(area), data = fig1_analysis)
lm_links_plant <- lm(log10(links.sp_plant) ~ log10(area), data = fig1_analysis)
lm_links_poll <- lm(log10(links.sp_insect) ~ log10(area), data = fig1_analysis)

cat("\nPanel A: Species Richness vs Area\n")
cat("Plants: Slope =", round(coef(lm_plant_rich)[2], 3), 
    ", R² =", round(summary(lm_plant_rich)$r.squared, 3),
    ", p =", round(summary(lm_plant_rich)$coefficients[2,4], 4), "\n")
cat("Pollinators: Slope =", round(coef(lm_poll_rich)[2], 3),
    ", R² =", round(summary(lm_poll_rich)$r.squared, 3),
    ", p =", round(summary(lm_poll_rich)$coefficients[2,4], 4), "\n")
cat("Total: Slope =", round(coef(lm_total_rich)[2], 3),
    ", R² =", round(summary(lm_total_rich)$r.squared, 3),
    ", p =", round(summary(lm_total_rich)$coefficients[2,4], 4), "\n")

cat("\nPanel B: Total Links vs Area\n")
cat("Slope =", round(coef(lm_total_links)[2], 3),
    ", R² =", round(summary(lm_total_links)$r.squared, 3),
    ", p =", round(summary(lm_total_links)$coefficients[2,4], 4), "\n")

cat("\nPanel C: Links/Species vs Area\n")
cat("Plants: Slope =", round(coef(lm_links_plant)[2], 3),
    ", R² =", round(summary(lm_links_plant)$r.squared, 3),
    ", p =", round(summary(lm_links_plant)$coefficients[2,4], 4), "\n")
cat("Pollinators: Slope =", round(coef(lm_links_poll)[2], 3),
    ", R² =", round(summary(lm_links_poll)$r.squared, 3),
    ", p =", round(summary(lm_links_poll)$coefficients[2,4], 4), "\n")

# ============================================================================
# Panel A: Species richness vs Area
# ============================================================================

p1a <- ggplot(fig1_data, aes(x = area)) +
  geom_point(data = fig1_data[!fig1_data$is_difuntito,], 
             aes(y = total_species), color = color_total, size = 3, alpha = 0.7) +
  geom_point(data = fig1_data[fig1_data$is_difuntito,], 
             aes(y = total_species), color = 'red', shape = 17, size = 4) +
  geom_smooth(data = fig1_analysis, aes(y = total_species),
              method = 'lm', color = color_total, se = TRUE, alpha = 0.2, linewidth = 1.2) +
  geom_point(data = fig1_data[!fig1_data$is_difuntito,], 
             aes(y = total_insect_richness), color = color_pollinator, size = 3, alpha = 0.7) +
  geom_point(data = fig1_data[fig1_data$is_difuntito,], 
             aes(y = total_insect_richness), color = 'red', shape = 17, size = 4) +
  geom_smooth(data = fig1_analysis, aes(y = total_insect_richness),
              method = 'lm', color = color_pollinator, se = TRUE, alpha = 0.2, linewidth = 1.2) +
  geom_point(data = fig1_data[!fig1_data$is_difuntito,], 
             aes(y = total_plant_richness), color = color_plant, size = 3, alpha = 0.7) +
  geom_point(data = fig1_data[fig1_data$is_difuntito,], 
             aes(y = total_plant_richness), color = 'red', shape = 17, size = 4) +
  geom_smooth(data = fig1_analysis, aes(y = total_plant_richness),
              method = 'lm', color = color_plant, se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(x = "Area", y = "Species Richness", title = "A. Species-Area Relationship") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)  # Made all margins equal
  )

# ============================================================================
# Panel B: Total links vs Area
# ============================================================================

p1b <- ggplot(fig1_data, aes(x = area)) +
  geom_point(data = fig1_data[!fig1_data$is_difuntito,], 
             aes(y = total_number_of_links), color = 'grey30', size = 3, alpha = 0.7) +
  geom_point(data = fig1_data[fig1_data$is_difuntito,], 
             aes(y = total_number_of_links), color = 'red', shape = 17, size = 4) +
  geom_smooth(data = fig1_analysis, aes(y = total_number_of_links),
              method = 'lm', color = 'grey30', se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(x = "Area", y = "Total Number of Links", title = "B. Links-Area Relationship") +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Slope = %.2f\nR² = %.2f\np = %.3f",
                           coef(lm_total_links)[2],
                           summary(lm_total_links)$r.squared,
                           summary(lm_total_links)$coefficients[2,4]),
           hjust = 1.1, vjust = -0.3, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)  # Made all margins equal
  )

# ============================================================================
# Panel C: Links/Species vs Area
# ============================================================================

p1c <- ggplot(fig1_data, aes(x = area)) +
  geom_point(data = fig1_data[!fig1_data$is_difuntito,], 
             aes(y = links.sp_insect), color = color_pollinator, size = 3, alpha = 0.7) +
  geom_point(data = fig1_data[fig1_data$is_difuntito,], 
             aes(y = links.sp_insect), color = 'red', shape = 17, size = 4) +
  geom_smooth(data = fig1_analysis, aes(y = links.sp_insect),
              method = 'lm', color = color_pollinator, se = TRUE, alpha = 0.2, linewidth = 1.2) +
  geom_point(data = fig1_data[!fig1_data$is_difuntito,], 
             aes(y = links.sp_plant), color = color_plant, size = 3, alpha = 0.7) +
  geom_point(data = fig1_data[fig1_data$is_difuntito,], 
             aes(y = links.sp_plant), color = 'red', shape = 17, size = 4) +
  geom_smooth(data = fig1_analysis, aes(y = links.sp_plant),
              method = 'lm', color = color_plant, se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(x = "Area", y = "Links per Species", title = "C. Links/Species-Area Relationship") +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Pollinators: slope = %.2f, R² = %.2f, p = %.3f\nPlants: slope = %.2f, R² = %.2f, p = %.3f",
                           coef(lm_links_poll)[2], summary(lm_links_poll)$r.squared, summary(lm_links_poll)$coefficients[2,4],
                           coef(lm_links_plant)[2], summary(lm_links_plant)$r.squared, summary(lm_links_plant)$coefficients[2,4]),
           hjust = 1.1, vjust = -0.3, size = 3.2, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)  # Made all margins equal - THIS WAS THE KEY FIX
  )

# ============================================================================
# Combine panels and save
# ============================================================================

fig1_combined <- grid.arrange(
  p1a, p1b, p1c,
  ncol = 3,
  top = textGrob("", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

ggsave("figure1_area_patterns.pdf", fig1_combined, 
       width = 10.5, height = 3.5, units = "in", dpi = 300)
ggsave("figure1_area_patterns.png", fig1_combined, 
       width = 10.5, height = 3.5, units = "in", dpi = 300)


cat("\nFigure 1 saved: figure1_area_patterns.pdf\n")
cat("  Dimensions: 14 × 5 inches\n")
cat("  3 panels: Species Richness, Total Links, Links/Species vs Area\n")

# ============================================================================
# TABLE FOR FIGURE 1: Statistical Results of Area-Dependent Patterns
# ============================================================================

cat("\n=== GENERATING TABLE FOR FIGURE 1 STATISTICS ===\n")

# ============================================================================
# Run all statistical models (assuming data is prepared)
# ============================================================================

# Models for Panel A: Species Richness vs Area
lm_plant_rich <- lm(log10(total_plant_richness) ~ log10(area), data = fig1_analysis)
lm_poll_rich <- lm(log10(total_insect_richness) ~ log10(area), data = fig1_analysis)
lm_total_rich <- lm(log10(total_species) ~ log10(area), data = fig1_analysis)

# Models for Panel B: Total Links vs Area
lm_total_links <- lm(log10(total_number_of_links) ~ log10(area), data = fig1_analysis)

# Models for Panel C: Links/Species vs Area
lm_links_plant <- lm(log10(links.sp_plant) ~ log10(area), data = fig1_analysis)
lm_links_poll <- lm(log10(links.sp_insect) ~ log10(area), data = fig1_analysis)

# ============================================================================
# Extract statistics function
# ============================================================================

extract_stats <- function(model, response_name, predictor_name = "log10(area)") {
  summary_model <- summary(model)
  coefs <- summary_model$coefficients
  
  # Get slope statistics
  slope <- coefs[2, 1]
  slope_se <- coefs[2, 2]
  t_value <- coefs[2, 3]
  p_value <- coefs[2, 4]
  
  # Get intercept
  intercept <- coefs[1, 1]
  intercept_se <- coefs[1, 2]
  
  # Get model fit statistics
  r_squared <- summary_model$r.squared
  adj_r_squared <- summary_model$adj.r.squared
  f_stat <- summary_model$fstatistic[1]
  f_df1 <- summary_model$fstatistic[2]
  f_df2 <- summary_model$fstatistic[3]
  f_p_value <- pf(f_stat, f_df1, f_df2, lower.tail = FALSE)
  
  # Sample size
  n <- nobs(model)
  
  # Format p-value
  p_formatted <- ifelse(p_value < 0.001, "< 0.001",
                        ifelse(p_value < 0.01, "< 0.01",
                               ifelse(p_value < 0.05, "< 0.05",
                                      sprintf("%.3f", p_value))))
  
  f_p_formatted <- ifelse(f_p_value < 0.001, "< 0.001",
                          ifelse(f_p_value < 0.01, "< 0.01",
                                 ifelse(f_p_value < 0.05, "< 0.05",
                                        sprintf("%.3f", f_p_value))))
  
  return(data.frame(
    Response = response_name,
    n = n,
    Intercept = sprintf("%.3f ± %.3f", intercept, intercept_se),
    Slope = sprintf("%.3f ± %.3f", slope, slope_se),
    t_value = sprintf("%.2f", t_value),
    p_value = p_formatted,
    R2 = sprintf("%.3f", r_squared),
    Adj_R2 = sprintf("%.3f", adj_r_squared),
    F_statistic = sprintf("%.2f", f_stat),
    F_p_value = f_p_formatted,
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# Create comprehensive statistics table
# ============================================================================

stats_table <- rbind(
  extract_stats(lm_plant_rich, "Plant Richness"),
  extract_stats(lm_poll_rich, "Pollinator Richness"),
  extract_stats(lm_total_rich, "Total Species Richness"),
  extract_stats(lm_total_links, "Total Links"),
  extract_stats(lm_links_plant, "Links per Plant Species"),
  extract_stats(lm_links_poll, "Links per Pollinator Species")
)

# Add panel information
stats_table$Panel <- c("A", "A", "A", "B", "C", "C")

# Reorder columns
stats_table <- stats_table[, c("Panel", "Response", "n", "Intercept", "Slope", 
                               "t_value", "p_value", "R2", "Adj_R2", 
                               "F_statistic", "F_p_value")]

# ============================================================================
# Print and save table
# ============================================================================

cat("\n=== TABLE 1: Area-Dependence Statistics ===\n\n")
print(stats_table, row.names = FALSE)

# Save as CSV
write.csv(stats_table, "table1_figure1_statistics.csv", row.names = FALSE)
cat("\nTable saved: table1_figure1_statistics.csv\n")

# ============================================================================
# Create publication-ready formatted table
# ============================================================================

# More compact version for publication
pub_table <- data.frame(
  Panel = stats_table$Panel,
  `Response Variable` = stats_table$Response,
  n = stats_table$n,
  `Slope (z)` = stats_table$Slope,
  `R²` = stats_table$R2,
  `p-value` = stats_table$p_value,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n=== PUBLICATION TABLE (Compact) ===\n\n")
print(pub_table, row.names = FALSE)

write.csv(pub_table, "table1_figure1_compact.csv", row.names = FALSE)
cat("\nCompact table saved: table1_figure1_compact.csv\n")

# ============================================================================
# Create detailed statistics summary
# ============================================================================

cat("\n\n=== DETAILED STATISTICS SUMMARY ===\n")

summary_stats <- data.frame(
  Panel = c("A", "A", "A", "B", "C", "C"),
  Variable = c("Plant Richness", "Pollinator Richness", "Total Species", 
               "Total Links", "Links/Plant", "Links/Pollinator"),
  Slope = c(coef(lm_plant_rich)[2], coef(lm_poll_rich)[2], coef(lm_total_rich)[2],
            coef(lm_total_links)[2], coef(lm_links_plant)[2], coef(lm_links_poll)[2]),
  SE = c(summary(lm_plant_rich)$coefficients[2,2],
         summary(lm_poll_rich)$coefficients[2,2],
         summary(lm_total_rich)$coefficients[2,2],
         summary(lm_total_links)$coefficients[2,2],
         summary(lm_links_plant)$coefficients[2,2],
         summary(lm_links_poll)$coefficients[2,2]),
  R2 = c(summary(lm_plant_rich)$r.squared,
         summary(lm_poll_rich)$r.squared,
         summary(lm_total_rich)$r.squared,
         summary(lm_total_links)$r.squared,
         summary(lm_links_plant)$r.squared,
         summary(lm_links_poll)$r.squared),
  p_value = c(summary(lm_plant_rich)$coefficients[2,4],
              summary(lm_poll_rich)$coefficients[2,4],
              summary(lm_total_rich)$coefficients[2,4],
              summary(lm_total_links)$coefficients[2,4],
              summary(lm_links_plant)$coefficients[2,4],
              summary(lm_links_poll)$coefficients[2,4])
)

# Add confidence intervals
for(i in 1:nrow(summary_stats)) {
  models <- list(lm_plant_rich, lm_poll_rich, lm_total_rich, 
                 lm_total_links, lm_links_plant, lm_links_poll)
  ci <- confint(models[[i]], level = 0.95)[2,]
  summary_stats$CI_lower[i] <- ci[1]
  summary_stats$CI_upper[i] <- ci[2]
}

print(summary_stats)

write.csv(summary_stats, "table1_figure1_detailed.csv", row.names = FALSE)
cat("\nDetailed table saved: table1_figure1_detailed.csv\n")

# ============================================================================
# FIGURE 2 ANALYSES: CORRECTED Degree Distribution Analysis
# ============================================================================

cat("\n=== FIGURE 2 ANALYSES: DEGREE DISTRIBUTIONS ===\n")

# Calculate degree distributions
# For Plants
degree_dist_data_plant <- overall_data %>%
  dplyr::group_by(season, Sierra, Plant_species) %>%
  dplyr::mutate(
    number_of_interactions_sp = length(unique(Insect_species)),
    number_of_interactions_sp_norm = length(unique(Insect_species))/(total_number_of_links/(total_plant_richness+total_insect_richness))
  )

# For Insects
degree_dist_data_insect <- overall_data %>%
  dplyr::group_by(season, Sierra, Insect_species) %>%
  dplyr::mutate(
    number_of_interactions_sp = length(unique(Plant_species)),
    number_of_interactions_sp_norm = length(unique(Plant_species))/(total_number_of_links/(total_plant_richness+total_insect_richness))
  )

# CORRECTED: Build cumulative distributions properly
output_degree_plants <- NULL
for(s in unique(degree_dist_data_plant$Sierra)){
  df_temp <- degree_dist_data_plant[degree_dist_data_plant$Sierra==s,]
  for(seas in unique(df_temp$season)){
    df_seas <- df_temp[df_temp$season==seas,]
    df <- unique(df_seas[,c('Plant_species','number_of_interactions_sp')])
    
    # Count occurrences of each degree
    degree_counts <- table(df$number_of_interactions_sp)
    
    # Calculate probability mass function (PMF)
    p <- as.vector(degree_counts) / sum(degree_counts)
    
    # Calculate cumulative distribution (CCDF = complementary cumulative)
    y <- rev(cumsum(rev(p)))
    
    # Degree values
    x <- as.numeric(names(degree_counts))
    
    cur_out <- data.frame(Sierra=s, Season=seas, links=x, prob_x=y)
    
    if(is.null(output_degree_plants)){
      output_degree_plants <- cur_out
    } else {
      output_degree_plants <- rbind(output_degree_plants, cur_out)
    }
  }
}

output_degree_insects <- NULL
for(s in unique(degree_dist_data_insect$Sierra)){
  df_temp <- degree_dist_data_insect[degree_dist_data_insect$Sierra==s,]
  for(seas in unique(df_temp$season)){
    df_seas <- df_temp[df_temp$season==seas,]
    df <- unique(df_seas[,c('Insect_species','number_of_interactions_sp')])
    
    # Count occurrences of each degree
    degree_counts <- table(df$number_of_interactions_sp)
    
    # Calculate probability mass function (PMF)
    p <- as.vector(degree_counts) / sum(degree_counts)
    
    # Calculate cumulative distribution (CCDF = complementary cumulative)
    y <- rev(cumsum(rev(p)))
    
    # Degree values
    x <- as.numeric(names(degree_counts))
    
    cur_out <- data.frame(Sierra=s, Season=seas, links=x, prob_x=y)
    
    if(is.null(output_degree_insects)){
      output_degree_insects <- cur_out
    } else {
      output_degree_insects <- rbind(output_degree_insects, cur_out)
    }
  }
}

# ============================================================================
# DEGREE DISTRIBUTION FITTING - SINGLE MODEL PER TROPHIC LEVEL
# Power-law for Pollinators, Truncated Power-law for Plants
# ============================================================================

cat("\n=== FITTING SINGLE MODEL PER TROPHIC LEVEL ===\n")
cat("Pollinators: Power-law (P(k >= x) = x^(-(a-1)))\n")
cat("Plants: Truncated Power-law (P(k >= x) = x^(-(a-1)) * exp(-x/b))\n\n")

# ============================================================================
# POLLINATORS: Fit Power-law only
# ============================================================================

degree_dist_params_insects <- NULL
fit_quality_insects <- NULL

for(seas in unique(output_degree_insects$Season)){
  df_seas <- output_degree_insects[output_degree_insects$Season==seas,]
  for(s in unique(df_seas$Sierra)){
    df_sierra <- df_seas[df_seas$Sierra==s,][c('links','prob_x')]
    
    # Filter out zero probabilities
    df_sierra <- df_sierra[df_sierra$prob_x > 0 & df_sierra$links > 0,]
    
    if(nrow(df_sierra) < 3) next
    
    cat(sprintf("\nFitting Pollinators (Power-law) - %s, %s (n=%d points)\n", seas, s, nrow(df_sierra)))
    
    # Fit Power-law model
    mod_powerlaw <- tryCatch({
      nls(prob_x ~ links^(-a+1), 
          data = df_sierra, 
          start = list(a = 2.0), 
          control = nls.control(maxiter = 1e3))
    }, error = function(e) { 
      cat("  Power-Law: FAILED -", e$message, "\n")
      NULL 
    })
    
    if(!is.null(mod_powerlaw)){
      # Extract model summary
      model_summary <- summary(mod_powerlaw)
      
      # Calculate R-squared
      residuals <- residuals(mod_powerlaw)
      ss_res <- sum(residuals^2)
      ss_tot <- sum((df_sierra$prob_x - mean(df_sierra$prob_x))^2)
      r_squared <- 1 - (ss_res / ss_tot)
      
      # Calculate RMSE
      rmse <- sqrt(mean(residuals^2))
      
      cat(sprintf("  SUCCESS: a = %.3f ± %.3f, R² = %.3f, RMSE = %.4f, AIC = %.2f\n", 
                  model_summary$coefficients[1,1], 
                  model_summary$coefficients[1,2],
                  r_squared, rmse, AIC(mod_powerlaw)))
      
      # Store parameters
      params_row <- data.frame(
        Trophic_level = 'Insects',
        Season = seas,
        Sierra = s,
        model = 'PowerLaw',
        a_ccdf = model_summary$coefficients[1,1] - 1,  # FIXED: a-1 for CCDF exponent
        a_pdf = model_summary$coefficients[1,1],  # FIXED: a for PDF exponent
        a.std.err = model_summary$coefficients[1,2],
        a.pval = model_summary$coefficients[1,4],
        b = NA,
        b.std.err = NA,
        b.pval = NA,
        mu = NA,
        sigma = NA,
        stringsAsFactors = FALSE
      )
      
      # Store fit quality
      quality_row <- data.frame(
        Trophic_level = 'Insects',
        Season = seas,
        Sierra = s,
        model = 'PowerLaw',
        n_points = nrow(df_sierra),
        AIC = AIC(mod_powerlaw),
        logLik = as.numeric(logLik(mod_powerlaw)),
        R_squared = r_squared,
        RMSE = rmse,
        residual_SE = model_summary$sigma,
        converged = TRUE,
        stringsAsFactors = FALSE
      )
      
      if(is.null(degree_dist_params_insects)){
        degree_dist_params_insects <- params_row
      } else {
        degree_dist_params_insects <- rbind(degree_dist_params_insects, params_row)
      }
      
      if(is.null(fit_quality_insects)){
        fit_quality_insects <- quality_row
      } else {
        fit_quality_insects <- rbind(fit_quality_insects, quality_row)
      }
      
    } else {
      # Store failed fit
      quality_row <- data.frame(
        Trophic_level = 'Insects',
        Season = seas,
        Sierra = s,
        model = 'PowerLaw',
        n_points = nrow(df_sierra),
        AIC = NA,
        logLik = NA,
        R_squared = NA,
        RMSE = NA,
        residual_SE = NA,
        converged = FALSE,
        stringsAsFactors = FALSE
      )
      
      if(is.null(fit_quality_insects)){
        fit_quality_insects <- quality_row
      } else {
        fit_quality_insects <- rbind(fit_quality_insects, quality_row)
      }
    }
  }
}


# ============================================================================
# PLANTS: Fit Truncated Power-law only
# ============================================================================

degree_dist_params_plants <- NULL
fit_quality_plants <- NULL

for(seas in unique(output_degree_plants$Season)){
  df_seas <- output_degree_plants[output_degree_plants$Season==seas,]
  for(s in unique(df_seas$Sierra)){
    df_sierra <- df_seas[df_seas$Sierra==s,][c('links','prob_x')]
    
    # Filter out zero probabilities
    df_sierra <- df_sierra[df_sierra$prob_x > 0 & df_sierra$links > 0,]
    
    if(nrow(df_sierra) < 3) next
    
    cat(sprintf("\nFitting Plants (Truncated Power-law) - %s, %s (n=%d points)\n", seas, s, nrow(df_sierra)))
    
    # Fit Truncated Power-law model
    mod_trunc <- tryCatch({
      nls(prob_x ~ (links^(-a+1)) * exp(-links/b), 
          data = df_sierra, 
          start = list(a = 2.0, b = 10), 
          control = nls.control(maxiter = 1e3))
    }, error = function(e) { 
      cat("  Truncated Power-Law: FAILED -", e$message, "\n")
      NULL 
    })
    
    if(!is.null(mod_trunc)){
      # Extract model summary
      model_summary <- summary(mod_trunc)
      
      # Calculate R-squared
      residuals <- residuals(mod_trunc)
      ss_res <- sum(residuals^2)
      ss_tot <- sum((df_sierra$prob_x - mean(df_sierra$prob_x))^2)
      r_squared <- 1 - (ss_res / ss_tot)
      
      # Calculate RMSE
      rmse <- sqrt(mean(residuals^2))
      
      cat(sprintf("  SUCCESS: a = %.3f ± %.3f, b = %.2f ± %.2f, R² = %.3f, RMSE = %.4f, AIC = %.2f\n",
                  model_summary$coefficients[1,1],
                  model_summary$coefficients[1,2],
                  model_summary$coefficients[2,1],
                  model_summary$coefficients[2,2],
                  r_squared, rmse, AIC(mod_trunc)))
      
      # Store parameters
      params_row <- data.frame(
        Trophic_level = 'Plants',
        Season = seas,
        Sierra = s,
        model = 'Truncated_PowerLaw',
        a_ccdf = model_summary$coefficients[1,1] - 1,  # FIXED: a-1 for CCDF exponent
        a_pdf = model_summary$coefficients[1,1],  # FIXED: a for PDF exponent
        a.std.err = model_summary$coefficients[1,2],
        a.pval = model_summary$coefficients[1,4],
        b = model_summary$coefficients[2,1],
        b.std.err = model_summary$coefficients[2,2],
        b.pval = model_summary$coefficients[2,4],
        mu = NA,
        sigma = NA,
        stringsAsFactors = FALSE
      )
      
      # Store fit quality
      quality_row <- data.frame(
        Trophic_level = 'Plants',
        Season = seas,
        Sierra = s,
        model = 'Truncated_PowerLaw',
        n_points = nrow(df_sierra),
        AIC = AIC(mod_trunc),
        logLik = as.numeric(logLik(mod_trunc)),
        R_squared = r_squared,
        RMSE = rmse,
        residual_SE = model_summary$sigma,
        converged = TRUE,
        stringsAsFactors = FALSE
      )
      
      if(is.null(degree_dist_params_plants)){
        degree_dist_params_plants <- params_row
      } else {
        degree_dist_params_plants <- rbind(degree_dist_params_plants, params_row)
      }
      
      if(is.null(fit_quality_plants)){
        fit_quality_plants <- quality_row
      } else {
        fit_quality_plants <- rbind(fit_quality_plants, quality_row)
      }
      
    } else {
      # Store failed fit
      quality_row <- data.frame(
        Trophic_level = 'Plants',
        Season = seas,
        Sierra = s,
        model = 'Truncated_PowerLaw',
        n_points = nrow(df_sierra),
        AIC = NA,
        logLik = NA,
        R_squared = NA,
        RMSE = NA,
        residual_SE = NA,
        converged = FALSE,
        stringsAsFactors = FALSE
      )
      
      if(is.null(fit_quality_plants)){
        fit_quality_plants <- quality_row
      } else {
        fit_quality_plants <- rbind(fit_quality_plants, quality_row)
      }
    }
  }
}


# ============================================================================
# COMBINE AND CREATE SUMMARY TABLES
# ============================================================================

cat("\n\n=== CREATING SUMMARY TABLES ===\n")

# Combine parameters
all_params <- rbind(degree_dist_params_insects, degree_dist_params_plants)
write.csv(all_params, "degree_distribution_parameters.csv", row.names = FALSE)
cat("\nParameters saved: degree_distribution_parameters.csv\n")

# Combine fit quality
all_fit_quality <- rbind(fit_quality_insects, fit_quality_plants)
write.csv(all_fit_quality, "degree_distribution_fit_quality.csv", row.names = FALSE)
cat("Fit quality saved: degree_distribution_fit_quality.csv\n")

# Create comprehensive table (parameters + fit quality)
comprehensive_table <- all_params %>%
  left_join(all_fit_quality, by = c("Trophic_level", "Season", "Sierra", "model"))

write.csv(comprehensive_table, "degree_distribution_comprehensive.csv", row.names = FALSE)
cat("Comprehensive table saved: degree_distribution_comprehensive.csv\n")

# Filter for 2007-2008 season
params_2008 <- all_params %>% filter(Season == "2007-2008")
fit_quality_2008 <- all_fit_quality %>% filter(Season == "2007-2008")
comprehensive_2008 <- comprehensive_table %>% filter(Season == "2007-2008")

write.csv(params_2008, "degree_distribution_parameters_2007-2008.csv", row.names = FALSE)
write.csv(fit_quality_2008, "degree_distribution_fit_quality_2007-2008.csv", row.names = FALSE)
write.csv(comprehensive_2008, "degree_distribution_comprehensive_2007-2008.csv", row.names = FALSE)

cat("2007-2008 season tables saved\n")

# ============================================================================
# SUPPLEMENTARY TABLE S1: Degree Distribution Model Fits
# ============================================================================

cat("\n=== GENERATING SUPPLEMENTARY TABLE S1: DEGREE DISTRIBUTION FITS ===\n")

# Create comprehensive supplementary table
supplementary_table_S1 <- comprehensive_2008 %>%
  select(Sierra, Trophic_level, model, 
         a_pdf, a.std.err, b, b.std.err,
         n_points, R_squared, RMSE, AIC) %>%
  arrange(Trophic_level, Sierra) %>%
  mutate(
    alpha_formatted = sprintf("%.3f ± %.3f", a_pdf, a.std.err),
    b_formatted = ifelse(!is.na(b), sprintf("%.2f ± %.2f", b, b.std.err), "—"),
    R2_formatted = sprintf("%.3f", R_squared),
    RMSE_formatted = sprintf("%.4f", RMSE),
    AIC_formatted = sprintf("%.2f", AIC)
  ) %>%
  select(
    Sierra,
    `Trophic Level` = Trophic_level,
    Model = model,
    `Data Points` = n_points,
    `α (±SE)` = alpha_formatted,
    `b (±SE)` = b_formatted,
    `R²` = R2_formatted,
    RMSE = RMSE_formatted,
    AIC = AIC_formatted
  )

# Save supplementary table
write.csv(supplementary_table_S1, 
          "supplementary_table_S1_degree_distributions.csv", 
          row.names = FALSE)

cat("\n")
print(supplementary_table_S1, row.names = FALSE)
cat("\n✓ Table S1 saved: supplementary_table_S1_degree_distributions.csv\n")

# ============================================================================
# Calculate summary statistics for main text citation
# ============================================================================

cat("\n=== SUMMARY STATISTICS FOR CITING IN MAIN TEXT ===\n\n")

# Pollinator summary
poll_summary <- comprehensive_2008 %>%
  filter(Trophic_level == "Insects") %>%
  summarize(
    n = n(),
    mean_R2 = mean(R_squared, na.rm = TRUE),
    sd_R2 = sd(R_squared, na.rm = TRUE),
    min_R2 = min(R_squared, na.rm = TRUE),
    max_R2 = max(R_squared, na.rm = TRUE)
  )

# Plant summary  
plant_summary <- comprehensive_2008 %>%
  filter(Trophic_level == "Plants") %>%
  summarize(
    n = n(),
    mean_R2 = mean(R_squared, na.rm = TRUE),
    sd_R2 = sd(R_squared, na.rm = TRUE),
    min_R2 = min(R_squared, na.rm = TRUE),
    max_R2 = max(R_squared, na.rm = TRUE)
  )

cat("POLLINATORS:\n")
cat(sprintf("  n = %d networks\n", poll_summary$n))
cat(sprintf("  Mean R² = %.2f (SD = %.2f)\n", poll_summary$mean_R2, poll_summary$sd_R2))
cat(sprintf("  Range: %.2f - %.2f\n\n", poll_summary$min_R2, poll_summary$max_R2))

cat("PLANTS:\n")
cat(sprintf("  n = %d networks\n", plant_summary$n))
cat(sprintf("  Mean R² = %.2f (SD = %.2f)\n", plant_summary$mean_R2, plant_summary$sd_R2))
cat(sprintf("  Range: %.2f - %.2f\n\n", plant_summary$min_R2, plant_summary$max_R2))

cat("=== SUGGESTED TEXT FOR RESULTS SECTION ===\n")
cat(sprintf("Model fits were consistently high across all sierras (Table S1): pollinator networks averaged R² = %.2f (range: %.2f-%.2f), and plant networks averaged R² = %.2f (range: %.2f-%.2f).\n",
            poll_summary$mean_R2, poll_summary$min_R2, poll_summary$max_R2,
            plant_summary$mean_R2, plant_summary$min_R2, plant_summary$max_R2))

# ============================================================================
# FIGURE 2: DEGREE DISTRIBUTION SCALING WITH AREA
# ============================================================================

cat("\n\n=== GENERATING FIGURE 2 ===\n")

# Color scheme
color_pollinator <- "#E69F00"  
color_plant <- "#009E73"

# Merge with sierra metadata (assuming sierras data frame exists)
degree_sierras <- rbind(
  degree_dist_params_insects %>% mutate(Trophic_level = 'Insects'),
  degree_dist_params_plants %>% mutate(Trophic_level = 'Plants')
)

# Merge with area data
if(exists("sierras")) {
  degree_sierras <- merge(degree_sierras, sierras, by='Sierra')
} else {
  cat("WARNING: 'sierras' data frame not found. Please load sierra metadata.\n")
}

# Add Difuntito flag
degree_sierras$is_difuntito <- degree_sierras$Sierra == 'Difuntito'

# Filter for 2007-2008 season and prepare data
fig2_poll <- degree_sierras %>%
  filter(Season == '2007-2008', Trophic_level == 'Insects') %>%
  mutate(log_area = log10(as.numeric(area)))

fig2_plant <- degree_sierras %>%
  filter(Season == '2007-2008', Trophic_level == 'Plants') %>%
  mutate(log_area = log10(as.numeric(area)))

# Separate Difuntito for analysis
fig2_poll_analysis <- fig2_poll %>% filter(Sierra != 'Difuntito')
fig2_plant_analysis <- fig2_plant %>% filter(Sierra != 'Difuntito')

# Linear models: exponent vs log(area)
lm_poll_exp <- lm(a_pdf ~ log_area, data = fig2_poll_analysis)
lm_plant_exp <- lm(a_pdf ~ log_area, data = fig2_plant_analysis)

cat("\n=== FIGURE 2 STATISTICS ===\n")
cat("\nPanel B: Pollinator Power-law Exponent vs log10(Area)\n")
cat("Slope =", round(coef(lm_poll_exp)[2], 4),
    ", R² =", round(summary(lm_poll_exp)$r.squared, 3),
    ", p =", round(summary(lm_poll_exp)$coefficients[2,4], 4), "\n")
print(summary(lm_poll_exp))

cat("\nPanel C: Plant Truncated Power-law Exponent vs log10(Area)\n")
cat("Slope =", round(coef(lm_plant_exp)[2], 4),
    ", R² =", round(summary(lm_plant_exp)$r.squared, 3),
    ", p =", round(summary(lm_plant_exp)$coefficients[2,4], 4), "\n")
print(summary(lm_plant_exp))

# ============================================================================
# NEW PANEL A: Example Degree Distributions for Amarante
# ============================================================================

cat("\n=== Creating Panel A: Example Degree Distributions ===\n")

# Example sierra
example_sierra <- "Amarante"

# Get CCDF data for both trophic levels
amarante_poll_data <- output_degree_insects %>%
  filter(Sierra == example_sierra, Season == '2007-2008')

amarante_plant_data <- output_degree_plants %>%
  filter(Sierra == example_sierra, Season == '2007-2008')

# Get fitted parameters
amarante_poll_params <- degree_dist_params_insects %>%
  filter(Sierra == example_sierra, Season == '2007-2008')

amarante_plant_params <- degree_dist_params_plants %>%
  filter(Sierra == example_sierra, Season == '2007-2008')

# Create fitted lines
x_poll_seq <- exp(seq(log(min(amarante_poll_data$links)), 
                      log(max(amarante_poll_data$links)), 
                      length.out = 200))
y_poll_fit <- x_poll_seq^(-amarante_poll_params$a_ccdf)
fit_poll_data <- data.frame(links = x_poll_seq, prob_x = y_poll_fit, 
                            trophic = "Pollinators")

x_plant_seq <- exp(seq(log(min(amarante_plant_data$links)), 
                       log(max(amarante_plant_data$links)), 
                       length.out = 200))
y_plant_fit <- (x_plant_seq^(-amarante_plant_params$a_ccdf)) * 
  exp(-x_plant_seq/amarante_plant_params$b)
fit_plant_data <- data.frame(links = x_plant_seq, prob_x = y_plant_fit, 
                             trophic = "Plants")

# Combine data for plotting
amarante_poll_data$trophic <- "Pollinators"
amarante_plant_data$trophic <- "Plants"
observed_data <- rbind(amarante_poll_data[, c("links", "prob_x", "trophic")],
                       amarante_plant_data[, c("links", "prob_x", "trophic")])

fit_data_combined <- rbind(fit_poll_data, fit_plant_data)

# Get R² values
r2_poll <- all_fit_quality %>%
  filter(Sierra == example_sierra, Season == '2007-2008', Trophic_level == 'Insects') %>%
  pull(R_squared)

r2_plant <- all_fit_quality %>%
  filter(Sierra == example_sierra, Season == '2007-2008', Trophic_level == 'Plants') %>%
  pull(R_squared)

# Create Panel A
p2a <- ggplot() +
  # Fitted lines
  geom_line(data = fit_data_combined, 
            aes(x = links, y = prob_x, color = trophic), 
            linewidth = 1.3, alpha = 0.9) +
  # Observed points
  geom_point(data = observed_data, 
             aes(x = links, y = prob_x, color = trophic, shape = trophic), 
             size = 3, alpha = 0.7) +
  scale_x_log10(
    breaks = c(1, 3, 10, 30),
    labels = c("1", "3", "10", "30")
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1),
    labels = c("0.01", "0.1", "1")
  ) +
  scale_color_manual(
    values = c("Pollinators" = color_pollinator, "Plants" = color_plant),
    name = ""
  ) +
  scale_shape_manual(
    values = c("Pollinators" = 16, "Plants" = 16),
    name = ""
  ) +
  labs(
    x = "Degree (k)", 
    y = "P(K ≥ k)",
    title = "A. Example: Amarante"
  ) +
  # Add model info
  annotate("text", x = 1, y = 0.012,
           label = sprintf("Pollinators: Power-law\na = %.2f, R^2 = %.2f",
                           amarante_poll_params$a_pdf, r2_poll),
           hjust = 0, vjust = 0, size = 3.2, 
           color = color_pollinator, fontface = "italic") +
  annotate("text", x = 2.7, y = 0.63,
           label = sprintf("Plants: Truncated power-law\na = %.2f, b = %.1f, R^2 = %.2f",
                           amarante_plant_params$a_pdf, 
                           amarante_plant_params$b, r2_plant),
           hjust = 0, vjust = 0, size = 3.2, 
           color = color_plant, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )


# Panel B: Pollinator exponent vs Area
p2b <- ggplot(fig2_poll, aes(x = area, y = a_pdf)) +
  geom_point(data = fig2_poll[!fig2_poll$is_difuntito,],
             size = 4, color = color_pollinator, alpha = 0.7) +
  geom_point(data = fig2_poll[fig2_poll$is_difuntito,],
             size = 5, color = 'red', shape = 17) +
  geom_smooth(data = fig2_poll_analysis, method = 'lm', 
              color = color_pollinator, fill = color_pollinator,
              se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10",
                     breaks = c(10, 100, 1000, 10000),
                     labels = c("10", "100", "1000", "10000")) +
  labs(x = "Area", 
       y = expression(bold("Power-law Exponent ("*alpha*")")),
       title = "B. Pollinators") +
  annotate("text", x = 12, y = 1.8 ,
           label = sprintf("Slope = %.3f\nR² = %.3f\np = %.3f",
                           coef(lm_poll_exp)[2],
                           summary(lm_poll_exp)$r.squared,
                           summary(lm_poll_exp)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel C: Plant exponent vs Area
p2c <- ggplot(fig2_plant, aes(x = area, y = a_pdf)) +
  geom_point(data = fig2_plant[!fig2_plant$is_difuntito,],
             size = 4, color = color_plant, alpha = 0.7) +
  geom_point(data = fig2_plant[fig2_plant$is_difuntito,],
             size = 5, color = 'red', shape = 17) +
  geom_smooth(data = fig2_plant_analysis, method = 'lm',
              color = color_plant, fill = color_plant,
              se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10",
                     breaks = c(10, 100, 1000, 10000),
                     labels = c("10", "100", "1000", "10000")) +
  labs(x = "Area", 
       y = expression(bold("Truncated Power-law Exponent ("*alpha*")")),
       title = "C. Plants") +
  annotate("text", x = 360, y = 0.92,
           label = sprintf("Slope = %.3f\nR² = %.3f\np = %.3f",
                           coef(lm_plant_exp)[2],
                           summary(lm_plant_exp)$r.squared,
                           summary(lm_plant_exp)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# Combine panels
fig2_combined <- grid.arrange(
  p2a, p2b, p2c,
  ncol = 3,
  top = textGrob("", 
                 gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save figure
ggsave("figure2_degree_distributions_updated.pdf", fig2_combined, 
       width = 10.5, height = 3.5, units = "in", dpi = 300)

ggsave("figure2_degree_distributions_updated.png", fig2_combined, 
       width = 10.4, height = 3.4, units = "in", dpi = 300)

cat("\nFigure 2 saved: figure2_degree_distributions_updated.pdf\n")

# Also save individual panels
ggsave("figure2a_pollinators.pdf", p2a, width = 5, height = 5, units = "in", dpi = 300)
ggsave("figure2b_plants.pdf", p2b, width = 5, height = 5, units = "in", dpi = 300)

cat("Individual panels saved: figure2a_pollinators.pdf, figure2b_plants.pdf\n")


# ============================================================================
# SUPPLEMENTARY FIGURE: All Degree Distribution Fits (2007-2008)
# Shows all sierras in a comprehensive grid
# ============================================================================

cat("\n=== GENERATING SUPPLEMENTARY FIGURE: ALL SIERRAS ===\n")


color_pollinator <- "#E69F00"  
color_plant <- "#009E73"

# Function to create individual fit plot
plot_degree_fit_clean <- function(sierra_name, trophic, output_degree_data, 
                                  params_data, fit_quality_data, season = "2007-2008") {
  
  # Get observed CCDF data
  df_data <- output_degree_data %>%
    filter(Sierra == sierra_name, Season == season)
  
  # Get fitted parameters
  params <- params_data %>%
    filter(Sierra == sierra_name, Season == season)
  
  if(nrow(df_data) == 0 | nrow(params) == 0) return(NULL)
  
  # Create fine-grained x sequence for smooth fitted line
  x_min <- min(df_data$links)
  x_max <- max(df_data$links)
  x_seq <- exp(seq(log(x_min), log(x_max), length.out = 200))
  
  # Calculate fitted CCDF values
  if(trophic == "Insects") {
    # Power-law CCDF: P(k >= x) = x^(-(a-1)) = x^(-a_ccdf)
    y_fit <- x_seq^(-params$a_ccdf)
    title_text <- sprintf("%s\nPollinators", sierra_name)
    param_text <- sprintf("α = %.2f", params$a_pdf)
  } else {
    # Truncated power-law CCDF: P(k >= x) = x^(-(a-1)) * exp(-x/b)
    y_fit <- (x_seq^(-params$a_ccdf)) * exp(-x_seq/params$b)
    title_text <- sprintf("%s\nPlants", sierra_name)
    param_text <- sprintf("α = %.2f\nb = %.1f", params$a_pdf, params$b)
  }
  
  fit_data <- data.frame(links = x_seq, prob_x = y_fit)
  
  # Get fit quality
  quality <- fit_quality_data %>%
    filter(Sierra == sierra_name, Season == season, Trophic_level == trophic)
  
  r2_text <- ifelse(nrow(quality) > 0, sprintf("R² = %.3f", quality$R_squared), "")
  
  # Set color
  color <- ifelse(trophic == "Insects", color_pollinator, color_plant)
  
  # Create plot
  p <- ggplot() +
    # Fitted line first (behind points)
    geom_line(data = fit_data, aes(x = links, y = prob_x), 
              color = color, linewidth = 1.2, alpha = 0.9) +
    # Observed data points on top
    geom_point(data = df_data, aes(x = links, y = prob_x),
               size = 2.5, alpha = 0.7, color = color, shape = 16) +
    scale_x_log10(
      breaks = c(1, 3, 10, 30),
      labels = c("1", "3", "10", "30")
    ) +
    scale_y_log10(
      breaks = c(0.001, 0.01, 0.1, 1),
      labels = c("0.001", "0.01", "0.1", "1")
    ) +
    labs(
      title = title_text,
      x = "Degree (k)",
      y = "P(K ≥ k)"
    ) +
    # Parameters in top-right
    annotate("text", x = Inf, y = Inf,
             label = param_text,
             hjust = 1.05, vjust = 1.2, size = 2.8, fontface = "bold",
             lineheight = 0.85) +
    # R² in bottom-left
    annotate("text", x = x_min * 1.1, y = min(df_data$prob_x) * 1.5,
             label = r2_text,
             hjust = 0, vjust = 0, size = 2.5, fontface = "italic") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5, 
                                lineheight = 0.9, margin = margin(b = 3)),
      axis.title = element_text(face = "bold", size = 7),
      axis.title.x = element_text(margin = margin(t = 3)),  # Reduce top margin
      axis.title.y = element_text(margin = margin(r = 3)),  # Reduce right margin
      axis.text = element_text(size = 6.5),
      panel.background = element_rect(fill = "white"),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  return(p)
}

# Get all sierras from 2007-2008 season
all_sierras_2008 <- unique(degree_dist_params_insects$Sierra[
  degree_dist_params_insects$Season == "2007-2008"
])

cat("Generating plots for", length(all_sierras_2008), "sierras\n")

# Generate all plots (Pollinators and Plants for each sierra)
all_plots <- list()

for(sierra in sort(all_sierras_2008)) {
  # Pollinators
  p_poll <- plot_degree_fit_clean(sierra, "Insects", output_degree_insects, 
                                  degree_dist_params_insects, all_fit_quality)
  if(!is.null(p_poll)) {
    all_plots[[length(all_plots) + 1]] <- p_poll
    cat("  ✓", sierra, "- Pollinators\n")
  }
  
  # Plants
  p_plant <- plot_degree_fit_clean(sierra, "Plants", output_degree_plants, 
                                   degree_dist_params_plants, all_fit_quality)
  if(!is.null(p_plant)) {
    all_plots[[length(all_plots) + 1]] <- p_plant
    cat("  ✓", sierra, "- Plants\n")
  }
}

if(length(all_plots) > 0) {
  n_plots <- length(all_plots)
  n_cols <- 4
  n_rows <- ceiling(n_plots / n_cols)
  
  cat("\nCreating grid:", n_plots, "plots in", n_rows, "rows ×", n_cols, "columns\n")
  
  # Create title
  main_title <- textGrob(
    "Degree Distribution Fits - All Sierras (2007-2008)",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  # Combine plots
  comprehensive_grid <- do.call(grid.arrange, 
                                c(all_plots, 
                                  list(
                                    ncol = n_cols,
                                    top = arrangeGrob(main_title, 
                                                      heights = c(0.6, 0.4))
                                  )))
  
  # Calculate height based on number of rows
  fig_height <- n_rows * 2.5 + 1  # Add 1 inch for title
  
  ggsave("supplementary_all_sierras_degree_fits.pdf", comprehensive_grid, 
         width = 12, height = fig_height, units = "in", dpi = 300, limitsize = FALSE)
  
  cat("\n✓ Supplementary figure saved: supplementary_all_sierras_degree_fits.pdf\n")
  cat("  Dimensions:", 12, "×", fig_height, "inches\n")
  cat("  Total plots:", n_plots, "\n")
} else {
  cat("\n✗ ERROR: No plots were generated. Check data availability.\n")
}


# ============================================================================
# EXTENDED DATA FIGURE: Degree Distributions with Area-based Color Gradient
# ============================================================================

cat("\n=== GENERATING EXTENDED DATA FIGURE: DEGREE DISTRIBUTIONS BY AREA ===\n")

# Prepare data - exclude Difuntito and merge with area
degree_insects_plot <- output_degree_insects %>%
  filter(Season == '2007-2008', Sierra != 'Difuntito') %>%
  left_join(select(sierras, Sierra, area), by = "Sierra")

degree_plants_plot <- output_degree_plants %>%
  filter(Season == '2007-2008', Sierra != 'Difuntito') %>%
  left_join(select(sierras, Sierra, area), by = "Sierra")

# Panel A: Pollinators (orange-yellow gradient)
p_s3a <- ggplot(degree_insects_plot, aes(x = links, y = prob_x, 
                                         color = log10(area), group = Sierra)) +
  geom_line(size = 1.3, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_x_log10(
    breaks = c(1, 3, 10, 30),
    labels = c("1", "3", "10", "30")
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  scale_color_gradientn(
    colors = c("#FFF4A3", "#FFD93D", "#FFA500", "#FF8C00", "#E65100"),
    name = "Area (ha)",
    breaks = log10(c(10, 100, 1000)),
    labels = c("10", "100", "1000"),
    guide = guide_colorbar(barwidth = 1, barheight = 10)
  ) +
  labs(
    x = "Degree (k)",
    y = "P(K ≥ k)",
    title = "A. Pollinators"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel B: Plants (green gradient)
p_s3b <- ggplot(degree_plants_plot, aes(x = links, y = prob_x, 
                                        color = log10(area), group = Sierra)) +
  geom_line(size = 1.3, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_x_log10(
    breaks = c(1, 3, 10, 30),
    labels = c("1", "3", "10", "30")
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  scale_color_gradientn(
    colors = c("#E8F5E9", "#A5D6A7", "#66BB6A", "#43A047", "#2E7D32"),
    name = "Area (ha)",
    breaks = log10(c(10, 100, 1000)),
    labels = c("10", "100", "1000"),
    guide = guide_colorbar(barwidth = 1, barheight = 10)
  ) +
  labs(
    x = "Degree (k)",
    y = "P(K ≥ k)",
    title = "B. Plants"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

# Combine panels
extended_fig_combined <- grid.arrange(
  p_s3a, p_s3b,
  ncol = 2,
  top = textGrob("", gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save figure
ggsave("Extended_data_figure_degree_distributions_gradient.pdf", extended_fig_combined,
       width = 9, height = 3.5, units = "in", dpi = 300)

ggsave("Extended_data_figure_degree_distributions_gradient.png", extended_fig_combined,
       width = 9, height = 3.5, units = "in", dpi = 300)

cat("\n✓ Extended data figure saved: Extended_data_figure_degree_distributions_gradient.pdf/png\n")
cat("  Shows degree distributions colored by area gradient\n")
cat("  Panel A: Pollinators (orange-yellow gradient)\n")
cat("  Panel B: Plants (green gradient)\n")
cat("  Excludes Difuntito from analysis\n")

# ============================================================================
# FIGURE 3 ANALYSES: Temporal Opportunity Mechanism
# ============================================================================

cat("\n=== FIGURE 3 ANALYSES: TEMPORAL MECHANISM ===\n")

# Calculate phenophase and overlap analysis
calculate_phenophase_analysis <- function(plant_data, insect_data, sierra_data) {
  plant_filtered <- plant_data %>% filter(season == "2007-2008")
  insect_filtered <- insect_data %>% filter(season == "2007-2008")
  
  # Phenophase duration for plants
  plant_phenophase <- plant_filtered %>%
    select(Sierra, Date, Plant_species, season) %>%
    distinct() %>%
    group_by(Sierra, Plant_species, season) %>%
    summarize(
      first_observed = min(Date),
      last_observed = max(Date),
      num_observations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      phenophase_duration = as.numeric(difftime(last_observed, first_observed, units = "days")) + 1,
      trophic_level = "Plants"
    )
  
  # Phenophase duration for pollinators
  pollinator_phenophase <- insect_filtered %>%
    select(Sierra, Date, Insect_species, season) %>%
    filter(Insect_species != "SV") %>%
    distinct() %>%
    group_by(Sierra, Insect_species, season) %>%
    summarize(
      first_observed = min(Date),
      last_observed = max(Date),
      num_observations = n(),
      .groups = "drop"
    ) %>%
    mutate(
      phenophase_duration = as.numeric(difftime(last_observed, first_observed, units = "days")) + 1,
      trophic_level = "Pollinators"
    )
  
  # Get degree data
  plant_degree_data <- plant_filtered %>%
    select(Sierra, Plant_species, season, number_of_interactions_sp) %>%
    distinct() %>%
    group_by(Sierra, Plant_species, season) %>%
    summarize(degree = first(number_of_interactions_sp), .groups = "drop")
  
  pollinator_degree_data <- insect_filtered %>%
    select(Sierra, Insect_species, season, number_of_interactions_sp) %>%
    distinct() %>%
    group_by(Sierra, Insect_species, season) %>%
    summarize(degree = first(number_of_interactions_sp), .groups = "drop")
  
  # Combine with degree
  plant_analysis <- plant_phenophase %>%
    left_join(plant_degree_data, by = c("Sierra", "Plant_species", "season")) %>%
    mutate(trophic_level = "Plant")
  
  pollinator_analysis <- pollinator_phenophase %>%
    left_join(pollinator_degree_data, by = c("Sierra", "Insect_species", "season")) %>%
    mutate(trophic_level = "Pollinator")
  
  # Calculate temporal overlap
  interactions <- plant_filtered %>%
    select(Sierra, Plant_species, Insect_species, season) %>%
    distinct()
  
  temporal_overlap <- interactions %>%
    left_join(plant_phenophase, by = c("Sierra", "Plant_species", "season")) %>%
    rename(
      plant_first = first_observed,
      plant_last = last_observed,
      plant_duration = phenophase_duration
    ) %>%
    left_join(pollinator_phenophase, by = c("Sierra", "Insect_species", "season")) %>%
    rename(
      pollinator_first = first_observed,
      pollinator_last = last_observed,
      pollinator_duration = phenophase_duration
    ) %>%
    mutate(
      overlap_start = pmax(plant_first, pollinator_first),
      overlap_end = pmin(plant_last, pollinator_last),
      overlap_days = pmax(0, as.numeric(difftime(overlap_end, overlap_start, units = "days")) + 1),
      overlap_proportion = overlap_days / ((plant_duration + pollinator_duration) / 2)
    )
  
  overlap_by_sierra <- temporal_overlap %>%
    group_by(Sierra, season) %>%
    summarize(
      mean_overlap_days = mean(overlap_days, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(select(sierra_data, Sierra, area), by = "Sierra")
  
  return(list(
    plant_analysis = plant_analysis,
    pollinator_analysis = pollinator_analysis,
    overlap_by_sierra = overlap_by_sierra,
    temporal_overlap = temporal_overlap
  ))
}

# Run phenophase analysis
results <- calculate_phenophase_analysis(
  degree_dist_data_plant, 
  degree_dist_data_insect,
  sierras
)

# Prepare phenophase data
phenophase_data_combined <- bind_rows(
  select(results$plant_analysis, Sierra, phenophase_duration, num_observations, trophic_level),
  select(results$pollinator_analysis, Sierra, phenophase_duration, num_observations, trophic_level)
) %>%
  filter(num_observations > 1)

# Panel A: Temporal overlap vs area
temporal_overlap_data <- results$overlap_by_sierra %>%
  filter(season == "2007-2008", Sierra != "Difuntito") %>%
  mutate(log_area = log10(area))

temporal_overlap_all <- results$overlap_by_sierra %>%
  filter(season == "2007-2008") %>%
  mutate(is_difuntito = Sierra == "Difuntito")

lm_temporal_overlap <- lm(mean_overlap_days ~ log_area, data = temporal_overlap_data)

cat("\nPanel A: Temporal Overlap vs Area\n")
cat("Slope =", round(coef(lm_temporal_overlap)[2], 3),
    ", R² =", round(summary(lm_temporal_overlap)$r.squared, 3),
    ", p =", round(summary(lm_temporal_overlap)$coefficients[2,4], 4), "\n")

# Panel B: Potential partners for BOTH trophic levels

# For pollinators
pollinator_phenology_for_plants <- results$pollinator_analysis %>%
  filter(num_observations > 1, season == "2007-2008") %>%
  select(Sierra, Insect_species, first_observed, last_observed)

plant_phenology_with_degree <- results$plant_analysis %>%
  filter(num_observations > 1, season == "2007-2008") %>%
  select(Sierra, Plant_species, first_observed, last_observed, degree)

# Calculate potential pollinator partners for each plant
plant_opportunities <- plant_phenology_with_degree %>%
  left_join(pollinator_phenology_for_plants, by = "Sierra", relationship = "many-to-many") %>%
  mutate(
    overlap_start = pmax(first_observed.x, first_observed.y),
    overlap_end = pmin(last_observed.x, last_observed.y),
    overlap_days = pmax(0, as.numeric(difftime(overlap_end, overlap_start, units = "days")) + 1),
    overlaps = overlap_days > 0
  ) %>%
  group_by(Sierra, Plant_species) %>%
  summarize(
    n_potential_partners = sum(overlaps),
    mean_overlap_with_pollinators = mean(overlap_days[overlaps]),
    degree = first(degree),
    .groups = "drop"
  )

plant_opportunities_by_sierra <- plant_opportunities %>%
  group_by(Sierra) %>%
  summarize(
    mean_potential_partners = mean(n_potential_partners),
    mean_overlap_duration = mean(mean_overlap_with_pollinators, na.rm = TRUE),
    n_plants = n(),
    .groups = "drop"
  ) %>%
  left_join(select(sierras, Sierra, area), by = "Sierra") %>%
  mutate(log_area = log10(area),
         is_difuntito = Sierra == 'Difuntito')

# Get plant phenology
plant_phenology <- results$plant_analysis %>%
  filter(num_observations > 1, season == "2007-2008") %>%
  select(Sierra, Plant_species, first_observed, last_observed)

# Get pollinator phenology  
pollinator_phenology <- results$pollinator_analysis %>%
  filter(num_observations > 1, season == "2007-2008") %>%
  select(Sierra, Insect_species, first_observed, last_observed, degree)

# For each pollinator, count how many plants overlapped temporally
pollinator_opportunities <- pollinator_phenology %>%
  left_join(plant_phenology, by = "Sierra", relationship = "many-to-many") %>%
  mutate(
    overlap_start = pmax(first_observed.x, first_observed.y),
    overlap_end = pmin(last_observed.x, last_observed.y),
    overlap_days = pmax(0, as.numeric(difftime(overlap_end, overlap_start, units = "days")) + 1),
    overlaps = overlap_days > 0
  ) %>%
  group_by(Sierra, Insect_species) %>%
  summarize(
    n_potential_partners = sum(overlaps),
    mean_overlap_with_plants = mean(overlap_days[overlaps]),
    degree = first(degree),
    .groups = "drop"
  )

opportunities_by_sierra <- pollinator_opportunities %>%
  group_by(Sierra) %>%
  summarize(
    mean_potential_partners = mean(n_potential_partners),
    mean_overlap_duration = mean(mean_overlap_with_plants, na.rm = TRUE),
    n_pollinators = n(),
    .groups = "drop"
  ) %>%
  left_join(select(sierras, Sierra, area), by = "Sierra") %>%
  mutate(log_area = log10(area),
         is_difuntito = Sierra == 'Difuntito')

# Models for Panel B (exclude Difuntito from statistical analyses)
lm_pollinator_area_partners <- lm(mean_potential_partners ~ log_area, 
                                  data = opportunities_by_sierra %>% filter(Sierra != 'Difuntito'))
lm_plant_area_partners <- lm(mean_potential_partners ~ log_area, 
                             data = plant_opportunities_by_sierra %>% filter(Sierra != 'Difuntito'))

cat("\nPanel B: Potential Partners vs Area\n")
cat("Pollinators: Slope =", round(coef(lm_pollinator_area_partners)[2], 3),
    ", R² =", round(summary(lm_pollinator_area_partners)$r.squared, 3),
    ", p =", round(summary(lm_pollinator_area_partners)$coefficients[2,4], 4), "\n")
cat("Plants: Slope =", round(coef(lm_plant_area_partners)[2], 3),
    ", R² =", round(summary(lm_plant_area_partners)$r.squared, 3),
    ", p =", round(summary(lm_plant_area_partners)$coefficients[2,4], 4), "\n")

# Panel C: Links/species vs potential partners
opportunities_links <- opportunities_by_sierra %>%
  left_join(
    overall_data_metadata %>%
      filter(season == '2007-2008') %>%
      select(Sierra, links.sp_insect) %>%
      distinct(),
    by = "Sierra"
  )

# For statistical models, exclude Difuntito
opportunities_links_analysis <- opportunities_links %>% filter(Sierra != 'Difuntito')

plant_opportunities_links <- plant_opportunities_by_sierra %>%
  left_join(
    overall_data_metadata %>%
      filter(season == '2007-2008') %>%
      select(Sierra, links.sp_plant) %>%
      distinct(),
    by = "Sierra"
  )

# For statistical models, exclude Difuntito
plant_opportunities_links_analysis <- plant_opportunities_links %>% filter(Sierra != 'Difuntito')

lm_pollinator_partners_links <- lm(links.sp_insect ~ mean_potential_partners, data = opportunities_links_analysis)
lm_plant_partners_links <- lm(links.sp_plant ~ mean_potential_partners, data = plant_opportunities_links_analysis)

cat("\nPanel C: Links/Species vs Potential Partners\n")
cat("Pollinators: Slope =", round(coef(lm_pollinator_partners_links)[2], 3),
    ", R² =", round(summary(lm_pollinator_partners_links)$r.squared, 3),
    ", p =", round(summary(lm_pollinator_partners_links)$coefficients[2,4], 4), "\n")
cat("Plants: Slope =", round(coef(lm_plant_partners_links)[2], 3),
    ", R² =", round(summary(lm_plant_partners_links)$r.squared, 3),
    ", p =", round(summary(lm_plant_partners_links)$coefficients[2,4], 4), "\n")

# Panel D: Species-level
lm_pollinator_species <- lm(degree ~ n_potential_partners, data = pollinator_opportunities)
lm_plant_species <- lm(degree ~ n_potential_partners, data = plant_opportunities)

cat("\nPanel D: Species-level (Degree vs Potential Partners)\n")
cat("Pollinators: Slope =", round(coef(lm_pollinator_species)[2], 3),
    ", R² =", round(summary(lm_pollinator_species)$r.squared, 3),
    ", p =", round(summary(lm_pollinator_species)$coefficients[2,4], 4),
    ", n =", nrow(pollinator_opportunities), "\n")
cat("Plants: Slope =", round(coef(lm_plant_species)[2], 3),
    ", R² =", round(summary(lm_plant_species)$r.squared, 3),
    ", p =", round(summary(lm_plant_species)$coefficients[2,4], 4),
    ", n =", nrow(plant_opportunities), "\n")

# Prepare comparison data for plotting
comparison_data <- data.frame(
  Sierra = opportunities_by_sierra$Sierra,
  log_area = opportunities_by_sierra$log_area,
  pollinator_potential_partners = opportunities_by_sierra$mean_potential_partners,
  plant_potential_partners = plant_opportunities_by_sierra$mean_potential_partners,
  links_per_pollinator = opportunities_links$links.sp_insect,
  links_per_plant = plant_opportunities_links$links.sp_plant,
  is_difuntito = opportunities_by_sierra$is_difuntito
)

# ============================================================================
# FIGURE 3: Temporal Mechanism (4 panels with BOTH trophic levels)
# ============================================================================

# Calculate temporal overlap vs area model (already done above)

# Panel A: Temporal overlap vs Area (both trophic levels - same data)
p3a <- ggplot(temporal_overlap_all, aes(x = area, y = mean_overlap_days)) +
  geom_point(data = temporal_overlap_all[!temporal_overlap_all$is_difuntito,],
             size = 4, color = "#2E86AB", alpha = 0.7) +
  geom_point(data = temporal_overlap_all[temporal_overlap_all$is_difuntito,],
             size = 4, color = "red", shape = 17) +
  geom_smooth(data = temporal_overlap_data, method = "lm",
              color = "#2E86AB", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area", y = "Mean Temporal Overlap (days)",title = "A. Temporal overlap") +
  annotate("text", x=11, y= 40, #x = 800, y = 20,
           label = sprintf("Slope = %.2f\nR² = %.2f\np = %.3f",
                           coef(lm_temporal_overlap)[2],
                           summary(lm_temporal_overlap)$r.squared,
                           summary(lm_temporal_overlap)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 0, 5)
  )

# Panel B: Mean potential partners vs Area (BOTH trophic levels)
p3b <- ggplot(opportunities_by_sierra, aes(x = area)) +
  # Pollinators (non-Difuntito)
  geom_point(data = opportunities_by_sierra[!opportunities_by_sierra$is_difuntito,],
             aes(y = mean_potential_partners), 
             size = 4, color = color_pollinator, alpha = 0.7) +
  # Pollinators (Difuntito)
  geom_point(data = opportunities_by_sierra[opportunities_by_sierra$is_difuntito,],
             aes(y = mean_potential_partners), 
             size = 5, color = 'red', shape = 17) +
  geom_smooth(data = opportunities_by_sierra[!opportunities_by_sierra$is_difuntito,],
              aes(y = mean_potential_partners), 
              method = "lm", color = color_pollinator, fill = color_pollinator,
              alpha = 0.2, linewidth = 1.2) +
  # Plants (non-Difuntito)
  geom_point(data = plant_opportunities_by_sierra[!plant_opportunities_by_sierra$is_difuntito,], 
             aes(y = mean_potential_partners), 
             size = 4, color = color_plant, alpha = 0.7) +
  # Plants (Difuntito)
  geom_point(data = plant_opportunities_by_sierra[plant_opportunities_by_sierra$is_difuntito,], 
             aes(y = mean_potential_partners), 
             size = 5, color = 'red', shape = 17) +
  geom_smooth(data = plant_opportunities_by_sierra[!plant_opportunities_by_sierra$is_difuntito,], 
              aes(y = mean_potential_partners),
              method = "lm", color = color_plant, fill = color_plant,
              alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area", y = "Mean Potential Partners", title = "B. Potential partners") +
  annotate("text", x = 10, y = 31.2, #x = 60, y = 5,
           label = sprintf("Slope = %.2f, R² = %.2f, p = %.3f",
                           coef(lm_pollinator_area_partners)[2],
                           summary(lm_pollinator_area_partners)$r.squared,
                           summary(lm_pollinator_area_partners)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic",color = color_pollinator) +
  annotate("text", x = 10, y = 33.2, #x = 60, y = 5,
           label = sprintf("Slope = %.2f, R² = %.2f, p = %.3f",
                           coef(lm_plant_area_partners)[2],
                           summary(lm_plant_area_partners)$r.squared,
                           summary(lm_plant_area_partners)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic",color = color_plant) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 0, 5)
  )

# Panel C: Links/species vs potential partners (BOTH trophic levels)
p3c <- ggplot(comparison_data, aes(x = pollinator_potential_partners)) +
  # Pollinators (non-Difuntito)
  geom_point(data = comparison_data[!comparison_data$is_difuntito,],
             aes(y = links_per_pollinator), 
             size = 4, color = color_pollinator, alpha = 0.7) +
  # Pollinators (Difuntito)
  geom_point(data = comparison_data[comparison_data$is_difuntito,],
             aes(y = links_per_pollinator), 
             size = 5, color = 'red', shape = 17) +
  geom_smooth(data = comparison_data[!comparison_data$is_difuntito,],
              aes(y = links_per_pollinator), 
              method = "lm", color = color_pollinator, fill = color_pollinator,
              alpha = 0.2, linewidth = 1.2) +
  # Plants (non-Difuntito - note: using plant_potential_partners for x-axis)
  geom_point(data = comparison_data[!comparison_data$is_difuntito,],
             aes(x = plant_potential_partners, y = links_per_plant), 
             size = 4, color = color_plant, alpha = 0.7) +
  # Plants (Difuntito)
  geom_point(data = comparison_data[comparison_data$is_difuntito,],
             aes(x = plant_potential_partners, y = links_per_plant), 
             size = 5, color = 'red', shape = 17) +
  geom_smooth(data = comparison_data[!comparison_data$is_difuntito,],
              aes(x = plant_potential_partners, y = links_per_plant), 
              method = "lm", color = color_plant, fill = color_plant,
              alpha = 0.2, linewidth = 1.2) +
  labs(x = "Mean Potential Partners", y = "Links per Species",title = "C. Realised interactions") +
  annotate("text", x = 5, y = 8,
           label = sprintf("Slope = %.3f, R² = %.2f, p = %.3f",
                           coef(lm_pollinator_partners_links)[2],
                           summary(lm_pollinator_partners_links)$r.squared,
                           summary(lm_pollinator_partners_links)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic",color=color_pollinator) +
  annotate("text", x = 5, y = 8.5,
           label = sprintf("Slope = %.3f, R² = %.2f, p = %.3f",
                           coef(lm_plant_partners_links)[2],
                           summary(lm_plant_partners_links)$r.squared,
                           summary(lm_plant_partners_links)$coefficients[2,4]),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic",color=color_plant) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(0, 5, 5, 5)
  )

# Panel D: Number of interactions vs potential partners (BOTH trophic levels)
set.seed(123)
sample_poll_fig3 <- pollinator_opportunities %>%
  sample_n(min(150, nrow(pollinator_opportunities)))

set.seed(456)
sample_plant_fig3 <- plant_opportunities %>%
  sample_n(min(150, nrow(plant_opportunities)))

p3d <- ggplot(sample_poll_fig3, aes(x = n_potential_partners, y = degree)) +
  # Pollinators
  geom_point(size = 2, alpha = 0.4, color = color_pollinator) +
  geom_smooth(method = "lm", color = color_pollinator, fill = color_pollinator,
              alpha = 0.2, linewidth = 1.2) +
  # Plants
  geom_point(data = sample_plant_fig3, 
             aes(x = n_potential_partners, y = degree),
             size = 2, alpha = 0.4, color = color_plant) +
  geom_smooth(data = sample_plant_fig3, 
              aes(x = n_potential_partners, y = degree),
              method = "lm", color = color_plant, fill = color_plant,
              alpha = 0.2, linewidth = 1.2) +
  labs(x = "Potential Partners", y = "Observed Interactions", title = "D. Species-level pattern") +
  annotate("text", x = 0 , y = 30,
           label = sprintf("Slope = %.2f, R² = %.2f, p < 0.001, n = %d\n",
                           coef(lm_pollinator_species)[2],
                           summary(lm_pollinator_species)$r.squared,
                           390),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic",color=color_pollinator) +
  annotate("text", x = 0 , y = 30,
           label = sprintf("Slope = %.2f, R² = %.2f, p < 0.001, n = %d",
                           coef(lm_plant_species)[2],
                           summary(lm_plant_species)$r.squared,
                           180),
           hjust = 0, vjust = 0, size = 3.2, fontface = "italic",color=color_plant) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(0, 5, 5, 5)
  )

# Combine Figure 3
fig3_combined <- grid.arrange(
  p3a, p3b,
  p3c, p3d,
  ncol = 2,
  top = textGrob("", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

ggsave("figure3_temporal_mechanism_both.pdf", fig3_combined, 
       width = 6.6, height = 6.4, units = "in", dpi = 300)
ggsave("figure3_temporal_mechanism_both.png", fig3_combined, 
       width = 6.6, height = 6.4, units = "in", dpi = 300)

cat("\n=== Figure 3 (with both trophic levels) saved ===\n")

# ============================================================================

# ============================================================================
# FIGURE 4 ANALYSES: Dissimilarity and Phenophase
# ============================================================================

cat("\n=== FIGURE 4 ANALYSES: BIOLOGICAL TRAITS ===\n")

# Temporal beta diversity
calculate_beta_diversity_temporal <- function(data_all) {
  data_2008 <- data_all %>% filter(season == "2007-2008")
  output_beta <- NULL
  
  for(sierra in unique(data_2008$Sierra)) {
    cur_sierra <- na.omit(data_2008[data_2008$Sierra == sierra, ])
    
    if(length(unique(cur_sierra$Date)) < 2) next
    
    # For plants
    data_beta_plants <- tapply(rep(1, nrow(cur_sierra)), 
                               list(cur_sierra$Date, cur_sierra$Plant_species), 
                               FUN = sum, default = 0)
    beta_plants_quanti <- NA
    tryCatch({
      beta_plants_quanti <- mean(vegdist(data_beta_plants, method = "bray"))
    }, error = function(e) {})
    
    beta_plants_quali <- NA
    tryCatch({
      data_beta_plants_quali <- data_beta_plants
      data_beta_plants_quali[data_beta_plants_quali > 0] <- 1
      beta_plants_quali <- beta.multi(data_beta_plants_quali, index.family = "sorensen")
    }, error = function(e) {})
    
    # For insects
    insect_data <- cur_sierra
    if("Insect_species" %in% colnames(insect_data)) {
      insect_data <- insect_data[insect_data$Insect_species != 'SV', ]
    }
    
    data_beta_insects <- tapply(rep(1, nrow(insect_data)), 
                                list(insect_data$Date, insect_data$Insect_species), 
                                FUN = sum, default = 0)
    beta_insects_quanti <- NA
    tryCatch({
      beta_insects_quanti <- mean(vegdist(data_beta_insects, method = "bray"))
    }, error = function(e) {})
    
    beta_insects_quali <- NA
    tryCatch({
      data_beta_insects_quali <- data_beta_insects
      data_beta_insects_quali[data_beta_insects_quali > 0] <- 1
      beta_insects_quali <- beta.multi(data_beta_insects_quali, index.family = "sorensen")
    }, error = function(e) {})
    
    # For interactions
    beta_interactions_quali <- NA
    if("interaction" %in% colnames(insect_data)) {
      data_beta_interactions <- table(insect_data$Date, insect_data$interaction)
      data_beta_interactions[data_beta_interactions > 0] <- 1
      
      tryCatch({
        beta_interactions_quali <- beta.multi(data_beta_interactions, index.family = "sorensen")
      }, error = function(e) {})
    }
    
    cur_out <- data.frame(
      Sierra = sierra,
      beta_plants_qualitative = ifelse(is.list(beta_plants_quali), beta_plants_quali$beta.SIM, NA),
      beta_plants_quantitative = beta_plants_quanti,
      beta_insects_qualitative = ifelse(is.list(beta_insects_quali), beta_insects_quali$beta.SIM, NA),
      beta_insects_quantitative = beta_insects_quanti,
      beta_interactions_qualitative = ifelse(is.list(beta_interactions_quali), beta_interactions_quali$beta.SIM, NA)
    )
    
    if(is.null(output_beta)) {
      output_beta <- cur_out
    } else {
      output_beta <- rbind(output_beta, cur_out)
    }
  }
  
  return(output_beta)
}

beta_results <- calculate_beta_diversity_temporal(degree_dist_data_plant)

df <- beta_results[, c('beta_plants_qualitative', 'beta_insects_qualitative', 'beta_interactions_qualitative')]
colnames(df) <- c('Plants', 'Pollinators', 'Interactions')
df_beta <- tidyr::gather(df, key = "Trophic_Level", value = "Beta_diversity")
df_beta$Trophic_Level <- factor(df_beta$Trophic_Level, levels = c("Pollinators", "Plants", "Interactions"))

cat("\nPanel B: Temporal Dissimilarity\n")
cat("Median - Pollinators:", round(median(df_beta$Beta_diversity[df_beta$Trophic_Level == "Pollinators"], na.rm = TRUE), 3), "\n")
cat("Median - Plants:", round(median(df_beta$Beta_diversity[df_beta$Trophic_Level == "Plants"], na.rm = TRUE), 3), "\n")
cat("Median - Interactions:", round(median(df_beta$Beta_diversity[df_beta$Trophic_Level == "Interactions"], na.rm = TRUE), 3), "\n")

# Spatial beta diversity
calculate_beta_diversity_spatial_simple <- function(data_all) {
  data_2008 <- data_all %>% filter(season == "2007-2008")
  sierras <- unique(data_2008$Sierra)
  n_sierras <- length(sierras)
  
  results_df <- data.frame(
    sierra1 = character(),
    sierra2 = character(),
    plant_dissimilarity = numeric(),
    pollinator_dissimilarity = numeric(),
    interaction_dissimilarity = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:(n_sierras-1)) {
    for(j in (i+1):n_sierras) {
      sierra1 <- sierras[i]
      sierra2 <- sierras[j]
      
      plants_sierra1 <- unique(data_2008$Plant_species[data_2008$Sierra == sierra1])
      plants_sierra2 <- unique(data_2008$Plant_species[data_2008$Sierra == sierra2])
      
      pollinators_sierra1 <- unique(data_2008$Insect_species[data_2008$Sierra == sierra1 & 
                                                               data_2008$Insect_species != 'SV'])
      pollinators_sierra2 <- unique(data_2008$Insect_species[data_2008$Sierra == sierra2 & 
                                                               data_2008$Insect_species != 'SV'])
      
      # Jaccard dissimilarity
      plant_intersection <- length(intersect(plants_sierra1, plants_sierra2))
      plant_union <- length(union(plants_sierra1, plants_sierra2))
      plant_jaccard <- 1 - (plant_intersection / plant_union)
      
      pollinator_intersection <- length(intersect(pollinators_sierra1, pollinators_sierra2))
      pollinator_union <- length(union(pollinators_sierra1, pollinators_sierra2))
      pollinator_jaccard <- 1 - (pollinator_intersection / pollinator_union)
      
      interaction_jaccard <- NA
      if("interaction" %in% colnames(data_2008)) {
        interactions_sierra1 <- unique(data_2008$interaction[data_2008$Sierra == sierra1 & 
                                                               data_2008$Insect_species != 'SV'])
        interactions_sierra2 <- unique(data_2008$interaction[data_2008$Sierra == sierra2 & 
                                                               data_2008$Insect_species != 'SV'])
        
        interaction_intersection <- length(intersect(interactions_sierra1, interactions_sierra2))
        interaction_union <- length(union(interactions_sierra1, interactions_sierra2))
        interaction_jaccard <- 1 - (interaction_intersection / interaction_union)
      }
      
      results_df <- rbind(results_df, data.frame(
        sierra1 = sierra1,
        sierra2 = sierra2,
        plant_dissimilarity = plant_jaccard,
        pollinator_dissimilarity = pollinator_jaccard,
        interaction_dissimilarity = interaction_jaccard,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results_df)
}

spatial_results_simple <- calculate_beta_diversity_spatial_simple(degree_dist_data_plant)

df_beta_spatial <- tidyr::gather(spatial_results_simple, 
                                 key = "trophic_level", 
                                 value = "beta_diversity", 
                                 plant_dissimilarity, 
                                 pollinator_dissimilarity, 
                                 interaction_dissimilarity)

df_beta_spatial$trophic_level <- gsub("_dissimilarity", "", df_beta_spatial$trophic_level)
df_beta_spatial$trophic_level <- factor(df_beta_spatial$trophic_level, 
                                        levels = c("pollinator", "plant", "interaction"),
                                        labels = c("Pollinators", "Plants", "Interactions"))

df_beta_spatial_clean <- df_beta_spatial %>% filter(!is.na(beta_diversity))

cat("\nPanel C: Spatial Dissimilarity\n")
cat("Median - Pollinators:", round(median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Pollinators"], na.rm = TRUE), 3), "\n")
cat("Median - Plants:", round(median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Plants"], na.rm = TRUE), 3), "\n")
cat("Median - Interactions:", round(median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Interactions"], na.rm = TRUE), 3), "\n")

cat("\nPanel A: Phenophase Duration\n")
cat("Median - Pollinators:", round(median(phenophase_data_combined$phenophase_duration[phenophase_data_combined$trophic_level == "Pollinator"], na.rm = TRUE), 1), "days\n")
cat("Median - Plants:", round(median(phenophase_data_combined$phenophase_duration[phenophase_data_combined$trophic_level == "Plant"], na.rm = TRUE), 1), "days\n")

# ============================================================================
# STATISTICAL TESTS FOR FIGURE 4: Phenological Patterns and Community Turnover
# SEASON 2007-2008
# ============================================================================

cat("\n=== STATISTICAL ANALYSES FOR FIGURE 4 (2007-2008 SEASON) ===\n\n")

# ============================================================================
# PANEL A: Phenophase Duration (Plants vs Pollinators)
# ============================================================================

cat("PANEL A: PHENOPHASE DURATION\n")
cat("=" , rep("=", 50), "\n", sep="")

# Prepare data (only species observed more than once, 2007-2008 season only)
phenophase_data_2008 <- phenophase_data_combined %>%
  filter(num_observations > 1)  # Already filtered by season in the main analysis

phenophase_plants <- phenophase_data_2008 %>%
  filter(trophic_level == "Plant") %>%
  pull(phenophase_duration)

phenophase_pollinators <- phenophase_data_2008 %>%
  filter(trophic_level == "Pollinator") %>%
  pull(phenophase_duration)

# Summary statistics
cat("\nSummary Statistics:\n")
cat("Plants: n =", length(phenophase_plants), 
    ", median =", median(phenophase_plants, na.rm = TRUE),
    ", mean =", round(mean(phenophase_plants, na.rm = TRUE), 1),
    ", SD =", round(sd(phenophase_plants, na.rm = TRUE), 1), "\n")
cat("Pollinators: n =", length(phenophase_pollinators),
    ", median =", median(phenophase_pollinators, na.rm = TRUE),
    ", mean =", round(mean(phenophase_pollinators, na.rm = TRUE), 1),
    ", SD =", round(sd(phenophase_pollinators, na.rm = TRUE), 1), "\n\n")

# Test for normality
cat("Normality tests (Shapiro-Wilk):\n")
shapiro_plants <- shapiro.test(phenophase_plants)
shapiro_pollinators <- shapiro.test(phenophase_pollinators)
cat("Plants: W =", round(shapiro_plants$statistic, 4), ", p =", 
    format.pval(shapiro_plants$p.value, digits = 3), "\n")
cat("Pollinators: W =", round(shapiro_pollinators$statistic, 4), ", p =", 
    format.pval(shapiro_pollinators$p.value, digits = 3), "\n\n")

# Since data is likely non-normal, use Mann-Whitney U test (Wilcoxon rank-sum)
wilcox_phenophase <- wilcox.test(phenophase_plants, phenophase_pollinators, 
                                 alternative = "two.sided")
cat("Mann-Whitney U test (Wilcoxon rank-sum):\n")
cat("W =", wilcox_phenophase$statistic, 
    ", p =", format.pval(wilcox_phenophase$p.value, digits = 4), "\n")

# Effect size (rank-biserial correlation)
r_effect <- abs(wilcox_phenophase$statistic / 
                  (length(phenophase_plants) * length(phenophase_pollinators)) - 0.5) * 2
cat("Effect size (rank-biserial r) =", round(r_effect, 3), "\n")

# Also perform t-test for comparison (if someone prefers parametric)
t_test_phenophase <- t.test(phenophase_plants, phenophase_pollinators)
cat("\nt-test (for comparison):\n")
cat("t =", round(t_test_phenophase$statistic, 3), 
    ", df =", round(t_test_phenophase$parameter, 1),
    ", p =", format.pval(t_test_phenophase$p.value, digits = 4), "\n\n")

# ============================================================================
# PANEL B: Temporal Dissimilarity (Pollinators vs Plants vs Interactions)
# Season 2007-2008 only (already filtered in df_beta from main analysis)
# ============================================================================

cat("\n\nPANEL B: TEMPORAL DISSIMILARITY\n")
cat("=" , rep("=", 50), "\n", sep="")

# Note: df_beta is already filtered for 2007-2008 in the main analysis
# Summary statistics for each group
temporal_summary <- df_beta %>%
  group_by(Trophic_Level) %>%
  summarize(
    n = n(),
    median = median(Beta_diversity, na.rm = TRUE),
    mean = mean(Beta_diversity, na.rm = TRUE),
    sd = sd(Beta_diversity, na.rm = TRUE),
    min = min(Beta_diversity, na.rm = TRUE),
    max = max(Beta_diversity, na.rm = TRUE)
  )

cat("\nSummary Statistics:\n")
print(temporal_summary)

# Kruskal-Wallis test (non-parametric ANOVA)
kw_temporal <- kruskal.test(Beta_diversity ~ Trophic_Level, data = df_beta)
cat("\nKruskal-Wallis test:\n")
cat("χ² =", round(kw_temporal$statistic, 3), 
    ", df =", kw_temporal$parameter,
    ", p =", format.pval(kw_temporal$p.value, digits = 4), "\n")

# Post-hoc pairwise comparisons (Dunn's test with Bonferroni correction)
cat("\nPost-hoc pairwise comparisons (Wilcoxon with Bonferroni correction):\n")
pairwise_temporal <- pairwise.wilcox.test(df_beta$Beta_diversity, 
                                          df_beta$Trophic_Level,
                                          p.adjust.method = "bonferroni")
print(pairwise_temporal)

# ============================================================================
# PANEL C: Spatial Dissimilarity (Pollinators vs Plants vs Interactions)
# Season 2007-2008 only (already filtered in df_beta_spatial_clean)
# ============================================================================

cat("\n\nPANEL C: SPATIAL DISSIMILARITY\n")
cat("=" , rep("=", 50), "\n", sep="")

# Note: df_beta_spatial_clean is already filtered for 2007-2008 in the main analysis
# Summary statistics for each group
spatial_summary <- df_beta_spatial_clean %>%
  group_by(trophic_level) %>%
  summarize(
    n = n(),
    median = median(beta_diversity, na.rm = TRUE),
    mean = mean(beta_diversity, na.rm = TRUE),
    sd = sd(beta_diversity, na.rm = TRUE),
    min = min(beta_diversity, na.rm = TRUE),
    max = max(beta_diversity, na.rm = TRUE)
  )

cat("\nSummary Statistics:\n")
print(spatial_summary)

# Kruskal-Wallis test
kw_spatial <- kruskal.test(beta_diversity ~ trophic_level, 
                           data = df_beta_spatial_clean)
cat("\nKruskal-Wallis test:\n")
cat("χ² =", round(kw_spatial$statistic, 3), 
    ", df =", kw_spatial$parameter,
    ", p =", format.pval(kw_spatial$p.value, digits = 4), "\n")

# Post-hoc pairwise comparisons
cat("\nPost-hoc pairwise comparisons (Wilcoxon with Bonferroni correction):\n")
pairwise_spatial <- pairwise.wilcox.test(df_beta_spatial_clean$beta_diversity, 
                                         df_beta_spatial_clean$trophic_level,
                                         p.adjust.method = "bonferroni")
print(pairwise_spatial)

# ============================================================================
# CREATE SUMMARY TABLE FOR SUPPLEMENTARY MATERIAL
# ============================================================================

cat("\n\n=== CREATING SUPPLEMENTARY TABLE: FIGURE 4 STATISTICS ===\n")

# Panel A statistics
panel_a_stats <- data.frame(
  Panel = "A",
  Comparison = "Plants vs Pollinators",
  Test = "Mann-Whitney U",
  Statistic = paste0("W = ", wilcox_phenophase$statistic),
  p_value = format.pval(wilcox_phenophase$p.value, digits = 4),
  Effect_size = round(r_effect, 3),
  n_group1 = length(phenophase_plants),
  n_group2 = length(phenophase_pollinators),
  median_group1 = median(phenophase_plants),
  median_group2 = median(phenophase_pollinators),
  stringsAsFactors = FALSE
)

# Panel B statistics
panel_b_stats <- data.frame(
  Panel = "B",
  Comparison = "Temporal dissimilarity",
  Test = "Kruskal-Wallis",
  Statistic = paste0("χ² = ", round(kw_temporal$statistic, 3)),
  p_value = format.pval(kw_temporal$p.value, digits = 4),
  Effect_size = NA,
  n_group1 = sum(df_beta$Trophic_Level == "Pollinators"),
  n_group2 = sum(df_beta$Trophic_Level == "Plants"),
  median_group1 = median(df_beta$Beta_diversity[df_beta$Trophic_Level == "Pollinators"], na.rm = TRUE),
  median_group2 = median(df_beta$Beta_diversity[df_beta$Trophic_Level == "Plants"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Panel C statistics
panel_c_stats <- data.frame(
  Panel = "C",
  Comparison = "Spatial dissimilarity",
  Test = "Kruskal-Wallis",
  Statistic = paste0("χ² = ", round(kw_spatial$statistic, 3)),
  p_value = format.pval(kw_spatial$p.value, digits = 4),
  Effect_size = NA,
  n_group1 = sum(df_beta_spatial_clean$trophic_level == "Pollinators"),
  n_group2 = sum(df_beta_spatial_clean$trophic_level == "Plants"),
  median_group1 = median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Pollinators"], na.rm = TRUE),
  median_group2 = median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Plants"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Combine all
figure4_stats_table <- rbind(panel_a_stats, panel_b_stats, panel_c_stats)

cat("\n")
print(figure4_stats_table, row.names = FALSE)

# Save table
write.csv(figure4_stats_table, 
          "supplementary_table_figure4_statistics.csv", 
          row.names = FALSE)

cat("\n✓ Supplementary table saved: supplementary_table_figure4_statistics.csv\n")

# ============================================================================
# SUGGESTED TEXT FOR RESULTS SECTION
# ============================================================================

cat("\n\n=== SUGGESTED TEXT FOR RESULTS ===\n\n")

cat("Panel A:\n")
cat(sprintf("Pollinator phenophases (median = %.0f days, n = %d) were significantly longer than plant flowering periods (median = %.0f days, n = %d; Mann-Whitney U: W = %.0f, p %s).\n\n",
            median(phenophase_pollinators),
            length(phenophase_pollinators),
            median(phenophase_plants),
            length(phenophase_plants),
            wilcox_phenophase$statistic,
            format.pval(wilcox_phenophase$p.value, digits = 3)))

cat("Panel B:\n")
cat(sprintf("Temporal dissimilarity differed significantly among groups (Kruskal-Wallis: χ² = %.2f, df = %d, p %s), with interactions showing the highest turnover.\n\n",
            kw_temporal$statistic,
            kw_temporal$parameter,
            format.pval(kw_temporal$p.value, digits = 3)))

cat("Panel C:\n")
cat(sprintf("Spatial dissimilarity also differed significantly among groups (Kruskal-Wallis: χ² = %.2f, df = %d, p %s), with interactions showing the highest turnover (median = %.3f) compared to plants (median = %.3f) and pollinators (median = %.3f).\n",
            kw_spatial$statistic,
            kw_spatial$parameter,
            format.pval(kw_spatial$p.value, digits = 3),
            median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Interactions"], na.rm = TRUE),
            median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Plants"], na.rm = TRUE),
            median(df_beta_spatial_clean$beta_diversity[df_beta_spatial_clean$trophic_level == "Pollinators"], na.rm = TRUE)))

cat("\n=== END OF STATISTICAL ANALYSES ===\n")

# ============================================================================
# FIGURE 4: Biological traits ------------------------------------------------
# ============================================================================

# Panel A: Phenophase duration
p4a <- ggplot(phenophase_data_combined, 
              aes(x = trophic_level, y = phenophase_duration, fill = trophic_level)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
  scale_fill_manual(values = c("Plant" = color_plant, "Pollinator" = color_pollinator)) +
  labs(x = "", y = "Phenophase Duration (days)", title = "A. Phenophase duration") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel B: Temporal dissimilarity
p4b <- ggplot(df_beta, aes(x = Trophic_Level, y = Beta_diversity, fill = Trophic_Level)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(size = 2, width = 0.2, height = 0, alpha = 0.5) +
  ylim(0.7, 1) +
  scale_fill_manual(values = c("Pollinators" = color_pollinator, 
                               "Plants" = color_plant, 
                               "Interactions" = "grey50")) +
  labs(x = "", y = "Dissimilarity in Time", title = "D. Temporal dissimilarity") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# Panel C: Spatial dissimilarity
p4c <- ggplot(df_beta_spatial_clean, 
              aes(x = trophic_level, y = beta_diversity, fill = trophic_level)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(size = 2, width = 0.2, height = 0, alpha = 0.5) +
  scale_fill_manual(values = c("Pollinators" = color_pollinator, 
                               "Plants" = color_plant, 
                               "Interactions" = "grey50")) +
  labs(x = "", y = "Dissimilarity in Space", title = "C. Spatial dissimilarity") +
  ylim(NA, 1) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

fig4_combined <- grid.arrange(
  p4a, p4b, p4c,
  ncol = 3,
  top = textGrob("", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

ggsave("figure4_biological_traits.pdf", fig4_combined, 
       width = 14, height = 5, units = "in", dpi = 300)

ggsave("figure4_biological_traits.png", fig4_combined, 
       width = 10.5, height = 3.5, units = "in", dpi = 300)


cat("\n=== ALL ANALYSES AND FIGURES COMPLETE ===\n")
cat("Figure 1: Area-dependent patterns saved\n")
cat("Figure 2: Degree distributions (CORRECTED) saved\n")
cat("Figure 3: Temporal mechanism saved\n")
cat("Figure 4: Biological traits saved\n")

# ============================================================================
# CONNECTANCE ANALYSIS: Realized vs Potential Interactions
# Similar to Galiana et al. analysis of co-occurrence vs interaction networks
# ============================================================================

cat("\n=== CONNECTANCE ANALYSIS ===\n")

# Get 2007-2008 season data
data_2008 <- overall_data %>% filter(season == "2007-2008")

# ============================================================================
# PART 1: METACOMMUNITY LEVEL (ALL SIERRAS COMBINED)
# ============================================================================

cat("\n=== METACOMMUNITY LEVEL ANALYSIS ===\n")

# REALIZED INTERACTIONS
realized_interactions <- data_2008 %>%
  select(Plant_species, Insect_species) %>%
  filter(Insect_species != 'SV') %>%
  distinct()

n_realized <- nrow(realized_interactions)

# POTENTIAL INTERACTIONS - Simple spatial co-occurrence
potential_simple <- NULL

for(sierra in unique(data_2008$Sierra)) {
  sierra_data <- data_2008 %>% filter(Sierra == sierra)
  
  plants <- unique(sierra_data$Plant_species)
  pollinators <- unique(sierra_data$Insect_species[sierra_data$Insect_species != 'SV'])
  
  sierra_pairs <- expand.grid(
    Plant_species = plants,
    Insect_species = pollinators,
    stringsAsFactors = FALSE
  )
  
  if(is.null(potential_simple)) {
    potential_simple <- sierra_pairs
  } else {
    potential_simple <- rbind(potential_simple, sierra_pairs)
  }
}

potential_simple <- potential_simple %>% distinct()
n_potential_simple <- nrow(potential_simple)

# POTENTIAL INTERACTIONS - Temporal co-occurrence
# Recalculate from phenology data (not from observed interactions)

# Get plant phenology
plant_pheno <- results$plant_analysis %>%
  filter(season == "2007-2008", num_observations > 1) %>%
  select(Sierra, Plant_species, first_observed, last_observed)

# Get pollinator phenology
poll_pheno <- results$pollinator_analysis %>%
  filter(season == "2007-2008", num_observations > 1) %>%
  select(Sierra, Insect_species, first_observed, last_observed)

# Calculate ALL potential temporal overlaps
potential_temporal_all <- NULL

for(sierra in unique(data_2008$Sierra)) {
  # Get species in this sierra
  plants_sierra <- unique(data_2008$Plant_species[data_2008$Sierra == sierra])
  polls_sierra <- unique(data_2008$Insect_species[data_2008$Sierra == sierra & 
                                                    data_2008$Insect_species != 'SV'])
  
  # Get phenology for these species
  plant_pheno_sierra <- plant_pheno %>% filter(Sierra == sierra, Plant_species %in% plants_sierra)
  poll_pheno_sierra <- poll_pheno %>% filter(Sierra == sierra, Insect_species %in% polls_sierra)
  
  # Calculate overlap for all combinations
  for(i in 1:nrow(plant_pheno_sierra)) {
    for(j in 1:nrow(poll_pheno_sierra)) {
      plant <- plant_pheno_sierra$Plant_species[i]
      poll <- poll_pheno_sierra$Insect_species[j]
      
      overlap_start <- max(plant_pheno_sierra$first_observed[i], 
                           poll_pheno_sierra$first_observed[j])
      overlap_end <- min(plant_pheno_sierra$last_observed[i], 
                         poll_pheno_sierra$last_observed[j])
      overlap_days <- as.numeric(difftime(overlap_end, overlap_start, units = "days")) + 1
      
      if(overlap_days > 0) {
        pair <- data.frame(
          Sierra = sierra,
          Plant_species = plant,
          Insect_species = poll,
          overlap_days = overlap_days,
          stringsAsFactors = FALSE
        )
        
        if(is.null(potential_temporal_all)) {
          potential_temporal_all <- pair
        } else {
          potential_temporal_all <- rbind(potential_temporal_all, pair)
        }
      }
    }
  }
}

# Get unique pairs across all sierras
potential_temporal <- potential_temporal_all %>%
  select(Plant_species, Insect_species) %>%
  distinct()

n_potential_temporal <- nrow(potential_temporal)

# CALCULATE METACOMMUNITY CONNECTANCE
connectance_simple_meta <- n_realized / n_potential_simple
connectance_temporal_meta <- n_realized / n_potential_temporal

cat("\nMetacommunity results:\n")
cat("Realized interactions: ", n_realized, "\n")
cat("Potential (simple): ", n_potential_simple, "\n")
cat("Potential (temporal): ", n_potential_temporal, "\n")
cat("Connectance (simple): ", round(connectance_simple_meta, 3), 
    " (", round(connectance_simple_meta * 100, 1), "%)\n", sep="")
cat("Connectance (temporal): ", round(connectance_temporal_meta, 3), 
    " (", round(connectance_temporal_meta * 100, 1), "%)\n", sep="")

# ============================================================================
# PART 2: SIERRA LEVEL (Connectance vs Area)
# ============================================================================

cat("\n=== SIERRA LEVEL ANALYSIS ===\n")

# Initialize results dataframe
connectance_by_sierra <- data.frame(
  Sierra = character(),
  area = numeric(),
  realized_interactions = numeric(),
  potential_simple = numeric(),
  potential_temporal = numeric(),
  connectance_simple = numeric(),
  connectance_temporal = numeric(),
  stringsAsFactors = FALSE
)

# Calculate for each sierra
for(sierra in unique(data_2008$Sierra)) {
  cat("\nProcessing:", sierra, "\n")
  
  sierra_data <- data_2008 %>% filter(Sierra == sierra)
  
  # Realized interactions
  realized <- sierra_data %>%
    select(Plant_species, Insect_species) %>%
    filter(Insect_species != 'SV') %>%
    distinct() %>%
    nrow()
  
  cat("  Realized interactions:", realized, "\n")
  
  # Number of species
  n_plants <- length(unique(sierra_data$Plant_species))
  n_pollinators <- length(unique(sierra_data$Insect_species[sierra_data$Insect_species != 'SV']))
  
  cat("  Plants:", n_plants, ", Pollinators:", n_pollinators, "\n")
  
  # Potential interactions - simple
  potential_simple <- n_plants * n_pollinators
  
  # Connectance - simple
  connectance_simple <- realized / potential_simple
  
  # Potential interactions - temporal (from our recalculated data)
  potential_temporal <- potential_temporal_all %>%
    filter(Sierra == sierra) %>%
    nrow()
  
  cat("  Potential (simple):", potential_simple, ", Potential (temporal):", potential_temporal, "\n")
  
  # Connectance - temporal
  connectance_temporal <- ifelse(potential_temporal > 0, realized / potential_temporal, NA)
  
  # Get area
  sierra_area <- sierras$area[sierras$Sierra == sierra]
  
  if(length(sierra_area) == 0) {
    cat("  WARNING: No area found for", sierra, "- skipping\n")
    next
  }
  
  cat("  Area:", sierra_area, "\n")
  cat("  Connectance (simple):", round(connectance_simple, 3), 
      ", Connectance (temporal):", round(connectance_temporal, 3), "\n")
  
  # Store results
  connectance_by_sierra <- rbind(connectance_by_sierra, data.frame(
    Sierra = sierra,
    area = sierra_area[1],
    realized_interactions = realized,
    potential_simple = potential_simple,
    potential_temporal = potential_temporal,
    connectance_simple = connectance_simple,
    connectance_temporal = connectance_temporal,
    stringsAsFactors = FALSE
  ))
}

# Display results
cat("\n=== Results by Sierra ===\n")
print(connectance_by_sierra)

# Add Difuntito flag
connectance_by_sierra$is_difuntito <- connectance_by_sierra$Sierra == 'Difuntito'

# Exclude Difuntito for statistical analysis
connectance_analysis <- connectance_by_sierra %>% filter(Sierra != 'Difuntito')

# Statistical models
lm_simple <- lm(connectance_simple ~ log10(area), data = connectance_analysis)
lm_temporal <- lm(connectance_temporal ~ log10(area), data = connectance_analysis)

cat("\n=== STATISTICS: Connectance vs Area ===\n")
cat("\nSimple co-occurrence:\n")
cat("Slope:", round(coef(lm_simple)[2], 4), "\n")
cat("R²:", round(summary(lm_simple)$r.squared, 3), "\n")
cat("p-value:", format.pval(summary(lm_simple)$coefficients[2,4], digits = 3), "\n")

cat("\nTemporal co-occurrence:\n")
cat("Slope:", round(coef(lm_temporal)[2], 4), "\n")
cat("R²:", round(summary(lm_temporal)$r.squared, 3), "\n")
cat("p-value:", format.pval(summary(lm_temporal)$coefficients[2,4], digits = 3), "\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\n=== GENERATING FIGURES ===\n")

# Figure: Metacommunity bar plot
metacommunity_data <- data.frame(
  Category = c("Realized", "Potential\n(Simple)", "Potential\n(Temporal)"),
  Value = c(n_realized, n_potential_simple, n_potential_temporal),
  Connectance = c(NA, connectance_simple_meta, connectance_temporal_meta)
)

metacommunity_data$Category <- factor(metacommunity_data$Category, 
                                      levels = c("Realized", "Potential\n(Simple)", "Potential\n(Temporal)"))

p_meta <- ggplot(metacommunity_data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = Value), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Realized" = "grey30", 
                               "Potential\n(Simple)" = "#A8DADC", 
                               "Potential\n(Temporal)" = "#457B9D")) +
  labs(x = "", 
       y = "Number of Interactions",
       title = "Metacommunity: Realized vs Potential Interactions") +
  annotate("text", x = 2, y = n_potential_simple * 0.5, 
           label = sprintf("Connectance = %.3f", connectance_simple_meta),
           size = 4, fontface = "italic") +
  annotate("text", x = 3, y = n_potential_temporal * 0.5,
           label = sprintf("Connectance = %.3f", connectance_temporal_meta),
           size = 4, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# Figure: Connectance vs Area (two panels)
p_area_simple <- ggplot(connectance_by_sierra, aes(x = area, y = connectance_simple)) +
  geom_point(data = connectance_by_sierra[!connectance_by_sierra$is_difuntito,],
             size = 4, color = "#A8DADC", alpha = 0.7) +
  geom_point(data = connectance_by_sierra[connectance_by_sierra$is_difuntito,],
             size = 5, color = "red", shape = 17) +
  geom_smooth(data = connectance_analysis, method = "lm",
              color = "#457B9D", fill = "#A8DADC", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area (ha)", 
       y = "Connectance",
       title = "A. Simple Co-occurrence") +
  annotate("text", x = min(connectance_analysis$area) * 1.2, 
           y = max(connectance_analysis$connectance_simple) * 0.95,
           label = sprintf("Slope = %.4f\nR² = %.3f\np = %s",
                           coef(lm_simple)[2],
                           summary(lm_simple)$r.squared,
                           format.pval(summary(lm_simple)$coefficients[2,4], digits = 3)),
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

p_area_temporal <- ggplot(connectance_by_sierra, aes(x = area, y = connectance_temporal)) +
  geom_point(data = connectance_by_sierra[!connectance_by_sierra$is_difuntito,],
             size = 4, color = "#A8DADC", alpha = 0.7) +
  geom_point(data = connectance_by_sierra[connectance_by_sierra$is_difuntito,],
             size = 5, color = "red", shape = 17) +
  geom_smooth(data = connectance_analysis, method = "lm",
              color = "#457B9D", fill = "#A8DADC", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area (ha)", 
       y = "Connectance",
       title = "B. Temporal Co-occurrence") +
  annotate("text", x = min(connectance_analysis$area) * 1.2, 
           y = max(connectance_analysis$connectance_temporal) * 0.95,
           label = sprintf("Slope = %.4f\nR² = %.3f\np = %s",
                           coef(lm_temporal)[2],
                           summary(lm_temporal)$r.squared,
                           format.pval(summary(lm_temporal)$coefficients[2,4], digits = 3)),
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# Combine area panels
fig_connectance_area <- grid.arrange(
  p_area_simple, p_area_temporal,
  ncol = 2,
  top = textGrob("Connectance vs Area", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save figures
ggsave("figure_connectance_metacommunity.pdf", p_meta, 
       width = 8, height = 6, units = "in", dpi = 300)

ggsave("figure_connectance_vs_area.pdf", fig_connectance_area, 
       width = 12, height = 5, units = "in", dpi = 300)

ggsave("figure_connectance_metacommunity.png", p_meta, 
       width = 8, height = 6, units = "in", dpi = 300)

ggsave("figure_connectance_vs_area.png", fig_connectance_area, 
       width = 12, height = 5, units = "in", dpi = 300)

cat("\n✓ Figures saved:\n")
cat("  - figure_connectance_metacommunity.pdf/png\n")
cat("  - figure_connectance_vs_area.pdf/png\n")

# Save data tables
write.csv(metacommunity_data, "connectance_metacommunity_summary.csv", row.names = FALSE)
write.csv(connectance_by_sierra, "connectance_by_sierra.csv", row.names = FALSE)

cat("✓ Data tables saved:\n")
cat("  - connectance_metacommunity_summary.csv\n")
cat("  - connectance_by_sierra.csv\n")

cat("\n=== CONNECTANCE ANALYSIS COMPLETE ===\n")


# ============================================================================
# ADDITIONAL ANALYSES: Temporal filtering effect across area gradient
# ============================================================================

cat("\n=== TEMPORAL FILTERING ANALYSIS ===\n")

# Add new metrics to the connectance_by_sierra dataframe
connectance_by_sierra <- connectance_by_sierra %>%
  mutate(
    # Option 2: Ratio of temporal to simple connectance
    connectance_ratio = connectance_temporal / connectance_simple,
    
    # Option 3: How much does temporal overlap constrain the interaction space?
    potential_ratio = potential_temporal / potential_simple,
    potential_difference = potential_simple - potential_temporal,
    potential_reduction_pct = (1 - potential_ratio) * 100
  )

# Exclude Difuntito for analysis
filtering_analysis <- connectance_by_sierra %>% filter(Sierra != 'Difuntito')

# Statistical models
# Option 2: Connectance ratio vs area
lm_connectance_ratio <- lm(connectance_ratio ~ log10(area), data = filtering_analysis)

# Option 3 models
lm_potential_ratio <- lm(potential_ratio ~ log10(area), data = filtering_analysis)
lm_potential_diff <- lm(potential_difference ~ log10(area), data = filtering_analysis)
lm_potential_reduction <- lm(potential_reduction_pct ~ log10(area), data = filtering_analysis)

cat("\n=== STATISTICS ===\n")

cat("\nOption 2: Connectance ratio (temporal/simple) vs Area\n")
cat("Slope:", round(coef(lm_connectance_ratio)[2], 4), "\n")
cat("R²:", round(summary(lm_connectance_ratio)$r.squared, 3), "\n")
cat("p-value:", format.pval(summary(lm_connectance_ratio)$coefficients[2,4], digits = 3), "\n")

cat("\nOption 3a: Potential ratio (temporal/simple) vs Area\n")
cat("Slope:", round(coef(lm_potential_ratio)[2], 4), "\n")
cat("R²:", round(summary(lm_potential_ratio)$r.squared, 3), "\n")
cat("p-value:", format.pval(summary(lm_potential_ratio)$coefficients[2,4], digits = 3), "\n")

cat("\nOption 3b: Potential difference (simple - temporal) vs Area\n")
cat("Slope:", round(coef(lm_potential_diff)[2], 4), "\n")
cat("R²:", round(summary(lm_potential_diff)$r.squared, 3), "\n")
cat("p-value:", format.pval(summary(lm_potential_diff)$coefficients[2,4], digits = 3), "\n")

cat("\nOption 3c: Potential reduction (%) vs Area\n")
cat("Slope:", round(coef(lm_potential_reduction)[2], 4), "\n")
cat("R²:", round(summary(lm_potential_reduction)$r.squared, 3), "\n")
cat("p-value:", format.pval(summary(lm_potential_reduction)$coefficients[2,4], digits = 3), "\n")

# Display summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
summary_stats <- filtering_analysis %>%
  summarize(
    mean_connectance_ratio = mean(connectance_ratio),
    sd_connectance_ratio = sd(connectance_ratio),
    mean_potential_ratio = mean(potential_ratio),
    sd_potential_ratio = sd(potential_ratio),
    mean_reduction_pct = mean(potential_reduction_pct),
    sd_reduction_pct = sd(potential_reduction_pct)
  )
print(summary_stats)

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("\n=== GENERATING FIGURES ===\n")

# Option 2: Connectance ratio vs Area
p_ratio_connectance <- ggplot(connectance_by_sierra, aes(x = area, y = connectance_ratio)) +
  geom_point(data = connectance_by_sierra[!connectance_by_sierra$is_difuntito,],
             size = 4, color = "#457B9D", alpha = 0.7) +
  geom_point(data = connectance_by_sierra[connectance_by_sierra$is_difuntito,],
             size = 5, color = "red", shape = 17) +
  geom_smooth(data = filtering_analysis, method = "lm",
              color = "#1D3557", fill = "#457B9D", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area (ha)", 
       y = "Connectance Ratio\n(Temporal / Simple)",
       title = "Temporal vs Simple Connectance") +
  annotate("text", x = min(filtering_analysis$area) * 1.2, 
           y = max(filtering_analysis$connectance_ratio) * 0.95,
           label = sprintf("Slope = %.4f\nR² = %.3f\np = %s",
                           coef(lm_connectance_ratio)[2],
                           summary(lm_connectance_ratio)$r.squared,
                           format.pval(summary(lm_connectance_ratio)$coefficients[2,4], digits = 3)),
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

# Option 3a: Potential ratio vs Area
p_ratio_potential <- ggplot(connectance_by_sierra, aes(x = area, y = potential_ratio)) +
  geom_point(data = connectance_by_sierra[!connectance_by_sierra$is_difuntito,],
             size = 4, color = "#457B9D", alpha = 0.7) +
  geom_point(data = connectance_by_sierra[connectance_by_sierra$is_difuntito,],
             size = 5, color = "red", shape = 17) +
  geom_smooth(data = filtering_analysis, method = "lm",
              color = "#1D3557", fill = "#457B9D", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area (ha)", 
       y = "Proportion of Potential Interactions\nwith Temporal Overlap",
       title = "A. Temporal Filtering Strength") +
  annotate("text", x = min(filtering_analysis$area) * 1.2, 
           y = min(filtering_analysis$potential_ratio) * 1.05,
           label = sprintf("Slope = %.4f\nR² = %.3f\np = %s",
                           coef(lm_potential_ratio)[2],
                           summary(lm_potential_ratio)$r.squared,
                           format.pval(summary(lm_potential_ratio)$coefficients[2,4], digits = 3)),
           hjust = 0, vjust = 0, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# Option 3b: Potential reduction percentage vs Area
p_reduction_pct <- ggplot(connectance_by_sierra, aes(x = area, y = potential_reduction_pct)) +
  geom_point(data = connectance_by_sierra[!connectance_by_sierra$is_difuntito,],
             size = 4, color = "#457B9D", alpha = 0.7) +
  geom_point(data = connectance_by_sierra[connectance_by_sierra$is_difuntito,],
             size = 5, color = "red", shape = 17) +
  geom_smooth(data = filtering_analysis, method = "lm",
              color = "#1D3557", fill = "#457B9D", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Area (ha)", 
       y = "Interaction Space Reduced by\nTemporal Constraints (%)",
       title = "B. Temporal Constraint Effect") +
  annotate("text", x = max(filtering_analysis$area) * 0.8, 
           y = max(filtering_analysis$potential_reduction_pct) * 0.95,
           label = sprintf("Slope = %.4f\nR² = %.3f\np = %s",
                           coef(lm_potential_reduction)[2],
                           summary(lm_potential_reduction)$r.squared,
                           format.pval(summary(lm_potential_reduction)$coefficients[2,4], digits = 3)),
           hjust = 1, vjust = 1, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# Combine Option 3 panels
fig_temporal_filtering <- grid.arrange(
  p_ratio_potential, p_reduction_pct,
  ncol = 2,
  top = textGrob("Temporal Filtering Across Area Gradient", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save figures
ggsave("figure_connectance_ratio.pdf", p_ratio_connectance, 
       width = 8, height = 6, units = "in", dpi = 300)

ggsave("figure_temporal_filtering.pdf", fig_temporal_filtering, 
       width = 12, height = 5, units = "in", dpi = 300)

ggsave("figure_connectance_ratio.png", p_ratio_connectance, 
       width = 8, height = 6, units = "in", dpi = 300)

ggsave("figure_temporal_filtering.png", fig_temporal_filtering, 
       width = 12, height = 5, units = "in", dpi = 300)

cat("\n✓ Additional figures saved:\n")
cat("  - figure_connectance_ratio.pdf/png (Option 2)\n")
cat("  - figure_temporal_filtering.pdf/png (Option 3)\n")

# Save updated data table
write.csv(connectance_by_sierra, "connectance_by_sierra_complete.csv", row.names = FALSE)
cat("✓ Updated data table saved: connectance_by_sierra_complete.csv\n")

cat("\n=== TEMPORAL FILTERING ANALYSIS COMPLETE ===\n")