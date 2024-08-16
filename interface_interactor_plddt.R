library(readr)
library(tidyverse)
library(ggplot2)

# set maximum n of elements to be printed (useful for large lists)
options(max.print = 100)

# theme for plots
theme_set(
  theme_bw() +
    theme(
      text = element_text(colour = 'grey10'),
      panel.grid = element_blank(),
      plot.tag = element_text(
        size = 14,
        face = 'bold',
        colour = 'grey10'
      )
    )
)


# read data
af2_fusions <- read_tsv('./alphafold2_fusions_residue_level.tsv')
af3_fusions <- read_tsv('./alphafold3_fusions_residue_level.tsv')
af3_fusions_interaction <- read_tsv('./alphafold3_fusion_interaction_residue_level.tsv')
af2_fusions_interaction <- read_tsv('./alphafold2_fusion_interaction_residue_level.tsv')
af2_breakpoint <- read_csv('./breakpoint_all_residues_af2.csv')
af3_breakpoint <- read_csv('./breakpoint_all_residues_af3.csv')

library(stringr)
# you can come up with your own fusion name cleaning function
clean_fusion_names <- function(x) {
  str_replace_all(x, '_', ' ')
}


# 确保 fusion 列存在并执行 mutate 操作
af2_fusions <- af2_fusions %>%
  mutate(fusion = if_else(str_detect(fusion, 'ossifying'), 'CREBBP-BCORL1_oft', fusion))
af3_fusions <- af3_fusions %>%
  mutate(fusion = if_else(str_detect(fusion, 'ossifying'), 'CREBBP-BCORL1_oft', fusion))

# comparing AF2 to AF3 pLDDT at the breakpoint
inner_join(
  # average pLDDT on the breakpoint
  summarise(af2_breakpoint, AF2 = mean(plddt), .by = fusion),
  summarise(af3_breakpoint, AF3 = mean(plddt), .by = fusion),
  by = join_by(fusion)
) |> 
  # change long name to short
  mutate(fusion = if_else(str_detect(fusion, 'ossifying'), 
                          'CREBBP-BCORL1_oft', fusion)) |> 
  # pivot to long-format data
  pivot_longer(cols = where(is.numeric)) |> 
  # plot
  mutate(fusion = clean_fusion_names(fusion)) |>
  ggplot(aes(y = fct_reorder(fusion, value), x = value, fill = name)) +
  geom_point(size = 2.5, shape = 21, alpha = 0.9, stroke = 0.2) +
  coord_cartesian(xlim = c(0, 100)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8, 0.05),
        legend.background = element_blank()) +
  labs(y = NULL, x = 'pLDDT', fill = NULL)+
  scale_color_manual(values = c("AF2" = "red", "AF3" = "blue")) +
  scale_fill_manual(values = c("AF2" = "red", "AF3" = "blue"))



# USEFUL PLOTS BELOW -----------------------------------------------------------

# 1)
# correlation between pLDDT values at the breakpoint-af2
plddt_correlation_breakpoint_af2<- af2_fusions |>
  # -10/+10 from breakpoint
  filter(res < breakpoint + 10 & res > breakpoint - 10) |> 
  # create a relative location column
  mutate(location = if_else(res < breakpoint, '-10 AA', '+10 AA')) |> 
  # keep columns 
  select(fusion, gene, res, location, plddt) |> 
  # take mean per location for each fusion
  summarise(mean_plddt = mean(plddt), .by = c(fusion, location))

# calculate correlation
plddt_correlation_af2 <- plddt_correlation_breakpoint_af2 |> 
  pivot_wider(names_from = location,
              values_from = mean_plddt) |> 
#   mutate(abs_diff = abs(`-10 AA` - `+10 AA`))
# 
# # Filter points with the largest absolute differences
# threshold <- quantile(plddt_correlation$abs_diff, 0.9)  # Change the quantile value as needed
# low_correlation_points <- filter(plddt_correlation, abs_diff > threshold)
# 
# print(low_correlation_points)
  with(cor.test(`-10 AA`, `+10 AA`, method = 'spearman'))

print(plddt_correlation_af2)

# correlation between pLDDT values at the breakpoint-af3
plddt_correlation_breakpoint_af3<- af3_fusions |>
  # -10/+10 from breakpoint
  filter(res < breakpoint + 10 & res > breakpoint - 10) |> 
  # create a relative location column
  mutate(location = if_else(res < breakpoint, '-10 AA', '+10 AA')) |> 
  # keep columns 
  select(fusion, gene, res, location, plddt) |> 
  # take mean per location for each fusion
  summarise(mean_plddt = mean(plddt), .by = c(fusion, location))

# calculate correlation
plddt_correlation_af3 <- plddt_correlation_breakpoint_af3 |> 
  pivot_wider(names_from = location,
              values_from = mean_plddt) |> 
  #   mutate(abs_diff = abs(`-10 AA` - `+10 AA`))
  # 
  # # Filter points with the largest absolute differences
  # threshold <- quantile(plddt_correlation$abs_diff, 0.9)  # Change the quantile value as needed
  # low_correlation_points <- filter(plddt_correlation, abs_diff > threshold)
  # 
  # print(low_correlation_points)
  with(cor.test(`-10 AA`, `+10 AA`, method = 'spearman'))

print(plddt_correlation_af3)

# you can come up with your own fusion name cleaning function
clean_fusion_names <- function(x) {
  str_replace_all(x, '_', ' ')
}

# visualise
a <- plddt_correlation_breakpoint_af2 |>
  mutate(fusion = clean_fusion_names(fusion)) |> 
  ggplot(mapping = aes(
    y = fct_reorder(fusion, -mean_plddt),
    x = mean_plddt,
    fill = location
  )) +
  geom_point(
    size = 2.5,
    shape = 21,
    alpha = 0.9,
    stroke = 0.2
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.9),
    legend.background = element_blank()
  ) +
  labs(y = NULL, x = 'pLDDT', fill = NULL, tag = 'A') +
  annotate(
    geom = 'text',
    colour = 'grey10',
    label = paste('atop(italic(rho) ==', round(plddt_correlation_af2$estimate, 2), ', italic(p) == 1.35e-07)'),
    parse = TRUE,
    x = 50,
    y = 50
  )

b <- plddt_correlation_breakpoint_af3 |>
  mutate(fusion = clean_fusion_names(fusion)) |> 
  ggplot(mapping = aes(
    y = fct_reorder(fusion, -mean_plddt),
    x = mean_plddt,
    fill = location
  )) +
  geom_point(
    size = 2.5,
    shape = 21,
    alpha = 0.9,
    stroke = 0.2
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.9),
    legend.background = element_blank()
  ) +
  labs(y = NULL, x = 'pLDDT', fill = NULL, tag = 'B') +
  annotate(
    geom = 'text',
    colour = 'grey10',
    label = paste('atop(italic(rho) ==', round(plddt_correlation_af3$estimate, 2), ', italic(p) < 2.2e-16 )'),
    parse = TRUE,
    x = 50,
    y = 50
  )
library(patchwork)
a | b


# 2)
# extent of within-fusion interface
extent_of_fusion_interface_af2 <- af2_fusions |>
  mutate(bsa = sasa_split - sasa_ch) |>
  summarise(
    fusion_interface = sum(bsa),
    total_surface = sum(sasa_ch),
    .by = fusion
  ) |> 
  mutate(fusion_interface_to_surface_ratio = fusion_interface / total_surface) 

extent_of_fusion_interface_plddt70_af2 <- af2_fusions |>
  mutate(bsa = sasa_split - sasa_ch) |>
  filter(plddt > 70) |> 
  summarise(
    fusion_interface = sum(bsa),
    total_surface = sum(sasa_ch),
    .by = fusion
  ) |> 
  mutate(fusion_interface_to_surface_ratio = fusion_interface / total_surface) 


# extent of within-fusion interface
extent_of_fusion_interface_af3 <- af3_fusions |>
  mutate(bsa = sasa_split - sasa_ch) |>
  summarise(
    fusion_interface = sum(bsa),
    total_surface = sum(sasa_ch),
    .by = fusion
  ) |> 
  mutate(fusion_interface_to_surface_ratio = fusion_interface / total_surface) 

extent_of_fusion_interface_plddt70_af3 <- af3_fusions |>
  mutate(bsa = sasa_split - sasa_ch) |>
  filter(plddt > 70) |> 
  summarise(
    fusion_interface = sum(bsa),
    total_surface = sum(sasa_ch),
    .by = fusion
  ) |> 
  mutate(fusion_interface_to_surface_ratio = fusion_interface / total_surface) 

# visualise
a <- extent_of_fusion_interface_af2 |>
  mutate(fusion = clean_fusion_names(fusion)) |> 
  ggplot(
    mapping = aes(
      x = fusion_interface_to_surface_ratio,
      y = fct_reorder(fusion, fusion_interface_to_surface_ratio),
      fill = fusion_interface
    )
  ) +
  geom_col(width = 0.5) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_blank()
  ) +
  labs(x = 'Intra-fusion interfaces-to-surface ratio',
       y = NULL,
       fill = 'Intra-fusion interfaces size (Å²)',
       tag = 'A')

b <- extent_of_fusion_interface_plddt70_af2 |>
  mutate(fusion = clean_fusion_names(fusion)) |> 
  ggplot(
    mapping = aes(
      x = fusion_interface_to_surface_ratio,
      y = fct_reorder(fusion, fusion_interface_to_surface_ratio),
      fill = fusion_interface
    )
  ) +
  geom_col(width = 0.5) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_blank()
  ) +
  labs(x = 'Intra-fusion interfaces-to-surface ratio',
       y = NULL,
       fill = 'Intra-fusion interface size (Å²)',
       tag = 'B')

c <- extent_of_fusion_interface_af3 |>
  mutate(fusion = clean_fusion_names(fusion)) |> 
  ggplot(
    mapping = aes(
      x = fusion_interface_to_surface_ratio,
      y = fct_reorder(fusion, fusion_interface_to_surface_ratio),
      fill = fusion_interface
    )
  ) +
  geom_col(width = 0.5) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_blank()
  ) +
  labs(x = 'Intra-fusion interfaces-to-surface ratio',
       y = NULL,
       fill = 'Intra-fusion interfaces size (Å²)',
       tag = 'C')

d <- extent_of_fusion_interface_plddt70_af3 |>
  mutate(fusion = clean_fusion_names(fusion)) |> 
  ggplot(
    mapping = aes(
      x = fusion_interface_to_surface_ratio,
      y = fct_reorder(fusion, fusion_interface_to_surface_ratio),
      fill = fusion_interface
    )
  ) +
  geom_col(width = 0.5) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_blank()
  ) +
  labs(x = 'Intra-fusion interfaces-to-surface ratio',
       y = NULL,
       fill = 'Intra-fusion interface size (Å²)',
       tag = 'D')

library(patchwork)
a | b
c | d


# add fusion interaction-based data
a <- extent_of_fusion_interaction_interface_af2 <- af2_fusions_interaction |> 
  mutate(bsa = sasa_ch - sasa_cx) |>
  summarise(
    fusion_interaction_interface = sum(bsa),
    total_surface_interaction = sum(sasa_ch),
    .by = fusion
  ) |> 
  mutate(fusion_interaction_interface_to_surface_ratio = 
           fusion_interaction_interface / total_surface_interaction) 

fusion_interface_ratio_comparsion_af2 <- extent_of_fusion_interface_af2 |> 
  inner_join(extent_of_fusion_interaction_interface_af2, by = join_by(fusion)) |> 
  select(fusion, fusion_interface, fusion_interaction_interface) 

# correlation between fusion interface (as in the fusion),
# and the fusion interaction interface  (as in the predicted interaction)
fusion_interface_ratio_comparsion_af2 |>
  with(
    cor.test(
      fusion_interface,
      fusion_interaction_interface
    ),
    method = 'spearman'
  )

# add fusion interaction-based data
b <- extent_of_fusion_interaction_interface_af3 <- af3_fusions_interaction |> 
  mutate(bsa = sasa_ch - sasa_cx) |>
  summarise(
    fusion_interaction_interface = sum(bsa),
    total_surface_interaction = sum(sasa_ch),
    .by = fusion
  ) |> 
  mutate(fusion_interaction_interface_to_surface_ratio = 
           fusion_interaction_interface / total_surface_interaction) 

fusion_interface_ratio_comparsion_af3 <- extent_of_fusion_interface_af3 |> 
  inner_join(extent_of_fusion_interaction_interface_af3, by = join_by(fusion)) |> 
  select(fusion, fusion_interface, fusion_interaction_interface) 

# correlation between fusion interface (as in the fusion),
# and the fusion interaction interface (as in the predicted interaction)
fusion_interface_ratio_comparsion_af3 |>
  with(
    cor.test(
      fusion_interface,
      fusion_interaction_interface
    ),
    method = 'spearman'
  )

# plot
a <- fusion_interface_ratio_comparsion_af2 |> 
  pivot_longer(cols = c(fusion_interface, fusion_interaction_interface ),
               names_to = 'type',
               values_to = 'interface_size') |> 
  ggplot(
    mapping = aes(
      x = interface_size,
      y = fct_reorder(fusion, interface_size),
      fill = type
    )
  ) +
  geom_col(width = 0.5,
           position = position_dodge(width = 0.9)) +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_blank()
  ) +
  labs(x = 'Interface size (Å²)',
       y = NULL,
       fill = NULL,
       tag = 'A')+

  annotate(
    geom = 'text',
    colour = 'grey10',
    label = paste('t = 5.1056\n', 'p-value = 5.38e-06\n', 'cor = 0.5892804'),
    parse = TRUE,
    x = 2000,
    y = 20
  )

b <- fusion_interface_ratio_comparsion_af3 |> 
  pivot_longer(cols = c(fusion_interface, fusion_interaction_interface),
               names_to = 'type',
               values_to = 'interface_size') |> 
  ggplot(
    mapping = aes(
      x = interface_size,
      y = fct_reorder(fusion, interface_size),
      fill = type
    )
  ) +
  geom_col(width = 0.5,
           position = position_dodge(width = 0.9)) +
    theme(
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.3),
    legend.background = element_blank()
  ) +
  labs(x = 'Interface size (Å²)',
       y = NULL,
       fill = NULL,
       tag = 'B')+
  annotate(
  geom = 'text',
  colour = 'grey10',
  label = paste('t = 1.3868\n', 'p-value = 0.1718\n', 'cor = 0.1943361'),
  parse = TRUE,
  x = 2000,
  y = 20
)
library(patchwork)
a | b

A <- fusion_interface_ratio_comparsion_af2 |> 
  pivot_longer(cols = c(fusion_interface, fusion_interaction_interface),
               names_to = 'type',
               values_to = 'interface_size') |> 
  mutate(type = if_else(condition = type == 'fusion_interface', 
                        true = 'Interface in fusion',
                        false = 'Interface in complex')) |> 
  ggplot(aes(type, interface_size, fill = type)) +
  geom_boxplot(outliers = FALSE, staplewidth = 0.2, show.legend = FALSE) +
  labs(x = NULL, y = 'Interface size (Å²)', fill = NULL, tag='A')

B <- fusion_interface_ratio_comparsion_af3 |> 
  pivot_longer(cols = c(fusion_interface, fusion_interaction_interface),
               names_to = 'type',
               values_to = 'interface_size') |> 
  mutate(type = if_else(condition = type == 'fusion_interface', 
                        true = 'Interface in fusion',
                        false = 'Interface in complex')) |> 
  ggplot(aes(type, interface_size, fill = type)) +
  geom_boxplot(outliers = FALSE, staplewidth = 0.2, show.legend = FALSE) +
  labs(x = NULL, y = 'Interface size (Å²)', fill = NULL, tag='B')

library(patchwork)
A | B


# 3)
# quantify the extent of overlap
af3_fusions_bsa <- af3_fusions |> 
  mutate(fusion_bsa = sasa_split - sasa_ch) |> 
  select(fusion, uniprot_id, gene, res, fusion_bsa)

af3_fusions_interaction_bsa <- af3_fusions_interaction |> 
  mutate(fusion_interaction_bsa = sasa_ch - sasa_cx) |> 
  select(fusion, uniprot_id, gene, res, fusion_interaction_bsa)

af3_fusions_both_bsa <- af3_fusions_bsa |> 
  inner_join(af3_fusions_interaction_bsa,
             by = join_by(fusion, uniprot_id, gene, res))

af3_fusions_both_bsa_cor <- af3_fusions_both_bsa |> 
  group_by(fusion) |> 
  rstatix::cor_test(fusion_bsa, fusion_interaction_bsa) |> 
  drop_na(cor)

af3_fusions_both_bsa_cor |>
  ggplot(mapping = aes(
    y = fct_reorder(fusion, -cor),
    x = cor,
    fill = p
  )) +
  geom_col(colour = 'grey10') +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.1,
                 colour = 'grey10') +
  scale_fill_viridis_c(option = 'C') +
  labs(y = NULL, x = 'r', fill = 'P-value',
       title = 'Pearson correlation between interface in fusion and in interaction')
