library(tidyverse)
library(magrittr)
library(conflicted)
library(BMATools)
library(ggpubr)
library(patchwork)
library(tidybox)
library(rpart)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(ggpmisc)
library(viridis)
library(patchwork)
library(caret)
library(ggpattern)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)

# initial data ----
excluded = read_csv("../tnbcLeap/experiments/exclusions_mono.csv") %>%
  filter(! Node %in% c('CellDeath', 'Proliferation')) %>%
  deframe

clinical = read_csv("../patient_data/parsed_clinical_data.csv") %>%
  mutate(treatment = ifelse(arm == "arm1", "chemo", "chemo+olaparib"), .keep = "unused")

perts = read_csv("../patient_data/simulations/results/all_perturbations_full_cnv.csv") %>%
  distil(id=="cnv")

cnv_result_list = read_rds("../patient_data/simulations/results/cnv_results_list_220124.rds")

data = cnv_result_list$full_data %>%
  pivot_wider(names_from = phenotype, values_from = c(control, diff)) %>%
  mutate(control = control_Proliferation - control_CellDeath,
         diff = diff_Proliferation - diff_CellDeath, .keep = "unused") %>%
  select(-c(control, arm)) %>%
  filter(str_starts(condition, "arm")) %>%
  rename(patient = "background", survival_diff = "diff",
         sim_treatment = "condition") %>%
  inner_join(perts)

# decision trees with rpart ----
input = data %>%
  select(-response) %>%
  complete(node, patient, fill = list(activity = 1)) %>%
  select(node, patient, activity) %>%
  left_join(data %>% select(patient, sim_treatment, response, survival_diff),
            by = join_by(patient)) %>%
  rowwise() %>%
  mutate(perturbation = paste(node, activity, sep = "_"))

tree = rpart(
  survival_diff ~ sim_treatment + perturbation,
  data = input,
  #method = "class",
  minsplit = 2,
  minbucket = 1
)

fancyRpartPlot(tree)


# power sets simulations ----
source("cnv_into_network.R")

unique_perts = perts %>%
  mutate(perturbation = paste(node, activity, sep = "@_")) %>%
  pull(perturbation) %>%
  unique()

# maximum size of background
max = perts %>% summarise(n = n(), .by = patient) %>% pull(n) %>% max()

power_set = rje::powerSet(unique_perts, m = max)

all_perts = power_set %>%
  discard(\(x) length(x) == 0) %>%
  map(as_tibble) %>%
  list_rbind(names_to = "id") %T>%
  write_csv("../patient_data/simulations/all_perts_combinations.csv")

all_perts = read_csv("../patient_data/simulations/all_perts_combinations.csv")

comb_background = all_perts %>%
  separate_wider_delim(value, delim = "@_", names = c("node", "value")) %>%
  group_by(id) %>%
  filter(anyDuplicated(node) == 0) %T>%  # filter out groups that have a repeated node i.e.: a combination of A 0.5 and A 0.
  write_csv("../patient_data/simulations/all_combinations_backgrounds.csv")

# generate networks
comb_background = read_csv("F:/patient_simulations/all_combinations_backgrounds.csv")
network_path = "../tnbcLeap/tnbcLeap.json"
bma_path = 'C:\\"Program Files (x86)"\\BMA\\BioCheckConsole.exe'

cnv_to_networks(comb_background, network_path, bma_path,
                results_path = "F:/patient_simulations/powerset_networks/",
                cores = 14)

# check network files
check_cnv_network = Vectorize(function(id, subdir){
  file.exists(paste0("F:/patient_simulations/powerset_networks/",
                     subdir, "/", id, "_network.json"))
})

cnv_ids = comb_background %>%
  distinct(id) %>%
  mutate(subdir = as.character(floor(id/1000))) %>%
  write_csv("F:/patient_simulations/cnv_ids.csv")

## computationally expensive tasks ----
cnv_ids = read_csv("F:/patient_simulations/cnv_ids.csv")

# run backgrounds
registerDoParallel(14)
foreach(cuSubdir = unique(cnv_ids$subdir),
        .packages = c('magrittr', 'dplyr', 'tidybox',
                      'jsonlite', 'BMATools', 'purrr', 'stringr', 'readr')) %dopar% {
                        ids = cnv_ids %>%
                          filter(cuSubdir == subdir) %>%
                          distinct(id) %>%
                          mutate(id = as.character(id)) %>%
                          pull(id)

                        subdir_path = paste0("F:/patient_simulations/powerset_results/", cuSubdir, "/")
                        dir.create(subdir_path)

                        run_cnv_networks(networks_path = paste0("F:/patient_simulations/powerset_networks/",
                                                                cuSubdir, "/"),
                                         ids = ids,
                                         results_path = subdir_path,
                                         bma_path = 'C:\\"Program Files (x86)"\\BMA\\BioCheckConsole.exe')
                      }

# merge results
registerDoParallel(14)
foreach(cuSubdir = unique(cnv_ids$subdir),
        .packages = c('dplyr', 'purrr', 'readr')) %dopar% {
          subdir_path = paste0("F:/patient_simulations/powerset_results/", cuSubdir, "/")
          list.files(path = subdir_path,
                     pattern = ".csv") %>%
            map(\(x) read_csv(paste0(subdir_path, x), show_col_types = FALSE)) %>%
            list_rbind() %>%
            write_csv(paste0("F:/patient_simulations/powerset_results/merged/", cuSubdir, ".csv"))
        }

# final merge
path = "F:/patient_simulations/powerset_results/merged/"
list.files(path = path,
           pattern = ".csv") %>%
  map(\(x) read_csv(paste0(path, x), show_col_types = FALSE)) %>%
  list_rbind() %>%
  write_csv("F:/patient_simulations/results.csv")


# cleaning up results ----
powerset_results = read_csv("../patient_data/simulations/powerset/results.csv")
powerset_backgrounds = read_csv("../patient_data/simulations/powerset/all_combinations_backgrounds.csv")

background_wide = powerset_backgrounds %>%
  pivot_wider(names_from = node, values_from = value, values_fill = 1)

background_oneliner = powerset_backgrounds %>%
  summarise(pert_line = paste(node, value, sep = "=", collapse = ";"), .by = id)

powerset_data = powerset_results %>%
  #slice_head(n = 10^5) %>% # for quick testing
  select(-c(lo, hi, ...6)) %>%
  pivot_wider(names_from = phenotype, values_from = mean) %>%
  mutate(survival = Proliferation - CellDeath, .keep = "unused") %>%
  separate_wider_delim(background, "@", names = c("treatment", "id"),
                       too_few = "align_end") %>%
  mutate(treatment = str_remove(treatment, "_$"),
         id = as.numeric(id),
         treatment = ifelse(is.na(treatment), "control", treatment)) %>%
  pivot_wider(names_from = treatment, values_from = survival) %>%
  mutate(across(3:7, ~ .x - control)) %>%
  pivot_longer(!c(id, control),
               names_to = "condition", values_to = "diff") %>%
  left_join(background_wide) %>%
  left_join(background_oneliner) %>%
  relocate(pert_line, .after = diff) %T>%
  write_csv("../patient_data/simulations/powerset/full_results.csv")

powerset_data = read_csv("../patient_data/simulations/powerset/full_results.csv")

# boxplot per single gene ----
powerset_gene = powerset_data %>%
  select(id, condition, diff, pert_line) %>%
  separate_longer_delim(pert_line, delim = ";") %>%
  mutate(n = n(), .by = c(id, condition, diff)) %T>%
  write_csv("../patient_data/simulations/powerset/results_per_gene.csv")

powerset_gene = read_csv("../patient_data/simulations/powerset/results_per_gene.csv")

selected_condition = "chemo_olaparib"
for(selected_condition in unique(powerset_gene$condition)){
  powerset_gene %>%
    filter(condition == selected_condition) %>%
    ggplot(aes(x=factor(n), y = diff)) +
    geom_boxplot() +
    stat_poly_line(se = FALSE, method = "lm", aes(group = 1)) +
    stat_poly_eq(method = "lm", use_label(c("eq", "R2"), aes(group = 1))) +
    facet_wrap(vars(pert_line)) +
    theme_pubr()

  ggsave(paste0("../patient_data/simulations/powerset/", selected_condition, "_boxplot.png"), height = 10, width = 16)
}


# k-modes clustering ----
library(klaR)

binary_background = read_csv("../patient_data/simulations/powerset/all_combinations_backgrounds.csv") %>%
  mutate(presence = TRUE) %>%
  mutate(cna = paste(node, value, sep = ":"), .keep = "unused") %>%
  pivot_wider(names_from = cna, values_from = presence, values_fill = FALSE)


cluster_input = read_csv("../patient_data/simulations/powerset/full_results.csv") %>%
  select(id, condition, diff) %>%
  left_join(binary_background) %T>%
  write_csv("../patient_data/simulations/powerset/cluster_binary_data.csv")

cluster_input = read_csv("../patient_data/simulations/powerset/cluster_binary_data.csv")
selected_condition = "chemo_olaparib"

cluster_input_condition = cluster_input %>%
  distil(condition == selected_condition)

cx = cluster_input_condition %>% select(-id)

iter_max = 100
kmin = 1
kmax = 30

compute_withinss = function(model, k = NULL) {
  clusters = assign_clusters(model = model, k = k)

  centers = clusters |>
    group_by(cluster) |>
    summarize_all(mean) |>
    pivot_longer(cols = -cluster, names_to = "question", values_to = "cluster_mean")

  withinss = clusters |>
    pivot_longer(cols = -cluster, names_to = "question", values_to = "response") |>
    left_join(centers, by = c("cluster", "question")) |>
    summarize(k = max(cluster),
              withinss = sum((response - cluster_mean)^2)) |>
    mutate(model = class(model)[1])

  return(withinss)
}

cluster_results = tibble(k = kmin:kmax) %>%
  mutate(kclust = map(k, ~ kmodes(cx, modes = ., iter.max = iter_max)))


saveRDS(cluster_results, "../patient_data/simulations/powerset/cluster_results.rds")

cluster_results = readRDS("../patient_data/simulations/powerset/cluster_results.rds")

clustered_data = cluster_input_condition %>%
  select(id, diff) %>%
  mutate(cluster = factor(cluster_results$cluster))

summary_cluster = clustered_data %>%
  count(diff, cluster) %>%
  mutate(n = n/sum(n), .by = diff)

ggplot(clustered_data, aes(x = cluster, y = diff)) +
  geom_bin_2d() +
  theme_pubr() +
  theme(legend.key.width= unit(2, 'cm'))

ggsave("../patient_data/simulations/powerset/cluster_bin2d.png", height = 6, width = 7)

ggplot(clustered_data, aes(x = diff, fill = cluster)) +
  geom_bar(position = "dodge") +
  theme_pubr()

ggplot(summary_cluster, aes(x = diff, y = n, fill = cluster)) +
  geom_col(position = "dodge") +
  theme_pubr()

# powerset decision trees with rpart ----
library(ggdendro)

tree_input = read_csv("../patient_data/simulations/powerset/cluster_binary_data.csv")

selected_condition = "chemo_olaparib"
selected_condition = "chemo"
selected_condition = "olaparib"
selected_condition = "CaPa"
selected_condition = "EpCy"

conditions = c("chemo_olaparib" = "Chemotherapy + Olaparib",
               "chemo" = "Chemotherapy",
               "olaparib" = "Olaparib",
               "CaPa" = "Carboplatin + Paclitaxel",
               "EpCy" = "Epirubicin + Cyclophosphamide")

for(selected_condition in unique(tree_input$condition)){
  tree = tree_input %>%
    distil(condition == selected_condition) %>%
    select(-id) %>%
    mutate(diff = factor(diff)) %>%
    rpart(
      diff ~ .,
      data = .,
      method = "class",
      minsplit = 2,
      minbucket = 1
    )

  saveRDS(tree, paste0("../patient_data/simulations/powerset/", selected_condition, "_tree.rds"))
}

# tree plot
tree = readRDS(paste0("../patient_data/simulations/powerset/", selected_condition, "_tree.rds"))

for(selected_condition in unique(tree_input$condition)){
  height_upset = 7
  tree = readRDS(paste0("../patient_data/simulations/powerset/", selected_condition, "_tree.rds"))
  selected_input = tree_input %>% distil(condition == selected_condition) %>%
    mutate(diff = factor(diff))

  prediction = as_tibble(predict(tree, selected_input, type = "class"), rownames = "id")

  binary_input = selected_input %>%
    select(id, diff) %>%
    mutate(diff = as.numeric(levels(diff))[diff], # conversion from factor to numeric
           response = factor(diff <= -5))

  binary_prediction = prediction %>%
    select(id, value) %>%
    mutate(value = as.numeric(levels(value))[value], # conversion from factor to numeric
           response = factor(value <= -5))

  cm = confusionMatrix(binary_prediction$response, binary_input$response,
                       positive = "TRUE", mode = "everything")

  f = cm$byClass["F1"]

  dendro_data(tree) %>%
    modify_in("labels", \(x) x %>%
                mutate(yes = str_detect(label, ">")) %>%
                mutate(left = ifelse(yes, "Yes", "No"),
                       right = ifelse(yes, "No", "Yes")) %>%
                mutate(clean = str_split_i(label, "[<>]", 1)) %>%
                mutate(clean = str_replace_all(clean, "BRCA1_full", "BRCA1"),
                       clean = str_replace_all(clean, "FA_core_complex", "FA_core")) %>%
                mutate(label = paste(left, "-", clean, "-", right)) %>%
                select(x, y, label)) %>%
    {ggplot() +
        geom_segment(data = .$segments,
                     aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_label(data = .$labels, aes(x = x, y = y+0.005, label = label), size = 5) +
        geom_label(data = .$leaf_labels,
                   aes(x = x, y = y, label = abs(as.numeric(label)),
                       fill = abs(as.numeric(label)),
                       colour = abs(as.numeric(label))),
                   hjust = "center", vjust = "center", size = 8) +
        labs(fill = "Cell response") +
        theme_dendro() +
        scale_fill_gradient(high = "#2E008B", low = "#00B6ED") +
        scale_colour_gradient2(high = "lightgrey",mid = "lightgrey", low = "black",
                               midpoint = 4) +
        guides(colour = "none") +
        coord_cartesian(clip = "off") +
        ggtitle(paste0(conditions[selected_condition], ". F-score: ", round(f, 2)))}

  ggsave(paste0("../patient_data/simulations/powerset/",selected_condition, "_tree.png"), width = 14, height = 10)
  ggsave(paste0("../patient_data/simulations/powerset/",selected_condition, "_tree.svg"), width = 14, height = 10)

  #upset plot
  if(selected_condition == "chemo"){
    selected_condition = "chemo_validated"
    height_upset = 10
  }

  df = read_csv(paste0("../patient_data/simulations/powerset/", selected_condition, "_tree.csv")) %>%
    mutate(response = abs(value), .keep = "unused") %>%
    arrange(response) %>%
    mutate(index = factor(row_number())) %>%
    rowwise() %>%
    filter(if_any(contains(":"))) # filter out column with all CNVs = FALSE. All CNV columns will have ':' in the title

  y_text_size = 18

  bars = ggplot(df, aes(x=index, y = response, fill = response)) +
    geom_col() +
    theme_pubr() +
    scale_fill_gradient(high = "#2E008B", low = "#00B6ED") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size = y_text_size)) +
    labs(fill = "Cell response", y = "Cell response")

  if(selected_condition == "chemo_validated"){
    dumbbells = df %>%
      select(-c(12:18)) %>%
      pivot_longer(!c(response, index), names_to = "CNA") %>%
      distil(value) %>%
      ggplot(aes(x=index, y = CNA, colour = response)) +
      geom_point(size = 5) +
      geom_line(aes(group = index), linewidth = 1) +
      scale_colour_gradient(high = "#2E008B", low = "#00B6ED") +
      theme_pubr() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.x = element_line( size=.1, color="lightgrey"),
            legend.position = "none",
            axis.title.y = element_text(size = y_text_size)) +
      labs(colour = "Cell response")

    validation = df %>%
      select(12:20) %>%
      pivot_longer(!c(response, index), names_to = "condition", values_transform = as.character)%>%
      mutate(condition = case_match(condition,
                                    "pre_chemo_olaparib" ~ "Pre-treatment (chemo + olaparib)",
                                    "pre_chemo" ~ "Pre-treatment (chemo)",
                                    "pre" ~ "Pre-treatment",
                                    "post_chemo_olaparib" ~ "Post-treatment (chemo + olaparib)",
                                    "post_chemo" ~ "Post-treatment (chemo)",
                                    "post" ~ "Post-treatment",
             "imc_pre" ~ "Pre-treatment (IMC)")) %>%
      mutate(pattern = ifelse(value == "P", "p", "none")) %>%
      mutate(`Differentially expressed` =
               case_match(value,
                          "TRUE" ~ "Spatial transcriptomics",
                          "FALSE" ~ "No",
                          "T" ~ "Spatial proteomics",
                          "F" ~ "No",
                          "P" ~ "Spatial proteomics")) %>%
      ggplot(aes(x = index, y = condition,
                 fill = `Differentially expressed`,
                 pattern = pattern)) +
      geom_tile_pattern(pattern_color = "#2E008B",
                        pattern_density = 0.1,
                        pattern_spacing = 0.05,
                        pattern_key_scale_factor = 0.6) +
      scale_fill_manual(values = c(
        "Spatial transcriptomics" = "#E40087",
        "Spatial proteomics" = "#00B6ED",
        "No" = "grey90"),
        breaks = c("Spatial transcriptomics",
                   "Spatial proteomics")) +
      scale_pattern_manual(values = c("p" = "stripe", "none" = "none")) +
      guides(fill = guide_legend(override.aes = list(pattern = "none")),
             pattern = "none") +
      theme_pubr() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "bottom",
            axis.title.y = element_text(size = y_text_size)) +
      labs(y = "Condition")

    bars/dumbbells/validation + plot_layout(heights = c(2,2,1))
  } else{
    dumbbells = df %>%
      pivot_longer(!c(response, index), names_to = "CNA") %>%
      distil(value) %>%
      ggplot(aes(x=index, y = CNA, colour = response)) +
      geom_point(size = 5) +
      geom_line(aes(group = index), linewidth = 1) +
      theme_pubr() +
      scale_colour_gradient(high = "#2E008B", low = "#00B6ED") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.x = element_line( size=.1, color="grey90"),
            legend.position = "none",
            axis.title.y = element_text(size = y_text_size)) +
      labs(colour = "Cell response")

    bars/dumbbells
  }

  ggsave(
    paste0("../patient_data/simulations/powerset/",
           selected_condition, "_upset.png"),
         width = 11, height = height_upset)
}
