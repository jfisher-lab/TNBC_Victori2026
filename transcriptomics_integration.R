library(tidyverse)
library(readxl)
library(viridis)
library(patchwork)
library(ggpubr)
library(ggvenn)
library(magrittr)
library(khroma)
conflicted::conflicts_prefer(dplyr::filter)

path = "../comparisons_rna-seq_mar24/"

network_nodes = read_csv("../tnbcLeap/model_gene_list.csv")
dict = read_csv("../tnbcLeap/model_node_dict.csv")

rna = list.files(path) %>%
  purrr::set_names() %>%
  map(\(x) read_excel(paste0(path, x))) %>%
  list_rbind(names_to = "comparison") %>%
  rename(n = 2) %>%
  mutate(comparison = str_remove(comparison, ".xlsx")) %>%
  mutate(abs_logFC = abs(avg_logFC)) %>%
  filter(! str_detect(comparison, "PBCP")) # exclude patient specific comparisons


rna_path = "../patient_data/transcriptomics/filtered_rna_v1.csv"

if(!file.exists(rna_path)){
  filtered_rna = rna %>%
    left_join(dict) %>%
    select(-c(max, id)) %>%
    filter(!is.na(node)) %T>%
    write_csv(rna_path)
}

filtered_rna = read_csv(rna_path)

# Add markers from visium
path = "../omics_analyses/tables/"

visium = list.files(path) %>%
  purrr::set_names() %>%
  map(\(x) read_csv(paste0(path, x))) %>%
  list_rbind(names_to = "comparison") %>%
  mutate(comparison = str_remove(comparison, ".csv")) %>%
  rename(avg_logFC = "avg_log2FC") %>%
  mutate(abs_logFC = abs(avg_logFC)) %>%
  filter(str_detect(comparison, "visium")) %>%
  filter(cluster == "complete") # no need to include both comparisons as they will be the same, there are only two groups

visium_path = "../patient_data/transcriptomics/filtered_visium_v1.csv"

if(!file.exists(visium_path)){
  filtered_visium = visium %>%
    left_join(dict) %>%
    select(-c(max, id)) %>%
    filter(!is.na(node)) %T>%
    write_csv(visium_path)
}

filtered_visium = read_csv(visium_path)

# Add markers from imc
imc_path = "../patient_data/transcriptomics/filtered_imc_v1.csv"

if(!file.exists(imc_path)){
ab_dict = read_csv("../omics_analyses/tables/ab_gene_dict_p4.csv")

imc = read_csv("../omics_analyses/tables/imc_p4_markers.csv") %>%
  rename(avg_logFC = "avg_log2FC",
         antibody = gene) %>%
  mutate(abs_logFC = abs(avg_logFC)) %>%
  filter(cluster == "complete") %>% # no need to include both comparisons as they will be the same, there are only two groups
  inner_join(ab_dict)

  filtered_imc = imc %>%
    left_join(dict) %>%
    select(-c(max, id)) %>%
    filter(!is.na(node)) %T>%
    write_csv(imc_path)
}

filtered_imc = read_csv(imc_path)

# Patient DE
# todo analyse
patient_de = read_xlsx("../Patient_Cancer_DEs_20240726.xlsx")

#for v1 we used corr_biomarkers_powerset.csv

corr_biomarkers = read_csv("../patient_data/simulations/results/corr_biomarkers_powerset.csv") %>%
  mutate(node = str_remove(node, "_full"),
         node = str_replace(node, "_", "\n"),
         node = str_replace(node, "BAKBAX", "BAK/BAX"),
         sim_treatment = case_match(sim_treatment,
                                    "CaPa" ~ "Carboplatin + Paclitaxel",
                                    "EpCy" ~ "Epirubicin + Cyclophosphamide",
                                    "olaparib" ~ "Olaparib")) %>%
  mutate(node = fct(node, levels = sort(unique(node)))) %>%
  mutate(pearson = -pearson, diff = -diff) # making diff represent intensity of response

corr_sign = corr_biomarkers %>%
  select(node, sim_treatment, pearson) %>%
  distinct()

shared = bind_rows(filtered_rna %>%
                     select(comparison, node, abs_logFC, avg_logFC),
                   filtered_visium %>%
                     select(comparison, node, abs_logFC, avg_logFC)) %>%
        bind_rows(filtered_imc %>%
                    mutate(comparison = "IMC Pre") %>%
                    select(comparison, node, abs_logFC, avg_logFC)) %>%
  mutate(node = str_remove(node, "_full"),
         node = str_replace(node, "_", "\n"),
         node = str_replace(node, "BAKBAX", "BAK/BAX")) %>%
  inner_join(corr_sign) %>%
  mutate(node = fct(node, levels = sort(unique(node)))) %>%
  mutate(comparison = case_match(comparison,
                                 "Retro_CancerCell_NonResponderGroup_Pre_vs_Post" ~
                                   "SN NR Pre vs. Post",
                                 "Retro_CancerCell_PreTreatmentGroup_CR_vs_NR" ~
                                   "SN CR vs. NR",
                                 "Retro_CompleteResponsePreTreatmentGroup_Cancer_vs_NonCancer" ~
                                   "SN CR cancer cells",
                                 "Retro_PartialResponsePreTreatmentGroup_Cancer_vs_NonCancer" ~
                                   "SN PR cancer cells",
                                 "visium_post_cancer_markers" ~
                                   "Sp. CR vs. NR post",
                                 "visium_post_chemo_cancer_markers" ~
                                   "Sp. CR vs. NR post (chemo)",
                                 "visium_post_chemo_olaparib_cancer_markers" ~
                                   "Sp. CR vs. NR post (chemo + ola.)",
                                 "visium_pre_cancer_markers" ~
                                   "Sp. CR vs. NR pre",
                                 "visium_pre_chemo_cancer_markers" ~
                                   "Sp. CR vs. NR pre (chemo)",
                                 "visium_pre_chemo_olaparib_cancer_markers" ~
                                   "Sp. CR vs. NR pre (chemo + ola.)",
                                 "IMC Pre" ~
                                   "IMC CR vs. NR pre")) %>%
  filter(abs_logFC > 0.5 | comparison == "IMC CR vs. NR pre", # IMC avg_logFC not big enough
         (avg_logFC * pearson > 0 & comparison != "NR Pre vs. Post") | # are they the same sign
           (avg_logFC * pearson < 0 & comparison == "NR Pre vs. Post"))

# horizontal biomarker plot
p_bio = ggplot(shared, aes(x = avg_logFC, y = comparison)) +
  geom_point(aes(colour = avg_logFC), size = 4) +
  geom_vline(xintercept = 0, linetype="dashed", color = "darkgrey") +
  facet_grid(cols = vars(node), space = "free") +
  scale_colour_sunset(reverse = FALSE) +
  coord_cartesian(clip = "off") +
  theme_pubr() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 24)) +
  labs(x = "", y = "Average logFC")

p_model = corr_biomarkers %>%
  filter(node %in% shared$node) %>%
  ggplot(aes(x = activation_diff, y = diff)) +
  facet_grid(cols = vars(node), rows = vars(sim_treatment)) +
  geom_point() +
  geom_smooth(aes(fill = pearson, colour = pearson), method = "lm") +
  scale_fill_gradient2(low = "#3f88c5", high = "#93032E", mid = "#ffba08") +
  scale_colour_gradient2(low = "#3f88c5", high = "#93032E", mid = "#ffba08") +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4)) +
  theme_pubr() +
  theme(legend.key.width= unit(2, 'cm'),
        axis.title = element_text(size = 24),
        panel.background = element_rect(fill = NA, color = "grey")) +
  labs(y = "Treatment efficacy",
       x = "Node activation change after treatment",
       fill = "Pearson correlation",
       colour = "Pearson correlation")

p_model / plot_spacer() / p_bio +
  plot_layout(heights = c(20, 1, 4))

ggsave("../patient_data/simulations/results/plots/biomarker_comparison_reduced.png",
       height = 14, width = 18)

# vertical biomarker plot
p_bio = ggplot(shared, aes(x = avg_logFC, y = comparison)) +
  geom_point(aes(fill = avg_logFC, colour = avg_logFC), size = 5) +
  geom_vline(xintercept = 0, linetype="dashed", color = "darkgrey") +
  facet_grid(rows = vars(node), space = "free") +
  scale_colour_sunset(reverse = FALSE) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title = element_text(size = 24)) +
  labs(x = "", y = "Average logFC")

p_model = ggplot(corr_biomarkers, aes(x = activation_diff, y = -diff)) +
  facet_grid(rows = vars(node), cols = vars(sim_treatment), scales = "free") +
  geom_point() +
  geom_smooth(aes(fill = pearson, colour = pearson), method = "lm") +
  scale_fill_gradient2(low = "#3f88c5", high = "#93032E", mid = "#ffba08") +
  scale_colour_gradient2(low = "#3f88c5", high = "#93032E", mid = "#ffba08") +
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4)) +
  theme_pubr() +
  theme(legend.key.width= unit(2, 'cm'),
        axis.title = element_text(size = 24),
        panel.background = element_rect(fill = NA, color = "grey")) +
  labs(y = "Tumour cell predicted net survival\n decrease after treatment",
       x = "Node activation change after treatment",
       fill = "Pearson correlation",
       colour = "Pearson correlation")

p_model + plot_spacer() + p_bio +
  plot_layout(widths = c(12, 1, 4))

ggsave("../patient_data/simulations/results/plots/biomarker_comparison_vert.png",
       height = 22, width = 12)

## venn ----
venn_corr = read_csv("../patient_data/biomarkers/results/corr_biomarkers.csv") %>%
  left_join(dict) %>%
  pull(gene) %>%
  unique()

venn_bio = biomarkers %>%
  filter(abs(avg_logFC) > 1) %>%
  filter(p_val_adj < 0.01)

venn_bio_filtered = filtered_biomarkers %>%
  filter(abs(avg_logFC) > 1) %>%
  filter(p_val_adj < 0.01)

venn_list = list(scRNAseq = unique(venn_bio$gene), model = unique(venn_corr))
ggvenn(venn_list, c("scRNAseq", "model"))
ggsave("../patient_data/biomarkers/results/plots/venn_all_vs_predicted.png")

venn_list = list(scRNAseq = unique(venn_bio_filtered$gene), model = unique(venn_corr))
ggvenn(venn_list, c("scRNAseq", "model"))
ggsave("../patient_data/biomarkers/results/plots/venn_present_vs_predicted.png")

venn_list = list(scRNAseq = unique(venn_bio$gene), model = unique(network_nodes$name))
ggvenn(venn_list, c("scRNAseq", "model"))
ggsave("../patient_data/biomarkers/results/plots/venn_all_vs_present.png")

## lists ----
model_biom = read_csv("../patient_data/biomarkers/results/corr_biomarkers.csv") %>%
  mutate(predicted_biomarker = TRUE) %>%
  full_join(dict) %>%
  replace_na(list(predicted_biomarker = FALSE)) %>%
  select(node, gene, predicted_biomarker) %>%
  distinct()

all_biom = biomarkers %>%
  select(comparison, cluster, gene, avg_logFC) %>%
  distinct() %>%
  full_join(model_biom) %>%
  mutate(present_in_model = ifelse(is.na(predicted_biomarker), FALSE, TRUE)) %>%
  distinct()

