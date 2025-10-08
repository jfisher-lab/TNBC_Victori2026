library(tidyverse)
library(magrittr)
library(patchwork)
library(ggpubr)
library(viridis)
library(khroma)
library(conflicted)
conflicts_prefer(dplyr::filter)

clinical = read_csv("../patient_data/clin_retro_data.csv")
druggable = read_csv("../tnbcLeap/druggable_nodes.csv") %>%
  select(node, gene, drug) %>%
  filter(!is.na(drug))

tag = "KO"
label = "Inhibition"
#tag = "KI"
#label = "Activation"

dir = "../screens/2024-06-14_mono_results_whole_arms_cnvmuts/"

excluded_genes = c("Unrepaired_teDSB", "sBRCAmut", "sBRCAcnv", "RepairNHEJ_teDSB",
                   "RepairNHEJ_seDSB", "PARTNER_olaparib", "PARTNER_Ep_Cy",
                   "PARTNER_Ca_Pa", "Paclitaxel", "Olaparib", "HR_teDSB",
                   "HR_seDSB", "gBRCAmut", "gBRCAcnv", "G2M_checkpoint",
                   "Epirubicin", "DNAreplicated", "Cyclophosphamide",
                   "CTL_killing", "CTL_activity", "BRCAmut", "ATMi")

patient_list = c("LEAP-008",
                    "LEAP-010",
                    "LEAP-014",
                    "LEAP-016",
                    "LEAP-017",
                    "LEAP-022",
                    "LEAP-023",
                    "LEAP-027",
                    "LEAP-032",
                    "LEAP-034",
                    "LEAP-037",
                    "LEAP-041",
                    "LEAP-045")
# Data processing ----
# separate background into patient and treatment
data = read_csv(paste0(dir, "node_integrated_results.csv")) %>%
  select(background, muta, leva, mean, pheno = node) %>%
  filter(!muta %in% excluded_genes) %>%
  mutate(patient = paste0("LEAP-", str_extract(background, "\\d{3}"))) %>%
  filter(patient %in% patient_list_v1) %>%
  mutate(treatment = str_remove(background, "leap_\\d{3}"),
         treatment = ifelse(treatment == "", "control", treatment),
         treatment = str_remove(treatment, "^_")) %>%
  #filter(patient %in% read_csv("../patient_data/simulations/patient_list_v2.3.csv")$background) %>%
  mutate(type = ifelse(muta == "baseline", "baseline", "mut")) %>%
  mutate(muta = ifelse(muta == "BRCA1_full", "BRCA1", muta)) %>%
  distil(((leva != 0 & tag == "KI")| leva == 0 & tag == "KO") |
           is.na(leva))

# compute difference between control and treatment
baseline_control = data %>%
  filter(type == "baseline",
         treatment == "control") %>%
  select(patient, baseline_control = mean, pheno)

data_diff = data %>%
  mutate(treatment = ifelse(treatment == "control", "control", "treatment")) %>%
  select(-background) %>%
  distinct() %>%
  distil(treatment != "control") %>%
  left_join(baseline_control) %>%
  # diff should be  (KO + chemo) - (no treatment)
  mutate(diff = mean - baseline_control, .keep = "unused")

# classify predictions from baseline and keep correct ones
predictions = data_diff %>%
  filter(type == "baseline") %>%
  left_join(clinical, by = join_by(patient == leap_id)) %>%
  select(pheno, patient, diff, pcr) %>%
  pivot_wider(names_from = pheno, values_from = diff) %>%
  mutate(prediction = CellDeath >=3 & Proliferation <=-1,
         correct_prediction = prediction == pcr)

predictions_filtered = predictions %>%
  filter(correct_prediction & ! pcr)

# prepare baseline data for plots
baseline = data_diff %>%
  filter(type == "baseline") %>%
  select(pheno, patient, diff) %>%
  pivot_wider(names_from = pheno, values_from = diff) %>%
  rename(baseline_death = "CellDeath", baseline_prolif = "Proliferation")

# Keep only correctly predicted patients
# and perturbations that produce a difference to the baseline
data_filtered = data_diff %>%
  filter(patient %in% predictions_filtered$patient) %>%
  left_join(baseline)

# Keep only targets that improve response
resensitising_targets = data_filtered %>%
  filter((pheno == "CellDeath" & diff >= 1) |
           (pheno == "Proliferation" & diff <= -1) |
           type == "baseline")

# classify into resensitising and improvement, add druggable information
truly_resensitising = resensitising_targets %>%
  pivot_wider(names_from = pheno, values_from = diff, values_fill = 0) %>% # assuming 0 diff won't meet any cutoff
  mutate(pcr = ifelse(((CellDeath >= 3) |
                         baseline_death >= 3 ) &
                        ((Proliferation <= -1) |
                           baseline_prolif <= -1), "Yes", "Improvement"),
         pcr = ifelse(pcr=="Improvement" & ((CellDeath - baseline_death < 1) &
                                              (Proliferation - baseline_prolif > -1)),
                      "No improvement", pcr)) %>%
  select(muta, patient, pcr, type) %>%
  distinct() %>%
  left_join(druggable,
            by = join_by(muta == node),
            relationship = "many-to-many") %>%
  select(muta, patient, pcr, drug) %>%
  mutate(Druggable = ifelse(is.na(drug), "Not druggable", "Druggable"), .keep = "unused") %>%
  mutate(Druggable = factor(Druggable, levels = c("Druggable", "Not druggable"))) %T>%
  write_csv(paste0(dir, '/truly_resensitising.csv'))

# Summary heatmap ----
if(tag == "KO"){ # only interested in producing summary plots for KO for the time being
  # Prepare labels for plot
  hm_targets_data = truly_resensitising %>%
    mutate(pcr = case_when(
      pcr == "Yes" ~ "Predicted pCR",
      pcr == "Improvement" ~ "Predicted improvement",
      .default = "No improvement")) %>%
    filter(any(pcr != "No improvement"), .by = muta) %>%
    # for v1
    mutate(combination_needed = ifelse((patient == "LEAP-010" & (muta == "mTORC1" | muta == "wee1")) |
                                         (patient == "LEAP-027" & (muta == "ATM" | muta == "4EBP1")),
                                       # for v2
                                       # mutate(combination_needed = ifelse(
                                       #         (patient %in% c("LEAP-027", "LEAP-032", "LEAP-034") & (muta == "ATM" | muta == "4EBP1")),
                                       "*", "")) %>% # highlighting patients that need a combination
    complete(muta, patient, fill = list(pcr = "No improvement")) %>%
    filter(muta != "baseline")

  order_drug = truly_resensitising %>%
    count(muta, pcr) %>%
    pivot_wider(names_from = pcr, values_from = n,values_fill = 0) %>%
    select(muta, `Predicted pCR` = "Yes",
           `Predicted improvement` = "Improvement")

  hm_targets = hm_targets_data %>%
    inner_join(order_drug) %>%
    ggplot(aes(x = patient, y = fct_reorder(muta, `Predicted pCR`,
                                            .desc = FALSE),
               fill = factor(
                 pcr, levels = c("Predicted pCR",
                                 "Predicted improvement",
                                 "No improvement")))) +
    geom_tile() +
    geom_text(aes(label = combination_needed), color = "white", size = 9, nudge_y = -0.4) +
    scale_fill_manual(values = c("#2E008B","#E40087", "#E3E3E3")) +
    theme_pubr() +
    rotate_x_axis_labels() +
    theme(legend.text = element_text(size = 12),
          legend.direction = "vertical",
          legend.position = "top") +
    labs(y = "Inhibition target", x = "Patient",
         fill = "")

  hm_drugs = truly_resensitising %>%
    filter(muta %in% hm_targets_data$muta) %>%
    mutate(t = "druggable") %>%
    filter(muta %in% hm_targets_data$muta) %>%
    inner_join(order_drug) %>%
    ggplot(aes(x = t, y = fct_reorder(muta, `Predicted pCR`,
                                      .desc = FALSE),
               fill = Druggable)) +
    geom_tile(colour = "white") +
    scale_fill_manual(values = c("grey20", "grey85")) +
    ggpubr::theme_pubr() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.direction = "vertical",
          legend.position = "top") +
    labs(x = "")

  hm_targets + hm_drugs +
    plot_layout(widths = c(20,1))

  ggsave(paste0(dir, "summary_pcr.png"), width = 12, height = 12, limitsize = FALSE)

  # with background
  perts = read_csv("../patient_data/funmuts/all_perts_nonresponders_v1.csv") %>%
    right_join(hm_targets_data %>% select(patient), relationship = "many-to-many") %>% # ensure all and only patients in the heatmap are present
    distinct() %>%
    complete(patient, node) %>%
    filter(!is.na(node)) %>% # this and previous lines are needed to keep patient with no background
    mutate(node = ifelse(node == "BRCA1_full", "BRCA1", node)) %>%
    rename(Patient = "patient") %>%
    mutate(effect = factor(effect, levels = c("Loss (both)",
                                              "Loss (CNV)",
                                              "Loss (Mut.)",
                                              "Gain (Mut.)",
                                              "Gain (CNV)",
                                              "Gain (both)"))) %>%
    ggplot(aes(x = Patient, y = node, fill = effect)) +
    geom_point(shape = 21, colour = "white", stroke = 0, size = 8) +
    labs(y = "Patient genetic background", fill = "") +
    scale_fill_manual(values = c("Loss (both)" = "#364B9A",
                                 "Loss (CNV)" = "#777DA7",
                                 "Loss (Mut.)" = "#A5D2E5",
                                 "Gain (Mut.)" = "#FDC072",
                                 "Gain (CNV)" = "#FF6700",
                                 "Gain (both)" = "#A50026"),
                      na.value = "white", na.translate = FALSE) +
    theme_pubr() +
    rotate_x_axis_labels() +
    theme(legend.text = element_text(size = 12),
          legend.direction = "vertical",
          legend.position = "top")

  perts

  ggsave(paste0(dir, "summary_background.png"), width = 10, height = 6, limitsize = FALSE)

  guide_area() + hm_targets + perts + hm_drugs +
    plot_layout(axes = "collect", guides = "collect",
                design = c(area(t = 1, l = 1, b = 2, r = 21),
                           area(t = 3, l = 1, b = 6, r = 20),
                           area(t = 7, l = 1, b = 8, r = 20),
                           area(t = 3, l = 21, b = 6, r = 21)))

  hm_targets + perts + hm_drugs +
    plot_layout(axes = "collect", guides = "collect",
                design = c(area(t = 1, l = 1, b = 4, r = 20),
                           area(t = 5, l = 1, b = 6, r = 20),
                           area(t = 1, l = 21, b = 4, r = 21)))

  ggsave(paste0(dir, "summary_pcr_background.png"), width = 10, height = 14, limitsize = FALSE)

  # only druggable ones, with backgrounds
  hm_targets_druggable = hm_targets_data %>%
    filter(muta %in% druggable$node) %>%
    inner_join(order_drug) %>%
    ggplot(aes(x = patient, y = fct_reorder(muta, `Predicted pCR`,
                                            .desc = FALSE),
               fill = factor(
                 pcr, levels = c("Predicted pCR",
                                 "Predicted improvement",
                                 "No improvement")))) +
    geom_tile() +
    geom_text(aes(label = combination_needed), color = "white", size = 5, nudge_y = -0.2) +
    scale_fill_manual(values = c("#2E008B","#E40087", "#E3E3E3")) +
    theme_pubr() +
    rotate_x_axis_labels() +
    theme(legend.text = element_text(size = 12),
          legend.direction = "vertical") +
    labs(y = "Inhibition target", x = "Patient",
         fill = "")

  hm_targets_druggable / perts + plot_layout(heights = c(3, 1),
                                             axes = "collect", guides = "collect")

  ggsave(paste0(dir, "summary_pcr_druggable_background.png"), width = 12, height = 12, limitsize = FALSE)
}

# Diff heatmaps ----
data_filtered %>%
  distil(pheno == "CellDeath") %>%
  select(muta, patient, type, diff) %>%
  ggplot(aes(x = patient, y = muta, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = diff, colour = diff), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_nightfall() +
  scale_colour_gradient(low = "black", high = "white",
                        na.value = NA, guide = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Patient",
       fill = "Change in cell death\nafter PARTNER treatment")

ggsave(paste0(dir, "reproduced_non_responders_celldeath_", tag, ".png"), width = 8, height = 16, limitsize = FALSE)

data_filtered %>%
  distil(pheno == "Proliferation") %>%
  select(muta, patient, type, diff) %>%
  complete(muta, patient, fill = list(diff = 0, type = "mut")) %>%
  ggplot(aes(x = patient, y = muta, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = diff, colour = diff), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), space = "free", scales = "free") +
  scale_fill_sunset() +
  scale_colour_gradient(low = "white", high = "black",
                        na.value = NA, guide = "none") +
  rotate_x_axis_labels() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Patient",
       fill = "Change in cell proliferation\nafter PARTNER treatment")

ggsave(paste0(dir, "reproduced_non_responders_proliferation_", tag, ".png"),
       width = 8, height = 10, limitsize = FALSE)

# Regular heatmaps ----
data %>%
  mutate(treatment = case_match(treatment,
                                "chemo" ~ "Chemotherapy",
                                "chemo_olaparib" ~ "Chemotherapy + Olabarib",
                                "control" ~ "Control"),
         treatment = factor(treatment, levels = c("Control",
                                                  "Chemotherapy",
                                                  "Chemotherapy + Olabarib"))) %>%
  distil(pheno == "CellDeath") %>%
  ggplot(aes(x = treatment, y = muta, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = mean, colour = mean), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), cols = vars(patient), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_nightfall() +
  scale_colour_gradient(low = "black", high = "white",
                        na.value = NA, guide = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Treatment",
       fill = "Predicted cell death")

ggsave(paste0(dir, "cell_death_complete_", tag, ".png"), width = 30, height = 30, limitsize = FALSE)

data %>%
  mutate(treatment = case_match(treatment,
                                "chemo" ~ "Chemotherapy",
                                "chemo_olaparib" ~ "Chemotherapy + Olabarib",
                                "control" ~ "Control"),
         treatment = factor(treatment, levels = c("Control",
                                                  "Chemotherapy",
                                                  "Chemotherapy + Olabarib"))) %>%
  distil(pheno == "Proliferation") %>%
  ggplot(aes(x = treatment, y = muta, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = mean), colour = "black", size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), cols = vars(patient), space = "free", scales = "free") +
  scale_fill_sunset(reverse = TRUE) +
  rotate_x_axis_labels() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Treatment",
       fill = "Predicted cell proliferation")

ggsave(paste0(dir, "proliferation_complete_", tag, ".png"), width = 30, height = 30, limitsize = FALSE)

# Boxplots ----
if(tag == "KO"){
  data_bp =
    data_diff %>%
    filter(patient %in% predictions_filtered$patient, # only correct predictions
           muta != "baseline") %>%
    left_join(baseline) %>%
    pivot_wider(names_from = pheno, values_from = diff) %>%
    mutate(`Cell death` = CellDeath - baseline_death,
           `Cell proliferation` = Proliferation - baseline_prolif) %>%
    filter(muta %in% druggable$node) %>%
    filter(any(`Cell death` > 0) | any(`Cell proliferation` < 0), .by = muta) %>%
    select(Target = "muta", patient,
           `Cell death`, `Cell proliferation`)

  order = data_bp %>%
    summarise(cd_mean = mean(`Cell death`), cp_mean = -mean(`Cell proliferation`), # change sign so positive = better response
              .by = Target)

  data_bp %>%
    pivot_longer(!c(Target, patient), names_to = "type",
                 values_to = "change",
                 values_drop_na = TRUE) %>%
    inner_join(order) %>%
    ggplot(aes(x = change,
               y = fct_reorder2(Target, cp_mean, cd_mean, .desc = FALSE), colour = type)) +
    geom_boxplot(position = position_identity()) +
    geom_jitter(width = 0, height = 0.25) +
    scale_colour_manual(values = c("#E40087", "#2E008B")) +
    facet_grid(cols = vars(type)) +
    theme_pubr() +
    theme(legend.direction = "horizontal") +
    labs(colour = "Phenotype", x = "Effect of additional target on cell phenotype",
         y = "Target")

  ggsave(paste0(dir, "drug_boxplot.png"), width = 8, height = 10, limitsize = FALSE)
}

# specific plots for jasmin ----
data_reduced = data %>%
  mutate(treatment = case_match(treatment,
                                "chemo" ~ "Chemotherapy",
                                "chemo_olaparib" ~ "Chemotherapy + Olabarib",
                                "control" ~ "Control"),
         treatment = factor(treatment, levels = c("Control",
                                                  "Chemotherapy",
                                                  "Chemotherapy + Olabarib"))) %>%
  distil(pheno == "CellDeath") %>%
  filter(muta %in% hm_targets_data$muta | muta == "baseline",
         patient %in% c("LEAP-014", "LEAP-022", "LEAP-034", "LEAP-037", "LEAP-045"))

hmr_cd = data_reduced %>%
  ggplot(aes(x = treatment, y = muta, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = mean), color = "black", size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), cols = vars(patient), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_gradient(low = "yellow", high = "#FF9FE5", na.value = NA) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 17),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 17)) +
  labs(y = "Inhibition target", x = "Treatment",
       fill = "Predicted cell death")

hmr_drugs = truly_resensitising %>%
  filter(muta %in% data_reduced$muta) %>%
  mutate(t = "druggable") %>%
  ggplot(aes(x = t, y = muta, fill = Druggable)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("grey20", "grey85")) +
  ggpubr::theme_pubr() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 17)) +
  labs(x = "")

hmr_cd + hmr_drugs +
  plot_layout(widths = c(20,1), guides = "collect")

ggsave(paste0(dir, "_cell_death_reduced.png"), width = 12, height = 12, limitsize = FALSE)

# reduced diff
cdr_data = resensitising_targets %>%
  distil(pheno == "CellDeath") %>%
  filter(muta %in% hm_targets_data$muta | muta == "baseline",
         patient %in% c("LEAP-014", "LEAP-027", "LEAP-034", "LEAP-041", "LEAP-045")) %>%
  mutate(diff = round(diff))

cdr = cdr_data %>%
  ggplot(aes(x = patient, y = muta, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = diff), color = "black", size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_gradient(low = "yellow", high = "#FF9FE5", na.value = NA) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = "Inhibition target", x = "Patient",
       fill = "Change in cell death\nafter PARTNER treatment")


hmr_drugs_diff = truly_resensitising %>%
  filter(muta %in% cdr_data$muta) %>%
  mutate(t = "druggable") %>%
  ggplot(aes(x = t, y = muta, fill = Druggable)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("grey20", "grey85")) +
  ggpubr::theme_pubr() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 17)) +
  labs(x = "")

cdr + hmr_drugs_diff +
  plot_layout(widths = c(20,1), guides = "collect")

ggsave(paste0(dir, "reduced_diff_celldeath.png"), width = 13, height = 11, limitsize = FALSE)

# 008 clones screen ----
dir = "../screens/2024-10-23_mono_results_clones_cnv/"
tag = "KO"
label = "Inhibition"
#tag = "KI"
#label = "Activation"

# separate background into patient and treatment
clonal_data = read_csv(paste0(dir, "node_integrated_results.csv")) %>%
  select(treatment = background, muta, leva, mean, pheno = node, clone) %>%
  filter(!muta %in% excluded_genes) %>%
  mutate(type = ifelse(muta == "baseline", "baseline", "mut")) %>%
  mutate(clone = case_match(clone,
                            "clone4" ~ "A",
                            "clone5" ~ "B",
                            "clone6" ~ "C",
                            "clonebulk" ~ "Bulk")) %>%
  distil(((leva != 0 & tag == "KI")| leva == 0 & tag == "KO") |
           is.na(leva)) %>%
  filter(treatment == "chemo_olaparib" |  # this is the treatment 008 got
           treatment == "control") %>%
  pivot_wider(names_from = "clone", values_from = "mean") %>%
  pivot_longer(A:Bulk, names_to = "clone", values_to = "mean") %>%
  mutate(clone = factor(clone, levels = c("A", "B", "C", "Bulk")))

# compute difference between control and treatment
clonal_data_diff = clonal_data %>%
  mutate(treatment = ifelse(treatment == "control", "control", "treatment")) %>%
  distinct() %>%
  pivot_wider(names_from = treatment, values_from = mean) %>%
  mutate(diff = treatment - control, .keep = "unused")

# prepare baseline data for plots
clonal_baseline = clonal_data_diff %>%
  filter(type == "baseline") %>%
  select(pheno, clone, baseline = diff) %>%
  pivot_wider(names_from = clone, values_from = baseline, names_prefix = "baseline")

# Keep only perturbations that produce a difference to the baseline in at least one clone
clonal_data_filtered = clonal_data_diff %>%
  left_join(clonal_baseline) %>%
  pivot_wider(names_from = clone, values_from = diff, names_prefix = "diff") %>%
  filter(!(diffA == baselineA & diffB == baselineB &
             diffC == baselineC & diffBulk == baselineBulk) | type == "baseline") %>%
  pivot_longer(starts_with("diff"), names_to = "clone", values_to = "diff",
               names_prefix = "diff") %>%
  select(-starts_with("baseline")) %>%
  mutate(clone = factor(clone, levels = c("A", "B", "C", "Bulk")))

## Clone diff heatmaps ----
clonal_data_filtered %>%
  distil(pheno == "CellDeath") %>%
  ggplot(aes(x = clone, y = muta, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = diff, colour = diff), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_nightfall() +
  scale_colour_gradient(low = "black", high = "white",
                        na.value = NA, guide = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Clone",
       fill = "Change in cell death\nafter PARTNER treatment")

ggsave(paste0(dir, "clones_celldeath_", tag, ".png"), width = 5, height = 13, limitsize = FALSE)

clonal_data_filtered %>%
  distil(pheno == "Proliferation") %>%
  ggplot(aes(x = clone, y = muta, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = diff, colour = diff), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_sunset() +
  scale_colour_gradient(low = "white", high = "black",
                        na.value = NA, guide = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Clone",
       fill = "Change in cell proliferation\nafter PARTNER treatment")

ggsave(paste0(dir, "clones_proliferation_", tag, ".png"), width = 5, height = 11, limitsize = FALSE)


## Clone regular heatmaps ----
clonal_data %>%
  mutate(treatment = case_match(treatment,
                                "chemo" ~ "Chemotherapy",
                                "chemo_olaparib" ~ "Chemotherapy + Olabarib",
                                "control" ~ "Control"),
         treatment = factor(treatment, levels = c("Control",
                                                  "Chemotherapy",
                                                  "Chemotherapy + Olabarib"))) %>%
  distil(pheno == "CellDeath") %>%
  ggplot(aes(x = treatment, y = muta, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = mean, , colour = mean), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), cols = vars(clone), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_nightfall() +
  scale_colour_gradient(low = "black", high = "white",
                        na.value = NA, guide = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Treatment",
       fill = "Predicted cell death")

ggsave(paste0(dir, "clones_cell_death_complete_", tag, ".png"), width = 8, height = 16, limitsize = FALSE)

clonal_data %>%
  mutate(treatment = case_match(treatment,
                                "chemo" ~ "Chemotherapy",
                                "chemo_olaparib" ~ "Chemotherapy + Olabarib",
                                "control" ~ "Control"),
         treatment = factor(treatment, levels = c("Control",
                                                  "Chemotherapy",
                                                  "Chemotherapy + Olabarib"))) %>%
  distil(pheno == "Proliferation") %>%
  ggplot(aes(x = treatment, y = muta, fill = mean)) +
  geom_tile() +
  geom_text(aes(label = mean, colour = mean), size = 4) +
  theme_pubr() +
  facet_grid(rows = vars(type), cols = vars(clone), space = "free", scales = "free") +
  rotate_x_axis_labels() +
  scale_fill_sunset(reverse = TRUE) +
  scale_colour_gradient(low = "black", high = "white",
                        na.value = NA, guide = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  ) +
  labs(y = paste0(label, " target"), x = "Treatment",
       fill = "Predicted cell proliferation")

ggsave(paste0(dir, "clones_proliferation_complete_", tag, ".png"), width = 8, height = 16, limitsize = FALSE)

## Clone summary heatmap ----
# baseline by pheno
clonal_pheno_baseline = clonal_data_diff %>%
  filter(type == "baseline") %>%
  select(pheno, clone, diff) %>%
  pivot_wider(names_from = pheno, values_from = diff) %>%
  rename(baseline_death = "CellDeath", baseline_prolif = "Proliferation")

# Keep only targets that resensitise to treatment
clonal_resensitising_targets = clonal_data_filtered %>%
  filter((pheno == "CellDeath" & diff >= 3) |
           (pheno == "Proliferation" & diff <= -1) |
           type == "baseline")

# classify into resensitising and improvement, add druggable information
clonal_truly_resensitising = clonal_resensitising_targets %>%
  left_join(clonal_pheno_baseline) %>%
  pivot_wider(names_from = pheno, values_from = diff, values_fill = 0) %>% # assuming 0 diff won't meet any cutoff
  mutate(pcr = ifelse(((CellDeath >= 3) |
                         baseline_death >= 3) &
                        ((Proliferation <= -1) |
                           baseline_prolif <= -1), "Yes", "Improvement"),
         pcr = ifelse(pcr=="Improvement" & ((CellDeath - baseline_death < 1) &
                                              (Proliferation - baseline_prolif > -1)),
                      "No improvement", pcr)) %>%
  select(muta, clone, pcr, type) %>%
  distinct() %>%
  left_join(druggable,
            by = join_by(muta == node),
            relationship = "many-to-many") %>%
  select(muta, clone, pcr, drug) %>%
  mutate(Druggable = ifelse(is.na(drug), "Not druggable", "Druggable"), .keep = "unused") %>%
  mutate(Druggable = factor(Druggable, levels = c("Druggable", "Not druggable"))) %>%
  mutate(clone = factor(clone, levels = c("A", "B", "C", "Bulk"))) %T>%
  write_csv(paste0(dir, '/truly_resensitising.csv'))

# Prepare labels for plot
hm_targets_data = clonal_truly_resensitising %>%
  mutate(pcr = case_when(
    pcr == "Yes" ~ "Predicted pCR",
    pcr == "Improvement" ~ "Predicted improvement",
    .default = "No improvement")) %>%
  filter(any(pcr != "No improvement"), .by = muta) %>%
  complete(muta, clone, fill = list(pcr = "No improvement")) %>%
  filter(muta != "baseline")

order_drug = clonal_truly_resensitising %>%
  count(muta, pcr) %>%
  pivot_wider(names_from = pcr, values_from = n,values_fill = 0) %>%
  select(muta, `Predicted pCR` = "Yes",
         `Predicted improvement` = "Improvement")

hm_targets = hm_targets_data %>%
  inner_join(order_drug) %>%
  ggplot(aes(x = clone, y = fct_reorder(muta, `Predicted pCR`,
                                        .desc = FALSE),
             fill = factor(
               pcr, levels = c("Predicted pCR",
                               "Predicted improvement",
                               "No improvement")))) +
  geom_tile() +
  scale_fill_manual(values = c("#2E008B","#E40087", "#E3E3E3")) +
  theme_pubr() +
  rotate_x_axis_labels() +
  theme(legend.text = element_text(size = 12),
        legend.direction = "vertical",
        legend.position = "top") +
  labs(y = "Inhibition target", x = "Clone",
       fill = "")

hm_drugs = clonal_truly_resensitising %>%
  filter(muta %in% hm_targets_data$muta) %>%
  mutate(t = "druggable") %>%
  filter(muta %in% hm_targets_data$muta) %>%
  inner_join(order_drug) %>%
  ggplot(aes(x = t, y = fct_reorder(muta, `Predicted pCR`,
                                    .desc = FALSE),
             fill = Druggable)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("grey20", "grey85")) +
  ggpubr::theme_pubr() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.direction = "vertical",
        legend.position = "top") +
  labs(x = "")

hm_targets + hm_drugs +
  plot_layout(widths = c(20,3))

ggsave(paste0(dir, "summary_pcr.png"), width = 5, height = 8, limitsize = FALSE)

# with background
perts = read_csv("../patient_data/clonal_data/cnv_background.csv") %>%
  mutate(Clone = case_match(Clone,
                            "4" ~ "A",
                            "5" ~ "B",
                            "6" ~ "C",
                            "bulk" ~ "Bulk"), .keep = "unused") %>%
  mutate(Clone = factor(Clone, levels = c("A", "B", "C", "Bulk"))) %>%
  ggplot(aes(x = Clone, y = node, fill = effect)) +
  geom_point(shape = 21, colour = "white", stroke = 0, size = 6) +
  labs(y = "Patient genetic background", fill = "") +
  scale_fill_manual(values = c(
    "CN loss" = "#777DA7",
    "CN gain" = "#FF6700"),
    na.value = "white", na.translate = FALSE) +
  theme_pubr() +
  rotate_x_axis_labels() +
  theme(legend.text = element_text(size = 12),
        legend.direction = "vertical",
        legend.position = "top")

perts

ggsave(paste0(dir, "summary_background.png"), width = 10, height = 6, limitsize = FALSE)

guide_area() + hm_targets + perts + hm_drugs +
  plot_layout(axes = "collect", guides = "collect",
              design = c(area(t = 1, l = 1, b = 2, r = 21),
                         area(t = 3, l = 1, b = 6, r = 20),
                         area(t = 7, l = 1, b = 8, r = 20),
                         area(t = 3, l = 21, b = 6, r = 21)))

hm_targets + perts + hm_drugs +
  plot_layout(axes = "collect", guides = "collect",
              design = c(area(t = 1, l = 1, b = 4, r = 20),
                         area(t = 5, l = 1, b = 8, r = 20),
                         area(t = 1, l = 21, b = 4, r = 21)))

ggsave(paste0(dir, "summary_pcr_background.png"), width = 8, height = 14, limitsize = FALSE)


## tmp perts ----
perts = read_csv("tmp_perts.csv")
