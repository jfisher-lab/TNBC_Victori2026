library(tidyverse)
library(magrittr)
library(conflicted)
library(NANSEN)
conflicts_prefer(dplyr::filter)

excluded = read_csv("../tnbcLeap/experiments/exclusions_mono.csv") %>%
  filter(! Node %in% c('CellDeath', 'Proliferation')) %>%
  deframe

clinical = read_csv("../patient_data/clin_retro_data.csv")

path = "../patient_data/simulations/results/cnvmut"

patient_regex = "LEAP-\\d{3}"

data = list.dirs(path = path,recursive=FALSE) %>%
  map(\(x) parse_biocheck_dir(x,
                              get_netw_variables("../tnbcLeap/tnbcLeap.json"))) %>%
  list_rbind() %T>%
  write_csv("../patient_data/simulations/results/parsed_cnvmut.csv")

proc = read_csv("../patient_data/simulations/results/parsed_cnvmut.csv") %>%
  filter(time == max(time), .by = filename)%>%
  rowwise() %>%
  mutate(mean = mean(c(lo, hi))) %>%
  ungroup() %>%
  mutate(leap_id = str_extract(filename, patient_regex),
         sim_treatment = str_remove(filename, ".json"),
         sim_treatment = str_remove(sim_treatment, paste0(patient_regex, "_")),
         sim_treatment = ifelse(str_detect(sim_treatment, patient_regex), "control", sim_treatment)) %>%
  select(name, mean, leap_id, sim_treatment) %>%
  filter(! name %in% excluded) %>%
  select(name, leap_id, sim_treatment, mean) %>%
  mutate(across(name:sim_treatment, as_factor)) %>%
  rename(node = "name")

prediction = readRDS("../patient_data/simulations/results/cnvmut.rds")$filtered_arms_data %>%
  select(phenotype, background, diff, pcr) %>%
  pivot_wider(names_from = phenotype, values_from = diff) %>%
  mutate(prediction = CellDeath >= 2 & Proliferation <= -1)

## differences by response ----
by_response = proc %>%
  pivot_wider(names_from = sim_treatment, values_from = mean) %>%
  mutate(across(chemo:chemo_olaparib, ~ .x - control)) %>%
  select(-control) %>%
  pivot_longer(chemo:chemo_olaparib, names_to = "sim_treatment", values_to = "activation_diff") %>%
  filter(!(is.na(activation_diff) | activation_diff == 0), .by = node) %>%
  mutate(node = fct_drop(node)) %>% # otherwise the 'complete' function completes factors not present in the data but still on the factor levels
  complete(leap_id, node, sim_treatment, fill = list(activation_diff = 0)) %>%
  filter(! node %in% c("CellDeath", "Proliferation")) %>%
  left_join(prediction, by = join_by(leap_id == background)) %>%
  rename(patient = "leap_id") %T>%
  write_csv("../patient_data/simulations/results/node_activity_by_response_cnvmut.csv")

by_response = read_csv("../patient_data/simulations/results/node_activity_by_response_cnvmut.csv")

# filter correct predictions
by_response_good = by_response %>%
  distil(prediction)

## correlation between node activation and survival diff, by node, in the correct predictions ----
corr = by_response_good %>%
  filter(! node %in% c("CellDeath", "Proliferation")) %>%
  select(node, activation_diff, CellDeath, Proliferation, sim_treatment) %>%
  distinct() %>%
  group_by(node, sim_treatment) %>%
  mutate(pvalue = cor.test(activation_diff, Proliferation - CellDeath)$p.value,
         pearson = cor.test(activation_diff, Proliferation - CellDeath)$estimate) %>%

  #mutate(across(c(CellDeath, Proliferation),
  #              list(pvalue = ~ cor.test(activation_diff, .x)$p.value,
  #                   pearson = ~ cor.test(activation_diff, .x)$estimate))) %>%
  ungroup() #%>%
  filter(pvalue <= 0.05) %T>%
  write_csv("../patient_data/simulations/results/corr_biomarkers_cnvmut.csv")