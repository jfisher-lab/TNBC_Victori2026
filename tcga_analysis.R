# TODO
# 4. Annotate deleterious mutations
#
# For each gene in network:
#
#   Count deleterious mutations where delet_sift is true
#
# If multiple mutations exist for the same gene:
#
#   Sum the count of deleterious mutations

# data loading and cleaning ----
library(tidyverse)
library(tidybox)
library(magrittr)
library(TCGAretriever)
conflicted::conflicts_prefer(dplyr::filter)
library(janitor)
library(survival)
library(scales)
library(survminer)
library(contsurvplot)
library(interactions)
library(jtools)

# from Matthew
surv_custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12)
    )
}

model_genes = read_csv("../tnbcLeap/model_gene_list.csv") %>%
  pull(name)

node_dict = read_csv("../tnbcLeap/model_node_dict.csv") %>%
  select(node, gene)

# downloaded from cbioportal on 6/8/24
mutations_raw = read_tsv("../patient_data/other_datasets/tcga_2018/raw_data/data_mutations.txt")

# get data from cbioportal
brca_case_lists = get_case_lists("brca_tcga")

brca_tcga_studies = get_cancer_studies() %>%
  filter(cancerTypeId == "brca" & str_detect(studyId, "tcga"))

# from there we want to look at the study with id brca_tcga_pan_can_atlas_2018,
# which is the same I downloaded above.
# But later check if we can use any other
study_id = "brca_tcga_pan_can_atlas_2018"

brca_profiles = get_genetic_profiles(study_id)

clinical_raw = read_tsv("../patient_data/other_datasets/tcga_2018/raw_data/data_clinical_patient.txt",
                        comment = "#")

# get relevant data
mutations = get_mutation_data(case_list_id = "brca_tcga_all",
                              gprofile_id = "brca_tcga_pan_can_atlas_2018_mutations",
                              glist = model_genes)

cna_raw = get_molecular_data(case_list_id = "brca_tcga_all",
                             gprofile_id = "brca_tcga_pan_can_atlas_2018_gistic",
                             glist = model_genes)

rna_raw = get_molecular_data(case_list_id = "brca_tcga_all",
                             gprofile_id = "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores",
                             glist = model_genes)

# for some reason the mutation data is not being retrieved, but it should be the
# same I have loaded above as mutations_raw

mutations = mutations_raw %>%
  filter(Hugo_Symbol %in% model_genes) %>%
  select(gene = "Hugo_Symbol", Consequence, Variant_Classification, Variant_Type,
         patient_id = "Tumor_Sample_Barcode", Matched_Norm_Sample_Barcode, BIOTYPE,
         COSMIC, IMPACT, polyphen = "PolyPhen", SIFT, all_effects) %>%
  clean_names() %>%
  mutate(patient_id = str_remove(patient_id, "-01$")) %>%
  mutate(delet_sift = str_detect(sift, "deleterious")) %>%
  filter(delet_sift)

cna = cna_raw %>%
  select(-c(entrezGeneId, type)) %>%
  rename(gene = "hugoGeneSymbol") %>%
  pivot_longer(!gene, names_to = "patient_id", values_to = "cn") %>%
  mutate(patient_id = str_remove(patient_id, "-01$")) %>%
  filter(cn!=0) %>%
  # convert cna to categorical
  mutate(cn_binary = case_when(
    cn > 0 ~ "amp",
    cn < 0 ~ "del"))

rna = rna_raw %>%
  select(-c(entrezGeneId, type)) %>%
  rename(gene = "hugoGeneSymbol") %>%
  pivot_longer(!gene, names_to = "patient_id", values_to = "mrna") %>%
  mutate(patient_id = str_remove(patient_id, "-01$"))

# Clinical selections ----
clinical = clinical_raw %>%
  clean_names() %>%
  select(patient_id, subtype, ajcc_pathologic_tumor_stage, radiation_therapy,
         os_status_str = os_status, os_months, dss_status, dss_months,
         dfs_status_str = dfs_status, dfs_months, pfs_status, pfs_months) %>%
  mutate(os_status = ifelse(str_starts(os_status_str, "0"), 0, 1)) %>%
  mutate(dfs_status = ifelse(str_starts(dfs_status_str, "0"), 0, 1),
         dfs_status = ifelse(is.na(dfs_status), os_status, dfs_status),
         dfs_months = ifelse(is.na(dfs_months), os_months, dfs_months))

treatment = read_tsv("../patient_data/other_datasets/tcga_2018/raw_data/data_timeline_treatment.txt",
                     comment = "#") %>%
  clean_names() %>%
  select(patient_id, treatment_type, measure_of_response)

clinical_relevant = inner_join(clinical, treatment) %>%
  filter(subtype == "BRCA_Basal",
         radiation_therapy == "No",
         treatment_type == "Chemotherapy")

clinical_all = inner_join(clinical, treatment) %>%
  filter(treatment_type == "Chemotherapy")

p53_mutants = mutations %>%
  filter(gene == "TP53", delet_sift) %>%
  pull(patient_id)

clinical_p53_mutants = clinical_all %>%
  filter(patient_id %in% p53_mutants)

clinical_basal = clinical_all %>%
  filter(subtype == "BRCA_Basal")

clinical_chemo = clinical_all %>%
  filter(radiation_therapy == "No")

clinical_sets = list("all_patients" = clinical_all, "chemo" = clinical_chemo,
                     "p53_mutants" = clinical_p53_mutants, "basal" = clinical_basal,
                     "chemo" = clinical_chemo, "relevant" = clinical_relevant)

length(unique(clinical_relevant$patient_id))


# Stratify by RNA biomarkers ----
#rna_biomarkers = read_csv("../patient_data/simulations/results/corr_biomarkers_powerset.csv") %>%
rna_biomarkers = read_csv("F:/OneDrive - University College London/patient_data/simulations/results/corr_biomarkers_powerset.csv") %>%
  filter(sim_treatment != "olaparib") %>%
  select(node, pearson) %>%
  distinct() %>%
  mutate(pearson = pearson * -1) %>%  # positive numbers should indicate better response
  mutate(type = ifelse(pearson >0, "response", "resistance")) %>%
  left_join(node_dict,relationship = "many-to-many") %>%
  summarise(pearson = max(abs(pearson)), .by=c(gene, type)) %>%
  mutate(pearson = ifelse(type == "resistance",
                          pearson * -1, pearson),
         .keep = "unused")

# survival curves
clinical_selection = clinical_relevant
#for(clinical_selection in clinical_sets){
  clinical_tag = names(keep(clinical_sets, ~ identical(.x, clinical_selection)))

   stratification_rna = rna %>%
    filter(abs(mrna) >= 1) %>%
    inner_join(rna_biomarkers, relationship = "many-to-many") %>%
    # filter good results
    filter((mrna * pearson) > 0) %>%  # positive if same sign
    select(patient_id) %>%
    group_by(patient_id) %>%
    mutate(n_good = n()) %>%
    ungroup() %>%
    right_join(clinical_selection, relationship = "many-to-many") %>%
    select(patient_id, n_good, dfs_months, dfs_status, os_months, os_status) %>%
    distinct() %>%
    mutate(n_good = ifelse(is.na(n_good), 0, n_good),
           predicted_response = ifelse(n_good > median(n_good), "good", "bad"),
           predicted_response = factor(predicted_response,
                                       levels = c("good", "bad")))


  survcurve = survfit(Surv(os_months) ~ predicted_response, data = stratification_rna)

  #png(filename = paste0("../patient_data/other_datasets/tcga_2018/plots/rna_tcga_surcurve_",
  #                      clinical_tag, ".png"), width = 500, height = 650)
  print(ggsurvplot(survcurve,
             pval = TRUE,
             risk.table = TRUE,
             risk.table.height = 0.15,
             surv.scale = "percent",
             ggtheme = surv_custom_theme(),
             title = "RNA",
             palette = c("#00B6ED","#E40087"),
             xlab = "Overall survival time (months)",
             legend.labs = c("Complete response", "No response"),
             legend.title = ""))
  #dev.off()
#}

# Stratify by mutations that resensitise patients to treatment ----
resensitising_muts = read_csv(
  "../screens/2024-06-14_mono_results_whole_arms_cnvmuts/truly_resensitising.csv") %>%
  distil(pcr == "Yes") %>%
  left_join(read_csv("../patient_data/clin_retro_data.csv"), by = join_by(patient == leap_id)) %>%
  distil(arm == "chemo") %>%
  select(node = muta) %>%
  distinct() %>%
  inner_join(node_dict) %>%
  select(gene)

# survival curves
for(clinical_selection in clinical_sets){
  clinical_tag = names(keep(clinical_sets, ~ identical(.x, clinical_selection)))

  stratification_muts = resensitising_muts %>%
    inner_join(mutations) %>%
    select(gene, patient_id) %>%
    distinct() %>%
    group_by(patient_id) %>%
    mutate(n_muts = n()) %>%
    ungroup() %>%
    right_join(clinical_selection, relationship = "many-to-many") %>%
    select(patient_id, n_muts, os_months, os_status, dfs_months, dfs_status) %>%
    distinct() %>%
    mutate(predicted_response = ifelse(n_muts >= 1 & !is.na(n_muts), "good", "bad"),
           n_muts = ifelse(is.na(n_muts), 0, n_muts),
           predicted_response = factor(predicted_response,
                                       levels = c("good", "bad")))

  survcurve = survfit(Surv(os_months) ~ predicted_response, data = stratification_muts)

  png(filename = paste0("../patient_data/other_datasets/tcga_2018/plots/mut_tcga_survcurve_",
                        clinical_tag, ".png"), width = 500, height = 650)
  print(ggsurvplot(survcurve,
                   pval = TRUE,
                   risk.table = TRUE,
                   risk.table.height = 0.15,
                   surv.scale = "percent",
                   ggtheme = surv_custom_theme(),
                   title = "Resensitising mutations",
                   palette = c("#00B6ED","#E40087"),
                   xlab = "Overall survival time (months)",
                   legend.labs = c("Complete response", "No response"),
                   legend.title = ""))
  dev.off()
}

# Stratify by CNA minimal sets ----

# cna data legend:
# profile_description: Putative copy-number from GISTIC 2.0.
# Values: -2 = homozygous deletion; -1 = hemizygous deletion;
# 0 = neutral / no change; 1 = gain; 2 = high level amplification.

# read sets from powerset analysis and convert to format above
chemo_sets = read_csv("../patient_data/simulations/powerset/chemo_tree.csv") %>%
  filter(value <= -5) %>%
  mutate(predicted_response = "good") %>%
  select(-value) %>%
  rownames_to_column(var = "set_id") %>%
  pivot_longer(!c(set_id, predicted_response), names_to = "cna") %>%
  distil(value) %>%
  mutate(node = str_split_i(cna, ":",1),
         cna = as.numeric(str_split_i(cna, ":",2))) %>%
  # convert cna to categorical
  mutate(cna = case_when(
    cna > 1 ~ "amp",
    cna < 1 ~ "del")) %>%
  left_join(node_dict,relationship = "many-to-many") %>%
  mutate(gene = ifelse(node == "BRCA1", node, gene)) %>%
  mutate(n = ifelse(set_id == 12, 1, n()), .by = set_id) %>%
  mutate(cna = paste(gene, cna, sep = ":")) %>%
  select(set_id, predicted_response, cna, n)

stratification_cna = full_data %>%
  select(patient_id, gene, cn_binary) %>%
  mutate(cna = paste(gene, cn_binary, sep = ":")) %>%
  inner_join(chemo_sets, relationship = "many-to-many") %>%
  group_by(set_id, patient_id) %>%
  distinct() %>%
  mutate(n_match = n()) %>%
  # if number of relevant cna in one patient is
  # equal to a set size, then that means that the full set is found in that patient
  # and the patient is classified as responding according to
  # the predicted response from the chemo set
  filter(n_match >= n) %>%
  ungroup() %>%
  select(patient_id, predicted_response, set_id) %>%
  distinct() %>%
  mutate(n_sets = n(), .by = patient_id) %>%  # count how many sets are present in each patient
  right_join(clinical_selection, relationship = "many-to-many") %>%
  select(patient_id, predicted_response, dfs_months, dfs_status, os_months, os_status, n_sets) %>%
  distinct() %>%
  mutate(predicted_response = ifelse(is.na(predicted_response), "bad", predicted_response),
         n_sets = ifelse(is.na(n_sets), 0, n_sets),
         predicted_response = factor(predicted_response, levels = c("good", "bad")))

survcurve = survfit(Surv(os_months) ~ predicted_response, data = stratification_cna)

png(filename = paste0("../patient_data/other_datasets/cna_tcga_surcurve_",
                      clinical_tag, ".png"), width = 500, height = 650)
ggsurvplot(survcurve,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.15,
           surv.scale = "percent",
           ggtheme = surv_custom_theme(),
           title = "CNA",
           palette = c("#00B6ED","#E40087"),
           legend.labs = c("Complete response", "No response"),
           legend.title = "")
dev.off()

# make coxph model based on continuous response prediction from partition trees
chemo_sets_cont_response = read_csv("../patient_data/simulations/powerset/chemo_tree.csv") %>%
  mutate(predicted_response = abs(value), .keep = "unused") %>%
  rownames_to_column(var = "set_id") %>%
  pivot_longer(!c(set_id, predicted_response), names_to = "cna")%>%
  distil(value) %>%
  mutate(node = str_split_i(cna, ":",1),
         cna = as.numeric(str_split_i(cna, ":",2))) %>%
  # convert cna to categorical
  mutate(cna = case_when(
    cna > 1 ~ "amp",
    cna < 1 ~ "del")) %>%
  left_join(node_dict,relationship = "many-to-many") %>%
  mutate(gene = ifelse(node == "BRCA1", node, gene)) %>%
  mutate(n = n(), .by = set_id) %>%
  mutate(cna = paste(gene, cna, sep = ":")) %>%
  select(set_id, predicted_response, cna, n)

stratification_cna_cont = full_data %>%
  select(patient_id, gene, cn_binary) %>%
  mutate(cna = paste(gene, cn_binary, sep = ":")) %>%
  right_join(chemo_sets_cont_response, relationship = "many-to-many") %>%
  group_by(set_id, patient_id) %>%
  distinct() %>%
  mutate(n_match = n()) %>%
  filter(n_match >= n) %>%
  ungroup() %>%
  select(patient_id, predicted_response, set_id) %>%
  distinct() %>%
  right_join(clinical_selection, relationship = "many-to-many") %>%
  select(patient_id, predicted_response, dfs_months, dfs_status, os_months, os_status) %>%
  distinct() %>%
  mutate(predicted_response = ifelse(is.na(predicted_response), 0, predicted_response))

# Stratify using CNA partition tree ----
chemo_tree = readRDS("../patient_data/simulations/powerset/chemo_tree.rds")
chemo_cnas =  read_csv("../patient_data/simulations/powerset/chemo_tree.csv") %>%
  select(-value) %>%
  colnames() %>%
  tibble(cna = .) %>%
  separate_wider_delim(cna, delim = ":", names = c("node", "cn")) %>%
  mutate(cn_binary = case_when(
    cn > 1 ~ "amp",
    cn < 1 ~ "del")) %>%
  mutate(node = ifelse(node == "BRCA1", "BRCA1_full", node)) %>%
  left_join(node_dict, relationship = "many-to-many")


stratification_tree = full_data %>%
  select(patient_id, gene, cn_binary) %>%
  inner_join(chemo_cnas, relationship = "many-to-many") %>%
  mutate(cna = paste(node, cn, sep = ":"),
         value = TRUE) %>%
  select(patient_id, cna, value) %>%
  distinct() %>%
  complete(patient_id, cna, fill = list(value = FALSE)) %>%
  pivot_wider(names_from = cna, values_from = value) %>%
  # these columns are missing from all cna sets, but need to be present
  # for the next step. They were completed manually as errors were thrown
  mutate(`p27:0.5` = FALSE,
         `PTEN:0` = FALSE,
         `PTEN:0.5` = FALSE,
         `pRB:0.5` = FALSE,
         `pRB:0` = FALSE,
         `p53:0.5` = FALSE,
         `p16:0.5` = FALSE,
         `p16:0` = FALSE,
         `CycD:2` = FALSE,
         `ARF:0.5` = FALSE,
         `CDKN2A:0.5` = FALSE,
         `EGFR:4` = FALSE,
         `PI3K:2` = FALSE,
         `cMyc:2` = FALSE,
         `cMyc:3` = FALSE)

prediction = tibble(patient_id = stratification_tree$patient_id,
                    prediction =
                      predict(chemo_tree, newdata = stratification_tree, type = "class")) %>%
  mutate(prediction = abs(as.numeric(as.character(prediction)))) %>% #weird but necessary to convert from factor
  left_join(clinical_selection) %>%
  select(patient_id, prediction, dfs_months, dfs_status, os_status, os_months) %>%
  distinct() %>%
  mutate(prediction_binary = factor(ifelse(prediction >= median(prediction),
                                           "good", "bad"),
                                    levels = c("good", "bad")))

survcurve = survfit(Surv(dfs_months) ~ prediction_binary, data = prediction)

png(filename = paste0("../patient_data/other_datasets/cnatree_tcga_surcurve_",
                      clinical_tag, ".png"), width = 500, height = 650)
ggsurvplot(survcurve,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.15,
           surv.scale = "percent",
           ggtheme = surv_custom_theme(),
           title = "CNA partition tree",
           palette = c("#00B6ED","#E40087"),
           legend.labs = c("Complete response", "No response"),
           legend.title = "")
dev.off()

# Stratify by both CNA and RNA ----
stratification_cnarna = inner_join(stratification_cna, stratification_rna,
                                   by = join_by(patient_id, os_months, dfs_months,
                                                os_status, dfs_status),
                                   suffix = c("_cna", "_rna"))

survcurve = survfit(Surv(os_months) ~ predicted_response_cna + predicted_response_rna,
                    data = stratification_cnarna)

png(filename = paste0("../patient_data/other_datasets/cnarna_tcga_survarea_",
                      clinical_tag, ".png"), width = 700, height = 1000)

ggsurvplot(survcurve,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.15,
           surv.scale = "percent",
           ggtheme = surv_custom_theme(),
           title = "CNA + RNA",
           legend.title = "")
dev.off()
