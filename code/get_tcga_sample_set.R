
get_tumor_samples <- TRUE
get_unknown_fusion_tumor_samples <- TRUE

if(get_tumor_samples){
  # selected tumor samples with published fusion
  # tab "Final fusion call set" from https://www.sciencedirect.com/science/article/pii/S2211124718303954
  read_truesets<-read_tsv("../input/raw/final_fusion_mmc2.txt",quote = "Final fusion call")
  
  # select 10 samples per Cancer group
  count <-read_truesets %>% dplyr::select("Cancer", "Sample") %>% 
    unique() %>%
    group_by(Cancer) %>%
    slice(tail(row_number(), 10))
  
  # keep only cancer types with atleast 10 samples
  cancer_groups_more_than10<-table(count$Cancer)[unname(table(count$Cancer)==10)]
  
  # keep samples with atleast 10 samples 
  count_more_than_10<-count[count$Cancer %in% 
                              names(cancer_groups_more_than10),]
  count_more_than_10 %>% write_tsv("../input/processed/TCGA_final_fusion_morethan_10_samples_driver_fusion_2018.tsv")
}

if(get_unknown_fusion_tumor_samples){
  # selected tumor samples with published fusion
  # tab "Final fusion call set" from https://www.sciencedirect.com/science/article/pii/S2211124718303954
  read_truesets<-read_tsv("../input/raw/final_fusion_mmc2.txt",quote = "Final fusion call")
  
  # tab "TCGA samples used in this study" from https://www.sciencedirect.com/science/article/pii/S2211124718303954
  read_allset<-read_tsv("../input/raw/all_tcga.txt")
  
  # select 10 samples per Cancer group
  select_cancer_typ<-read_truesets %>% 
    dplyr::group_by(Cancer) %>% 
    tally() %>% 
    arrange(n) %>% 
    slice(tail(row_number(),10))
  
  read_unknown_fusion_set<-read_allset %>% 
    # select Cancer which have samples in trueset
    dplyr::filter(Cancer %in% select_cancer_typ$Cancer) %>%
    # select samples with no published fusion
    dplyr::filter(!Sample %in% read_truesets$Sample) %>%
    group_by(Cancer) %>% 
    # select 5 
    slice(tail(row_number(), 5)) %>% 
    write_tsv("../input/processed/TCGA_unknown_fusion_morethan_10_samples_driver_fusion_2018.tsv")
  
}
