# get access to cavatica projects
library("sevenbridges")
library("tidyr")
(a <- Auth(profile_name = "cavatica", from = "file"))

# selected tumor samples with published fusion
count_more_than_10 <- read_tsv("../input/processed/TCGA_final_fusion_morethan_10_samples_driver_fusion_2018.tsv")
# select samples with no published fusion
unknown_fusion_set <- read_tsv("../input/processed/TCGA_unknown_fusion_morethan_10_samples_driver_fusion_2018.tsv")

# get cavatica project 
p<-a$project(id="gaonkark/d3b-0007-annofuse-dev")

# generate rnaseq tasks
genertate_rnaseq_task<-function(x){
  input_file <- p$file(name="bam",metadata = list(aliquot_id = x))
  if (!is.null(input_file)){
  input_sample <- input_file$meta()$sample_id 
  input_aliquot_id <- x
  input_platform <- input_file$meta()$platform
  star_out_sam_att <- paste0("ID:",input_sample," LB:", input_aliquot_id," PL:", input_platform," SM:",input_sample)
  tsk_name=paste0("RNAseq_wf_",input_file$meta()$aliquot_id)
  tsk <- p$task_add(name=tsk_name,
                    app = "gaonkark/d3b-0007-annofuse-dev/kids-first-drc-rna-seq-workflow/0",
                    inputs = list(input_type ="BAM",
                                  sample_name=input_sample,
                                  wf_strand_param="default",
                                  runThread=36,
                                  reads1=input_file,
                                  STAR_outSAMattrRGline=star_out_sam_att),
                    input_check = F)
  }
  #Run task; better to validate correct task generation for a subset of samples in your list. Then continue with other samples uncommenting the following tsk_TF_Wf$run()
  #Verify correct generation of task on website and run
  
  tsk$run()
}

# run for 2 sets of samples
lapply(count_more_than_10$Sample, function(x) genertate_rnaseq_task(x))
lapply(unknown_fusion_set$Sample, function(x) genertate_rnaseq_task(x))


