library("sevenbridges")
library("tidyverse")
(a <- Auth(profile_name = "cavatica", from = "file"))
a$user(username = "gaonkark")

# read fusion files manifest
manifest<-read_csv("../input/processed/1597413367787-manifest.csv") %>%
  dplyr::filter(!grepl("pdf",name))
manifest$sample_id<-unlist(lapply(strsplit(manifest$name,"[.]"), function(x) x[1]))

# read rsem file manifest
manifest_rsem<-read_csv("../input/processed/1597415171566-manifest.csv") %>%
  dplyr::select("name","id") 
manifest_rsem$sample_id<-unlist(lapply(strsplit(manifest_rsem$name,"[.]"), function(x) x[1]))
  
# select arriba output files
arriba_files<-manifest %>%
  dplyr::filter(grepl("arriba",name,perl = TRUE) ) %>%
  dplyr::select("sample_id","name","id")  

# select starfusion output files
starfusion_files<-manifest %>%
  dplyr::filter(grepl("STAR.fusion",name)) %>%
  dplyr::select("sample_id","name","id")  

# merge by sample_id to match arriba, starfusion and rsem files
fusion_files<-full_join(arriba_files,starfusion_files,
                        by="sample_id",
                        suffix = c(".arriba", ".starfusion")) %>%
  full_join(manifest_rsem,by="sample_id") %>%
  as.data.frame() %>% na.exclude()


run_format_fusion <-function(i=1){
  # run annoFuse 
  p<-a$project("d3b-0007-annofuse-dev")
  arriba_file<-p$file(id=fusion_files[i,"id.arriba"])
  starfusion_file<- p$file(id=fusion_files[i,"id.starfusion"]) 
  rsem_file<- p$file(id=fusion_files[i,"id"]) 
  app_id<-"gaonkark/d3b-0007-annofuse-dev/kfdrc-annofuse-wf"
  tsk<-a$project("d3b-0007-annofuse-dev")$task_add(
    name = paste0("annoFuse runs",fusion_files[i,"sample_id"]),
    app = app_id,
    inputs = list(arriba_output_file=arriba_file,
                  star_fusion_output_file=starfusion_file,
                  rsem_expr_file = rsem_file, 
                  output_basename=fusion_files[i,"sample_id"],
                  col_num=25,
                  sample_name=fusion_files[i,"sample_id"])
  )
}

lapply(1:nrow(fusion_files),function(x) run_format_fusion(i=x))


run_format_arriba <-function(i=1){
  # run formatting
  p<-a$project("d3b-0007-annofuse-dev")
  arriba_file<-p$file(id=fusion_files[i,"id.arriba"])
  app_id<-"gaonkark/d3b-0007-annofuse-dev/formatfusioncalls-annofuse"
  tsk<-a$project("d3b-0007-annofuse-dev")$task_add(
    name = paste0("Format fusion",fusion_files[i,"sample_id"]),
    app = app_id,
    inputs = list(inputfile = arriba_file,
                  outputfile = paste0(fusion_files[i,"sample_id"],".arriba_formatted.tsv"),
                  tumor_id = fusion_files[i,"sample_id"],
                  caller = "arriba"),
    input_check = F
  )
  tsk$run()
}

lapply(1:nrow(fusion_files),function(x) run_format_arriba(i=x))


arriba_formatted<-read_csv("../input/processed/1598278296271-manifest.csv")
arriba_formatted$sample_id<-unlist(lapply(strsplit(arriba_formatted$name,"[.]"), function(x) x[1]))


run_annotate_arriba <-function(i=1){
  # run annotation
  p<-a$project("d3b-0007-annofuse-dev")
  arriba_file<-p$file(id=arriba_formatted[i,"id"])
  # searching for id from name takes too lonng GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
  reference <-  p$file(id="5f2991b7e4b0e4c5ab17f0e1")
  app_id<-"gaonkark/d3b-0007-annofuse-dev/fusionannotator"
  tsk<-a$project("d3b-0007-annofuse-dev")$task_add(
    name = paste0("Annotate arriba fusion",arriba_formatted[i,"sample_id"]),
    app = app_id,
    inputs = list(input_fusion_file = arriba_file,
                  reference = reference,
                  colnum = 25,
                  outputname = paste0(arriba_formatted[i,"sample_id"],".arriba_formatted_annotated.tsv")
                  ),
    input_check = F
  )
  tsk$run()
}

lapply(1:nrow(arriba_formatted),function(x) run_annotate_arriba(i=x))




