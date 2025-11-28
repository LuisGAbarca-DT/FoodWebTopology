#  
#  
#               Food Web Topology Analysis
#                   Master Script
#                       V.1.1
#                   26 NOVEMBER 2025
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#               IIB Universidad Veracruzana
#   
#   

# INICIO ------------------------------------------------------------------
#setwd("D:/PROYECTO_FW_AZAR_GITHUB/FW_SCRIPTS/FWTopo_V2")
#      INSTALL DEPENDENCIES
source("install-dependencies.R")


#REMOVE-DATA -------------------------------------------------------------
####-----------ERASE ALL DATA AND FUNCTIONS
####--------USE IT IF YOU ARE SURE TO REMOVE ALL DATA
rm(list=ls(all.names = TRUE))

#   READS AND RUNS ALL FUNCTIONS
#   
source("FWTopo_functions.R")

    #CLEAR SCREEN (CONSOLE)
    #
    cat("\014")

    # cat("  _________________________________________________\n")
    # cat("               Food Web Topology Analysis\n")
    # cat("  _________________________________________________\n")
    
#   READ DATA
#   the adjacency matrix as an .csv type
#  file_name -------------------------------------------------------------
    cat("\014")
    cat("  _________________________________________________\n")
    cat("               Food Web Topology Analysis\n")
    cat("                         V 1.4 \n")
    cat("  _________________________________________________\n")
    cat("\n")
    
    cat("Choose the data file \n")
file_address <- file.choose()
dat <- read.table(file_address,sep = ",", header = T)
file_name <- basename(file_address)  

head(dat,3)
tail(dat,3)

   
names_1 <- as.matrix(dat[,1])
#head(names_1)

dat<-as.matrix(dat[,-1])

rownames(dat) <- paste("SPS_", names_1, sep = "")
colnames(dat) <- paste("SPS_", names_1, sep = "")

# ARRENGING DATA FOR cheddar --------------------------------------------------------

NODE<-rownames(dat)
cama<- Community(nodes=data.frame(node=NODE),
                 trophic.links=PredationMatrixToLinks(dat),
                 properties=list(title="Community"))

#ARRENGING DATA FOR iraph
dat_mat <- as.matrix(dat)
gr<-graph_from_adjacency_matrix(dat_mat, weighted = FALSE, mode = c("directed"))

if(igraph::is_connected(gr) == "FALSE"){
    stop("Can not proceed with the analysis. 
         The graph is not completely connected")
}else {
    cat("All good")
}

# CHECK THAT THERE IS AT LEAST 1 BASAL AND 1 TOP
# 
#Especies BSALES
basal <- BasalNodes(cama)
b <- length(basal)

#Especies TOPE
top <- TopLevelNodes(cama)
tope <- length(top)

if (b  == 0 | tope == 0) {
    stop("There are no basal or top nodes. Can not proceed with the analysis")
} else {
    cat("All good\n")
}

# CREATE_DIRS -------------------------------------------------------------

# Create a new directory in the current working directory
# all results will be saved in there
# 
fwt_results_dir <- gsub(".csv", "", file_name)
fwt_results_dir <- str_c(fwt_results_dir, "_", "FWTA")
dir.create(fwt_results_dir, recursive = TRUE) 

# Create a directory with a specific path
# This will create original with the results of the original fw
 
fwt_results_dir_orig_res <- str_c(fwt_results_dir, "/RESULTS")
dir.create(fwt_results_dir_orig_res, recursive = TRUE) 

fwt_results_dir_orig <- str_c(fwt_results_dir_orig_res, "/ORIGINAL")
dir.create(fwt_results_dir_orig, recursive = TRUE) 

#create the random directory eith the reults of the random fw
fwt_results_dir_rand <- str_c(fwt_results_dir_orig_res, "/RANDOM")
dir.create(fwt_results_dir_rand, recursive = TRUE) 

# SELECT_RAND_ALGO --------------------------------------------------------

#    ALGORITHM TO USE TO GENERATE RANDOM FOOD WEBS

# FWTopo Algorithm Selection System
# Version: 1.0

# Algorithm mapping (number to name)
algorithm_map <- c(
    "1" = "erdos-renyi",
    "2" = "cascade",
    "3" = "niche-Williams-Cohen",
    "4" = "niche-Allesina-Alonso-Pascal",  
    "5" = "none"
)

    algo <- get_algorithm_choice()

    algo <- as.data.frame(algo)
    algo_num <- as.numeric(rownames(algo))
    
if(algo_num == 1){
    random_model_type <- "erdos-renyi"
} else{
    if (algo_num == 2) {
        random_model_type <- "cascade"
    } else {
        if(algo_num == 3) {
            random_model_type <- "niche-model"
        } else {
            if (algo_num == 4){
                random_model_type <- "niche_allesina"
            }
            else {
                if (algo_num == 5){
                    random_model_type <- "none"
                    num_rand_webs <- 0
                    tiempo <- 0
                }
            }
        }
    }
}

if (algo_num != 5) {
#   DEFINE NUMBER OF RANDOM FOOD WEBS TO BE GENERATED
    #   
    algo <- get_num_fw_choice()
    algo
    num_rand_webs <- as.numeric(algo)
    
    tl_rnd <- get_TL_choice()
    
    if (tl_rnd == "1") {
        tl_rnd_y <- "YES"
    } else if (tl_rnd == "2") {
        tl_rnd_y <- "NO"
    }
}
    
# DEFINE RESOLUTION FOR THE LEIDEN ALGORITHM
    algo <- get_resolution_choice()
    algo
    resolucion <- as.numeric(algo)

#DEFINE if you want a not very nice plot of the food web
    #
    figure <- get_figura_choice()
    if (figure == "1") {
        figure <- "YES"
    } else if (figure == "2"){
        figure <- "NO"
    }
    figure
    

# COMPUTE THE STRUCTURE FOR EACH MODULE FOR EACH RANDOM MATRIX?
# NOTE: BIG FOOD WEBS WILL TAKE LONGER TIME AND A LOT OF RAM
# YES
# NO
compute_module_str <- "NO"

# 
# LIST OF THE MODULES STRUCTURE FOR THE ORIGINAL WEB
orig_module_str <- list()

#cheddar option to continue computting trophic level without producing an error
#for big food webs
    options(cheddarMaxQueue = 0)

        #CLEAR SCREEN (CONSOLE)
        #
        cat("\014")
        cat("\n")
        cat("\n")
        cat("\n")
        cat("      WORKING ON THE ORIGINAL FOOD WEB\n")
        cat("................................................\n")

if (figure == "YES") {
    draw_1(gr, cama)
}

# COMPUTE STRUCTURE ORIGINAL -----------------------------------------------------------

original_str <- fw_struct_2 (gr, "YES", names_1)
    #TROPHIC LEVEL FOR EACH SPECIES
    niv_trof_sps <- as.data.frame(original_str$nivel_trof_sps)
    colnames(niv_trof_sps)[1] <- "TL"
    
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_ORIGINAL_", "SPSsTL.csv")
        write.table(niv_trof_sps, file = write_file_1, append = F, sep = ",")
        
    original_str$nivel_trof_sps <- NULL
        
#   COMPUTE TRANSITIVITY
tt <- transi(gr)
original_str$trans <- tt 
#remove tt object
rm(tt)

# KEYSTONE ----------------------------------------------------------

kb<-k.parameter(dat)
kt<-k.parameter(t(dat))

    keystone<-data.frame(Kbu=kb[,3], Ktd=kt[,3], Kdir=kb[,1]+kt[,1], 
                         Kindir=kb[,2]+kt[,2], K=kb[,3]+kt[,3])

    keystone <- as.data.frame(keystone)
    keystone$SPS <- paste("SPS", names_1, sep = "")
#       make nice table
#       move last col to first place
    keystone <- keystone %>%
        relocate(SPS)
    gt_tabla <- gt(keystone)
    gt_tabla
    # gt_tabla %>%
    #     cols_move_to_start(columns = SPS)
    # 
    tabla_K <- as.data.frame(gt_tabla)

xx <- gsub(".csv", "", file_name)
write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_KEYSTONE.csv")
write.table(tabla_K, file = write_file_1, append = F, sep = ",")

#remove keystone object
rm(keystone)

# CENTRALIDAD -------------------------------------------------------------

deg <-  grados(gr)

    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_CENTRALITIES.csv")
    write.table(deg, file = write_file_1, append = F, sep = ",")
#remove deg object
rm(deg)

#------------ COMPUTE TOPOLOGICAL IMPORTANCE ----------
#
    pasos <- 1
        TopoImp_1 <- TopologicalImportance(dat, pasos)
    
    pasos <- 3
        TopoImp_3 <- TopologicalImportance(dat, pasos)
    
    pasos <- 5
        TopoImp_5 <- TopologicalImportance(dat, pasos)
    
            TopoImpor.1.5 <- data.frame(TI1 = TopoImp_1$TI1, 
                                        TI3 = TopoImp_3$TI3,
                                        TI5 = TopoImp_5$TI5)
    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "TOPO_IMPORT.csv")
    write.table(TopoImpor.1.5, file = write_file_1, append = F, sep = ",")
#remove pasos object
rm(pasos)

#-------COMPUTE STATUS Y CONTRASTATUS-------------

Status <- StatusContrastatus(gr)
xx <- gsub(".csv", "", file_name)
write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "STATUS.csv")
write.table(Status, file = write_file_1, append = F, sep = ",")

#remove Status object
rm(Status)

#-------------   COMPUTE MODULARITY FOR THE ORIGINAL FOOD WEB
# MODULARITY --------------------------------------------------------------

    resulta7_a <- modul_l(gr, resolucion)
    
    members_leid <- resulta7_a$membership  
    members_leid <- as.numeric(members_leid)
    modularity_leid <- modularity(gr, members_leid, directed = TRUE) 
                                  
    module_numb_leid <- max(as.numeric(members_leid))
    modularidad <- modularity_leid
    
    original_str$modularity <- modularidad
    original_str$No_Modules <- module_numb_leid
    
    #ASSIGN VALUE TO THE WHOLE DATA FRAME
        estruct_original_1 <- do.call(rbind, lapply(original_str, as.data.frame))
    
            xx <- gsub(".csv", "", file_name)
            write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_STRC.csv")
            write.table(estruct_original_1, file = write_file_1, append = F, sep = ",")
           
             if (figure == "YES") {

                draw_modules(gr, module_numb_leid, resulta7_a)
                
            }

# PRTITION ----------------------------------------------------------------
#  to which module the nodes pertain to 
parti <- NULL
            
    parti <- as.data.frame(resulta7_a$membership)
            parti[,2] <- row.names(parti)
            parti[,3] <- as.numeric(parti$'resulta7_a$membership')
            parti <- parti %>% 
                select(-'resulta7_a$membership')
            rownames(parti) <- NULL

    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_LEIDEN_PART.csv")
    write.table(parti, file = write_file_1, append = F, sep = ",")

#remove parti object
rm(parti)

# MODULE_STRUCTURE --------------------------------------------------------
#STRUCTURE OF EACH MODULE FOR THE ORIGINAL WEB

orig_module_str <- list(prop_modules(groups = as.numeric(resulta7_a$membership  ), g_rand = gr, graphic ="no"))
    my_dataframe <- bind_rows(orig_module_str)

        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_MODULE_STRC.csv")
        write.table(my_dataframe, file = write_file_1, append = F, sep = ",")
#remove orig_module_str
rm(orig_module_str)

# CROSSED_NODES -----------------------------------------------------------
        
shared_nodes <- intersected(gr, resulta7_a)
    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_SHARED_NODES.csv")
    write.table(shared_nodes, file = write_file_1, append = F, sep = ",")
#remove shared_nodes object
rm(shared_nodes)

# NODE_ROLE ---------------------------------------------------------------

roles <- rol_nodos(gr)
    roles[[2]] <- NULL
# roles
    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "ORIGINAL_NODES_ROL.csv")
    write.table(roles, file = write_file_1, append = F, sep = ",")
#remove roles object
rm(roles)

# GEN_RND_FW --------------------------------------------------------------
 #if we want to generate random food webs go ahead, whwn algo_num less 
 #than 5 if not jump this section

 if (algo_num != 5) {
         
    cat("GENERATING RANDOM FOOD WEBS...\n")
    #GENERATE RANDOM FOOD WEBS ACCORDING TO THE CHOSEN ALGORITHM 
        rand_matrix <- list()
        rand_matrix <- gen_rnd_fw(gr, names_1, num_rand_webs, random_model_type)
        
    # RND_FW_STRUCT -----------------------------------------------------------
    # this routine computes the topology of the random food webs one by one
    # this will avoid loosing the information if there is a crash due to 
    # RAM overflow
        estructura_azar <- NULL
        start_time <- Sys.time()
        
        #INITIAL MEMORY
        mem <- get_memory_usage()
        
        
    for (i in 1: num_rand_webs) {
        cat("\014")
        cat("  _________________________________________________\n")
        cat("               Food Web Topology Analysis\n")
        cat("                         V 1.4 \n")
        cat("  _________________________________________________\n")
        cat("\n")
        
        cat("\n")
        cat("\n")
        print(". . . w o r k i n g. . . ")
        cat("\n")
        cat("Random Food Web being analyzed Number: \n")
        print(i)
        
        #STRUCTURE...
            res_estr_rnd <- NULL
            res_estr_rnd <- fw_struct_rnd(rand_matrix[[i]], tl_rnd_y, names_1)
            
            # Convert results to a data frame
            estructura_azar <- rbind(estructura_azar, res_estr_rnd)
            #estructura_azar <- do.call(rbind, lapply(res_estr_rnd, as.data.frame))
            
            # print("Number of random webs analyzed:")
            # print(i)
            
            #MEMORY USE
                cat("Memory used: \n")
                mem_2 <- get_memory_usage()
                print(mem$total - mem_2$total)
            
    }
        
        tiempo<-difftime(Sys.time(), start_time, units = "secs")    
        cat("\n")
        cat("Time spent :\n")
        
        cat(tiempo/60, "MINUTES\n")
        
        
        
        #MODULARITY...       
            
            #MODULARITY FOR EACH RANDOM FOOD WEB
            #
            modularity_leid_rnd <- list()
            modularity_leid_rnd <- lapply(1:num_rand_webs, function(i) modul_l(rand_matrix[[i]], 1))
            
            #MODULARITY OF RANDOM FOOD WEBS AS LIST
            mod_leid_azar <- lapply(1:num_rand_webs, function(i) modularity(rand_matrix[[i]], modularity_leid_rnd[[i]]$membership, directed = TRUE))
            #MODULARITY OF RANDOM FOOD WEBS AS DATA FRAME
            modul_leiden_rnd <- do.call(rbind, lapply(mod_leid_azar, as.data.frame))
            
            #NUMBER OF MODULES FOR EACH RANDOM FOOD WEB
            num_modul_leid_rnd <- lapply(1:num_rand_webs, function(i) max(as.numeric(modularity_leid_rnd[[i]]$membership)))
            num_modul_leid_rnd <- do.call(rbind, lapply(num_modul_leid_rnd, as.data.frame))
            
            estructura_azar[,14] <- num_modul_leid_rnd
            estructura_azar[,15] <- modul_leiden_rnd
            
            colnames(estructura_azar)[14] <- "No_Modules"
            colnames(estructura_azar)[15] <- "Modularity"
    
            
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_rand, "/", xx, "_", 
                              "STRUCT_AZAR_", random_model_type, ".csv")
            # Instead of overwriting, auto-archive:
            if (file.exists(write_file_1)) {
                xx2 <- str_c(fwt_results_dir_rand, "/", xx, "_", random_model_type, 
                             format(Sys.time(), "%Y%m%d_%H%M%S_"),".csv")
                #archive_name <- paste0(write_file_1, format(Sys.time(), "%Y%m%d_%H%M%S_"), xx2)
                file.rename(write_file_1, xx2)
            } else if(!file.exists(write_file_1)) {
                write.table(estructura_azar, file = write_file_1, 
                    append = F, sep = ",")
            }
        
}
        
# VALIDATION REPORT -----------------------
    # produces a log file of the analyses
    # 
    validation <- generate_validation_report(dat, file_name, random_model_type,
                                             num_rand_webs, resolucion, tiempo, 
                                             fwt_results_dir_orig_res)
                
    #CLEAN MEMORY
    gc()
    #CLEAR SCREEN (CONSOLE)
    #
    cat("\014")
    validation <- as.data.frame(do.call(rbind, validation))
    validation$V2 <- NULL
    print(validation)
    
# END ---------------------------------------------------------------------
     
#                       END           
   
    cat("=====================\n")
    cat("  ...LISTO...READY...\n")
    cat("    ... ALL DONE ...\n")
    cat("=====================\n")



