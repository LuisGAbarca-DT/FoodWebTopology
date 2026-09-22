#  
#  
#               Food Web Topology Analysis
#                   Master Script
#                       V.1.2
#                   12 SEPTEMBER 2026
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#               IIB Universidad Veracruzana
#   
#   



#   ____________________________________________________________________________
#   1. Inicio                                                               ####
####-----------ERASE ALL DATA AND FUNCTIONS
####--------USE IT IF YOU ARE SURE TO REMOVE ALL DATA
rm(list=ls(all.names = TRUE))

##  ............................................................................
##  1.1. Install Dependencies                                               ####

#setwd("D:/PROYECTO_FW_AZAR_GITHUB/FW_SCRIPTS/FWTopo_V2")
#      INSTALL DEPENDENCIES
source("fwt_dependencies.R")

##  ............................................................................
##  1.2. Load Functions                                                     ####

#   READS AND RUNS ALL FUNCTIONS
#   
source("fwt_functions.R")

# start counting the time used
    start_time <- Sys.time()

    #CLEAR SCREEN (CONSOLE)
    cat("\014")

    #cheddar option to run regardless of the number of chains
    options(cheddarMaxQueue = 0)
    
    
#   ____________________________________________________________________________
#   2. Data & Directories                                                   ####
##  ............................................................................
##  2.1. Read data                                                        ####
    
####***LOAD DATA AND ARRABGE IT TO BE USED BY cheddar and igraph
####*
    datos <- read_data()
    file_name <- datos$file_name
    dat <- datos$dat
    cama<- datos$cama
    gr <- datos$gr
    
    names_1 <- cama$node
        rm(datos)
        #FLAG FOR THE FIRST OR NOT ROUND OF ANALYSES
        FIRST_ROUND <- FALSE

##  ...............................................................
##  2.2. Create directories and various options                 ####
        
# Create a new directory in the current working directory
# all results will be saved in there
# 
fwt_results_dir <- gsub(".csv", "", file_name)
fwt_results_dir <- str_c(fwt_results_dir, "_", "FWTA")

#           CHECAR LA EXISTENCIA DEL DIRECTORIO
# Check for original results
original_dir <- file.path(fwt_results_dir)    #, "RESULTS/ORIGINAL")
original_done <- dir.exists(original_dir)

if (!original_done) {
    
    #FLAG TO SHOW THAT IT IS THE FIRST TIME ANALYSIS
    FIRST_ROUND <- TRUE
    
    # Create ORIGINAL directory and compute
    dir.create(original_dir, recursive = TRUE, showWarnings = FALSE)
    fwt_results_dir_orig_res <- str_c(fwt_results_dir, "/RESULTS")
    dir.create(fwt_results_dir_orig_res, showWarnings = FALSE, recursive = TRUE) 
    
    fwt_results_dir_orig <- str_c(fwt_results_dir_orig_res, "/ORIGINAL")
    dir.create(fwt_results_dir_orig, showWarnings = FALSE, recursive = TRUE) 

    #DEFINE if you want a not very nice plot of the food web
    #
    figure <- get_figura_choice()
    if (figure == "1") {
        figure <- "YES"
    } else if (figure == "2"){
        figure <- "NO"
    }

    # DEFINE RESOLUTION FOR THE LEIDEN ALGORITHM
    algo <- get_resolution_choice()
        resolucion <- as.numeric(algo)

# OPTION: COMPUTE THE STRUCTURE FOR EACH MODULE FOR EACH RANDOM MATRIX?
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


#   ____________________________________________________________________________
#   3. Compute structure and topology                                 ####

        #CLEAR SCREEN (CONSOLE)
        #
        cat("\014")
        cat("\n")
        cat("\n")
        cat("\n")
        cat("      WORKING ON THE ORIGINAL FOOD WEB\n")
        cat("................................................\n")


        
##  ............................................................................
##  3.1. Compute structural components                                      ####
        
# Compute the structural components of the food web
# passing: gr -> the food web in igraph format
#           "YES" to calculate trophic levels
#           names_1 -> the names of the columns (nodes)

        cat("\n")
        cat("\n")
        cat("------------STRUCTURE AND TOPOLOGY---------------")
        
        
    original_str <- fw_struct_2 (gr, "YES", names_1$node)
    
        #TROPHIC LEVEL FOR EACH SPECIES
    
    niv_trof_sps <- as.data.frame(original_str$nivel_trof_sps)
    rownames(niv_trof_sps) <- names_1$node
    niv_trof_sps$SPS <- rownames(niv_trof_sps)
    colnames(niv_trof_sps)[1] <- "TL"
    
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_ORIGINAL_", 
                              "SPSsTL.csv")
        write.table(niv_trof_sps, file = write_file_1, append = F, sep = ",", row.names = FALSE)

    if (figure == "YES") {
        draw_1(gr, niv_trof_sps, fwt_results_dir)
    }

    
##  ............................................................................
##  3.2. Compute Efficiency                                                 ####
    cat("\n")
        cat("\n")
        cat("--------------------EFICIENCY--------------------")
    
#COMPUTE EFFICIENCY
    efic <- eficiency(cama)
    original_str$eficiency <- efic


##  ............................................................................
##  3.3. Compute Transitivity                                               ####

    cat("\n")
    cat("\n")
    cat("------------------TRANSITIVITY-------------------")
    
    tt <- transi(gr)
    original_str$trans <- tt 
    #remove tt object

##  ............................................................................
##  3.4. Compute Keystone                                                   ####
    
    cat("\n")
    cat("\n")
    cat("---------------------KEYSTONE--------------------")
    
    kb<-k.parameter(dat)
    kt<-k.parameter(t(dat))
    
        keystone<-data.frame(Kbu=kb[,3], Ktd=kt[,3], Kdir=kb[,1]+kt[,1], 
                             Kindir=kb[,2]+kt[,2], K=kb[,3]+kt[,3])
    #       make nice table
    #       move last col to first place
        keystone$SPS <- names_1$node
        
        keystone <- keystone %>%
            relocate(SPS)
    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", 
                          "ORIGINAL_KEYSTONE.csv")
    write.table(keystone, file = write_file_1, append = F, sep = ",", row.names = FALSE)

##  ............................................................................
##  3.5. Compute Centralities                                               ####
    
    cat("\n")
    cat("\n")
    cat("--------------------CENTRALITY-----------------")
    
    deg <-  grados(gr)
    
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", 
                              "ORIGINAL_CENTRALITIES.csv")
        write.table(deg, file = write_file_1, row.names = F, 
                    append = F, sep = ",")

##  ............................................................................
##  3.6. Compute Topological Importance                                     ####
    
    cat("\n")
    cat("\n")
    cat("--------------TOPOLOGICAL IMPORTANCE--------------")
    
    pasos <- 1
        TopoImp_1 <- TopologicalImportance(dat, pasos)
    
    pasos <- 3
        TopoImp_3 <- TopologicalImportance(dat, pasos)
    
    pasos <- 5
        TopoImp_5 <- TopologicalImportance(dat, pasos)
    
            TopoImpor.1.5 <- data.frame(TI1 = TopoImp_1$TI1, 
                                        TI3 = TopoImp_3$TI3,
                                        TI5 = TopoImp_5$TI5)
            TopoImpor.1.5$SPS <- names_1$node
            TopoImpor.1.5 <- TopoImpor.1.5 %>% 
                relocate(SPS)
    
    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "TOPO_IMPORT.csv")
    write.table(TopoImpor.1.5, file = write_file_1, row.names = F, 
                append = F, sep = ",")


##  ............................................................................
##  3.7. Compute Status                                                     ####

cat("\n")
cat("\n")
cat("--------------------STATUS-----------------------")

    Status <- StatusContrastatus(gr)
    
    Status$SPS <- names_1$node
    Status <- Status %>% 
        relocate(SPS)

    xx <- gsub(".csv", "", file_name)
    write_file_1 <- str_c(fwt_results_dir_orig, "/", xx, "_", "STATUS.csv")
    write.table(Status, file = write_file_1, row.names = F, 
                append = F, sep = ",")
    
##  ............................................................................
##  3.8. Compute Modularity                                                 ####
    
    
    cat("\n")
    cat("\n")
    cat("------------------MODULARITY---------------------")
    
    resulta7_a <- modul_l(gr, resolucion)
    
    members_leid <- resulta7_a$membership  
    members_leid <- as.numeric(members_leid)
    modularity_leid <- igraph::modularity(gr, members_leid, directed = TRUE) 
                                  
    module_numb_leid <- max(as.numeric(members_leid))
    modularidad <- modularity_leid
    
    original_str$modularity <- modularidad
    original_str$No_Modules <- module_numb_leid
    
    #ASSIGN VALUE TO THE WHOLE DATA FRAME
        estruct_original_1 <- do.call(rbind, lapply(original_str, 
                                                    as.data.frame))
    
            xx <- gsub(".csv", "", file_name)
            write_file_1 <- str_c(fwt_results_dir_orig, "/", 
                                  xx, "_", "ORIGINAL_STRC.csv")
            write.table(estruct_original_1, file = write_file_1, 
                        append = F, sep = ",")
           
             if (figure == "YES") {

                draw_modules(gr, module_numb_leid, resulta7_a)
                
            }


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### 3.8.1. Compute Module Partition                                         ####
            
#  to which module the nodes pertain to 
    parti <- NULL
                
        parti <- as.data.frame(resulta7_a$membership)
                parti[,2] <- row.names(parti)
                parti[,3] <- as.numeric(parti$'resulta7_a$membership')
                parti <- parti %>% 
                    select(-'resulta7_a$membership')
                rownames(parti) <- NULL
    
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", 
                              xx, "_", "ORIGINAL_LEIDEN_PART.csv")
        write.table(parti, file = write_file_1, row.names = F, 
                    append = F, sep = ",")
    

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### 3.8.2. Compute Module Structure                                         ####
    
#STRUCTURE OF EACH MODULE FOR THE ORIGINAL WEB

    orig_module_str <- list(prop_modules(groups = as.numeric(resulta7_a$membership), 
                                         g_rand = gr, graphic ="no"))
        my_dataframe <- bind_rows(orig_module_str)
    
            xx <- gsub(".csv", "", file_name)
            write_file_1 <- str_c(fwt_results_dir_orig, "/", 
                                  xx, "_", "ORIGINAL_MODULE_STRC.csv")
            write.table(my_dataframe, file = write_file_1,
                        append = F, sep = ",")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### 3.8.3. Compute Crossed Nodes                                            ####
    
    shared_nodes <- intersected(gr, resulta7_a)
    
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", 
                              xx, "_", "ORIGINAL_SHARED_NODES.csv")
        write.table(shared_nodes, file = write_file_1, row.names = T, 
                    append = F, sep = ",")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### 3.8.4. Compute Node Roles                                               ####
    
    roles <- rol_nodos(gr, members_leid)
    
        xx <- gsub(".csv", "", file_name)
        write_file_1 <- str_c(fwt_results_dir_orig, "/", 
                              xx, "_", "ORIGINAL_NODES_ROL.csv")
        write.table(roles, file = write_file_1, row.names = F, 
                    append = F, sep = ",")
    
        

}   #END OF THE QUESTION IF THERE IS A PREVIOUS ANALYSIS FOR THAT FILE


#   ____________________________________________________________________________
#   4. Randomization                                                        ####

#COMPUTE RANDOM WEBS REGARDLESS IF THERE IS OR NOT A PREVIOUS ANALYSIS
#IF THERE IS A PREVIOUS ANALYSIS THEN ASK IF A RANDOMIZED WEB WILL BE PRODUCED
#AND ANALYZED


##  ............................................................................
##  4.1. Create Directory for Random                                                  ####

#create the random directory eith the reults of the random fw
fwt_results_dir_orig_res <- str_c(fwt_results_dir, "/RESULTS")
fwt_results_dir_rand <- str_c(fwt_results_dir_orig_res, "/RANDOM")
dir.create(fwt_results_dir_rand, showWarnings = FALSE, recursive = TRUE) 

#create directory which contains the results of the random fw
fwt_results_dir_rand_webs <- str_c(fwt_results_dir_orig_res, "/RANDOM/RANDOMIZED_WEBS")
dir.create(fwt_results_dir_rand_webs, showWarnings = FALSE, recursive = TRUE) 


##  ............................................................................
##  4.2. Choose Algorithm                                                   ####

#    ALGORITHM TO USE TO GENERATE RANDOM FOOD WEBS

# Algorithm mapping (number to name)
algorithm_map <- c(
    "1" = "erdos-renyi",
    "2" = "cascade",
    "3" = "niche-Williams-Cohen",
    "4" = "niche-Allesina-Alonso-Pascal",
    "5" = "Link randomization",
    "6" = "none"
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
                    random_model_type <- "random_links"
                }
                else {
                    if (algo_num == 6) {
                        random_model_type <- "none"
                        num_rand_webs <- 0
                        tiempo <- 0
                    }
                }
            }
        }
    }
}


    
if (algo_num != 6) {
#   DEFINE NUMBER OF RANDOM FOOD WEBS TO BE GENERATED
    #
    algo <- get_num_fw_choice()
    #algo
    num_rand_webs <- as.numeric(algo)

    tl_rnd <- get_TL_choice()

    if (tl_rnd == "1") {
        tl_rnd_y <- "YES"
    } else if (tl_rnd == "2") {
        tl_rnd_y <- "NO"
    }
}

    # chose if you want a tl density plot of the mean random tl
    # 
    #DEFINE if you want a not very nice plot of the density of TL
    #
    figure_tldens <- get_figura_tldensity_choice()
        if (figure_tldens == "1") {
            figure_tldens <- "YES"
        } else if (figure_tldens == "2"){
            figure_tldens <- "NO"
        }
    
##  ............................................................................
##  4.3. Generate Random Web                                                ####
    
 #if we want to generate random food webs go ahead, when algo_num less 
 #than 5 if not jump this section

    if (algo_num != 6) {
    #Generates and analyzes the random food webs, and saves the results
        tl_mean <- gen_analiza_rand_fw(gr, names_1$node, num_rand_webs, 
                                       random_model_type, 
                                       fwt_results_dir_rand_webs)

        if (figure_tldens == "YES") {
            
            final_results_2 <- TL_Compute(gr, cama)
                mean_tl_obs <- mean(final_results_2$ShortWeightedTL)
            tl_mean_rnd <- as.data.frame(tl_mean)
                draw_rnd_dens(tl_mean_rnd, mean_tl_obs)
 
         }

    }
    

##  ............................................................................
##  4..4. Log file                                                          ####
    
        #       CHECAR PORQUE A LA SEGUNDA VUELTA LAS VARIABLES DE 
        #       RESULTADOS ORIGINALES YA NO EXISTEN
        #       

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### 4.4.1. First Time                                                       ####
    
if (FIRST_ROUND == TRUE) {
    
    tiempo<-difftime(Sys.time(), start_time, units = "secs")    
    cat("\n")
    cat("Time spent :\n")
    
    cat(tiempo/60, "MINUTES\n")
    
        # produces a log file of the analyses
        # 
        validation <- generate_validation_report(dat, file_name, 
                                                 random_model_type,
                                                 num_rand_webs, resolucion, 
                                                 tiempo, 
                                                 fwt_results_dir_orig_res)
    
    #CLEAR SCREEN (CONSOLE)
    #
    cat("\014")
    validation <- as.data.frame(do.call(rbind, validation))
    validation$V2 <- NULL
    validation[4,1] <- format(Sys.time(), "%Y.%m.%d___%H.%M.%S_")
    print(validation)
    #CLEAN MEMORY
    gc()


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### 4.4.2. Second,,,Time                                                    ####
    
    #else generate a validation report or log file for the new random analysis
}   else{
    
    tiempo<-difftime(Sys.time(), start_time, units = "secs")    
    cat("\n")
    cat("Time spent :\n")
    
    cat(tiempo/60, "MINUTES\n")
    
        validation <- gen_valid_for_rnd(dat, file_name, random_model_type, 
                                        num_rand_webs,
                                        fwt_results_dir_rand,
                                        tiempo)
        #CLEAR SCREEN (CONSOLE)
        #
        cat("\014")
        validation <- as.data.frame(do.call(rbind, validation))
        validation$V2 <- NULL
        validation[4,1] <- format(Sys.time(), "%Y.%m.%d___%H.%M.%S_")
        print(validation)
        #CLEAN MEMORY
        gc()
}
    

#   ____________________________________________________________________________
#   5. END                                                                  ####
    

    cat("               =====================\n")
    cat("                 ...LISTO...READY...\n")
    cat("                    ... ALL DONE ...\n")
    cat("               =====================\n")



