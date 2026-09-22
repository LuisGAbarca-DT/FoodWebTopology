#  
#  
#               Food Web Topology Analysis
#                   FUNCTIONS
#                     V 1.2
#                     
#               12 SEPTEMBER 2026
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#               IIB Universidad Veracruzana
#   
#   
#   
#   FUNCTION TO DRAW THE DENSITY OF TROPHIC LEVELS OF A RANDOM FOOD WEB AND
#  A NOTATION AS AN ARROW OF THE EMPIRICAL (OBSERVED) TROPHIC LEVEL OF 
#  THE ORIGINAL FOOD WEB  
#   

#FUNCTION TO CALCULATE THE TROPHIC LEVEL OF EACH NODE OF THE FOOD WEB
#RECEIVES   gr_tempo THE FOOD WEB IN IGRAPH FORMAT
#           comty THE FOOD WEB AS A COMMUNITY OF cheddar PCKG
#RETURNS final_results WITH NIMBER OF NODES AND THE TROPHIC LEVEL OF EACH NODE
# ACCORDING TO THREE DIFFERENT ALGORITHMS

TL_Compute <- function(g_rand, commty) {
    
    # chain.stats <- TrophicChainsStats(commty)
    # chains_nom <- length(chain.stats$chain.lengths)
    # rm(chain.stats) #remove this object to get RAM
    # browser()
    # 3. METRIC 1: PREY-AVERAGED TL (The native cheddar version)
    prey_avg_tl <- PreyAveragedTrophicLevel(commty)
    
    # 4. METRIC 2: SHORTEST TL (The 1st version via igraph)
    links <- TLPS(commty)
    basal_nodes <- BasalNodes(commty)
    
    path_lengths <- distances(g_rand, to = basal_nodes, mode = "in")
    
    shortest_tl <- apply(path_lengths, 1, min) + 1
    shortest_tl[is.infinite(shortest_tl)] <- NA
    
    # 5. ALIGN AND INTERSECT BY SPECIES NAMES
    # We find only the species names present in both outputs to prevent alignment NAs
    common_species <- intersect(names(shortest_tl), names(prey_avg_tl))
    
    # Filter and align both vectors to match perfectly
    shortest_tl_aligned <- shortest_tl[common_species]
    prey_avg_tl_aligned <- prey_avg_tl[common_species]
    
    
    # 6. COMPUTE THE SHORT-WEIGHTED TROPHIC LEVEL
    short_weighted_tl <- (shortest_tl_aligned + prey_avg_tl_aligned) / 2
    
    # 7. COMBINE RESULTS
    final_results <- data.frame(
        Node = common_species,
        ShortestTL = shortest_tl_aligned,
        PreyAveragedTL = prey_avg_tl_aligned,
        ShortWeightedTL = short_weighted_tl,
        row.names = NULL
    )    
    
    return(final_results)
}

#FUNCTION TO DRAW THE DENSITY OF TROPHIC LEVELS OF A RANDOM FOOD WEB AND
# A NOTATION AS AN ARROW OF THE EMPIRICAL (OBSERVED) TROPHIC LEVEL OF 
# THE ORIGINAL FOOD WEB
# 
draw_rnd_dens <- function(tl_azar, tl_observ) {
    
    #tl_observ = the trophic level of the original food web
    #rnd_tl = the trophic level of the random food web
    observed_df_1 <- data.frame(tl_observ)
    value <-as.data.frame(tl_azar)

    df <- data.frame(
        value = value$tl_mean,
        web = rep("RndWeb"),
            times = c(length(value$tl_mean))
        
    )
    obs <- data.frame(
        web = c("RndWeb"),
    x = c(tl_observ)
    )
    
    xxx<-split(df$value,df$web)
    
    media_null_1 <- mean(xxx$RndWeb)
    
    medias_rand <- data.frame(
        web = c("RndWeb"),
        x = c(media_null_1)
    )
    
    my_colors <- c(
        "RndWeb" = "coral"
    )
    # --- Compute max density height for arrow scaling ---
    dens_max <- max(sapply(split(df$value, df$web),
                           function(v) max(density(v)$y)))
    # --- Plot ---
    p = ggplot(df, aes(x = value, fill = web, color = web)) +
        geom_density(alpha = 0.5, linewidth = 0.4) +
        geom_segment(
            data = obs,
            aes(x = x, xend = x, y = 0, yend = dens_max * 0.0, color = web),
            arrow = arrow(length = unit(0.5, "cm")),
            linewidth = 0.9,
            inherit.aes = FALSE,
            show.legend = FALSE  
        ) +
        geom_segment(
            data = medias_rand,
            aes(x = x, xend = x, y = 0, yend = dens_max, color = web),
            linewidth = 0.5,
            inherit.aes = FALSE,
            show.legend = FALSE
        ) +
        scale_fill_manual(values = my_colors, name = "Web") +
        scale_color_manual(values = my_colors, name = "Web") +
        labs(
            x = "Trophic Level",
            y = "Density"
        ) +
        theme_few(base_size = 15)
    
    print (p)

return()
    
}

#   
#   FUNCTION TO GENERATE A FOOD WEB RANDOMIZING THE LINKS BETWEEN NODES
#   AND PRESERVING THE IN- AND OUT-DEGREES

random_links <- function(gr) {
    # g <- gr
    iteraciones <- 100*ecount(gr)
    # generate the random web using igraph function
    random_gr <- rewire(
        gr,
        with = keeping_degseq(loops = FALSE,
                              niter = iteraciones)
    )
    
    cat("Random web is connected?  \n",igraph::is_connected(random_gr))
    cat("In-degree for each node: \n", igraph::degree(random_gr, mode = "in"))
    cat("Out-degree for each node: \n", igraph::degree(random_gr, mode = "out"))
    
    return(random_gr)
}

#FUNCTION TO GENERATE RANDOM FOOD WEBS AND ANALYZE THEIR TOPOLOGY
#
gen_analiza_rand_fw <- function(gr, names_1, num_rand_webs, random_model_type, 
                                fwt_results_dir_rand_webs) {
    
    cat("GENERATING RANDOM FOOD WEBS...\n")
    #GENERATE RANDOM FOOD WEBS ACCORDING TO THE CHOSEN ALGORITHM 
    rand_matrix <- list()
    rand_matrix <- gen_rnd_fw(gr, names_1, num_rand_webs, random_model_type, 
                              fwt_results_dir_rand_webs)

        # RND_FW_STRUCT 
    # this routine computes the topology of the random food webs one by one
    # this will avoid loosing the information if there is a crash due to 
    # RAM overflow
    estructura_azar <- NULL

        cat("\014")
        cat("  _________________________________________________\n")
        cat("               Food Web Topology Analysis\n")
        cat("                         V 1.2\n")
        cat("  _________________________________________________\n")
        cat("\n")
        
        cat("\n")
        cat("\n")
        cat(". . . w o r k i n g. . . ")
        cat("\n")
        cat("Random Food Webs being analyzed \n")
        #cat("\n              ", i)    

        #STRUCTURE...
        res_estr_rnd <- NULL
        #res_estr_rnd <- fw_struct_rnd(rand_matrix[[i]], tl_rnd_y, names_1)
        
            res_estr_rnd <- lapply(1:num_rand_webs, function(i) 
                fw_struct_rnd(rand_matrix[[i]], tl_rnd_y, names_1))
        
        #browser()
        # Convert results to a data frame
            estructura_azar <- do.call(rbind, lapply(res_estr_rnd, as.data.frame))

    #MODULARITY...       
        #MODULES FOR EACH RANDOM FOOD WEB

    modularity_leid_rnd <- list()
    modularity_leid_rnd <- lapply(1:num_rand_webs, function(i) 
        modul_l(rand_matrix[[i]], 1))
    
    #MODULARITY OF RANDOM FOOD WEBS AS LIST
    mod_leid_azar <- lapply(1:num_rand_webs, function(i) 
        igraph::modularity(rand_matrix[[i]], modularity_leid_rnd[[i]]$membership,
                           resolition = 1, directed = TRUE))
    #MODULARITY OF RANDOM FOOD WEBS AS DATA FRAME
    modul_leiden_rnd <- do.call(rbind, lapply(mod_leid_azar, as.data.frame))
    
    #OBTAIN THE NUMBER OF MODULES FOR EACH RANDOM FOOD WEB
    num_modul_leid_rnd <- lapply(1:num_rand_webs, function(i) 
        max(as.numeric(modularity_leid_rnd[[i]]$membership)))
    num_modul_leid_rnd <- do.call(rbind, lapply(num_modul_leid_rnd, as.data.frame))
    
    estructura_azar[,15] <- num_modul_leid_rnd
    estructura_azar[,16] <- modul_leiden_rnd
    
    colnames(estructura_azar)[15] <- "No_Modules"
    colnames(estructura_azar)[16] <- "Modularity"
    
    xx <- gsub(".csv", "", file_name)
        xx2 <- str_c(fwt_results_dir_rand, "/", xx, "_", random_model_type, "_", 
                     format(Sys.time(), "%Y.%m.%d___%H.%M.%S_"),".csv")
    write.table(estructura_azar, file = xx2, 
                append = F, sep = ",")

 return(estructura_azar$TLMean)   

}

eficiency <- function(cama) {
    #computes efficiency of food web
    cam_corto_all <- ShortestPaths(cama)
    nodos <- NumberOfNodes(cama)
    # VALOR DE E 1/(n(n-1)) suma i-j (1/(d ij))
    # 
    inv_camino_corto_all <- 1 / cam_corto_all
    diag(inv_camino_corto_all) <- 0
    
    suma_caminos <- sum(rowSums(inv_camino_corto_all))
    max_l <- nodos * (nodos - 1)
    
    E <- (1 / max_l) * suma_caminos
    
    return(E)
}

read_data <- function () {
    #receives nothing
    #Returns a list with dat as the dayacency matrix
    #                    cama as a chedar community list
    #                    gr as an igraph list
    
    #   READ DATA
    #   the adjacency matrix as an .csv type
    
    cat("\014")
    cat("  _________________________________________________\n")
    cat("               Food Web Topology Analysis\n")
    cat("                         V 1.2 \n")
    cat("  _________________________________________________\n")
    cat("\n")
    
    cat("Choose the data file \n")
    file_address <- file.choose()
    dat <- read.table(file_address,sep = ",", header = T)
    file_name <- basename(file_address)  
    
    names_1 <- as.matrix(dat[,1])
    
    # #head(names_1)
    # 
    dat<-as.matrix(dat[,-1])
    # 
    rownames(dat) <- names_1  
    colnames(dat) <- rownames(dat)    
    head(dat,2)
    
    cat("Validating format for:", file_name, "\n")
    
    # 1. CHECK: Matrix is square
    if (nrow(dat) != ncol(dat)) {
        stop("ERROR: The adjacency matrix is not square.\n",
             "  -> Found ", nrow(dat), " rows and ", ncol(dat), " columns.\n",
             "  -> Please ensure the number of species (rows) matches the number of predators (columns).")
    }
    
    # 2. CHECK: Row names match column names IN THE SAME ORDER
    if (!identical(rownames(dat), colnames(dat))) {
        # Give a helpful hint about the first mismatch
        mismatch_idx <- which(rownames(dat) != colnames(dat))[1]
        stop("ERROR: Species order mismatch.\n",
             "  -> At position ", mismatch_idx, ", row name is '", 
             rownames(dat)[mismatch_idx],
             "' but column name is '", colnames(dat)[mismatch_idx], "'.\n",
             "  -> Please ensure the species list is identical and in the *same order* for rows and columns.")
    }
    
    # 3. CHECK: All values are 0 or 1 (Binary Adjacency)
    if (!all(as.matrix(dat) %in% c(0, 1))) {
        invalid_vals <- unique(as.vector(as.matrix(dat)[!as.matrix(dat) %in% c(0, 1)]))
        stop("ERROR: Matrix contains non-binary values.\n",
             "  -> Found values: ", paste(invalid_vals, collapse = ", "), "\n",
             "  -> Please ensure all interactions are coded as 0 (absent) or 1 (present).")
    }
    
    # 5. NEW CHECK: Names must not be pure numbers
    
    all_names <- c(rownames(dat), colnames(dat))
    
    # Identify names that are purely numeric (can be coerced to a number without error)
    numeric_names <- all_names[!is.na(suppressWarnings(as.numeric(all_names)))]
    
                    #browser()
    
    if (length(numeric_names) > 0) {
        cat("ATTENTION: Invalid species names detected.\n",
             "  -> The following names are purely numeric: ", 
             paste(unique(numeric_names), collapse = ", "), "\n",
             "  -> Package dependencies (like 'cheddar') require proper alphanumeric names.\n",
             "  -> Please rename them (e.g., '1' -> 'Sp1', 'Species_1', or 'S01').\n",
            "-> The program will add the prefix SPS_ to the number\n\n")
        
        rownames(dat) <- paste("SPS_", names_1, sep = "")
        colnames(dat) <- paste("SPS_", names_1, sep = "")
        cat("   Modified names use 'SPS_' prefix to ensure compatibility with analysis packages.\n\n")
        
    }
    
    # Optional but good: Also warn about names with spaces or special characters
    problematic_names <- all_names[grepl("[^[:alnum:]_]", all_names)] # Matches non-alphanumeric/underscore
    if(length(problematic_names) > 0) {
        warning("  NOTE: Some names contain spaces or special characters: ", 
                #paste(unique(problematic_names), collapse = ", "), "\n",
                "  -> Consider using only letters, numbers, and underscores for compatibility.")
        
    }
    
    cat("  ✓ Species name validation passed.\n")
    cat("  ✓ Format validation passed for", nrow(dat), "species.\n")

# ARRANGE DATA FOR cheddar
    NODE <- colnames(dat)
    
    cama<- Community(nodes=data.frame(node=NODE),
                     trophic.links=PredationMatrixToLinks(dat),
                     properties=list(title="Community"))
    
# ARRENGING DATA FOR iraph
    dat_mat <- as.matrix(dat)
    gr<-graph_from_adjacency_matrix(dat_mat, weighted = FALSE, 
                                    mode = c("directed"))
    
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
        stop("There are no basal or top nodes. Can not proceed with 
             the analysis")
    } else {
        cat("All good\n")
    }
    
        return(list(file_name = file_name,
                    dat = dat,
                    cama = cama,
                    gr = gr,
                    names_1 = names_1)
               )
}

#draw modules one by one if the option is chosen
draw_modules <- function(grafica, num_mod, modularidad) {
    
    comps <- num_mod
    colbar <- rainbow(max(comps)+1, start = 0, end = max(1, comps-1)/comps)
    #colbar <- rainbow(max(comps)+1, start = 0, end = max(comps+1))
    V(grafica)$color <- colbar[comps+1]
    plot(grafica, 
         vertex.shape="sphere", 
         layout = layout_with_lgl(grafica), 
         vertex.size=5, 
         edge.arrow.size =0.1, 
         edge.curved = 0.4, 
         vertex.label.cex = 1)
    
    # loop over each community Y graphic CADA MODULO
    for (i in unique(membership(modularidad))) {
        
        # extract the nodes in this community
        nodes <- which(membership(modularidad) == i)
        
        # create a subgraph for this community
        subgraph <- induced_subgraph(grafica, nodes)
        
        #graphic de cada modulo
        
        plot(subgraph,
             vertex.color = "gray",
             vertex.shape="sphere", 
             layout = layout_nicely, 
             vertex.size=8, 
             edge.arrow.size =0.1, 
             edge.curved = 0.1, 
             vertex.label.cex = 1.2,
             vertex.label.degree = 3.1416)
        
        readline(prompt = "enter para continuar...")
        
    }
}

# DEGREE AND CENTRALITY VALUES FOR THE FOOD WEB
# 
grados <- function(gr_deg) {
    
    #RECEIVES THE FOOD WEB AS AN igraph and returns a 
    #data frame with the results
    #
    #gr_deg an igraph type object
    #file_name the file name 
    #
    #   DEGREES
    
    degree_out <- igraph::degree(gr_deg, mode = "out")
    degree_out_norm <- igraph::degree(gr_deg, mode = "out", normalized = TRUE)
    degree_out_std <- (degree_out - mean(degree_out))/sd(degree_out)
    
    centr_grado_in <- igraph::degree(gr_deg, mode = "in")
    centr_grado_in_norm <- igraph::degree(gr_deg, mode = "in", normalized = TRUE)
    centr_grados_in_std <- (centr_grado_in - mean(centr_grado_in))/sd(centr_grado_in)
    
    centr_grado_all <- igraph::degree(gr_deg, mode = "all")
    centr_grado_all_norm <- igraph::degree(gr_deg, mode = "all", normalized = TRUE)
    centr_grados_all_std <- (centr_grado_all - mean(centr_grado_all))/sd(centr_grado_all)
    
    #                   BETWEENNESS
    
    centr_bet <- igraph::betweenness(gr_deg, directed = TRUE)
    centr_bet_norm <- igraph::betweenness(gr_deg, directed = TRUE, normalized = TRUE)
    centr_bet_std <- (centr_bet - mean(centr_bet))/sd(centr_bet)
    
    
    #                   CLOSENESS
    
    centr_cerca_out <- igraph::closeness(gr_deg, mode = "out")
    centr_cerca_out_norm <- igraph::closeness(gr_deg, mode = "out", normalized = TRUE)
    centr_cerca_out_std <- (centr_cerca_out - mean(centr_cerca_out))/sd(centr_cerca_out)
    
    centr_cerca_in <- igraph::closeness(gr_deg, mode = "in")
    centr_cerca_in_norm <- igraph::closeness(gr_deg, mode = "in", normalized = TRUE)
    centr_cerca_in_std <- (centr_cerca_in - mean(centr_cerca_in))/sd(centr_cerca_in)
    
    centr_cerca_all <- igraph::closeness(gr_deg, mode = "all", normalized = FALSE)
    centr_cerca_all_norm <- igraph::closeness(gr_deg, mode = "all", normalized = TRUE)
    centr_cerca_all_std <- (centr_cerca_all - mean(centr_cerca_all))/sd(centr_cerca_all)
    
    centralidades <- data.frame(SPS = names_1, 
                                DO = degree_out,
                                DO_N = degree_out_norm,
                                DO_STD = degree_out_std,
                                DI = centr_grado_in,
                                DI_N = centr_grado_in_norm,
                                DI_STD = centr_grados_in_std,
                                D_ALL = centr_grado_all,
                                D_ALL_N = centr_grado_all_norm,
                                D_ALL_STD = centr_grados_all_std,
                                BET = centr_bet,
                                BET_N = centr_bet_norm,
                                BET_STD = centr_bet_std,
                                CLO_I = centr_cerca_in,
                                CLO_I_N = centr_cerca_in_norm,
                                CLO_I_STD = centr_cerca_in_std,
                                CLO_O = centr_cerca_out,
                                CLO_O_N = centr_cerca_out_norm,
                                CLO_O_STD = centr_cerca_out_std,
                                CLO_ALL = centr_cerca_all,
                                CLO_ALL_N = centr_cerca_all_norm,
                                CLO_ALL_STD =centr_cerca_all_std)
    
    #replace NA fpor zeros
    #
        centralidades<- replace(centralidades, is.na(centralidades), 0)

    return(centralidades)
    
}

# --- 4. CÁLCULO DE MÉTRICAS TOPOLÓGICAS POR NODO -----------------------------
# ACCORDING TO A  KIMI(?) IDEA
# CONSIDERS IN AND OUT DEGREE
# 
# 
calcular_roles <- function(grafo, membresia) {
    
    nodos <- V(grafo)
    n_nodos <- vcount(grafo)
    ids_nodos <- 1:n_nodos
    comunidades_unicas <- unique(membresia)
    n_coms <- length(comunidades_unicas)
    
    # Grados de la red dirigida (total = in + out)
    grado_in <- igraph::degree(grafo, mode = "in")
    grado_out <- igraph::degree(grafo, mode = "out")
    grado_total <- igraph::degree(grafo, mode = "all")
    
    resultados <- data.frame(
        nodo = ids_nodos,
        nombre = if(is.null(V(grafo)$name)) as.character(ids_nodos) else 
            V(grafo)$name,
        comunidad = membresia,
        grado_in = grado_in,
        grado_out = grado_out,
        grado_total = grado_total,
        stringsAsFactors = FALSE
    )
    
    # --- Within-Module Degree (z-score) ---
    z_score <- numeric(n_nodos)
    
    for (com in comunidades_unicas) {
        nodos_en_com <- which(membresia == com)
        n_nodos_com <- length(nodos_en_com)
        
        if (n_nodos_com <= 1) {
            z_score[nodos_en_com] <- 0
            next
        }
        
        k_interno <- numeric(n_nodos_com)
        
        for (i in seq_along(nodos_en_com)) {
            nodo_i <- nodos_en_com[i]
            
            vecinos_out <- neighbors(grafo, nodo_i, mode = "out")
            vecinos_out_com <- vecinos_out[membresia[vecinos_out] == com]
            
            vecinos_in <- neighbors(grafo, nodo_i, mode = "in")
            vecinos_in_com <- vecinos_in[membresia[vecinos_in] == com]
            
            k_interno[i] <- length(unique(c(as.numeric(vecinos_out_com), 
                                            as.numeric(vecinos_in_com))))
        }
        
        media_k <- mean(k_interno)
        sd_k <- sd(k_interno)
        
        if (sd_k == 0) {
            z_score[nodos_en_com] <- 0
        } else {
            z_score[nodos_en_com] <- (k_interno - media_k) / sd_k
        }
    }
    
    resultados$z_score <- z_score
    
    # --- Participation Coefficient (P) ---
    P <- numeric(n_nodos)
    
    for (i in ids_nodos) {
        if (grado_total[i] == 0) {
            P[i] <- 0
            next
        }
        
        k_por_comunidad <- numeric(n_coms)
        
        vecinos_out <- neighbors(grafo, i, mode = "out")
        if (length(vecinos_out) > 0) {
            coms_out <- membresia[vecinos_out]
            for (j in seq_along(coms_out)) {
                idx_com <- which(comunidades_unicas == coms_out[j])
                k_por_comunidad[idx_com] <- k_por_comunidad[idx_com] + 1
            }
        }
        
        vecinos_in <- neighbors(grafo, i, mode = "in")
        if (length(vecinos_in) > 0) {
            coms_in <- membresia[vecinos_in]
            for (j in seq_along(coms_in)) {
                idx_com <- which(comunidades_unicas == coms_in[j])
                k_por_comunidad[idx_com] <- k_por_comunidad[idx_com] + 1
            }
        }
        
        P[i] <- 1 - sum((k_por_comunidad / grado_total[i])^2)
    }
    
    resultados$P <- P
    
    return(resultados)
}


# --- 5. CLASIFICACIÓN DE ROLES (Guimerà & Amaral, 2005) ----------------------
# ACCORDING TO A  KIMI(?) IDEA

clasificar_rol <- function(z, P) {
    if (z < 2.5) {
        if (P < 0.05) {
            return("R1: Ultra-peripheral")
        } else if (P < 0.62) {
            return("R2: Peripheral")
        } else if (P < 0.80) {
            return("R3: Satellite connector")
        } else {
            return("R4: Kinless")
        }
    } else {
        if (P < 0.30) {
            return("R5: Provincial hub")
        } else if (P < 0.75) {
            return("R6: Connector hub")
        } else {
            return("R7: Kinless hub")
        }
    }
}


# calculo del rol de nodes and the figure of them
# 
rol_nodos <- function(red_igraph, leiden_dat) {
    # 
    # Asignar membresía a los nodos
    
    #membresia_leiden <- membership(leiden_dat)
    V(red_igraph)$comunidad <- leiden_dat
    n_comunidades <- length(unique(leiden_dat))
    
    cat("Calculando métricas de roles topológicos...\n")
    roles_df <- calcular_roles(red_igraph, leiden_dat)
    roles_df$rol <- mapply(clasificar_rol, roles_df$z_score, roles_df$P)
    roles_df$rol_simple <- gsub("R[0-9]: ", "", roles_df$rol)
    
    cat("\n--- DISTRIBUCIÓN DE ROLES ---\n")
    tabla_roles <- table(roles_df$rol)
    print(tabla_roles)
    
    # --- 6. AÑADIR ATRIBUTOS AL GRAFO -------------------------------------------
    
    V(red_igraph)$z_score <- roles_df$z_score
    V(red_igraph)$P <- roles_df$P
    V(red_igraph)$rol <- roles_df$rol_simple
    
    # --- 7. VISUALIZACIONES ------------------------------------------------------
    
    colores_roles <- brewer.pal(n = 7, name = "Set1")
    names(colores_roles) <- c("Ultra-peripheral", "Peripheral", 
                              "Satellite connector",
                              "Kinless", "Provincial hub", "Connector hub", 
                              "Kinless hub")
    
    roles_presentes <- unique(roles_df$rol_simple)
    colores_usar <- colores_roles[roles_presentes]
    
    p1 <- ggplot(roles_df, aes(x = P, y = z_score, color = rol_simple)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_hline(yintercept = 2.5, linetype = "dashed", color = "black", 
                   linewidth = 0.7) +
        geom_vline(xintercept = c(0.05, 0.62, 0.80), linetype = "dotted", 
                   color = "black", linewidth = 0.7) +
        scale_color_manual(values = colores_usar, name = "Roll") +
        labs(
            title = "Nodes Topological Rolls",
            x = "Participation Coefficient (P)",
            y = "Within-Module Degree (z-score)"
        ) +
        theme_few(base_size = 15) + #  theme_minimal(base_size = 15) +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold", size = 14),
            panel.grid.minor = element_blank(),
            #panel.grid.major = element_blank()
        )
    
    print(p1)
    
    return(roles_df)
}

# FUNCTION TO GENERATE A FOOD WEB ACCORDING TO THE CASCADE MODEL
# Cohen, J. E., F. Briand, and C. M. Newman. 1990a. Community Food Webs:
#Data and Theory. New York: Springer-Verlag.
#
#   S = number of species
#   C = connectance
#   tol = tolerance of the connectance
#   N = number of food webs to be generated
#   
#   returns a list of food webs

Web.CascadeModel <- function(S, C, tol, names_1, ...) {

    # Connectance interval
    C.min <- C - tol
    C.max <- C + tol
    
    WebCascadeModel <- list()
        cat(". . . WORKING CASCADE MODEL . . .\n")
        
        iter <- 0
        
        while (iter < 10000) {
            iter <- iter + 1
            
            Web.Cascade <- matrix(0, S, S, dimnames = list(paste("S", 1:S, sep = ""), paste("S", 1:S, sep = "")))
                                                                 
            P <- 2 * C * S/(S - 1)
            r.value <- runif((S^2 - S)/2)
            
            # Adjacency matrix
            Web.Cascade[upper.tri(Web.Cascade)][r.value < P] = 1
            
            # Links
            L <- sum(Web.Cascade)
            
            # Connectance
            C.value = L/S^2
            
            # Connected
            connected.nodes <- sna::is.connected(Web.Cascade, connected = "weak") #Packages:sna
            
            # Loops
            net <- network(Web.Cascade) #Packages:network
            loop.value <- has.loops(net)
            
            #quitar loops
            #
            
            if (loop.value == FALSE) {
                #LEVINES METHOD FOR TROPHIC LEVEL
                #checar si la matriz es singular
                
                A <- t(Web.Cascade)
                
                # Assuming you have your food web matrix 'M' (the normalized flow matrix)
                # or adjacency matrix 'A' from which you'll create M
                
                # Create the matrix needed for Levine's method
                n <- nrow(A)  # number of species/size of your matrix
                
                # First, create matrix M (normalized by row sums)
                row_sums <- rowSums(A)
                # Handle division by zero for basal species (rows with no prey)
                row_sums[row_sums == 0] <- 1  
                M <- A / row_sums
                
                # Create the identity matrix
                I <- diag(n)
                
                # The key matrix for Levine's method
                levine_matrix <- I - M
                
                # Check if it's singular by looking at the determinant
                det_value <- det(levine_matrix)
                
                if(abs(det_value) < 1e-10) {  # Using a small tolerance for numerical precision
                    #cat("Matrix is singular or nearly singular - trophic levels will be NA\n")
                    Singular <- "TRUE"
                } else {
                    #cat("Matrix is invertible - trophic levels should be computable\n")
                    Singular <- "FALSE"
                    cat("Determinant of (I - M):", det_value, "\n")
                    
                }
                
            }
            if (C.value > C.min & C.value < C.max & loop.value == FALSE &
                connected.nodes == TRUE & Singular =="FALSE") {
                iter = 10000
            }
        }
        # rownames(Web.Cascade) <- paste("SPS", names_1, sep = "")
        # colnames(Web.Cascade) <- paste("SPS", names_1, sep = "")
        rownames(Web.Cascade) <- names_1
        colnames(Web.Cascade) <- names_1
        WebCascadeModel<-graph_from_adjacency_matrix(Web.Cascade, weighted = FALSE, mode = c("directed"))

    return(WebCascadeModel)
}

#   FUNCTION TO GENERATE FOOD WEBS ACCORDING TO THE NICHE MODEL
#   
#   S = number of species
#   C = connectance
#   tol = tolerance of the connectance
#   N = number of food webs to be generated
#   
#   returns a list of food webs

Web.NicheModel <- function(S, C, tol, names_1, ...) {
    #number of singular matrices generated
    #no_invertibles <- 0
    # Connectance interval
    C.min <- C - tol
    C.max <- C + tol
    
    WebNicheModel <- list()

        iter <- 0
        
        cat("WORKING...NICHE MODEL...\n")
        while (iter < 10000) {
            iter <- iter + 1
            
            Web.Niche <- matrix(0, S, S, dimnames = list(paste("S", 1:S, sep = ""), paste("S", 1:S, sep = "")))
                                                               
            n <- sort(runif(S))
            b <- 1/(2 * C) - 1
            x <- 1 - (1 - runif(S))^(1/b)
            r <- n * x
            r[1] <- 0
            center <- numeric()
            for (i in 1:S) {
                center[i] <- runif(1, r[i]/2, min(n[i], 1 - r[i]/2))
            }
            getMin <- center - r/2
            getMax <- center + r/2
            
            # Adjacency matrix
            for (i in 1:S) {
                Web.Niche[c(1:S)[n > getMin[i] & n < getMax[i]], i] <- 1
            }
            
            # Links
            L <- sum(Web.Niche)
            
            # Connectance
            C.value = L/S^2
            
            # Connected
            connected.nodes <- sna::is.connected(Web.Niche, connected = "weak")
            
            # Loops
            net <- network(Web.Niche)
            loop.value <- has.loops(net)
            
            #quitar loops
            #
            
            if (loop.value == FALSE) {
                #LEVINES METHOD FOR TROPHIC LEVEL
                #checar si la matriz es singular
                
                A <- t(Web.Niche)
                
                # Assuming you have your food web matrix 'M' (the normalized flow matrix)
                # or adjacency matrix 'A' from which you'll create M
                
                # Create the matrix needed for Levine's method
                n <- nrow(A)  # number of species/size of your matrix
                
                # First, create matrix M (normalized by row sums)
                row_sums <- rowSums(A)
                # Handle division by zero for basal species (rows with no prey)
                row_sums[row_sums == 0] <- 1  
                M <- A / row_sums
                
                # Create the identity matrix
                I <- diag(n)
                
                # The key matrix for Levine's method
                levine_matrix <- I - M
                
                # Check if it's singular by looking at the determinant
                det_value <- det(levine_matrix)
                
                if(abs(det_value) < 1e-10) {  # Using a small tolerance for numerical precision
                    #cat("Matrix is singular or nearly singular - trophic levels will be NA\n")
                    Singular <- "TRUE"
                    #no_invertibles <- no_invertibles + 1
                } else {
                    #cat("Matrix is invertible - trophic levels should be computable\n")
                    Singular <- "FALSE"
                    cat("Determinant of (I - M):", det_value, "\n")
                    
                }
            
            }
            
            if (C.value > C.min & C.value < C.max & loop.value == FALSE &
                connected.nodes == TRUE & Singular == "FALSE") {
                iter = 10000
            }
        }
        # rownames(Web.Niche) <- paste("SPS", names_1, sep = "")
        # colnames(Web.Niche) <- paste("SPS", names_1, sep = "")
        rownames(Web.Niche) <- names_1
        colnames(Web.Niche) <- names_1
        WebNicheModel <- graph_from_adjacency_matrix(Web.Niche, 
                                                     weighted = FALSE, 
                                                     mode = c("directed"))
        
    return(WebNicheModel)
}

#       FUNCTION TO COMPUTE FOOD WEB STRUCTURE
#       
#       g_rand = an igraph type object of the food web
#       tl_y_or_no = compute "si" on not "no" the trophic level for each node
#           If the food web is very large (S>150 and L > 2000) chedar sends a warning
#           that it will be imposible to compute the trophic level. 
#           The option  chedarMaxQueue = 0 circunvents this problem, never the less
#           if the size of the food web is large and not enough RAM memory is
#           available the program will crash.
#       
#       returns a list with the food web parameters
#       

fw_struct_2 <- function(g_rand, tl_y_or_no, names_1) {
    
    #browser()
    cat("\n")
    cat("\n")
    cat(". . . w o r k i n g. . . ")
    cat("\n")
    num_sps <- vcount(g_rand)
    num_links <- ecount(g_rand)
    conectance <- num_links / (num_sps ^2)
    
    DAG<-is_dag(g_rand)
    
    #Eliminar ciclos del grafo
    if(DAG ==FALSE){
        edges_cycle<-feedback_arc_set(g_rand, algo="approx_eades")
        g_rand<-delete_edges(g_rand, edges_cycle)
    }
    #   calcula el average path length
    
    mean_path_len <- mean_distance(g_rand, directed = TRUE)
    g_tempo <- as.matrix(as_adjacency_matrix(g_rand))
    rownames(g_tempo) <- names_1
    colnames(g_tempo) <- names_1
    
    
    NODE <- rownames(g_tempo)
    
    commty <- Community(nodes = data.frame(node=NODE),
                        trophic.links=PredationMatrixToLinks(g_tempo),
                        properties=list(title="Community"))
    #CLEAR g_tempo...RAM amount IS IMPORTANT AT THIS STAGE
    #rm(g_tempo) 
    
    if (tl_y_or_no == "YES") {
        
        final_results <- TL_Compute(g_rand, commty)
        chain.stats <- 0   #TrophicChainsStats(commty)
        chains_nom <- 0   #length(chain.stats$chain.lengths)
        # rm(chain.stats) #remove this object to get RAM
        # 
        # # 3. METRIC 1: PREY-AVERAGED TL (The native cheddar version)

        tl_mean <- mean(final_results$ShortWeightedTL)
        maxTL <- max(final_results$ShortWeightedTL)
        
    }
    
    #Especies BSALES
    b <- length(BasalNodes(commty))

    #Especies INTERMEDIAS
    int <- length(IntermediateNodes(commty))

    #Especies TOPE
    tope <- length(TopLevelNodes(commty))

    #"Especies OMNIVORAS
    omni <- length(Omnivores(commty))

    #vulnerability and generality normalized by L. use std for comparissons
    
    vul <- NormalisedTrophicVulnerability(commty)
    
    gen <- NormalisedTrophicGenerality(commty)
    
    return(list(N = num_sps, 
                L = num_links, 
                con = conectance, 
                B = b, 
                I = int, 
                T = tope, 
                O = omni, 
                TLMean = tl_mean,
                TLmax = maxTL, 
                NCh = chains_nom, 
                mean_path_length = mean_path_len,
                nivel_trof_sps = final_results$ShortWeightedTL,
                vul_std = sd(vul),
                gen_std = sd(gen)))
}

#   FUNCION draw_1 RED
#   FUNCTION TO DRAW A FOOD WEB
#   graph_data = an igraph object of the food web
#   
draw_1 <- function(gr, niv_trof_sps, fwt_results_dir) {
    #browser()
    # ==========================================
    # 2. PREPARE THE GRAPH WITH TROPHIC LEVELS
    # ==========================================
    # Add the calculated ShortWeightedTL as an attribute to our igraph object nodes
    V(gr)$TrophicLevel <- niv_trof_sps$TL[match(V(gr)$name, niv_trof_sps$SPS)]
    
    # Create a clean layout: 
    # X = random spread so nodes don't overlap horizontally
    # Y = the exact Trophic Level so energy flows upwards
    custom_layout <- matrix(NA, nrow = vcount(gr), ncol = 2)
    custom_layout[, 1] <- runif(vcount(gr), min = 0, max = 10)  # Horizontal spread
    custom_layout[, 2] <- V(gr)$TrophicLevel                  # Vertical Trophic Level
    
    # ==========================================
    # 3. CONVERT TO GGNETWORK FORMAT
    # ==========================================
    n <- ggnetwork::ggnetwork(gr, layout = custom_layout)
    
    # ==========================================
    # 4. PLOT THE FOOD WEB
    # ==========================================
    p <- ggplot2::ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
        # Draw the feeding links (light gray with transparency so it's not overwhelming)
        geom_edges(color = "grey70", alpha = 0.4, 
                   arrow = arrow(length = unit(6, "pt"), type = "closed")) +
        # Draw the species nodes (colored by their Trophic Level)
        geom_nodes(aes(color = TrophicLevel), size = 4) +
        # Styling the color palette (Darker colors for basal, lighter for apex predators)
        scale_color_viridis_c( option = "plasma") +
        # Clean up the background layout
        theme_blank() +
        labs(
            title = paste("Food Web Structure ", fwt_results_dir),
            subtitle = paste("Nodes:", vcount(gr), " | Links:", ecount(gr)),
            y = "Short-Weighted Trophic Level"
        ) +
        # Ensure the Y-axis acts as a true scale for your Trophic Levels
        # name = "Short Weighted TL",
        theme(
            axis.title.y = element_text(angle = 90, vjust = 2, size = 12),
            plot.title = element_text(size = 14, face = "bold")
        )
    
    print(p)

    return()
}

#FUNCTION TO CALCULATE MODULARITY ACCORDING TO THE Leiden ALORITHM
#matriz = igraph object of the food web
#resolucion = resolution to be used by the algorithm

modul_l <- function(matriz, resolucion) {
    
    resulta_l <- leiden.community(matriz, resolution = resolucion, 
                                   n.iterations = 1000)
    
    membership.fakeCommunities(resulta_l)
    miembros_azar <- resulta_l$membership
    miembros_azar <- as.numeric(miembros_azar)
    
    no_modulos_azar <- max(miembros_azar)

    return(resulta_l)
}

#  FUNTCION TO CALCULATE THE transitivity or clustering of a food web
#  g_trans = an igraph object of the food web
#  returns the transitivity value
transi <- function(g_trans) {
    trans <- transitivity(g_trans, type = "global")
    return(trans)
}

#FUNCTION TO GENERATE A FOOD WEB ACCORDING TO THE Erdös-Renyi model
# gr = an igraph object of the food web
# graphic = ?
# names_1 = node names
# 
# returns an igraph object

er <- function(gr, names_1) {
    
    v <- gsize(gr)
    n <- vcount(gr)
    # browser()
    #   DEFINE VARIABLES
    q_2 <- 1
    basal <- 0
    top <- 0
    out_loop <- 0  #FALSO
    conected <- "FALSE"
    
    #graph_data DUMMY to begin with the while loop 
    
    g1 <- sample_gnm(5, 6, directed = TRUE)

    while (conected == "FALSE") {
        
        basal <- 0
        top <- 0
        out_loop <- 0
        
        g_rand_erdos <- sample_gnm(n, v, directed = TRUE)
        
        # 
        #prueba ciclos
        DAG<-is_dag(g_rand_erdos)
        #Eliminar ciclos del grafo
        if(DAG == FALSE){
            edges_cycle<-feedback_arc_set(g_rand_erdos, algo="approx_eades")
            g_rand_erdos<-delete_edges(g_rand_erdos, edges_cycle)
            # Ex<-as.matrix(as_adjacency_matrix(gr))
        }

        if (conected == "TRUE") {
            #suma columnas para checar basal si es diferente de zero
            pred_mat_rand <- as_adjacency_matrix(g_rand_erdos)
            pred_mat_rand <- as.matrix(pred_mat_rand)
            q_1 <- apply(pred_mat_rand, 2, sum)
            basal <- sum(q_1 == 0)
            print(basal)
            q_3 <- apply(pred_mat_rand, 1, sum)
            top <- sum(q_3 == 0)
            print(top)
            
            if (basal == 0 && top == 0){
                print("basal Y top = 0")
            } else {
                print("ni top ni basal = 0")
                out_loop <- 1  #verdadero
            }
            if (out_loop == 1) {
                
                #   sale de todo el loop porque ya es conected y las top 
                #   y basal si existen
                
                break
            }
        }
        
        g_tempo <- as.matrix(as_adjacency_matrix(g_rand_erdos))
        # rownames(g_tempo) <- paste("SPS", names_1, sep = "")
        # colnames(g_tempo) <- paste("SPS", names_1, sep = "")
        
        rownames(g_tempo) <- names_1
        colnames(g_tempo) <- names_1
        g_rand_erdos<-graph_from_adjacency_matrix(g_tempo, 
                                                  weighted = FALSE, 
                                                  mode = c("directed"))
        
        conected <- igraph::is_connected(g_rand_erdos)
        
    }
    
    return(g_rand_erdos)
}

# 
# FUNCTION TO OBTAIN NODES THAT PERTAIN TO TWO OR MODE MODULES
# 
# g_cruz = THE IGRAPH OBJECT OF THE FOOD WEB
# modulos_leiden = the resulting modules from a modularity analysis
# 
# returns a data frame with the pairs of nodes that are present in the modules
# 
intersected <- function(g_cruz, modulos_leiden) {
    
    gr_nodir <- as_undirected(g_cruz)   #conviere a no direccional
    
    cruza <- igraph::crossing(modulos_leiden, gr_nodir) #nodos pares en mas de dos modulos
    cruza
    cruza_df <- as.data.frame(cruza)
    cruza_df
    
    return(cruza_df)
    
}

#FUNCTION TO CALCULATE THE VARIOUS KEYSTONE NUMBERS AS PRESENTED BY

# Cite: Jordán, F., Takács-Sánta, A., & Molnár, I. (1999). A reliability theoretical 
# quest for keystones. Oikos, 453-462.

# Ex: Adjacent matrix
# RETURNS A DATA FRAME WITH THE VALUES FOR EACH VARIABLE OF KEYSTONE

k.parameter<-function(Ex, ...){
    
    #Matriz binaria y eliminar diagonal
    diag(Ex)<-0
    Ex[Ex>0]<-1
    
    #Crear un grafo (igraph)
    require(igraph)
    gr<-graph_from_adjacency_matrix(as.matrix(Ex))
    DAG<-is_dag(gr)
    
    #Eliminar ciclos del grafo
    if(DAG ==FALSE){
        edges_cycle<-feedback_arc_set(gr, algo="approx_eades")
        gr<-delete_edges(gr, edges_cycle)
        Ex<-as.matrix(as_adjacency_matrix(gr))
    }
    
    #crear graph_data para guardar salidas
    K<-matrix(-1,ncol=1, nrow=nrow(Ex))
    rownames(K)<-rownames(K)<-rownames(Ex)
    Kdir<-matrix(0,ncol=1, nrow=nrow(Ex))
    rownames(Kdir)<-rownames(Kdir)<-rownames(Ex)
    Kindir<-matrix(0,ncol=1, nrow=nrow(Ex))
    rownames(Kindir)<-rownames(Kindir)<-rownames(Ex)  
    
    #Identificar especies tope                                      
    row.sum<-apply(Ex,1,sum)
    tope<-names(row.sum[row.sum==0])
    K[tope,]<-0
    
    #Iniciar el cálculo con especies tope
    Nodes<-names(K[,1][K[,1]<0])
    
    #Calcular el indicador con especies faltantes
    cond<-TRUE
    while (cond==TRUE) {
        #Buscar si el indicador de los depredadores fue calculado
        for(l in Nodes){
            if(K[l,1]<0){
                top<-names(Ex[l,])[Ex[l,]==1]
                cond2<-TRUE
                for(m in top){
                    if(K[m,1]<0){
                        cond2<-FALSE
                    }
                }
                #Calcular el indicador de una especie cuando exista la informacion de los depreadores
                if(cond2==TRUE){
                    Acum2<-0
                    Acum3<-0
                    for(n in top){
                        Acum2<-Acum2+(1/sum(Ex[,n]))
                        Acum3<-Acum3+(K[n,1]/sum(Ex[,n]))
                    }
                    K[l,1]<-Acum2+Acum3
                    Kdir[l,1]<-Acum2
                    Kindir[l,1]<-Acum3
                }
            }
        }
        #Terminar cuando el indicador se calcule para todas las especies
        if(length(K[K<0])==0){
            cond<-FALSE
        }
        Nodes<-names(K[,1])[K[,1]<0]
    }
    
    return(data.frame(Kdir,Kindir,K))
}

#   FUNCTION TO COMPUTE STRUCTURAL PROPERTIES FOR EACH MODULE OF THE FOOD WEB
#   groups = modules and nodes
#   g_rand = an igraph object of the food web
#   graphic = if a plot for each module is perfomed
#
prop_modules <- function(groups, g_rand, graphic) {
    
    module_str <- list()
    Ex <- NULL

    # loop over each community and compute some graph properties
    for (i in unique(groups)) {
        
        print(i)
        #browser()
        
        # extract the nodes in this community
        nodulos <- which(groups == i)
        
        # create a subgraph for this community
        subgraph <- induced_subgraph(g_rand, nodulos)
        
        #quitar todos los ciclos de la red para poder realizar los alculos 
        #con cheddar como chains_nom y nivel trofico
        #
        DAG<-is_dag(subgraph)
        # 
        # #Eliminar ciclos del grafo
        if(DAG ==FALSE){
            edges_cycle<-feedback_arc_set(subgraph, algo="approx_eades")
            subgraph <- delete_edges(subgraph, edges_cycle)
            Ex_adj<-as.matrix(as_adjacency_matrix(subgraph))
            NODE <- rownames(Ex_adj)
            Ex<-graph_from_adjacency_matrix(Ex_adj)
            
        } else {
            
            Ex_adj<-as.matrix(as_adjacency_matrix(subgraph))
            #Ex<-graph_from_adjacency_matrix(Ex_adj)
            NODE <- rownames(Ex_adj)
            Ex<-graph_from_adjacency_matrix(Ex_adj)
            
        }
        
        sps <- vcount(Ex)
        num_links <- ecount(Ex)
        conec_modulo <- num_links / (sps ^2)
        
        commty <- Community(nodes = data.frame(node=NODE),
                            trophic.links=PredationMatrixToLinks(Ex_adj),
                            properties=list(title="Community"))

        #opcion paraa que cheddar compute sin importar el tamaño de 
        #la red los niveles troficos
        options(cheddarMaxQueue = 0)
        
        
        chain.stats_modulo <- TrophicChainsStats(commty)
        
        #pone cero para muy grandes
        #chain.stats_modulo <- 0
        
        NivTrophModulo <-ShortWeightedTrophicLevel(commty)
        NivTrophModulo <- as.data.frame(NivTrophModulo)
        NivTrophModulo
        tl_mean_modulo <- mean(NivTrophModulo$NivTrophModulo)
        maxTL <- max(NivTrophModulo$NivTrophModulo)
        
        #print("Numero de chains_nom:")
        chains_nom <- length(chain.stats_modulo$chain.lengths)
        #print(chains_nom)
        
        #print("Especies BSALES")
        basal_modulo <- BasalNodes(commty)
        b <- length(basal_modulo)
        #print(b)
        
        #print("Especies INTERMEDIAS")
        interm_modulo <- IntermediateNodes(commty)
        int <- length(interm_modulo)
        # print(int)
        
        #print("Especies TOPE")
        top_modulo <- TopLevelNodes(commty)
        tope <- length(top_modulo)
        #print(tope)
        
        #print("Especies OMNIVORAS")
        omni_modulo <- Omnivores(commty)
        omni <- length(omni_modulo)
        #print(omni)
        
        modulo <- data.frame(N = sps, L = num_links, 
                             con = conec_modulo, B = b, I = int, 
                             T = tope, O = omni, TLm = tl_mean_modulo,
                             TLmax = maxTL, NCh = chains_nom)
        
        module_str[[paste0("modulo", i)]] <- modulo
        
        if (graphic == "si")  { 
            
            # 
            # #graphic de cada modulo
            # 
            plot(subgraph,
                 vertex.color = "gray",
                 vertex.shape="sphere",
                 layout = layout_nicely,
                 vertex.size=8,
                 edge.arrow.size =0.1,
                 edge.curved = 0.1,
                 vertex.label.cex = 1.2,
                 vertex.label.degree = 3.1416)
            
            readline(prompt = "enter to continue...")
            
        }
        
    }
    
    return(module_str)
}

#GENERATE RANDOM FOOD WEBS ACCORDING TO THE GIVEN ALGORITHM
#gr = original food web as an igraph object
#names_1 = names of the nodes
#num_rand_webs 0 number of random food webs to be generated
#random_model_type = name of the random model algorithm
#RETURN A LIST WITH THE RANDOM FOOD WEBS
#
gen_rnd_fw <- function(gr, names_1, num_rand_webs, random_model_type, 
                       fwt_results_dir_rand_webs) {

#browser()
#

    # =====================================================
    # SET UP RUN FOLDER AND LOG FILE
    # =====================================================
    
    # Build a unique run ID from algorithm + timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run_id <- paste0("run_", timestamp, "_", random_model_type)
    
    # Create the run folder inside a "runs" directory
    runs_root <- fwt_results_dir_rand_webs #"/RESULTS/RANDOM/RANDOMIZED_WEBS"
    
    run_dir <- file.path(runs_root, run_id)
    if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)
    
    cat("📁 Run folder created:", run_dir, "\n")
    

        if (random_model_type == "niche-model") {
            
            #modelo de nicho
            sps <- vcount(gr)
            num_links <- ecount(gr)
            conec <- num_links / (sps ^2)
            #tolerancia es conectividad * porcentaje de esa conectividad, 
            #Williams & Martinez (2000) mencionan este intervalo, + - el 3% de la C observada
            tolerancia <- conec * 0.03
            #   LLAMA A LA FUNCION PARA GENERAR vueltas NÚMERO DE REDES CON
            #   CIRTA TOLERANCIA PARA CONECTIVIDAD Y CIERTO NUMERO DE ESPSECIES
            #   EL RESULTADO ES UNA LISTA CON TODAS LAS graph_data GENERADAS

            #rand_matrix <- Web.NicheModel(sps, conec, tolerancia, num_rand_webs, names_1)
            rand_matrix <- lapply(1:num_rand_webs, function(i) 
                Web.NicheModel(sps, conec, tolerancia, names_1)) 
        }
        
        else {
            
            if (random_model_type == "cascade") {
                
                #modelo de cascada
                
                sps <- vcount(gr)
                num_links <- ecount(gr)
                conec <- num_links / (sps ^2)
                #tolerancia es conectividad * porcentaje de esa conectividad, 
                #Williams & Martinez (2000) mencionan este intervalo, + - el 3% de la C observada
                tolerancia <- conec * 0.03
                
                rand_matrix <- lapply(1:num_rand_webs, function(i) 
                    Web.CascadeModel(sps, conec, tolerancia, names_1))
                
            }
        else {
            
            if (random_model_type == "niche_allesina") {
                
                #niche model modified by Allesina et al

                sps <- vcount(gr)
                num_links <- ecount(gr)
                conec <- num_links / (sps ^2)
                tolerancia <- 0
                
                #rand_matrix <- niche_model_2(sps, conec, num_rand_webs)
                #browser()
                rand_matrix <- lapply(1:num_rand_webs, function(i) 
                    niche_model_2(sps, conec, names_1))
            }
               
        else {
            if (random_model_type == "erdos-renyi") {
                
                sps <- vcount(gr)
                num_links <- ecount(gr)
                conec <- num_links / (sps ^2)
                #tolerancia es conectividad * porcentaje de esa conectividad, 
                #Williams & Martinez (2000) mencionan este intervalo, + - el 3% de la C observada
                tolerancia <- 0
            #browser()
            #erdos-renyi
                rand_matrix <- lapply(1:num_rand_webs, function(i) er(gr, names_1))
            
            }
        else {
            
            
            if(random_model_type == "random_links") {
                
                sps <- vcount(gr)
                num_links <- ecount(gr)
                conec <- num_links / (sps ^2)
                tolerancia <- 0
                
                rand_matrix <- lapply(1:num_rand_webs, function(i) 
                                      random_links(gr))
            }
        }
        } 
        }
        }
    
        # Save if save_dir is provided
        
        if (!is.null(run_dir)) {
            # Create folder if it doesn't exist
            if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)
            
            # Save each web
            for (i in seq_along(rand_matrix)) {
                saveRDS(
                    rand_matrix[[i]], 
                    file = file.path(run_dir, paste0("_", random_model_type, i, ".rds"))
                )
            }
            cat("✅ Saved", length(rand_matrix), "random webs to:", run_dir, "\n")
        }

    # =====================================================
    # WRITE LOG FILE
    # =====================================================
    
    log_lines <- c(
        "==============================",
        "FWTopo Random Web Generation Log",
        "==============================",
        paste0("Run ID:          ", run_id),
        paste0("Date:            ", Sys.Date()),
        paste0("Time:            ", format(Sys.time(), "%H:%M:%S")),
        paste0("Algorithm:       ", random_model_type),
        paste0("Number of webs:  ", num_rand_webs),
        paste0("Species (S):     ", sps),
        paste0("Connectance (C): ", conec),
        paste0("Tolerance:       ", tolerancia),
        paste0("Output folder:   ", run_dir),
        "=============================="
    )
    
    writeLines(log_lines, con = file.path(run_dir, "run_log.txt"))
    cat("📝 Log file written to:", file.path(run_dir, "run_log.txt"), "\n")

    return(rand_matrix)
    
}

#COMPUTE THE MODIFIED NICHE MODEL (ALLESINA, ALONSO, PASCUAL 2008)
niche_model_2 <- function(sps, conec, names_1) {
    
        tempo <- create_niche_model(sps, conec) #FUNCTION FROM PACKAHE ATNr

        cat("FW-NICHE-MODEL-ALLESINA\n")

        # rownames(tempo) <- paste("SPS", names_1, sep = "")
        # colnames(tempo) <- paste("SPS", names_1, sep = "")
        rownames(tempo) <- names_1
        colnames(tempo) <- names_1
        g_rand_niche_2<-graph_from_adjacency_matrix(tempo, 
                                                    weighted = FALSE, 
                                                    mode = c("directed"))

        return(g_rand_niche_2)
        
}

# Display available algorithms
display_algorithm_menu <- function() {
    cat("\n")
    cat("  _________________________________________________\n")
    cat("               Food Web Topology Analysis\n")
    cat("  _________________________________________________\n")
    
    cat("=================================\n")
    cat("Available Algorithms\n")
    cat("=================================\n")
    for (i in 1:length(algorithm_map)) {
        cat(sprintf("%d. %s\n", i, algorithm_map[i]))
    }
    cat("=================================\n")
    
}

# Get user input with validation
get_algorithm_choice <- function() {
    display_algorithm_menu()
    
    while (TRUE) {
        choice <- readline(prompt = "Please select an algorithm (1-5): ")
        
        # Validate input
        if (choice %in% names(algorithm_map)) {
            selected_algorithm <- algorithm_map[choice]
            cat(sprintf("\n✅ Selected: %s\n", selected_algorithm))
            return(selected_algorithm)
        } else {
            cat("❌ Invalid choice. Please enter a number between 1 and 4.\n")
        }
    }
}

#VALIDATION REPORT
#
# Add this right after reading your adjacency matrix
generate_validation_report <- function(adj_matrix, filename, algo_rnd, 
                                       numb_fws, resol, tiempo, directorio) {
    cat("🔬 Generating FWTopo Validation Report...\n")
    
    # Basic validation
    report <- list(
        title = "....FOOD WEB TOPOLOGY ANALYSIS...",
        title_2 = "         V 1.2",
        input_file = filename,
        analysis_date <- Sys.Date(),
        dimensions = dim(adj_matrix),
        is_square = nrow(adj_matrix) == ncol(adj_matrix),
        total_interactions = sum(adj_matrix),
        connectance = round(sum(adj_matrix) / nrow(adj_matrix)^2, 3),
        top_species = sum(rowSums(adj_matrix) == 0),
        basal_predators = sum(colSums(adj_matrix) == 0),
        is_connected = igraph::is_connected(igraph::graph_from_adjacency_matrix(adj_matrix)),
        algorithm = algo_rnd,
        number_of_fws = numb_fws,
        Leiden_resolution = resol,
        Time_spent = tiempo
    )
    #browser()
    # Save validation report
    validation_file <- paste0(tools::file_path_sans_ext(filename), "_LOG.txt")
    validation_file <- paste0(directorio, "/", format(Sys.time(), "%Y.%m.%d_%H.%M.%S_"), validation_file)
    capture.output(print(report), file = validation_file)
    
    cat("✅ Validation report saved as:", validation_file, "\n")
    return(report)
}

fw_struct_rnd <- function(g_rand, tl_y_or_no, names_1) {
    #function to compute topology of the random food web one by one
    #g_rand = random food web as an igraph object
    #tl_y_or_no do we compute trophic level
    #names_1 = names of nodes
    
   # browser()
    
   # 
   #cheddar option to run regardless of the number of chains
    options(cheddarMaxQueue = 0)
    
    num_sps <- vcount(g_rand)
    num_links <- ecount(g_rand)
    conectance <- num_links / (num_sps ^2)
        cat("Number of node: \n")
        print(num_sps)
        cat("Number of links: \n")
        print(num_links)
    #   calcula el average path length
    mean_path_len <- mean_distance(g_rand, directed = TRUE)
    g_tempo <- as.matrix(as_adjacency_matrix(g_rand))
    rownames(g_tempo) <- names_1
    colnames(g_tempo) <- names_1
    
    NODE <- rownames(g_tempo)
    
    commty <- Community(nodes = data.frame(node=NODE),
                        trophic.links=PredationMatrixToLinks(g_tempo),
                        properties=list(title="Community"))
    
    if (tl_y_or_no == "YES") {
        # chain.stats <- TrophicChainsStats(commty)
        # chains_nom <- length(chain.stats$chain.lengths)
        #     cat("Number of chains: \n")
        #     print(chains_nom)
        
        cat("COMPUTING TROPHIC LEVELS OF TE RANDOM WEB \n\n")
        
        final_results <- TL_Compute(g_rand, commty)
        
        tl_mean <- mean(final_results$ShortWeightedTL)
        maxTL <- max(final_results$ShortWeightedTL)
    } else {
        tl_mean <- 0
        maxTL <- 0
        chains_nom <- 0
        
    }
    #Especies BSALES
    b <- length(BasalNodes(commty))
    #Especies INTERMEDIAS
    int <- length(IntermediateNodes(commty))
    #Especies TOPE
    tope <- length(TopLevelNodes(commty))
    #"Especies OMNIVORAS
    omni <- length(Omnivores(commty))

    #vulnerability and generality normalized by L. use std for comparissons
    
    vul <- NormalisedTrophicVulnerability(commty)
    
    gen <- NormalisedTrophicGenerality(commty)
    
    efic <- eficiency(commty)
    
    res_est_rnd_1 <- data.frame(N = num_sps, 
                L = num_links,
                con = conectance, 
                B = b, 
                I = int, 
                T = tope, 
                O = omni, 
                TLMean = tl_mean,
                TLmax = maxTL, 
                NCh = 0, #chains_nom, 
                mean_path_length = mean_path_len,
                vul_std = sd(vul),
                gen_std = sd(gen),
                effic = efic)
    
    return(res_est_rnd_1)
    
}

#FUNCTION FOR THE PLOTTING OPTION
#
display_figura_menu <- function() {
    cat("\n")
    cat("=====================================\n")
    cat("FWTopo: Plotting the Food Web Option\n")
    cat("=====================================\n")
    
    cat("=====================================\n")
    
}

display_figura_tldensity_menu <- function() {
    cat("\n")
    cat("=====================================\n")
    cat("FWTopo: Plotting the TL Random Density Option\n")
    cat("=====================================\n")
    
    cat("=====================================\n")
    
}



# Get user input with validation FOR THE PLOTTING OPTION
get_figura_choice <- function() {
    display_figura_menu()
    
    while (TRUE) {
        choice <- readline(prompt = "Please enter 1 = YES or 2 = NO : ")
        
        # Validate input
        if(choice == "1" | choice == "2" ) {
            return(choice)
        } else {
            
            cat("Enter YES or NO \n")
        }
    }
}


# Get user input with validation FOR THE PLOTTING OPTION
get_figura_tldensity_choice <- function() {
    display_figura_tldensity_menu()
    
    while (TRUE) {
        choice <- readline(prompt = "Please enter 1 = YES or 2 = NO : ")
        
        # Validate input
        if(choice == "1" | choice == "2" ) {
            return(choice)
        } else {
            
            cat("Enter YES or NO \n")
        }
    }
}



#MENU FOR THE OPTION OF THE LEIDEN ALGORITHM
#VALUE HAS TO BE GREATER THAN 0
#
display_resolution_menu <- function() {
    cat("\n")
    cat("=====================================\n")
    cat("FWTopo: Leiden Algorithm Resolution\n")
    cat("=====================================\n")
    
    cat("=====================================\n")
    
}

# Get user input with validation
# FOR THE LEIDEN RESOLUTION OPTION
get_resolution_choice <- function() {
    display_resolution_menu()
    
    while (TRUE) {
        choice <- readline(prompt = "Please enter a value : ")
        
        # Validate input
        if(choice > 0) {
            return(choice)
        } else {
            
            cat("Enter a value greater than 0\n")
        }
    }
}

#MENU FOR THE OPTION OF THE number of random webs to be geerated
#VALUE HAS TO BE GREATER THAN 0
#
display_numb_webs_menu <- function() {
    cat("\n")
    cat("===========================================\n")
    cat("FWTopo: Number of Random Webs to Generate\n")
    cat("===========================================\n")
    
    cat("===========================================\n")
    
}

# Get user input with validation
# FOR THE NUMBER OF WEBS TO BE GENERATED OPTION
get_num_fw_choice <- function() {
    display_numb_webs_menu()
    
    while (TRUE) {
        choice <- readline(prompt = "Enter the number of random food webs : \n ")
        
        # Validate input
        if(choice > 0) {
            return(choice)
        } else {
            
            cat("Enter a value greater than 0\n")
        }
    }
}

#Status and Contrastatus, Harary F. (1959)

StatusContrastatus<-function(gr, ...){
    #Nodes
    Nodes<-V(gr)
    Si<-numeric()
    Si.prime<-numeric()
    for(i in Nodes){
        Si.sum<-0
        Si.prime.sum<-0
        for(j in Nodes){
            Si.paths<-all_shortest_paths(gr, from = i, to = j, mode= "in")
            Si.prime.paths<-all_shortest_paths(gr, from = i, to = j, mode= "out")
            Si.sum<-ifelse(length(Si.paths$vpaths)>0,
                           Si.sum + (length(Si.paths$vpaths[[1]])-1),Si.sum)
            Si.prime.sum<-ifelse(length(Si.prime.paths$vpaths)>0,
                                 Si.prime.sum + (length(Si.prime.paths$vpaths[[1]])-1),Si.prime.sum)
        }
        Si[i]<-Si.sum
        Si.prime[i]<-Si.prime.sum
    }
    return(data.frame(Nodes=Nodes, Si=Si, Si.prime=Si.prime, Si.delta=Si-Si.prime))
}

#Positional importance based on indirect chain effects
#Function:
# Cite: Jordán, F., Liu, W., & van Veen, J.F. (2003). Quantifying 
# the importance of # species and their interactions in a host-parasitoid 
# community. Community # Ecology, 4(1), 79-88.

#Argument:
# x: Adjacent matrix
# n: number of steps

TopologicalImportance<-function(x, n, ...){
    DegreeTotal<-rowSums(x)+colSums(x)
    a.matrix<-as.matrix(x)
    for(i in 1:nrow(a.matrix)){
        ID.row<-which(a.matrix[i,]>0)
        ID.col<-which(a.matrix[,i]>0)
        a.matrix[i,ID.row]<-1/DegreeTotal[ID.row]
        a.matrix[i,ID.col]<-1/DegreeTotal[ID.col]
    }
    TI.list<-list()
    mult.matrix<-as.matrix(diag(nrow(a.matrix)))
    for(i in 1:n){
        mult.matrix<-a.matrix%*%mult.matrix
        TI.list[[i]]<-mult.matrix
    }
    a.sum<-matrix(0, ncol = ncol(a.matrix),nrow(a.matrix), 
                  dimnames=list(rownames(a.matrix), colnames(a.matrix)))
    for (i in 1:n) {
        a.sum<-a.sum+TI.list[[i]]
    }
    TIn<-data.frame(rowSums(a.sum)/n)
    colnames(TIn)<-paste("TI",n, sep="")
    return(TIn)
}

## CAPTURES AMOUNT OF RAM MEMORY
get_memory_usage <- function() {
    if (.Platform$OS.type == "windows") {
        mem <- system("wmic OS get FreePhysicalMemory,TotalVisibleMemorySize /Value", intern = TRUE)
        mem <- mem[grepl("FreePhysicalMemory|TotalVisibleMemorySize", mem)]
        free <- as.numeric(gsub("\\D", "", mem[1])) / 1024
        total <- as.numeric(gsub("\\D", "", mem[2])) / 1024
    } else {
        mem <- system('free -m | grep Mem:', intern = TRUE)
        mem <- strsplit(mem, "\\s+")[[1]]
        total <- as.numeric(mem[2])
        free <- as.numeric(mem[4]) + as.numeric(mem[6]) + as.numeric(mem[7])  # free + buffers + cache
    }
    return(list(total = total, free = free, used_pct = (total - free) / total * 100))
}

#FUNCTION FOR THE COMPUTING TROPHIC LEVEL FOR RANDOM WEBS OPTION
#
display_TL_menu <- function() {
    cat("\n")
    cat("=================================================\n")
    cat("FWTopo: Trophic Level for Random Webs Option\n")
    cat("=================================================\n")
    
    cat("=================================================\n")
    
}

# Get user input with validation FOR THE PLOTTING OPTION
get_TL_choice <- function() {
    display_TL_menu()
    
    while (TRUE) {
        choice <- readline(prompt = "Please enter 1 = YES or 2 = NO : ")
        
        # Validate input
        if(choice == "1" | choice == "2" ) {
            return(choice)
        } else {
            
            cat("Enter YES or NO \n")
        }
    }
}

#function to generate validation file when only the random web is produced and analyzed
#
gen_valid_for_rnd <- function(adj_matrix, filename, algo_rnd, 
                              numb_fws, directorio, tiempo) {
    cat("🔬 Generating FWTopo Validation Report...\n")
    
    # Basic validation
    report <- list(
        title = "....FOOD WEB TOPOLOGY ANALYSIS...",
        title_2 = "         V 1.2",
        input_file = filename,
        analysis_date <- Sys.Date(),
        dimensions = dim(adj_matrix),
        is_square = nrow(adj_matrix) == ncol(adj_matrix),
        total_interactions = sum(adj_matrix),
        top_species = sum(rowSums(adj_matrix) == 0),
        basal_predators = sum(colSums(adj_matrix) == 0),
        algorithm = algo_rnd,
        number_of_fws = numb_fws,
        Time_spent = tiempo
    )
    #browser()
    # Save validation report
    validation_file <- paste0(tools::file_path_sans_ext(filename), "_LOG.txt")
    validation_file <- paste0(directorio, "/", format(Sys.time(), "%Y.%m.%d_%H.%M.%S_"), validation_file)
    capture.output(print(report), file = validation_file)
    
    cat("✅ Validation report saved as:", validation_file, "\n")
    return(report)

}