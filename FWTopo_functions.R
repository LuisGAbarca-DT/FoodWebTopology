#  
#  
#               Food Web Topology Analysis
#                   FUNCTIONS
#                   V 1.4
#                   26 NOVEMBER 2025
#   
#   Luis Gerardo Abarca     gabarca@uv.mx   luisgaa@gmail.com
#   Israel Huesca Domínguez ihuesca@uv.mx
#   
#               IIB Universidad Veracruzana
#   
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

rol_nodos <- function(mat) {
    #simetrizar 
    matriz <- igraph::as_adjacency_matrix(mat, sparse=FALSE)
    matriz_simet <- sna::symmetrize(matriz, rule = "weak")
    
    roles <- netcarto(matriz_simet)
    
    return(roles)
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
            
            Web.Cascade <- matrix(0, S, S, dimnames = list(paste("S",
                                                                 1:S, sep = ""), paste("S", 1:S, sep = "")))
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
        rownames(Web.Cascade) <- paste("SPS", names_1, sep = "")
        colnames(Web.Cascade) <- paste("SPS", names_1, sep = "")
        
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
            
            Web.Niche <- matrix(0, S, S, dimnames = list(paste("S",
                                                               1:S, sep = ""), paste("S", 1:S, sep = "")))
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
        rownames(Web.Niche) <- paste("SPS", names_1, sep = "")
        colnames(Web.Niche) <- paste("SPS", names_1, sep = "")
        WebNicheModel <- graph_from_adjacency_matrix(Web.Niche, weighted = FALSE, mode = c("directed"))
        
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
    print(". . . w o r k i n g. . . ")
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
    rownames(g_tempo) <- paste("SPS", names_1, sep = "")
    colnames(g_tempo) <- paste("SPS", names_1, sep = "")
    NODE <- rownames(g_tempo)
    
    commty <- Community(nodes = data.frame(node=NODE),
                        trophic.links=PredationMatrixToLinks(g_tempo),
                        properties=list(title="Community"))
    #CLEAR g_tempo...RAM amount IS IMPORTANT AT THIS STAGE
    rm(g_tempo) 
    
    if (tl_y_or_no == "YES") {
        chain.stats <- TrophicChainsStats(commty)
        troph_level <- ShortWeightedTrophicLevel(commty)
        tl_mean <- mean(troph_level)
        maxTL <- max(troph_level)
        chains_nom <- length(chain.stats$chain.lengths)
        
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
                nivel_trof_sps = troph_level,
                vul_std = sd(vul),
                gen_std = sd(gen)))
    
}

#   FUNCION draw_1 RED
#   FUNCTION TO DRAW A FOOD WEB
#   graph_data = an igraph object of the food web
#   
draw_1 <- function(graph_data, comun) {
    
    plot(graph_data, 
         vertex.shape="sphere", 
         layout = layout_nicely(graph_data),   #layout_as_tree(gr,root = c(40,41,42,43), circular = TRUE, flip.y = TRUE),  #layout_with_kk(gr, dim=3),      #layout_nicely, 
         vertex.size=5, 
         edge.arrow.size =0.2, 
         edge.curved = 0.4, 
         #edge.width = E(gr)$weight * 10,
         vertex.label.cex = 1)
    #BY TROPHIC LEVEL USING SHORT WEIGHRTD TROPHIC LEVEL OPTION FROM cheddar
    #
    PlotWebByLevel(comun, level = "ShortWeightedTrophicLevel")
    
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
    # miembros_azar$sps <- row.names(miembros_azar)
    
    modularidad_azar <- modularity(matriz, membership(resulta_l),
                                   directed = TRUE)

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
        rownames(g_tempo) <- paste("SPS", names_1, sep = "")
        colnames(g_tempo) <- paste("SPS", names_1, sep = "")
        g_rand_erdos<-graph_from_adjacency_matrix(g_tempo, weighted = FALSE, mode = c("directed"))
        
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
    # i <- 1
    #browser()
    
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
            
            readline(prompt = "enter para continuar...")
            
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
gen_rnd_fw <- function(gr, names_1, num_rand_webs, random_model_type) {

#browser()
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
            rand_matrix <- lapply(1:num_rand_webs, function(i) Web.NicheModel(sps, conec, tolerancia, names_1)) 
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
                
                rand_matrix <- lapply(1:num_rand_webs, function(i) Web.CascadeModel(sps, conec, tolerancia, names_1))
                
            }
        else {
            
            if (random_model_type == "niche_allesina") {
                
                #niche model modified by Allesina et al

                sps <- vcount(gr)
                num_links <- ecount(gr)
                conec <- num_links / (sps ^2)
                
                #rand_matrix <- niche_model_2(sps, conec, num_rand_webs)
                #browser()
                rand_matrix <- lapply(1:num_rand_webs, function(i) 
                    niche_model_2(sps, conec, names_1))
            }
               
        else {
            if (random_model_type == "erdos-renyi") {
                
            #browser()
            #erdos-renyi
            rand_matrix <- lapply(1:num_rand_webs, function(i) er(gr, names_1))
            }
            
        } 
        }
        }

    return(rand_matrix)
    
}

#COMPUTE THE MODIFIED NICHE MODEL (ALLESINA, ALONSO, PASCUAL 2008)
niche_model_2 <- function(sps, conec, names_1) {
    
        tempo <- create_niche_model(sps, conec)

        cat("FW-NICHE-MODEL-ALLESINA\n")

        rownames(tempo) <- paste("SPS", names_1, sep = "")
        colnames(tempo) <- paste("SPS", names_1, sep = "")
        
        g_rand_niche_2<-graph_from_adjacency_matrix(tempo, weighted = FALSE, mode = c("directed"))

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
        title_2 = "         V 1.4",
        input_file = filename,
        analysis_date = Sys.Date(),
        dimensions = dim(adj_matrix),
        is_square = nrow(adj_matrix) == ncol(adj_matrix),
        total_interactions = sum(adj_matrix),
        connectance = round(sum(adj_matrix) / nrow(adj_matrix)^2, 3),
        basal_species = sum(rowSums(adj_matrix) == 0),
        top_predators = sum(colSums(adj_matrix) == 0),
        is_connected = igraph::is_connected(igraph::graph_from_adjacency_matrix(adj_matrix)),
        algorithm = algo_rnd,
        number_of_fws = numb_fws,
        Leiden_resolution = resol,
        Time_spent = tiempo
    )
    #browser()
    # Save validation report
    validation_file <- paste0(tools::file_path_sans_ext(filename), "_LOG.txt")
    validation_file <- paste0(directorio, "/", format(Sys.time(), "%Y%m%d_%H%M%S_"), validation_file)
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
    rownames(g_tempo) <- paste("SPS", names_1, sep = "")
    colnames(g_tempo) <- paste("SPS", names_1, sep = "")
    NODE <- rownames(g_tempo)
    
    commty <- Community(nodes = data.frame(node=NODE),
                        trophic.links=PredationMatrixToLinks(g_tempo),
                        properties=list(title="Community"))
    
    if (tl_y_or_no == "YES") {
        chain.stats <- TrophicChainsStats(commty)
        chains_nom <- length(chain.stats$chain.lengths)
            cat("Number of chains: \n")
            print(chains_nom)
        
        troph_level <- ShortWeightedTrophicLevel(commty)
        tl_mean <- mean(troph_level)
        maxTL <- max(troph_level)
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
    
    res_est_rnd_1 <- data.frame(N = num_sps, 
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
                vul_std = sd(vul),
                gen_std = sd(gen))

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
