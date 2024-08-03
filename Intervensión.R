library(igraph)
library(bnlearn)
library(causaleffect)
library(gRain)

datos <-read.csv("file", header = TRUE)

datos_ <- datos
red_carac <-read.csv("file",header = TRUE)

red_carac<-as.matrix.noquote(red_carac)
datos_ <- as.data.frame(datos_[,unique(red_carac)])  
datos <- datos_
for(i in 1:length(datos_[1,])){
  datos[,i] <- as.factor(datos_[,i])
  print(typeof(datos[,i]))
  if(length(unique(datos_[,i]))==1){
    datos <-datos_[,-i]
  }
}
base <- datos

causa <-  "variable cause" 
efecto <- "variable effect"
intervencion <- unique(base[causa])[1,]


RB1 <- empty.graph(unique(red_carac))
RB1$arcs <- red_carac
RB <- pdag2dag(RB1,ordering = names(base))
tablas_RB <- bn.fit(RB, base, method = "mle")
modeloRB <- bnlearn::modelstring(RB)#Modelo de la red


#Imprime RB
red_arcos <- arcs(RB)
formula <- graph.data.frame(red_arcos)

#Imprime RB-mutilada
causa_red <- which(red_arcos[, 2] == causa)
if (length(causa_red) == 0) {
  red_arcos_causal <- red_arcos
} else{
  red_arcos_causal <- red_arcos[-causa_red, ]
}

#formula_2 <- graph.data.frame(red_arcos_causal)

#Grafica las redes
roriginal <- graph_from_data_frame(red_arcos)
rcausal <- graph_from_data_frame(red_arcos_causal)
pal2 <- rainbow(5, alpha = .2)

grafo <- graph.data.frame(red_arcos,directed = T)
grafo_sc <- simplify(grafo,remove.multiple = T, remove.loops = T)
g <- as.directed(grafo_sc)
conexiones<- as_data_frame(g)
model <- bnlearn::modelstring(as.bn(grafo_sc))
bl.model <- bnlearn::model2network(model)
plot(bl.model,main="RW")
graphviz.plot(bl.model, highlight = list(nodes = c(causa,efecto),
                                         col = c("limegreen"), 
                                         fill = c("lightblue")),layout = "dot",shape = "ellipse")

# Expresion causal a partir de la libreria
ruta_causal <-
  causal.effect(efecto,
                causa,
                G = formula,
                expr = F,
                steps = F,simp = F)
RC <- causal.effect(efecto, causa, G = formula,simp = F)

# Separa la expresión, busca el corchete y elimina primera parte de la expresion
RC_split <- as.character(unlist(strsplit(RC, split = "")))
corchete_cierra <- which(RC_split == "}")
corchete_abre <- which(RC_split == "{")
if (length(corchete_cierra) == 0){
  RC_limpia <- RC_split
  RC_suma <- RC_limpia
} else{
  RC_limpia <- RC_split[-1:-corchete_cierra]
  RC_suma <- RC_split[+(corchete_abre + 1):+(corchete_cierra - 1)]
  
}

# Crea lista de expresiones para calcular probabilidades
inicia_proba <- c(which(RC_limpia == "("))
termina_proba <- c(which(RC_limpia == ")"))
i = 1
prob_cond <- data.frame()

while (i <= length(termina_proba)) {
  extrae_de <- inicia_proba[i] + 1
  extrae_a <- termina_proba[i] - 1
  prob_cond[i, 1] <-
    paste(RC_limpia[extrae_de:extrae_a], collapse = "")
  i = i + 1
}

### Se usa para extraer las variables que formaran parte de la ruta causal
vector_sum <- c()
extrae_de <- 1
vector_comas <- c(which(RC_suma == ","))
i = 1
if (length(vector_comas) >= 1 ) {
  while (i <= length(vector_comas)){
    if (i == length(vector_comas)) {
      if(length(vector_comas) != 1){
        j = i-1
        extrae_a <- vector_comas[i] - 1
        extrae_de <- vector_comas[j] + 1
        vector_sum[i] <- paste(RC_suma[extrae_de:extrae_a], collapse = "")
      }
      extrae_dee <- vector_comas[i] + 1
      extrae_aa <- length(RC_suma)
      vector_sum[i+1] <- paste(RC_suma[extrae_dee:extrae_aa], collapse = "")
    }
    extrae_a <- vector_comas[i] - 1
    vector_sum[i] <- paste(RC_suma[extrae_de:extrae_a], collapse = "")
    extrae_de <- vector_comas[i] + 1
    i = i + 1
  }
} else{
  vector_sum <- paste(RC_suma, collapse = "")
}


#################################################
###########Solo para graficar rura causal########
#################################################
p = 1
nodos_causales <- c()
for (p in 1:length(prob_cond[, 1])) {
  p_colapsar <- unlist(strsplit(prob_cond[p, 1], split = ""))
  
  #Busca el | que en la expresion esta el condicionante y la condicion
  busca_condicion <- which(p_colapsar == "|")
  if (length(busca_condicion) == 0) {
    busca_condicion <- 0
    condicion <- c(paste(p_colapsar, collapse = ""))
  } else{
    condicion_hasta <- busca_condicion - 1
    condicion <-
      c(paste(p_colapsar[1:condicion_hasta], collapse = ""))
  }
  nodos_causales <-  c(nodos_causales, condicion)
  
  if (busca_condicion != 0) {
    #Separa y busca variables condicionantes
    cond_coma <-
      c(paste(p_colapsar[(busca_condicion + 1):length(p_colapsar)], collapse = ""))
    busca_condicionante <-
      which(unlist(strsplit(cond_coma, split = "")) == ",")
    condicionantes <- c()
    if (length(busca_condicionante) == 0) {
      busca_condicionante <- 0
      condicionantes <- cond_coma
      
    } else{
      vector_condicionantes <- unlist(strsplit(cond_coma, split = ""))
      i = 1
      for (j in 1:length(busca_condicionante)) {
        j2 <- busca_condicionante[j] - 1
        condicionantes[length(condicionantes) + 1] <-
          paste(vector_condicionantes[i:j2], collapse = "")
        i = j2 + 2
        if (j == length(busca_condicionante)) {
          j2 <- busca_condicionante[j] + 1
          condicionantes[length(condicionantes) + 1] <-
            paste(vector_condicionantes[j2:length(vector_condicionantes)], collapse = "")
          
        }
      }
      
    }
    nodos_causales <- c(nodos_causales, condicionantes)
  }
  
}

nodos_a_colorear <- unique(nodos_causales)

V(rcausal)$causa <- c(rep("no_causa", length(V(rcausal)$name)))
causa_pos <- c()
for (i in 1:length(nodos_a_colorear)) {
  causa_pos <-
    c(causa_pos, which(V(rcausal)$name == nodos_a_colorear[i]))
  V(rcausal)$causa[causa_pos] = "causa"
}
pal <- rainbow(2, alpha = .2)
plot(
 rcausal,
 vertex.color = c(pal)[1 + (V(rcausal)$causa == "no_causa")],
 vertex.size = 20,
 edge.arrow.size = .3,
 main = "Red Causal"
)

rcausal <- simplify(rcausal,remove.multiple = TRUE, remove.loops = FALSE)
g <- as.directed(rcausal,"arbitrary")
model <-bnlearn::modelstring(as.bn(rcausal))
bl.model <- model2network(model)
graphviz.plot(bl.model, highlight = list(nodes = nodos_a_colorear,col = "green", fill = "lightblue"),
              layout = "dot",shape = "ellipse",groups = nodos_a_colorear)


############################################################################################
###############  Calcula las tablas de probabilidad para todas las expresiones #############
############################################################################################
if(length(vector_sum) >= 1){
  tablas_de_prob <- list()
  tablas_probas_1 <- list()
  p = 1
  for (p in 1:length(prob_cond[, 1])) {
    p_colapsar <- unlist(strsplit(prob_cond[p, 1], split = ""))

    #Busca el dado que en la expresion
    busca_condicion <- which(p_colapsar == "|")
    if (length(busca_condicion) == 0) {
      busca_condicion <- 0
      condicion <- c(paste(p_colapsar, collapse = ""))
      condicionantes <- c()
    } else{
      # Busca variable condicion
      condicion_hasta <- busca_condicion - 1
      #Separa y busca variables condicionantes
      cond_coma <-
        c(paste(p_colapsar[(busca_condicion + 1):length(p_colapsar)], collapse = ""))
      condicion <-
        c(paste(p_colapsar[1:condicion_hasta], collapse = ""))
      busca_condicionante <-
        which(unlist(strsplit(cond_coma, split = "")) == ",")

      ########### busca condicionantes de la expresion############
      condicionantes <- c()
      i = 1
      if (length(busca_condicionante) == 0) {
        condicionantes <- cond_coma
      } else{
        vector_condicionantes <- unlist(strsplit(cond_coma, split = ""))
        for (j in 1:length(busca_condicionante)) {
          j2 <- busca_condicionante[j] - 1
          condicionantes[length(condicionantes) + 1] <-
            paste(vector_condicionantes[i:j2], collapse = "")
          i = j2 + 2
          j2 <- busca_condicionante[j + 1] - 1
          if (j == length(busca_condicionante)) {
            condicionantes[length(condicionantes) + 1] <-
              paste(vector_condicionantes[i:length(vector_condicionantes)], collapse = "")
          }
        }
      }
    }

    #Combinatoria, busca los valores de las variables para calcular las probabilidades
    pos_niveles <- c(which(names(base) == condicion))
    i = 1
    if (length(condicionantes) != 0) {
      while (i <= length(condicionantes)) {
        pos_niveles[length(pos_niveles) + 1] <-
          which(names(base) == condicionantes[i])
        i = i + 1
      }
    } else{
      pos_niveles <- pos_niveles
    }

    # Busca los niveles de las variables a calcular
    niveles <- list()
    i = 1
    while (i <= length(pos_niveles)) {
      niveles[[i]] <- levels(base[, pos_niveles[i]])
      i = i + 1
    }
    # Calcula las tablas de probabilidad para el conjunto de expresiones de la ruta causal
    tablas_probas <- expand.grid(niveles)
    tablas_probas <- as.matrix(tablas_probas)
    nombres <- names(base[pos_niveles])
    colnames(tablas_probas) <- nombres

    tablas_probas_1 <- expand.grid(niveles)
    tablas_probas_1 <- as.matrix(tablas_probas_1)
    nombres <- names(base[pos_niveles])
    colnames(tablas_probas_1) <- nombres

    if (length(pos_niveles) != 1) {
      probas <- data.frame()
      for (j in 1:length(tablas_probas[, 1])) {
        numerador <- rep(0, length(base[, 1]))
        for (i in 1:length(names(tablas_probas[1, ]))) {
          numerador[base[names(tablas_probas[1, i])] == tablas_probas[j, i]] = numerador[base[names(tablas_probas[1, i])] ==
                                                                                           tablas_probas[j, i]] + 1
        }
        numerador <-
          length(which(numerador == length(names(
            tablas_probas[1, ]
          ))))
        denominador <- rep(0, length(base[, 1]))
        for (i in 2:length(names(tablas_probas[1, ]))) {
          denominador[base[names(tablas_probas[1, i])] == tablas_probas[j, i]] = denominador[base[names(tablas_probas[1, i])] ==
                                                                                               tablas_probas[j, i]] + 1
        }
        denominador <-
          length(which(denominador == length(names(
            tablas_probas[1, -1]
          ))))
        probas[j, 1] <- numerador / denominador
      }
      tablas_de_prob[[p]] <- cbind(tablas_probas, probas)
    } else{
      tab_prob <- t(table(base[, prob_cond[p, 1]]) / length(base[, 1]))
      tab_prob <- t(tab_prob)
      tablas_probas_1 <- data.frame(tablas_probas_1, tab_prob[, 1])
      tablas_probas_1 <- type.convert(tablas_probas_1, "list")
      colnames(tablas_probas_1) <- c(nombres, "V1")
      tablas_de_prob[[p]] <- tablas_probas_1
    }
  }

  tablas <- list()
  for (tb in 1:length(tablas_de_prob)) {
    pos.NaN <- which(tablas_de_prob[[tb]][, "V1"] != "NaN")
    tablas[[tb]] <- tablas_de_prob[[tb]][pos.NaN, ]
  }

  ## Crea lista de tablas solo con variables y valores intervenidos
  tabla_intervenida <- tablas
  tcon_intervencion <- c()
  for (i in 1:length(tablas)) {
    if ((length(tablas[[i]][1,]) == 2) || sum(names(tablas[[i]]) == causa)==0){
      tabla_intervenida[[i]] <- tablas[[i]]
    } else{
      tcon_intervencion <- which(tablas[[i]][, causa] == intervencion)
      tabla_intervenida[[i]] <- tabla_intervenida[[i]][tcon_intervencion, ]

    }
  }

  vector_sum_1 <- vector_sum

  niveles.sum <- list()
  for (s in 1:length(vector_sum_1)) {
    niveles.sum[[s]] <- levels(base[, vector_sum_1[s]])
  }
  niveles.sum <- expand.grid(niveles.sum)
  colnames(niveles.sum) <- vector_sum_1

  categorias_efecto <- levels(base[, efecto])
  suma.probas.condicionales <- c()
  guarda.productos <- c()
  prob_categoria <- data.frame()
  coincidentes.p.suma = c()
  g = 1
  for (c in 1:length(categorias_efecto)) {
    for (s in 1:length(niveles.sum[, 1])) {
      for (t in 1:length(tabla_intervenida)){
        if (length(niveles.sum[1, ]) == 1) {
          coincidentes.p.suma = tabla_intervenida[[t]][tabla_intervenida[[t]][, vector_sum_1] == niveles.sum[s, ], ]
        } else{
          coincidentes.p.suma <- merge(niveles.sum[s, ], tabla_intervenida[[t]])
        }

        factor.normalizacion = sum(coincidentes.p.suma$V1)


        if (length(which(names(coincidentes.p.suma) == efecto)) != 0) {
          concidente.efecto <-which(coincidentes.p.suma[, efecto] == categorias_efecto[c])
          coincidentes.p.suma <- coincidentes.p.suma[concidente.efecto, ]
          suma.para.normalizar <-sum(coincidentes.p.suma[, "V1"]) / factor.normalizacion
          suma.probas.condicionales[t] <- suma.para.normalizar
        } else{
          if (length(names(tabla_intervenida[[t]])) > 2) {
            factor.normalizacion_2 = sum(merge(niveles.sum[s, ], tabla_intervenida[[t]][,-1])$V1)
            suma.probas.condicionales[t] <- sum(coincidentes.p.suma[, "V1"]) / factor.normalizacion_2
          } else{
            suma.probas.condicionales[t] <- sum(coincidentes.p.suma[, "V1"])
          }

          if (length(suma.probas.condicionales) == length(tablas_de_prob)) {
            if(is.na(prod(suma.probas.condicionales))){
              guarda.productos[s] <- 0
            }else{
              guarda.productos[s] <- prod(suma.probas.condicionales)
            }
            suma.probas.condicionales <- c()
            g = g + 1
          }

        }# es el de las 2das tablas
      }# es el de la T
    }
    prob_categoria[c, 1] <- categorias_efecto[c]
    prob_categoria[c, 2] <- sum(guarda.productos)
  }
  colnames(prob_categoria) <- c(efecto, "P.Causal")

  if(length(tabla_intervenida)==1){
    prob_categoria <- data.frame()
    for(i in 1:length(categorias_efecto)){
      pos <- which(tabla_intervenida[[1]][,efecto]==categorias_efecto[i])
      sum <- sum(tabla_intervenida[[1]][pos,"V1"]) / sum(tabla_intervenida[[1]][,"V1"])
      prob_categoria[i, 1] <- categorias_efecto[i]
      prob_categoria[i, 2] <- sum
    }
  }
  colnames(prob_categoria) <- c(efecto, "P.Causal")

  funcion = compile(as.grain(tablas_RB))
  causa_consulta = setEvidence(funcion, nodes = causa, states = intervencion)
  propagacion <- querygrain(causa_consulta,nodes = c(efecto, causa),type = "conditional")
  propagacion.observacion <- as.data.frame(propagacion[,intervencion])
  colnames(propagacion.observacion) <- "P.Observacion"
  prob_categoria

}

prob_causal <- prob_categoria
prob_causal[,2] <- prob_categoria[,2]/sum(prob_categoria[,2])
prob_causal

propagacion.observacion["P.Observacion"]
RC

which(prob_categoria[, 2] != propagacion.observacion["P.Observacion"])
intervencion

probabilidades <- list("Efecto_" = prob_causal,propagacion.observacion,"Causa"= causa,"Valor"=intervencion)

write.csv(probabilidades,"file)

