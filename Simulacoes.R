
if(!require(meta)){install.packages("meta")}
if(!require(metafor)){install.packages("metafor")}
if(!require(psych)){install.packages("psych")}

#####################
##### Simulação #####
#####################

K = c(3, 5, 10, 15, 20, 30, 40) #Número de estudos presentes na metanálise
p = 0.01 #Prevalência escolhida
sigma = 0.05 #Variação entre os estudos escolhida
nsim = 1000 #Número de simulações
set.seed(220154)

comparações = list()

for(k in 1:length(K)){
  
  aux = matrix(NA, nrow = nsim, ncol = 40) 
  
  N = as.integer(sample(20:30, K[k], replace = T)) #Alterar aqui o tamanho amostral dos estudos individuais
  
  for(i in 1:nsim){
    
    prob = rnorm(K[k], mean = p, sd = sigma*p)
    
    event = as.integer(rbinom(n = length(N), size = N, prob = prob))           
    
    tryCatch({
      logitoa = metaprop(event = event, n = N, method.ci = "NAsm", sm = "PLOGIT", comb.fixed = F, comb.random = T, method = "Inverse", backtransf = T, method.tau = "DL")
      
      aux[i, 1] = logitoa$TE.random
      aux[i, 2] = logitoa$lower.random
      aux[i, 3] = logitoa$upper.random
      
      if((p >= transf.ilogit(logitoa$lower.random)) && (p <= transf.ilogit(logitoa$upper.random))){aux[i, 4] = 1} else{aux[i, 4] = 0}
      
      aux[i, 5] = logitoa$Q
      aux[i, 6] = logitoa$I2
      aux[i, 7] = logitoa$tau
      
    }, error=function(e){cat()})
    
    tryCatch({
      arcosenoa = metaprop(event = event, n = N, method.ci = "NAsm", sm = "PAS", comb.fixed = F, comb.random = T, method = "Inverse", backtransf = T, method.tau = "DL")
      
      aux[i, 8] = arcosenoa$TE.random
      aux[i, 9] = arcosenoa$lower.random
      aux[i, 10] = arcosenoa$upper.random
      
      if((p >= transf.iarcsin(arcosenoa$lower.random)) && (p <= transf.iarcsin(arcosenoa$upper.random))){aux[i, 11] = 1} else{aux[i, 11] = 0}
      
      aux[i, 12] = arcosenoa$Q
      aux[i, 13] = arcosenoa$I2
      aux[i, 14] = arcosenoa$tau
      
    }, error=function(e){cat()})
    
    tryCatch({
      arcosenoduploa = metaprop(event = event, n = N, method.ci = "NAsm", sm = "PFT", comb.fixed = F, comb.random = T, method = "Inverse", backtransf = T, method.tau = "DL")
      
      n1 = harmonic.mean(N)
      n2 = mean(N)
      n3 = geometric.mean(N)
      n4 = sum(arcosenoduploa$w.random)
      
      aux[i, 15] = transf.ipft(arcosenoduploa$TE.random, ni = n1) 
      aux[i, 16] = transf.ipft(arcosenoduploa$lower.random, ni = n1)
      aux[i, 17] = transf.ipft(arcosenoduploa$upper.random, ni = n1) 
      if((p >= transf.ipft(arcosenoduploa$lower.random, ni = n1)) && (p <= transf.ipft(arcosenoduploa$upper.random, ni = n1))){aux[i, 18] = 1} else{aux[i, 18] = 0}
      
      aux[i, 19] = transf.ipft(arcosenoduploa$TE.random, ni = n2) 
      aux[i, 20] = transf.ipft(arcosenoduploa$lower.random, ni = n2)
      aux[i, 21] = transf.ipft(arcosenoduploa$upper.random, ni = n2) 
      if((p >= transf.ipft(arcosenoduploa$lower.random, ni = n2)) && (p <= transf.ipft(arcosenoduploa$upper.random, ni = n2))){aux[i, 22] = 1} else{aux[i, 22] = 0}
      
      
      aux[i, 23] = transf.ipft(arcosenoduploa$TE.random, ni = n3) 
      aux[i, 24] = transf.ipft(arcosenoduploa$lower.random, ni = n3)
      aux[i, 25] = transf.ipft(arcosenoduploa$upper.random, ni = n3) 
      if((p >= transf.ipft(arcosenoduploa$lower.random, ni = n3)) && (p <= transf.ipft(arcosenoduploa$upper.random, ni = n3))){aux[i, 26] = 1} else{aux[i, 26] = 0}
      
      
      aux[i, 27] = transf.ipft(arcosenoduploa$TE.random, ni = n4) 
      aux[i, 28] = transf.ipft(arcosenoduploa$lower.random, ni = n4)
      aux[i, 29] = transf.ipft(arcosenoduploa$upper.random, ni = n4) 
      if((p >= transf.ipft(arcosenoduploa$lower.random, ni = n4)) && (p <= transf.ipft(arcosenoduploa$upper.random, ni = n4))){aux[i, 30] = 1} else{aux[i, 30] = 0}
      
      aux[i, 31] = arcosenoduploa$Q
      aux[i, 32] = arcosenoduploa$I2
      aux[i, 33] = arcosenoduploa$tau
      
    }, error=function(e){cat()})
    
    
    tryCatch({
      glmma = metaprop(event = event, n = N, method.ci = "CP", comb.fixed = F, comb.random = T, method.tau = "ML")
      
      aux[i, 34] = glmma$TE.random
      aux[i, 35] = glmma$lower.random
      aux[i, 36] = glmma$upper.random
      
      if((p >= transf.ilogit(glmma$lower.random)) && (p <= transf.ilogit(glmma$upper.random))){aux[i, 37] = 1} else{aux[i, 37] = 0}
      
      aux[i, 38] = glmma$Q
      aux[i, 39] = glmma$I2
      aux[i, 40] = glmma$tau
      
    }, error=function(e){cat()})
    
  }
  
  comparações[[k]] = aux
  
}

for(k in 1:length(K)){
  
  comparações[[k]] = data.frame(comparações[[k]])
  
  names(comparações[[k]]) = c("logito_aleat", "Li",  "Ls", "coverage", "Q", "I2", "tau", "arcoseno_aleat", "Li", "Ls", "coverage", "Q", "I2", "tau", "arcosenod_aleat harm", "Li", "Ls", "coverage", "arcosenod_aleat arit", "Li", "Ls", "coverage", "arcosenod_aleat geom", "Li", "Ls", "coverage", "arcosenod_aleat poolv", "Li", "Ls", "Q", "I2", "tau", "glmm_aleat", "Li", "Ls", "coverage", "Q", "I2", "tau")
  
  comparações[[k]][, 1] =  transf.ilogit(comparações[[k]][, 1])
  comparações[[k]][, 2] =  transf.ilogit(comparações[[k]][, 2])
  comparações[[k]][, 3] =  transf.ilogit(comparações[[k]][, 3])
  comparações[[k]][, 34] =  transf.ilogit(comparações[[k]][, 34])
  comparações[[k]][, 35] =  transf.ilogit(comparações[[k]][, 35])
  comparações[[k]][, 36] =  transf.ilogit(comparações[[k]][, 36])
  
  comparações[[k]][, 8] =  transf.iarcsin(comparações[[k]][, 8])
  comparações[[k]][, 9] =  transf.iarcsin(comparações[[k]][, 9])
  comparações[[k]][, 10] =  transf.iarcsin(comparações[[k]][, 10])
  
}

resultados1 = list()

for(k in 1:length(K)){
  
  aux = data.frame(rbind(c(apply(comparações[[k]][, c(1:3)], 2, mean, na.rm = T), sum(comparações[[k]][, 4], na.rm = T), sum(is.na(comparações[[k]][, 1]) == "TRUE"), sum((p-comparações[[k]][, 1])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(5:7)], 2, mean, na.rm = T)),
                         
                         c(apply(comparações[[k]][, c(8:10)], 2, mean, na.rm = T), sum(comparações[[k]][, 11], na.rm = T), sum(is.na(comparações[[k]][, 8]) == "TRUE"), sum((p-comparações[[k]][, 8])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(12:14)], 2, mean, na.rm = T)),
                         
                         c(apply(comparações[[k]][, c(15:17)], 2, mean, na.rm = T), sum(comparações[[k]][, 18], na.rm = T), sum(is.na(comparações[[k]][, 15]) == "TRUE"), sum((p-comparações[[k]][, 15])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(31:33)], 2, mean, na.rm = T)),
                         
                         c(apply(comparações[[k]][, c(19:21)], 2, mean, na.rm = T), sum(comparações[[k]][, 22], na.rm = T), sum(is.na(comparações[[k]][, 19]) == "TRUE"), sum((p-comparações[[k]][, 19])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(31:33)], 2, mean, na.rm = T)),
                         
                         c(apply(comparações[[k]][, c(23:25)], 2, mean, na.rm = T), sum(comparações[[k]][, 26], na.rm = T), sum(is.na(comparações[[k]][, 23]) == "TRUE"), sum((p-comparações[[k]][, 24])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(31:33)], 2, mean, na.rm = T)),
                         
                         c(apply(comparações[[k]][, c(27:29)], 2, mean, na.rm = T), sum(comparações[[k]][, 30], na.rm = T), sum(is.na(comparações[[k]][, 27]) == "TRUE"), sum((p-comparações[[k]][, 27])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(31:33)], 2, mean, na.rm = T)),
                         
                         c(apply(comparações[[k]][, c(34:36)], 2, mean, na.rm = T), sum(comparações[[k]][, 37], na.rm = T), sum(is.na(comparações[[k]][, 34]) == "TRUE"), sum((p-comparações[[k]][, 34])^2, na.rm = T)/nsim, apply(comparações[[k]][, c(38:40)], 2, mean, na.rm = T))))
  
  names(aux) = c("me", "li", "ls", "coverage", "countna", "erro_quadr", "Q", "I2", "tau")
  aux$coverage = aux$coverage/nsim
  aux$trans = c("Logito", "Arco seno", "Arco seno duplo MH", "Arco seno duplo MA", "Arco seno duplo MG", "Arco seno duplo VC", "MLGM")
  aux$p = rep(p, 7)
  aux$sigma = rep(sigma, 7)
  aux$K = rep(K[k], 7)
  aux$amplitude = aux$ls-aux$li
  aux$erro = p-aux$me
  aux$amostra = rep("pequena", 7)
  
  resultados1[[k]] = aux
  
}

princ(resultados1)

###################
##### Figuras #####
###################

if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(gridExtra)){install.packages("gridExtra")}

dados1 = subset(dados, ) #Dados que serão utilizados para o primeiro gráfico da figura
dados2 = subset(dados, ) #Dados que serão utilizados para o segundo gráfico da figura
dados3 = subset(dados, ) #Dados que serão utilizados para o terceiro gráfico da figura
dados4 = subset(dados, ) #Dados que serão utilizados para o quarto gráfico da figura

#Calculando o limite mínimo e o máximo da variável que será plotada no eixo Y
minimo = min(dados1$variavel, dados2$variavel, dados3$variavel, dados4$variavel)
maximo = max(dados1$variavel, dados2$variavel, dados3$variavel, dados4$variavel)

p1 = ggplot(dados1, aes(x = K, y = variavel, group = forma_de_estimacao))  + geom_line(aes(color = forma_de_estimacao), alpha = 0.5, size = 1) + geom_point(aes(shape = forma_de_estimacao, color = forma_de_estimacao), alpha = 0.75, size = 4) + scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 30, 40)) + theme_bw() + labs(title = "Título do gráfico") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + ylim(c(minimo, maximo)) + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) 

p2 = ggplot(dados2, aes(x = K, y = variavel, group = forma_de_estimacao))  + geom_line(aes(color = forma_de_estimacao), alpha = 0.5, size = 1) + geom_point(aes(shape = forma_de_estimacao, color = forma_de_estimacao), alpha = 0.75, size = 4) + scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 30, 40)) + theme_bw() + labs(title = "Título do gráfico") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + ylim(c(minimo, maximo)) + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) + theme(legend.position='none')

p3 = ggplot(dados3, aes(x = K, y = variavel, group = forma_de_estimacao))  + geom_line(aes(color = forma_de_estimacao), alpha = 0.5, size = 1) + geom_point(aes(shape = forma_de_estimacao, color = forma_de_estimacao), alpha = 0.75, size = 4) + scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 30, 40)) + theme_bw() + labs(title = "Título do gráfico") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + ylim(c(minimo, maximo)) + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) + theme(legend.position='none')

p4 = ggplot(dados4, aes(x = K, y = variavel, group = forma_de_estimacao))  + geom_line(aes(color = forma_de_estimacao), alpha = 0.5, size = 1) + geom_point(aes(shape = forma_de_estimacao, color = forma_de_estimacao), alpha = 0.75, size = 4) + scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 30, 40)) + theme_bw() + labs(title = "Título do gráfico") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = 0, linetype = "dashed", color = "red") + ylim(c(minimo, maximo)) + theme(axis.title.y = element_blank()) + theme(axis.title.x = element_blank()) + theme(legend.position='none')

legend = get_legend(p1)
p1 = p1 + theme(legend.position='none')

grid.arrange(p1, p2, p3, p4, legend, nrow = 2, ncol = 3, left = text_grob("Variável", rot = 90, vjust = 1), bottom = text_grob("Número de estudos", hjust = 1.3), layout_matrix = rbind(c(1,1,2,2,5), c(3,3,4,4,5)))
