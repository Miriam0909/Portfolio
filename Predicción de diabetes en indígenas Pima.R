
#Librerias

library(readr)
library(tidyr)
library(dplyr)
library(car)
library(reshape2)
library(ggplot2)
library(CalvinBayes)
library(loo)
library(R2jags)
library(coda)
library(brms)

#Carga de datos
Diabetes <- read_delim("Diabetes.csv", delim = ";", 
                       escape_double = FALSE, trim_ws = TRUE)
View(Diabetes)
str(Diabetes)

#se descarta una variables para efectos del estudio no interesa
base = Diabetes[, -7]
base <- base %>%
  mutate(resu = ifelse(resu == "tested_positive", 1, 0)) # 1 es positivo a signos de diabetes 
str(base)

##Analisis descriptivo
par(mfrow = c(2, 2))
hist(base$preg)
hist(base$plas)
hist(base$pres) 
hist(base$skin)
hist(base$insu) 
hist(base$mass)
hist(base$age)
hist(base$resu) 

# Matriz de Correlación
cor_matrix <- cor(base[,sapply(base, is.numeric)], use = "complete.obs")

# Mapa de Calor de la Matriz de Correlación

cor_melt <- melt(cor_matrix)
ggplot(data = cor_melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlación") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()

cor(base)


#Supuestos de Multicolinealidad

# Calcular el VIF
vif_model <- lm(resu ~ ., data = base)
vif_values <- vif(vif_model)

print(vif_values)


#Verificacion de linealidad 

# Ajustar el modelo logit sin interacciones para evaluar la linealidad
logit_model <- glm(resu ~ preg + plas + pres + skin + insu + mass + age, 
                   family = binomial, data = base)

# Crear los component plus residual plots
crPlots(logit_model)

#Interaccion entre predictors

# Ajustar un modelo de regresión logística con interacciones
logit_model_interaction <- glm(resu ~ (preg + plas + pres + skin + insu + mass + age)^2, 
                               family = binomial, data = base)

# Resumen del modelo
summary(logit_model_interaction)

#Valores extremos

predictors = c('preg', 'plas', 'pres', 'skin', 'insu', 'mass', 'age')

# Crear boxplots para cada variable predictora
par(mfrow = c(3, 3))  # Configurar el gráfico en una cuadrícula de 3x3
for (pred in predictors) {
  boxplot(base[[pred]], main = pred, ylab = pred)
}

# Calcular la distancia de Cook para el modelo logit
cooksd <- cooks.distance(logit_model)
plot(cooksd, type = "h", main = "Cook's Distance", ylab = "Cook's Distance")
abline(h = 4/(nrow(base)-length(logit_model$coefficients)), col = "red") 

summary(base$resu)
#1:No , 2:Si
boxplot(base$resu,base$preg)
boxplot(base$resu,base$pres)
boxplot(base$resu,base$plas)
boxplot(base$resu,base$age)

##Analisis de modelos
# Define el modelo en R2JAGS
diabetes_model <- function() {
  for (i in 1:N) {
    # Modelo logístico
    logit(p[i]) <- beta0 + 
      beta1 * x[i, 1] +  # preg
      beta2 * x[i, 2] +  # plas
      beta3 * x[i, 3] +  # pres
      beta4 * x[i, 4] +  # skin
      beta5 * x[i, 5] +  # insu
      beta6 * x[i, 6] +  # mass
      beta7 * x[i, 7] +  # age
      beta8 * x[i, 3] * x[i, 7] +  # pres:age
      beta9 * x[i, 1] * x[i, 7] +  # preg:age
      beta10 * x[i, 2] * x[i, 7]    # plas:age
    y[i] ~ dbern(p[i])
  }
  
  # Priors
  beta0 ~ dnorm(0, 10)  # Prior para el intercepto
  beta1 ~ dnorm(0, 10)  # Prior para preg
  beta2 ~ dnorm(0, 10)  # Prior para plas
  beta3 ~ dnorm(0, 10)  # Prior para pres
  beta4 ~ dnorm(0, 10)  # Prior para skin
  beta5 ~ dnorm(0, 10)  # Prior para insu
  beta6 ~ dnorm(0, 10)  # Prior para mass
  beta7 ~ dnorm(0, 10)  # Prior para age
  beta8 ~ dnorm(0, 10)  # Prior para pres:age
  beta9 ~ dnorm(0, 10)  # Prior para preg:age
  beta10 ~ dnorm(0, 10)  # Prior para plas:age
}


# Crear la lista de datos para JAGS 
JAGS_data <- list(
  N = nrow(base),                             # Número de observaciones
  y = base$resu,                              # Variable dependiente binaria: 0 y 1
  x = as.matrix(base[, c('preg', 'plas', 'pres', 'skin', 'insu', 'mass', 'age')])  # Variables predictoras
)

# Ajuste del modelo en JAGS (con los datos actualizados)
set.seed(12)
jags_fit <- jags(data = JAGS_data,
                 inits = NULL,
                 parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "beta10"),
                 model.file = diabetes_model,
                 n.chains = 3,
                 n.iter = 5000,
                 n.burnin = 1000,
                 n.thin = 1,
                 jags.seed = 123)

# Resumen del ajuste del modelo
print(jags_fit)

#Diagnosticos
mcmc_samples <- as.mcmc(jags_fit)

diag_mcmc(mcmc_samples, par = "beta0")

# Realiza la validación cruzada LOO
(loonorm = loo(as.array(mcmc_samples)))
(waicnorm = waic(as.array(mcmc_samples)))



#JAGS con cauchy 
#hay una forma más directa de usar la distribución de Cauchy en JAGS utilizando la distribución de Cauchy estándar dt.

# Definir el modelo en JAGS con previas de Cauchy
diabetes_model <- function() {
  for (i in 1:N) {
    # Modelo logístico
    logit(p[i]) <- beta0 + 
      beta1 * x[i, 1] +  # preg
      beta2 * x[i, 2] +  # plas
      beta3 * x[i, 3] +  # pres
      beta4 * x[i, 4] +  # skin
      beta5 * x[i, 5] +  # insu
      beta6 * x[i, 6] +  # mass
      beta7 * x[i, 7] +  # age
      beta8 * x[i, 3] * x[i, 7] +  # pres:age
      beta9 * x[i, 1] * x[i, 7] +  # preg:age
      beta10 * x[i, 2] * x[i, 7]    # plas:age
    y[i] ~ dbern(p[i])
  }
  
  # Priors de Cauchy
  beta0 ~ dt(0, 2.5, 1) # Cauchy(0, 2.5)
  beta1 ~ dt(0, 2.5, 1)
  beta2 ~ dt(0, 2.5, 1)
  beta3 ~ dt(0, 2.5, 1)
  beta4 ~ dt(0, 2.5, 1)
  beta5 ~ dt(0, 2.5, 1)
  beta6 ~ dt(0, 2.5, 1)
  beta7 ~ dt(0, 2.5, 1)
  beta8 ~ dt(0, 2.5, 1)
  beta9 ~ dt(0, 2.5, 1)
  beta10 ~ dt(0, 2.5, 1)
}

library(R2jags)
library(coda)

# Crear la lista de datos para JAGS 
JAGS_data <- list(
  N = nrow(base),                             # Número de observaciones
  y = base$resu,                              # Variable dependiente binaria: 0 y 1
  x = as.matrix(base[, c('preg', 'plas', 'pres', 'skin', 'insu', 'mass', 'age')])  # Variables predictoras
)

# Ajuste del modelo en JAGS (con los datos actualizados)
set.seed(12)
jags_fit <- jags(data = JAGS_data,
                 inits = NULL,
                 parameters.to.save = c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "beta10"),
                 model.file = diabetes_model,
                 n.chains = 3,
                 n.iter = 5000,
                 n.burnin = 1000,
                 n.thin = 1,
                 jags.seed = 121)

# Resumen del ajuste del modelo
print(jags_fit)


#Diagnosticos
mcmc_samples <- as.mcmc(jags_fit)

diag_mcmc(mcmc_samples, par = "beta0")

# Realiza la validación cruzada LOO
(loocauchy = loo(as.array(mcmc_samples)))
(waicauchy = waic(as.array(mcmc_samples)))



# Datos para brms
brms_data <- list(
  N = nrow(base),                                      # Número de observaciones
  y = as.integer(base$resu),                           # Variable dependiente binaria: 0 y 1
  preg = base$preg,                                    # Variable predictora: preg
  plas = base$plas,                                    # Variable predictora: plas
  pres = base$pres,                                    # Variable predictora: pres
  skin = base$skin,                                    # Variable predictora: skin
  insu = base$insu,                                    # Variable predictora: insu
  mass = base$mass,                                    # Variable predictora: mass
  age = base$age                                       # Variable predictora: age
)

# Definir el modelo en brms
diabetes_model <- brm(
  formula = y ~ preg + plas + pres + skin + insu + mass + age + pres:age + preg:age + plas:age,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(normal(0,10), class = Intercept),
    prior(normal(0,10), class = b, coef = "preg"),
    prior(normal(0,10), class = b, coef = "plas"),
    prior(normal(0,10), class = b, coef = "pres"),
    prior(normal(0,10), class = b, coef = "skin"),
    prior(normal(0,10), class = b, coef = "insu"),
    prior(normal(0,10), class = b, coef = "mass"),
    prior(normal(0,10), class = b, coef = "age"),
    prior(normal(0,10), class = b, coef = "pres:age"),
    prior(normal(0,10), class = b, coef = "preg:age"),
    prior(normal(0,10), class = b, coef = "plas:age")
  ),
  data = brms_data,
  chains = 3,
  iter = 5000,
  warmup = 1000,
  thin = 1,
  seed = 123
)

# Resumen del ajuste del modelo
summary(diabetes_model)

#Diagnósticos

mcmc_results <- as.mcmc.list(as.mcmc(diabetes_model))
mcmc_combo(mcmc_results)

# Gráfico de trazas (trace plot)
mcmc_trace(mcmc_results, regex_pars = "^b_", inc_warmup = TRUE)

# Gráfico de densidad (density plot)
mcmc_dens(mcmc_results, regex_pars = "^b_")

# Gráfico de autocorrelación (autocorrelation plot)
mcmc_acf(mcmc_results, regex_pars = "^b_")

# Gráfico de intervalos de HPD (Highest Posterior Density)
mcmc_intervals(mcmc_results, regex_pars = "^b_", prob = 0.95)

#Validacion cruzada
(bwaicnorm = waic(diabetes_model))
(bloonorm = loo(diabetes_model))

# Modelo en brms con previas de Cauchy
diabetes_model <- brm(
  formula = y ~ preg + plas + pres + skin + insu + mass + age + pres:age + preg:age + plas:age,
  family = bernoulli(link = "logit"),
  prior = c(
    prior(cauchy(0, 2.5), class = Intercept),
    prior(cauchy(0, 2.5), class = b, coef = "preg"),
    prior(cauchy(0, 2.5), class = b, coef = "plas"),
    prior(cauchy(0, 2.5), class = b, coef = "pres"),
    prior(cauchy(0, 2.5), class = b, coef = "skin"),
    prior(cauchy(0, 2.5), class = b, coef = "insu"),
    prior(cauchy(0, 2.5), class = b, coef = "mass"),
    prior(cauchy(0, 2.5), class = b, coef = "age"),
    prior(cauchy(0, 2.5), class = b, coef = "pres:age"),
    prior(cauchy(0, 2.5), class = b, coef = "preg:age"),
    prior(cauchy(0, 2.5), class = b, coef = "plas:age")
  ),
  data = brms_data,
  chains = 3,
  iter = 5000,
  warmup = 1000,
  thin = 1,
  seed = 123
)

# Resumen del ajuste del modelo
summary(diabetes_model)



# Extraer la muestra de la cadena para diagnósticos

mcmc_results <- as.mcmc.list(as.mcmc(diabetes_model))
mcmc_combo(mcmc_results)

# Gráfico de trazas (trace plot)
mcmc_trace(mcmc_results, regex_pars = "^b_", inc_warmup = TRUE)

# Gráfico de densidad (density plot)
mcmc_dens(mcmc_results, regex_pars = "^b_")

# Gráfico de autocorrelación (autocorrelation plot)
mcmc_acf(mcmc_results, regex_pars = "^b_")

# Gráfico de intervalos de HPD (Highest Posterior Density)
mcmc_intervals(mcmc_results, regex_pars = "^b_", prob = 0.95)

#Validacion cruzada
(bwaiccauchy = waic(diabetes_model))
(bloocauchy = loo(diabetes_model))

#Comparacion de modelos

#brms
loo_compare(bloonorm,bloocauchy)

#jags
loo_compare(loonorm,loocauchy)


#Se deseaba realizar una comparacion entre los dos modelos finales de cada herramienta (jags vs brms) pero no se logro,
#indicaba que no presentan el mismo tamño los dos modelos.