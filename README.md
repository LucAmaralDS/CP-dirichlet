# Predição Conforme e Modelos Dirichlet

Este repositório contém códigos em R para análise de dados composicionais utilizando **Regressão Dirichlet** e **Predição Conforme (Conformal Prediction)**, com foco em cobertura preditiva e diagnóstico de modelos.

## Estrutura do Repositório
├── R6data.csv # Dados experimentais com 3 componentes (LMF, SMF, RMF)

├── aplicacao-biomass.R # Script principal para análise de biomassa (3 componentes)

├── dados-boot.xlsx # Dados de sono com 4 componentes (N1, N2, N3, REM)

├── sleep-data.R # Script de análise de sono e predição conforme (4 componentes)

├── Figuras/ # Figuras HDR 3D e ternário para visualização

├── README.md # Este arquivo


## Dependências

```r
install.packages(c(
  "DirichletReg", "progress", "tidyverse", 
  "ggtern", "gridExtra", "grid", "viridis", "readxl", "dplyr"
))

Como Executar:
Abra R ou RStudio no diretório do projeto.
Instale os pacotes listados acima.
Execute os scripts:

source("aplicacao-biomass.R")  # Análise de biomassa
source("sleep-data.R")         # Análise de sono


Os resultados e diagnósticos serão exibidos no console.

## **Referências**
Hyndman, R.J. (1996). Computing and Graphing Highest Density Regions.
Vovk, V., Gammerman, A., & Shafer, G. (2005). Algorithmic Learning in a Random World.
Lei, J., et al. (2018). Distribution-Free Predictive Inference for Regression.
