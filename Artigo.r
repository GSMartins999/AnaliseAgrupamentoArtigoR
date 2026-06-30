# ==============================================================================
# Artigo.r — Pipeline de Análise: Dengue × Saneamento e Clima
# ==============================================================================

# --- 1. SETUP E CARREGAMENTO DE PACOTES ---
pacotes <- c(
  "httr", "jsonlite", "sidrar", "dplyr", "tidyr", "readr", "purrr", 
  "stringr", "lubridate", "cluster", "factoextra", "NbClust", 
  "ggplot2", "corrplot", "geobr", "sf"
)

# Configuração de biblioteca local para evitar problemas de permissão
lib_usuario <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib_usuario)) {
  lib_usuario <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library",
                           paste(R.version$major, sub("\\..*", "", R.version$minor), sep = "."))
}
if (!dir.exists(lib_usuario)) dir.create(lib_usuario, recursive = TRUE)
.libPaths(c(lib_usuario, .libPaths()))

# Instalar pacotes faltantes e carregar
nao_instalados <- pacotes[!pacotes %in% installed.packages()[, "Package"]]
if (length(nao_instalados) > 0) {
  install.packages(nao_instalados, lib = lib_usuario, dependencies = TRUE, repos = "https://cloud.r-project.org")
}
invisible(lapply(pacotes, library, character.only = TRUE))

theme_set(theme_minimal(base_size = 12))
options(scipen = 999)

dir_dados <- file.path(getwd(), "dados")
if (!dir.exists(dir_dados)) dir.create(dir_dados)

# --- 2. COLETA DE DADOS (InfoDengue e IBGE/SIDRA) ---
cat("\n📡 Identificando estado com maior número de casos de dengue em 2025...\n")
capitais <- read_capitals() %>%
  st_drop_geometry() %>%
  transmute(uf_sigla = abbrev_state, geocode = as.character(code_muni), municipio = name_muni)

resultados_capitais <- list()
for (i in seq_len(nrow(capitais))) {
  uf <- capitais$uf_sigla[i]
  gc <- capitais$geocode[i]
  nom <- capitais$municipio[i]
  
  url <- paste0("https://info.dengue.mat.br/api/alertcity?geocode=", gc,
                "&disease=dengue&format=json&ew_start=1&ew_end=52&ey_start=2025&ey_end=2025")
  
  resp <- tryCatch(httr::GET(url, httr::timeout(15)), error = function(e) NULL)
  if (!is.null(resp) && httr::status_code(resp) == 200) {
    dados_api <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
    if (!is.null(dados_api) && nrow(dados_api) > 0) {
      resultados_capitais[[uf]] <- data.frame(
        uf_sigla = uf, municipio = nom, geocode = gc,
        total_casos = sum(dados_api$casos, na.rm = TRUE)
      )
    }
  }
  Sys.sleep(0.1)
}

capital_max <- bind_rows(resultados_capitais) %>% arrange(desc(total_casos)) %>% slice(1)
estado.alvo <- capital_max$uf_sigla
cat("-> Estado selecionado:", estado.alvo, "(Capital:", capital_max$municipio, ")\n")

# Obter municípios do estado
municipios_br <- jsonlite::fromJSON(httr::content(httr::GET("https://servicodados.ibge.gov.br/api/v1/localidades/municipios?view=nivelado"), "text", encoding = "UTF-8"))
municipios.estado <- municipios_br %>%
  filter(`UF-sigla` == estado.alvo) %>%
  transmute(geocode = as.character(`municipio-id`), municipio = `municipio-nome`)

# Coleta paralela de dados de dengue (InfoDengue)
cat("📡 Coletando dados de dengue para os", nrow(municipios.estado), "municípios de", estado.alvo, "...\n")
pool <- curl::new_pool(total_con = 100, host_con = 15)
resultados.lista <- list()

for (i in seq_len(nrow(municipios.estado))) {
  gc <- municipios.estado$geocode[i]
  nom <- municipios.estado$municipio[i]
  url <- paste0("https://info.dengue.mat.br/api/alertcity?geocode=", gc,
                "&disease=dengue&format=json&ew_start=1&ew_end=52&ey_start=2025&ey_end=2025")
  
  local({
    geocode_local <- gc
    nome_local <- nom
    curl::curl_fetch_multi(url, pool = pool, done = function(res) {
      tryCatch({
        dados_api <- jsonlite::fromJSON(rawToChar(res$content))
        if (!is.null(dados_api) && nrow(dados_api) > 0) {
          dados_api$geocode <- geocode_local
          dados_api$municipio <- nome_local
          resultados.lista[[geocode_local]] <<- dados_api
        }
      }, error = function(e) NULL)
    })
  })
}
curl::multi_run(pool = pool)

dados.mun <- bind_rows(resultados.lista) %>%
  mutate(across(c(casos, casos_est, tempmin, tempmed, tempmax, umidmin, umidmed, umidmax), as.numeric)) %>%
  group_by(geocode, municipio) %>%
  summarise(
    casos_total = sum(casos, na.rm = TRUE),
    casos_est = sum(casos_est, na.rm = TRUE),
    temp_min_media = mean(tempmin, na.rm = TRUE),
    temp_med_media = mean(tempmed, na.rm = TRUE),
    temp_max_media = mean(tempmax, na.rm = TRUE),
    umid_min_media = mean(umidmin, na.rm = TRUE),
    umid_med_media = mean(umidmed, na.rm = TRUE),
    umid_max_media = mean(umidmax, na.rm = TRUE),
    .groups = "drop"
  )

# Coleta de dados socioeconômicos do Censo (SIDRA)
codigo.uf <- substr(municipios.estado$geocode[1], 1, 2)
cat("📡 Coletando dados socioeconômicos do Censo para a UF", codigo.uf, "...\n")

# População (Tabela 6579)
pop <- get_sidra(x = 6579, api = paste0("/t/6579/n6/in%20n3%20", codigo.uf, "/v/9324/p/last%201")) %>%
  transmute(geocode = substr(`Município (Código)`, 1, 7), populacao = Valor)

# PIB (Tabela 5938) - Nota: A tabela do IBGE fornece o PIB Bruto (em milhares de reais). 
# Calcularemos o PIB per capita dividindo pela população na etapa de consolidação.
pib <- get_sidra(x = 5938, api = paste0("/t/5938/n6/in%20n3%20", codigo.uf, "/v/37/p/last%201")) %>%
  transmute(geocode = substr(`Município (Código)`, 1, 7), pib_bruto_mil = Valor)

# Esgotamento sanitário (Tabela 6805)
esgoto_bruto <- get_sidra(x = 6805, api = paste0("/t/6805/n6/in%20n3%20", codigo.uf, "/v/381/p/2022/c11558/all")) %>%
  transmute(geocode = substr(`Município (Código)`, 1, 7), tipo = `Tipo de esgotamento sanitário`, valor = Valor)
esgoto_total <- esgoto_bruto %>% filter(tipo == "Total") %>% select(geocode, dom_total = valor)
esgoto_rede <- esgoto_bruto %>% filter(grepl("^Rede geral, rede pluvial", tipo, ignore.case = TRUE)) %>%
  group_by(geocode) %>% summarise(dom_esgoto_rede = sum(valor, na.rm = TRUE), .groups = "drop")
esgoto <- esgoto_total %>% left_join(esgoto_rede, by = "geocode") %>%
  transmute(geocode, pct_esgoto = (coalesce(dom_esgoto_rede, 0) / dom_total) * 100)

# Abastecimento de água (Tabela 6803)
agua_bruto <- get_sidra(x = 6803, api = paste0("/t/6803/n6/in%20n3%20", codigo.uf, "/v/381/p/2022/c1821/all")) %>%
  transmute(geocode = substr(`Município (Código)`, 1, 7), tipo = `Existência de ligação à rede geral de distribuição de água e principal forma de abastecimento de água`, valor = Valor)
agua_total <- agua_bruto %>% filter(tipo == "Total") %>% select(geocode, dom_total = valor)
agua_rede <- agua_bruto %>% filter(grepl("^Possui liga.*rede geral e a utiliza", tipo, ignore.case = TRUE)) %>%
  group_by(geocode) %>% summarise(dom_agua_rede = sum(valor, na.rm = TRUE), .groups = "drop")
agua <- agua_total %>% left_join(agua_rede, by = "geocode") %>%
  transmute(geocode, pct_agua = (coalesce(dom_agua_rede, 0) / dom_total) * 100)

# Coleta de lixo (Tabela 6892)
lixo_bruto <- get_sidra(x = 6892, api = paste0("/t/6892/n6/in%20n3%20", codigo.uf, "/v/381/p/2022/c67/all")) %>%
  transmute(geocode = substr(`Município (Código)`, 1, 7), tipo = `Destino do lixo`, valor = Valor)
lixo_total <- lixo_bruto %>% filter(tipo == "Total") %>% select(geocode, dom_total = valor)
lixo_coleta <- lixo_bruto %>% filter(tipo == "Coletado") %>%
  group_by(geocode) %>% summarise(dom_lixo_coletado = sum(valor, na.rm = TRUE), .groups = "drop")
lixo <- lixo_total %>% left_join(lixo_coleta, by = "geocode") %>%
  transmute(geocode, pct_lixo_coletado = (coalesce(dom_lixo_coletado, 0) / dom_total) * 100)

# Área territorial e malha geográfica (geobr)
cat("🗺️ Carregando malha geográfica via geobr...\n")
mapa <- read_municipality(code_muni = estado.alvo, year = 2022)
geom_col <- if ("geom" %in% names(mapa)) "geom" else "geometry"
area.mun <- mapa %>% st_drop_geometry() %>%
  transmute(geocode = as.character(code_muni), area_km2 = as.numeric(st_area(mapa[[geom_col]])) / 1e6)

# Consolidação final do dataset com cálculo real do PIB per capita
dataset_final <- dados.mun %>%
  left_join(pop, by = "geocode") %>%
  left_join(pib, by = "geocode") %>%
  left_join(esgoto, by = "geocode") %>%
  left_join(agua, by = "geocode") %>%
  left_join(lixo, by = "geocode") %>%
  left_join(area.mun, by = "geocode") %>%
  filter(complete.cases(.)) %>%
  mutate(
    pib_percapita = (pib_bruto_mil * 1000) / populacao, # Correção de total PIB para per capita
    incidencia = (casos_total / populacao) * 100000,
    dens_demog = populacao / area_km2
  )

write_csv(dataset_final, "dados/dataset_final.csv")
saveRDS(estado.alvo, "dados/estado_alvo.rds")

# --- 3. ANÁLISE DE AGRUPAMENTO (CLUSTER) ---
cat("\n📊 Executando análise de agrupamento hierárquico (Ward.D2)...\n")

variaveis <- dataset_final %>%
  select(incidencia, pib_percapita, dens_demog, pct_esgoto, pct_agua, pct_lixo_coletado, temp_med_media, umid_med_media)

variaveis_mat <- as.data.frame(variaveis)
rownames(variaveis_mat) <- dataset_final$municipio
dados.p <- scale(variaveis_mat)

# Matriz de Distância e Clusters
d.eucl <- dist(dados.p, method = "euclidean")
metod.ward <- hclust(d.eucl, method = "ward.D2")
cat(sprintf("  - Coeficiente de Correlação Cofenética: %.4f\n", cor(d.eucl, cophenetic(metod.ward))))

# Classificação em 3 grupos
k <- 3
grupo <- cutree(metod.ward, k = k)
dataset_final$grupo <- grupo

# Perfil Médio dos Clusters
cat("\n📊 Perfil descritivo médio por Cluster:\n")
perfil <- dataset_final %>%
  group_by(grupo) %>%
  summarise(
    n_municipios      = n(),
    casos_total       = sum(casos_total, na.rm = TRUE),
    incidencia        = mean(incidencia, na.rm = TRUE),
    pib_percapita     = mean(pib_percapita, na.rm = TRUE),
    dens_demog        = mean(dens_demog, na.rm = TRUE),
    pct_esgoto        = mean(pct_esgoto, na.rm = TRUE),
    pct_agua          = mean(pct_agua, na.rm = TRUE),
    pct_lixo_coletado = mean(pct_lixo_coletado, na.rm = TRUE),
    temp_media        = mean(temp_med_media, na.rm = TRUE),
    umid_media        = mean(umid_med_media, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(perfil), row.names = FALSE)

# --- 4. EXPORTAÇÃO DOS GRÁFICOS (PNG) ---
cat("\n💾 Exportando gráficos e figuras...\n")

# 4.1. Método do Cotovelo (Elbow)
png("figura1_cotovelo.png", width = 800, height = 600, res = 120)
print(fviz_nbclust(dados.p, kmeans, method = "wss") + 
        labs(title = "Método do Cotovelo (Soma dos Quadrados)") + 
        theme_minimal(base_size = 12))
dev.off()

# 4.2. Método da Silhueta (Silhouette)
png("figura2_silhueta.png", width = 800, height = 600, res = 120)
print(fviz_nbclust(dados.p, kmeans, method = "silhouette") + 
        labs(title = "Método da Silhueta (Largura Média)") + 
        theme_minimal(base_size = 12))
dev.off()

# 4.3. Dendrograma Ward
png("figura3_dendrograma.png", width = 1200, height = 800, res = 120)
print(fviz_dend(metod.ward, k = k, cex = 0.35,
                    k_colors = c("#2E9FDF", "#00BB0C", "#E7B800"),
                    color_labels_by_k = TRUE, rect = TRUE,
                    main = "Dendrograma de Ward - Clusters de Municípios"))
dev.off()

# 4.4. Projeção PCA
km.res <- kmeans(dados.p, centers = k, nstart = 25)
png("figura4_pca.png", width = 800, height = 600, res = 120)
print(fviz_cluster(km.res, data = dados.p,
                      palette = c("#2E9FDF", "#00BB0C", "#E7B800"),
                      main = "Clusters de Municípios (Projeção PCA)",
                      ggtheme = theme_minimal(base_size = 12)))
dev.off()

# 4.5. Matriz de Correlação
png("figura5_correlacao.png", width = 800, height = 800, res = 120)
corrplot(cor(variaveis), method = "number", type = "upper",
         title = "Matriz de Correlação Linear",
         mar = c(1, 1, 3, 1), tl.col = "black")
dev.off()

# 4.6. Mapeamento Espacial dos Clusters
mapa$geocode <- as.character(mapa$code_muni)
mapa.dados <- mapa %>% left_join(dataset_final %>% select(geocode, grupo), by = "geocode")

p_mapa <- ggplot(mapa.dados) +
  geom_sf(aes(fill = as.factor(grupo)), color = "white", size = 0.1) +
  scale_fill_manual(
    values = c("#2E9FDF", "#00BB0C", "#E7B800"),
    name = "Cluster/Grupo",
    na.value = "grey80"
  ) +
  labs(
    title = paste("Mapeamento dos Clusters em", estado.alvo),
    subtitle = "Agrupamento com base em incidência de dengue, clima e saneamento"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text     = element_blank(),
    axis.ticks    = element_blank()
  )

ggsave("figura6_mapa.png", plot = p_mapa, width = 8, height = 6, dpi = 300)

cat("\n✅ Execução completa concluída com sucesso! Todos os arquivos foram unificados e gráficos salvos.\n")
