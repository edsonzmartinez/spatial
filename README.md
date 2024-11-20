# spatial

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/r/r-original.svg" width="40" height="40"/>

## Arquivos em formato shapefile

SP_RG_Intermediarias_2022: Regiões Intermediárias do Estado de Sâo Paulo

SP_Municipios_2022: Municípios do Estado de Sâo Paulo

BR_UF_2021: Unidades da Federação

Para ler um arquivo em formato shapefile no R:

```
library(sf)
shp <- st_read("/vsizip//vsicurl/https://raw.githubusercontent.com/edsonzmartinez/spatial/main/SP_RG_Intermediarias_2022.zip")
```

Para visualizar o mapa, usando o pacote mapview:

```
library(mapview)
mapview::mapview(shp)
```

Para ler um arquivo em formato shapefile no R diretamente do site do IBGE:

```
shp <- st_read("/vsizip//vsicurl/https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2022/UFs/SP/SP_UF_2022.zip")
```
