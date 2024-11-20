# spatial

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/r/r-original.svg" width="40" height="40"/>

## Bancos de dados

SP_RG_Intermediarias_2022: Regiões Intermediárias do Estado de Sâo Paulo

Para ler o arquivo em formato shapefile no R:

```
library(sf)
shp <- st_read("/vsizip//vsicurl/https://raw.githubusercontent.com/edsonzmartinez/spatial/main/SP_RG_Intermediarias_2022.zip")
```
