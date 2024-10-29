#! /usr/bin/Rscript
library(httpuv)
setwd("~/src/R/biostat3/inst/doc")
runServer(
  host = "127.0.0.1", port = 8080,
  app = list(
    staticPaths = list(
      "/" = staticPath(
        ".",
        headers = list(
          "Cross-Origin-Opener-Policy" = "same-origin",
          "Cross-Origin-Embedder-Policy" = "require-corp"
        )
      )
    )
  )
)
