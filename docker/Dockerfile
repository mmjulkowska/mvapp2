# Local + deployable image for MVApp 2.0.
# Build context is the repo root:  docker compose -f docker/docker-compose.yml up
FROM rocker/shiny:4.4.1

# App + backend runtime dependencies (add packages here as the app grows).
RUN install2.r --error --skipinstalled \
    remotes DT ggplot2 corrplot plotly RColorBrewer ggsci \
    pheatmap dendextend multcompView cluster \
 && rm -rf /tmp/downloaded_packages

# Install the mvapp backend package from source.
COPY mvapp /srv/mvapp
RUN R -e "remotes::install_local('/srv/mvapp', dependencies = TRUE, upgrade = 'never')"

# Deploy the Shiny app.
RUN rm -rf /srv/shiny-server/*
COPY app /srv/shiny-server/

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
