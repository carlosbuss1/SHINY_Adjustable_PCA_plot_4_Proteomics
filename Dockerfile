# Base image
FROM rocker/shiny:latest

# Set the working directory inside the container
WORKDIR /srv/shiny-server/app

# Copy the application files into the container
COPY proteomics_counts.csv /srv/shiny-server/app/proteomics_counts.csv
COPY shiny_pca_plot.R /srv/shiny-server/app/app.R

# Install required R packages
RUN R -e "install.packages(c(\"shiny\", \"ggplot2\", \"ggrepel\", \"limma\", \"edgeR\"), repos = 'https://cran.rstudio.com/')"

# Expose the default Shiny server port
EXPOSE 3838

# Set up permissions for Shiny to access the app directory
RUN chown -R shiny:shiny /srv/shiny-server

# Switch to the Shiny user
USER shiny

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app', host = '0.0.0.0', port = 3838)"]

