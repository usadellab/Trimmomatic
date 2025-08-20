# --- Stage 1: Unpack and Prepare ---
FROM eclipse-temurin:8-jre-focal as preparer
# Re-declare ARG for this stage
ARG TRIMMOMATIC_VERSION
RUN apt-get update && apt-get install -y unzip
WORKDIR /app
# Copy the zip file
COPY target/Trimmomatic-${TRIMMOMATIC_VERSION}.zip trimmomatic.zip
RUN unzip trimmomatic.zip

# --- Stage 2: Build the Final Image ---
FROM eclipse-temurin:8-jre-focal
# Re-declare ARG for this stage
ARG TRIMMOMATIC_VERSION
WORKDIR /opt/trimmomatic
# Copy the jar file
COPY --from=preparer /app/trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar
COPY --from=preparer /app/adapters ./adapters
# MODIFIED: Specify the main class explicitly
ENTRYPOINT ["java", "-cp", "trimmomatic.jar", "org.usadellab.trimmomatic.Trimmomatic"]
CMD ["-version"]