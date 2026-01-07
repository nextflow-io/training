# Makefile for nf-greeting plugin

.PHONY: build test install clean

# Build the plugin
build:
	./gradlew build

# Run tests
test:
	./gradlew test

# Install to local Nextflow plugins directory
install:
	./gradlew installPlugin

# Clean build artifacts
clean:
	./gradlew clean

# Assemble without tests
assemble:
	./gradlew assemble

# Show help
help:
	@echo "Available targets:"
	@echo "  build    - Build the plugin"
	@echo "  test     - Run unit tests"
	@echo "  install  - Install to local Nextflow plugins"
	@echo "  clean    - Clean build artifacts"
	@echo "  assemble - Assemble without running tests"
