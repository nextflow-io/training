---
description: Stop the running MkDocs preview server
---

Stop the MkDocs preview server that is running in the background.

This will:

1. Find any running Docker containers with the training-mkdocs image
2. Stop the container(s)
3. Confirm the server has been stopped

If no preview server is running, you'll be notified.

Use the following steps:

1. Run `docker ps --filter "ancestor=ghcr.io/nextflow-io/training-mkdocs:latest" --format "{{.ID}}"` to find running containers
2. If any container IDs are found, stop them with `docker stop <container_id>`
3. Confirm they have been stopped
