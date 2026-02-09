---
title: Opzioni per l'ambiente
description: Opzioni per configurare l'ambiente per i corsi di formazione Nextflow
hide:
  - toc
  - footer
---

# Opzioni per l'ambiente

Il nostro obiettivo è fornire un ambiente coerente e accuratamente testato che permetta agli studenti di concentrarsi sull'apprendimento di Nextflow senza dover dedicare tempo e sforzi alla gestione del software.
A tal fine, abbiamo sviluppato un ambiente containerizzato che contiene tutto il software necessario, i file di codice e i dati di esempio per lavorare attraverso tutti i nostri corsi.

Questo ambiente containerizzato può essere eseguito immediatamente su Github Codespaces o localmente in VS Code con l'estensione Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces è un servizio basato sul web che ci permette di fornire un ambiente pre-configurato per la formazione, con tutti gli strumenti e i dati inclusi, supportato da macchine virtuali nel cloud. È accessibile gratuitamente a chiunque abbia un account Github.

    [Usa Github Codespaces:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Devcontainers locali__

    ---

    VS Code con Devcontainers fornisce un ambiente di sviluppo containerizzato eseguito localmente con tutti gli strumenti di formazione pre-configurati. Offre lo stesso ambiente pre-configurato di Codespaces ma eseguito interamente sul vostro hardware locale.

    [Usa Devcontainers localmente :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Istruzioni per l'installazione manuale

Se nessuna delle opzioni sopra soddisfa le vostre esigenze, potete replicare questo ambiente sul vostro sistema locale installando manualmente le dipendenze software e clonando il repository di formazione.

[Installazione manuale :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Deprecazione di Gitpod"

    Nextflow Training utilizzava [Gitpod](https://gitpod.io) fino a febbraio 2025.
    Tuttavia, i creatori di Gitpod hanno deciso di ritirare la funzionalità gratuita a favore del sistema [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Per questo motivo, siamo passati all'utilizzo di GitHub Codespaces, che offre anche un ambiente di sviluppo con un solo clic senza configurazione preliminare.

    A seconda di quando vi siete registrati a Gitpod e di quando esattamente ritireranno il servizio, potreste ancora essere in grado di avviare la formazione nel loro vecchio IDE cloud, anche se non possiamo garantire un accesso affidabile in futuro:
    [Apri in Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
