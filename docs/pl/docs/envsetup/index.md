---
title: Opcje środowiska
description: Opcje konfiguracji środowiska do szkoleń Nextflow
hide:
  - toc
  - footer
---

# Opcje środowiska

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Naszym celem jest zapewnienie spójnego i dokładnie przetestowanego środowiska pracy, które pozwala uczącym się skupić się na nauce Nextflow'a bez konieczności poświęcania czasu i wysiłku na zarządzanie oprogramowaniem.
W tym celu opracowaliśmy skonteneryzowane środowisko zawierające wszystkie niezbędne narzędzia, pliki kodu i przykładowe dane do pracy ze wszystkimi naszymi kursami.

To skonteneryzowane środowisko można uruchomić bezpośrednio w GitHub Codespaces lub lokalnie w VS Code z rozszerzeniem Devcontainers.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces to usługa internetowa, która pozwala nam dostarczyć prekonfigurowane środowisko do szkoleń ze wszystkimi narzędziami i danymi, wspierane przez maszyny wirtualne w chmurze. Jest dostępna bezpłatnie dla każdego z kontem GitHub.

    [Użyj GitHub Codespaces :material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Lokalne Devcontainers__

    ---

    VS Code z Devcontainers zapewnia lokalnie uruchamiane skonteneryzowane środowisko programistyczne ze wszystkimi narzędziami szkoleniowymi prekonfigurowanymi. Oferuje to samo prekonfigurowane środowisko co Codespaces, ale działa całkowicie na Twoim lokalnym sprzęcie.

    [Użyj Devcontainers lokalnie :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Instrukcje ręcznej instalacji

Jeśli żadna z powyższych opcji nie odpowiada Twoim potrzebom, możesz odtworzyć to środowisko na własnym systemie lokalnym, ręcznie instalując zależności oprogramowania i klonując repozytorium szkoleniowe.

[Ręczna instalacja :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Wycofanie Gitpod"

    Szkolenie Nextflow korzystało z [Gitpod](https://gitpod.io) do lutego 2025 roku.
    Jednak twórcy Gitpod zdecydowali się wycofać darmową funkcjonalność na rzecz systemu [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex).
    Z tego powodu przeszliśmy na GitHub Codespaces, które również oferują jednoklikowe środowisko deweloperskie bez wcześniejszej konfiguracji.

    W zależności od daty rejestracji w Gitpod i terminu wycofania usługi, możesz nadal mieć dostęp do ich starego IDE chmurowego. Nie możemy jednak zagwarantować niezawodnego dostępu w przyszłości:
    [Otwórz w Gitpod](https://gitpod.io/#https://github.com/nextflow-io/training).
