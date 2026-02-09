# Lokalne Devcontainery

Jeśli masz lokalną instalację Dockera lub chcesz ją zainstalować, najłatwiejszym sposobem pracy z tymi materiałami lokalnie jest użycie funkcji devcontainer w Visual Studio Code. To podejście zapewnia wszystkie niezbędne narzędzia i zależności bez konieczności ręcznej instalacji.

## Wymagania

Aby korzystać z lokalnej konfiguracji devcontainera, będziesz potrzebować:

- [Visual Studio Code](https://code.visualstudio.com/)
- Lokalnej instalacji Dockera, na przykład:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (dla Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (dla Linuksa)
  - [Colima](https://github.com/abiosoft/colima) (alternatywa dla macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (dołączony do Docker Desktop, ale może wymagać osobnej instalacji w innych konfiguracjach Dockera)
- [Rozszerzenia Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) dla VS Code

Twoja instalacja Dockera musi być uruchomiona, zanim spróbujesz otworzyć devcontainer.

Aby sprawdzić, czy Docker buildx jest dostępny, uruchom:

```bash
docker buildx version
```

Jeśli to polecenie zakończy się niepowodzeniem, musisz zainstalować rozszerzenie buildx przed kontynuowaniem.

## Instrukcje konfiguracji

Wykonaj poniższe kroki, aby skonfigurować swoje lokalne środowisko przy użyciu devcontainerów VS Code:

### Zainstaluj rozszerzenie "Dev Containers" w VS Code

- Otwórz VS Code
- Przejdź do Rozszerzeń (Ctrl+Shift+X lub Cmd+Shift+X na macOS)
- Wyszukaj "Dev Containers"
- Kliknij "Install"

![Instalowanie rozszerzenia Dev Containers w VS Code](img/install_extension.png)

### Sklonuj repozytorium:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Otwórz repozytorium w VS Code:

- Uruchom VS Code
- Wybierz **File -> Open Folder** z menu
- Przejdź do folderu repozytorium szkoleniowego, który właśnie sklonowałeś, i wybierz go
- Kliknij **Open**

### Otwórz ponownie w kontenerze

Jeśli VS Code wyświetli monit "Reopen in Container", kliknij go. Alternatywnie:

- Naciśnij F1 (lub Ctrl+Shift+P / Cmd+Shift+P na macOS)
- Wpisz "Dev Containers: Reopen in Container"
- **Ważne**: Gdy zostaniesz poproszony o wybranie konfiguracji, wybierz konfigurację devcontainera **local-dev**

![Monit o ponowne otwarcie w kontenerze](img/reopen_prompt.png)

![Wybieranie konfiguracji lokalnej](img/select_local_config.png)

Poczekaj, aż kontener zostanie zbudowany. Może to zająć kilka minut przy pierwszym uruchomieniu, ponieważ pobiera i konfiguruje wszystkie niezbędne komponenty.

Gdy kontener zostanie zbudowany i uruchomiony, będziesz mieć w pełni skonfigurowane środowisko ze wszystkimi niezbędnymi narzędziami zainstalowanymi, w tym:

- Java
- Nextflow
- Docker
- Git
- Oraz wszystkimi innymi zależnościami wymaganymi do szkolenia

![VS Code z uruchomionym devcontainerem](img/running_container.png)

## Korzyści z używania Devcontainerów

Korzystanie z podejścia devcontainerowego oferuje kilka zalet:

- **Spójność**: Zapewnia spójne środowisko programistyczne na różnych maszynach
- **Prostota**: Wszystkie zależności są wstępnie zainstalowane i skonfigurowane
- **Izolacja**: Środowisko programistyczne jest odizolowane od Twojego systemu lokalnego
- **Powtarzalność**: Każdy korzystający z devcontainera otrzymuje tę samą konfigurację
- **Brak ręcznej instalacji**: Nie ma potrzeby ręcznego instalowania Javy, Nextflow'a i innych narzędzi

## Sprawdzanie środowiska

Gdy Twój devcontainer jest uruchomiony, możesz sprawdzić, czy wszystko jest poprawnie skonfigurowane, uruchamiając:

```bash
nextflow info
```

Powinno to wyświetlić wersję Nextflow'a i informacje o środowisku uruchomieniowym, potwierdzając, że Twoje środowisko jest prawidłowo skonfigurowane.

## Rozwiązywanie problemów

Jeśli napotkasz problemy z konfiguracją devcontainera:

1. Upewnij się, że Twoja instalacja Dockera (Docker Desktop, Colima, Docker Engine itp.) jest uruchomiona przed otwarciem devcontainera
2. Sprawdź, czy wybrałeś konfigurację **local-dev**, gdy zostałeś o to poproszony
3. Zweryfikuj, że Docker buildx jest zainstalowany i działa, uruchamiając `docker buildx version`
4. Jeśli kontener nie może zostać zbudowany, spróbuj go przebudować, uruchamiając polecenie "Dev Containers: Rebuild Container"
5. W przypadku uporczywych problemów zapoznaj się z [przewodnikiem rozwiązywania problemów VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
