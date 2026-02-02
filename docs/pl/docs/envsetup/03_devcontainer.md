# Lokalne Devcontainers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Jeśli masz lokalną instalację Docker lub chętnie ją skonfigurujesz, najłatwiejszym sposobem pracy lokalnej z tymi materiałami jest użycie funkcji devcontainer w Visual Studio Code. To podejście zapewnia wszystkie niezbędne narzędzia i zależności bez konieczności ręcznej konfiguracji.

## Wymagania

Aby użyć lokalnej konfiguracji devcontainer, będziesz potrzebować:

- [Visual Studio Code](https://code.visualstudio.com/)
- Lokalnej instalacji Docker, na przykład:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (dla Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (dla Linux)
  - [Colima](https://github.com/abiosoft/colima) (alternatywa dla macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (dołączony do Docker Desktop, ale może wymagać osobnej instalacji z innymi konfiguracjami Docker)
- [Rozszerzenie Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) dla VS Code

Twoja instalacja Docker musi być uruchomiona przed próbą otwarcia devcontainer.

Aby sprawdzić, czy Docker buildx jest dostępny, uruchom:

```bash
docker buildx version
```

Jeśli to polecenie nie powiedzie się, będziesz musiał zainstalować rozszerzenie buildx przed kontynuowaniem.

## Instrukcje konfiguracji

Wykonaj następujące kroki, aby skonfigurować lokalne środowisko za pomocą devcontainers VS Code:

### Zainstaluj rozszerzenie "Dev Containers" w VS Code

- Otwórz VS Code
- Przejdź do Extensions (Ctrl+Shift+X lub Cmd+Shift+X na macOS)
- Wyszukaj "Dev Containers"
- Kliknij "Install"

![Installing Dev Containers extension in VS Code](img/install_extension.png)

### Sklonuj repozytorium:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Otwórz repozytorium w VS Code:

- Uruchom VS Code
- Wybierz **File -> Open Folder** z menu
- Przejdź do i wybierz folder repozytorium szkoleniowego, który właśnie sklonowałeś
- Kliknij **Open**

### Otwórz ponownie w kontenerze

Jeśli VS Code wyświetli monit "Reopen in Container", kliknij go. Alternatywnie:

- Naciśnij F1 (lub Ctrl+Shift+P / Cmd+Shift+P na macOS)
- Wpisz "Dev Containers: Reopen in Container"
- **Ważne**: Gdy pojawi się monit o wybór konfiguracji, wybierz konfigurację devcontainer **local-dev**

![Reopen in Container prompt](img/reopen_prompt.png)

![Selecting local configuration](img/select_local_config.png)

Poczekaj, aż kontener się zbuduje. Za pierwszym razem może to zająć kilka minut, ponieważ pobiera i konfiguruje wszystkie niezbędne komponenty.

Po zbudowaniu i uruchomieniu kontenera będziesz mieć w pełni skonfigurowane środowisko ze wszystkimi niezbędnymi narzędziami, w tym:

- Java
- Nextflow
- Docker
- Git
- I wszystkie inne zależności wymagane do szkolenia

![VS Code with devcontainer running](img/running_container.png)

## Zalety korzystania z Devcontainers

Korzystanie z podejścia devcontainer oferuje kilka zalet:

- **Spójność**: Zapewnia spójne środowisko programistyczne na różnych maszynach
- **Prostota**: Wszystkie zależności są preinstalowane i skonfigurowane
- **Izolacja**: Środowisko programistyczne jest odizolowane od Twojego lokalnego systemu
- **Powtarzalność**: Każdy korzystający z devcontainer otrzymuje taką samą konfigurację
- **Brak ręcznej konfiguracji**: Nie ma potrzeby samodzielnego instalowania Java, Nextflow i innych narzędzi

## Sprawdzanie środowiska

Po uruchomieniu devcontainer możesz sprawdzić, czy wszystko jest poprawnie skonfigurowane, uruchamiając:

```bash
nextflow info
```

Powinno to wyświetlić wersję Nextflow i informacje o środowisku wykonawczym, potwierdzając, że Twoje środowisko jest poprawnie skonfigurowane.

## Rozwiązywanie problemów

Jeśli napotkasz problemy z konfiguracją devcontainer:

1. Upewnij się, że Twoja instalacja Docker (Docker Desktop, Colima, Docker Engine itp.) jest uruchomiona przed otwarciem devcontainer
2. Sprawdź, czy wybrałeś konfigurację **local-dev**, gdy pojawił się monit
3. Sprawdź, czy Docker buildx jest zainstalowany i działa, uruchamiając `docker buildx version`
4. Jeśli kontener nie zbuduje się, spróbuj go przebudować, uruchamiając polecenie "Dev Containers: Rebuild Container"
5. W przypadku trwałych problemów zapoznaj się z [przewodnikiem rozwiązywania problemów VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
