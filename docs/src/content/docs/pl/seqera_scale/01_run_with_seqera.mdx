# Część 1: Uruchamianie pipeline'ów z interfejsu webowego

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej części kursu Scale with Seqera skonfigurujesz dostęp do Seqera Platform i uruchomisz pipeline produkcyjny z interfejsu webowego.

Upewnij się, że Twój katalog roboczy jest ustawiony na `seqera-scale/`, zgodnie z instrukcjami na stronie [Pierwsze kroki](./00_orientation.md).

---

## 1. Pierwsze kroki z Seqera

Seqera oferuje kompleksową platformę do uruchamiania, monitorowania i zarządzania pipeline'ami Nextflow.
Ta sekcja przeprowadzi Cię przez proces rejestracji i zapoznania się z platformą przed uruchomieniem pierwszego pipeline'u.

### 1.1. Załóż darmowe konto

Przejdź na stronę [cloud.seqera.io](https://cloud.seqera.io) i utwórz darmowe konto.
Możesz się zarejestrować przy użyciu adresu e-mail, konta GitHub lub Google.

Darmowe konto daje Ci dostęp do:

- **Osobistego workspace'u**: własnej przestrzeni do dodawania pipeline'ów, konfigurowania środowisk obliczeniowych i zarządzania uruchomieniami
- **Community Showcase**: starannie dobranej kolekcji pipeline'ów nf-core i społecznościowych, ze wstępnie skonfigurowanymi ustawieniami i przykładowymi danymi

Pełny przegląd poziomów kont i dostępnych funkcji znajdziesz w [dokumentacji Seqera](https://docs.seqera.io).

### 1.2. Poznaj Community Showcase

Zanim uruchomisz własne pipeline'y, poświęć chwilę na zapoznanie się z Community Showcase.
Daje ono realistyczny podgląd tego, jak wygląda platforma z prawdziwymi pipeline'ami i danymi.

1. Zaloguj się na [cloud.seqera.io](https://cloud.seqera.io).
2. Na lewym pasku bocznym kliknij **Showcase**.
3. Przejrzyj dostępne pipeline'y — rozpoznasz kilka pipeline'ów nf-core z kursu Use nf-core.
4. Kliknij wybrany pipeline, aby zobaczyć jego konfigurację i ustawienia uruchomienia.
5. Kliknij **Runs**, aby przeglądać przykładowe historię uruchomień, w tym szczegóły na poziomie zadań i raporty z poprzednich wykonań.

Jest to widok tylko do odczytu, ale pozwala zapoznać się z interfejsem przed samodzielnym uruchomieniem czegokolwiek.

### 1.3. Uzyskaj dostęp do workspace'u z zasobami obliczeniowymi

Uruchamianie pipeline'ów wymaga workspace'u ze skonfigurowanym środowiskiem obliczeniowym.

Seqera obsługuje dwa sposoby zapewnienia zasobów obliczeniowych:

- **Podłączenie własnej infrastruktury**: AWS, Azure, Google Cloud oraz planisty HPC (SLURM, LSF, PBS i inne).
  Przewodniki konfiguracji znajdziesz w [dokumentacji środowisk obliczeniowych](https://docs.seqera.io).
- **Seqera Compute**: zarządzana usługa oferująca wstępnie skonfigurowane środowiska obliczeniowe na AWS, odpłatnie, bez konieczności konfigurowania własnego konta w chmurze.
  Możesz ją aktywować bezpośrednio z ustawień workspace'u.

**Szkolenie grupowe:**
Jeśli uczestniczysz w grupowej sesji szkoleniowej, możliwe, że zostałeś(-aś) dodany(-a) do organizacji i workspace'u z już skonfigurowanymi zasobami obliczeniowymi.
Prowadzący poda Ci nazwę organizacji, workspace'u i wszelkie inne potrzebne informacje.

**Praca samodzielna:**
Jeśli realizujesz to szkolenie we własnym zakresie, musisz skonfigurować środowisko obliczeniowe w swoim osobistym workspace'ie, korzystając z jednej z powyższych opcji.
Darmowe kredyty na wypróbowanie Seqera Compute są [dostępne na życzenie](https://seqera.io/platform/compute/).

!!! note "Uwaga"

    Dalsza część kursu zakłada, że masz dostęp do workspace'u ze skonfigurowanym środowiskiem obliczeniowym.
    Jeśli uczestniczysz w grupowej sesji szkoleniowej, prowadzący potwierdzi, którego workspace'u i środowiska obliczeniowego należy używać.

### Podsumowanie

Masz konto Seqera, zapoznałeś(-aś) się z Community Showcase i masz dostęp do workspace'u z zasobami obliczeniowymi.

### Co dalej?

Uruchom pipeline RNA-seq w skali produkcyjnej z interfejsu webowego Seqera Cloud.

---

## 2. Uruchamianie nf-core/rnaseq z interfejsu webowego

Jak omówiono w kursie Use nf-core, pipeline nf-core/rnaseq to opracowany przez społeczność pipeline do analizy danych z sekwencjonowania RNA (bulk RNA-seq).

W tej sekcji dodasz pipeline do swojego workspace'u, uruchomisz go i będziesz monitorować jego wykonanie.

### 2.1. Dodaj pipeline do workspace'u

Wygodnie, nf-core/rnaseq jest częścią starannie dobranej kolekcji pipeline'ów, które można dodać do workspace'u kilkoma kliknięciami za pośrednictwem usługi Seqera Pipelines.

_Pokażemy Ci, jak dodawać własne pipeline'y w dalszej części kursu._

1. Przejdź do [**Seqera Pipelines**](https://seqera.io/pipelines), aby przeglądać kolekcję społecznościową.
2. Wyszukaj `rnaseq` i wybierz **nf-core/rnaseq**.
3. Kliknij **Launch Pipeline** lub przewiń na dół strony do sekcji **Launch Pipeline**.
4. Upewnij się, że jesteś zalogowany(-a), i wybierz odpowiednie wartości z menu rozwijanych **Organizations**, **Workspace** i **Compute Environment**.
   **Wskazówka dla grup:** Jeśli korzystasz ze współdzielonego workspace'u, dodaj unikalny identyfikator (np. swoją nazwę użytkownika) do nazwy pipeline'u.
5. Kliknij **Add pipeline to your Seqera account**.

Pojawi się okno z komunikatem: **Pipeline added: View Pipeline**.
Kliknięcie linku przeniesie Cię do wpisu pipeline'u w Twoim launchpadzie.

Pipeline jest teraz widoczny w panelu **Launchpad** Twojego workspace'u i gotowy do uruchomienia.

### 2.2. Uruchom pipeline

Kliknij przycisk **Launch** przy pipeline'ie — w panelu **Launchpad** lub na stronie szczegółów pipeline'u.
Otworzy się interfejs konfiguracji.

Pipeline jest już skonfigurowany z profilem `test`, więc dane wejściowe, katalog wyjściowy i referencja genomu są wstępnie wypełnione.
Na razie możesz zignorować pozostałe parametry i ustawienia zaawansowane.

Kliknij niebieski przycisk **Launch**, aby faktycznie rozpocząć uruchomienie.

### 2.3. Monitoruj wykonanie

Po uruchomieniu zostaniesz przeniesiony(-a) do panelu **Runs** Twojego pipeline'u.

Widok uruchomienia pokazuje:

- **Status**: aktualny stan uruchomienia (submitted, running, succeeded, failed)
- **Command line**: dokładne polecenie `nextflow run`, które platforma skonstruowała i wysłała
- **Parameters**: wszystkie wartości parametrów użyte w tym uruchomieniu
- **Tasks**: tabelę wszystkich wywołań procesów ze statusem, czasem trwania i zużyciem zasobów

Kliknij dowolny wiersz zadania, aby sprawdzić szczegóły jego wykonania, w tym:

- Skrypt `.command.sh`, który został uruchomiony
- Logi stdout i stderr
- Metryki CPU, pamięci i I/O

Zakładka **Reports** wyświetli raport MultiQC po zakończeniu uruchomienia, agregując metryki kontroli jakości dla wszystkich próbek.

Wykonanie zajmie trochę czasu, więc na razie przejdziemy dalej i wrócimy później, żeby przyjrzeć się wynikom.

### Podsumowanie

Wiesz już, jak dodać pipeline do workspace'u Seqera, skonfigurować i uruchomić go oraz monitorować wykonanie na dużą skalę.

### Co dalej?

Przejdź do [Części 2](./02_launch_from_cli.md), gdzie nauczysz się robić to wszystko z wiersza poleceń przy użyciu `tw` CLI.

---

## Podsumowanie

W tej części nauczyłeś(-aś) się:

- Zakładać konto Seqera i poznawać Community Showcase
- Dodawać pipeline z katalogu, uruchamiać go w skali produkcyjnej i monitorować wykonanie
