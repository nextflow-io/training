# Część 3: Organizowanie danych wejściowych

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W Części 2 uruchomiliśmy molkart z wieloma parametrami w wierszu poleceń.
Teraz poznamy dwa lepsze podejścia do zarządzania wejściami: **pliki parametrów** oraz **arkusze próbek**.

## 1. Używanie plików parametrów

### 1.1. Problem z długimi wierszami poleceń

Przypomnijmy sobie nasze polecenie z Części 2:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

To działa, ale trudno to odtworzyć, udostępnić lub zmodyfikować.
Co jeśli musisz ponownie uruchomić tę samą analizę za miesiąc?
Co jeśli współpracownik chce użyć dokładnie Twoich ustawień?

### 1.2. Rozwiązanie: Użyj pliku parametrów

Utwórz plik o nazwie `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Teraz Twoje polecenie staje się:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

To wszystko! Plik parametrów dokumentuje dokładną konfigurację i ułatwia ponowne uruchomienie lub udostępnienie.

### 1.3. Nadpisywanie parametrów

Nadal możesz nadpisać konkretne parametry z wiersza poleceń:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

Powyższa linia zmienia `segmentation_method` na `stardist` oraz nazwę `--outdir` na `stardist_results` zamiast parametrów z pliku `params.yaml`.
Dodatkowo możesz zauważyć, że flaga `-resume` pozwoliła nam ponownie użyć wyników przetwarzania wstępnego z poprzedniego uruchomienia, oszczędzając czas.
Możesz użyć tego wzorca, aby szybko testować różne warianty pipeline'u.

### Podsumowanie

Pliki parametrów sprawiają, że Twoje analizy są odtwarzalne i łatwe do udostępnienia.
Używaj ich do każdej rzeczywistej pracy analitycznej.

### Co dalej?

Dowiedz się, jak arkusze próbek organizują informacje o wielu próbkach.

---

## 2. Wzorzec arkusza próbek

### 2.1. Czym jest arkusz próbek?

Arkusz próbek to plik CSV, który opisuje Twoje próbki wejściowe.
Każdy wiersz to próbka, a kolumny określają pliki i metadane dla tej próbki.

Spójrzmy na arkusz próbek, którego używaliśmy:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

Kolumny to:

- `sample`: Unikalny identyfikator próbki
- `nuclear_image`: Obraz barwienia jądrowego (TIFF)
- `spot_table`: Punkty transkryptów (TXT)
- `membrane_image`: Obraz barwienia błonowego (TIFF, opcjonalny)

### 2.2. Ścieżki do plików

Arkusze próbek akceptują wiele typów ścieżek:

- **URL**: Nextflow pobiera automatycznie (jak pokazano powyżej)
- **Ścieżki lokalne**: `data/nuclear.tiff` lub `/absolute/path/to/nuclear.tiff`
- **Przechowywanie w chmurze**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

Możesz mieszać typy ścieżek w tym samym arkuszu próbek.

### 2.3. Tworzenie własnego arkusza próbek

Najpierw pobierzmy pliki danych testowych lokalnie:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Teraz zmodyfikujmy arkusz próbek, aby odwoływał się do tych lokalnych plików:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! warning "Ostrzeżenie"

    Zauważ, że ścieżki w arkuszu próbek są względne względem miejsca, **z którego uruchamiasz** Nextflow'a, a nie względem miejsca, w którym znajduje się arkusz próbek.

Na koniec wykonajmy nf-core/molkart jeszcze raz z arkuszem próbek zawierającym lokalne ścieżki do plików:

`nextflow run ./molkart -params-file params.yaml -resume`

Jak widzisz, Nextflow wykonuje to uruchomienie podobnie jak wtedy, gdy pliki były pobierane z Github.
To jedna ze świetnych funkcji Nextflow'a - prawidłowo przygotowuje dane dla Ciebie, niezależnie od tego, gdzie się znajdują.

### Podsumowanie

Arkusze próbek organizują zestawy danych z wieloma próbkami w sposób umożliwiający jawne zdefiniowanie metadanych wraz ze ścieżkami do plików.
Większość pipeline'ów nf-core używa tego wzorca.

### Co dalej?

Teraz, gdy omówiliśmy dane wejściowe, zbadajmy, jak skonfigurować pipeline'y Nextflow dla różnych środowisk obliczeniowych.
