# Część 4: Dodawanie testów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W pierwszej części tego kursu zbudowałeś pipeline do wykrywania wariantów, który był całkowicie liniowy i przetwarzał dane każdej próbki niezależnie od innych.

W drugiej części pokazaliśmy Ci, jak używać kanałów i operatorów kanałów do implementacji wspólnego wykrywania wariantów za pomocą GATK.

W trzeciej części zmodularyzowaliśmy pipeline.

W tej części szkolenia pokażemy Ci, jak używać [**nf-test**](https://www.nf-test.com/), frameworka testowego, który dobrze integruje się z Nextflow i ułatwia dodawanie testów zarówno na poziomie modułów, jak i workflow'ów do Twojego pipeline'u. Aby podążać za tą częścią szkolenia, musisz ukończyć Część 1, Część 2 i Część 3, a także [zadanie dodatkowe nf-test](../../side_quests/nf-test.md), które omawia podstawy nf-test i wyjaśnia, dlaczego testowanie jest ważne.

---

## 0. Rozgrzewka

!!! note "Uwaga"

    Upewnij się, że jesteś w poprawnym katalogu roboczym:
    `cd /workspaces/training/nf4-science/genomics`

Jeśli przepracowałeś poprzednie części tego kursu szkoleniowego, masz działającą wersję pipeline genomiki z odpowiednią strukturą katalogów modułów.

??? abstract "Zawartość katalogów"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

Ten katalog modułów można znaleźć w katalogu `solutions`, jeśli jest potrzebny.

Rozpoczniemy od tego samego workflow'u co w Części 3, który przygotowaliśmy dla Ciebie w pliku `genomics-4.nf`. Dokładnie tak jak w [zadaniu dodatkowym nf-test](../../side_quests/nf-test.md), dodamy kilka różnych typów testów do trzech procesów w tym pipeline'ie, a także test na poziomie workflow'u.

### 0.1. Sprawdź, czy workflow działa

Zanim zaczniemy dodawać testy, upewnij się, że workflow działa zgodnie z oczekiwaniami.

```bash
nextflow run genomics-4.nf -resume
```

Powinno to wyglądać bardzo znajomo, jeśli pracowałeś nad tym kursem szkoleniowym od początku.

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Jak poprzednio, w Twoim katalogu projektu znajdzie się teraz katalog `work` i katalog `results_genomics`. Faktycznie wykorzystamy te wyniki później w naszych testach. Ale od teraz będziemy używać pakietu `nf-test` do testowania pipeline'u.

### 0.2. Zainicjuj `nf-test`

Tak jak w [zadaniu dodatkowym nf-test](../../side_quests/nf-test.md), musimy zainicjować pakiet `nf-test`.

```bash
nf-test init
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "Zawartość nf-test.config"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Tworzy to również katalog `tests` zawierający szkielet pliku konfiguracyjnego.

### Podsumowanie

Teraz jesteśmy gotowi, aby zacząć pisać testy dla naszego pipeline genomiki.

### Co dalej?

Napisz podstawowe testy, które oceniają, czy wywołania procesów zakończyły się sukcesem i wygenerowały poprawne wyjścia.

---

## 1. Testowanie procesu pod kątem sukcesu i zgodności wyjść

Zaczniemy od testowania procesu `SAMTOOLS_INDEX`, który tworzy pliki indeksów dla plików BAM, aby umożliwić efektywny losowy dostęp. Jest to dobry pierwszy przypadek testowy, ponieważ:

1. Ma pojedyncze, dobrze zdefiniowane wejście (plik BAM)
2. Produkuje przewidywalne wyjście (plik indeksu BAI)
3. Wyjście powinno być identyczne dla identycznych wejść

### 1.1. Wygeneruj szkielet pliku testowego

Najpierw wygeneruj szkielet pliku testowego:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Wynik polecenia"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

To tworzy plik w tym samym katalogu co `main.nf`.
Możesz przejść do katalogu w eksploratorze plików i otworzyć plik, który powinien zawierać następujący kod:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Początkowe asercje powinny być znajome z [zadania dodatkowego nf-test](../../side_quests/nf-test.md):

- `assert process.success` stwierdza, że oczekujemy, że proces uruchomi się pomyślnie i zakończy bez żadnych błędów.
- `snapshot(process.out).match()` stwierdza, że oczekujemy, że wynik uruchomienia będzie identyczny z wynikiem uzyskanym w poprzednim uruchomieniu (jeśli dotyczy).
  Omawiamy to bardziej szczegółowo później.

Używając tego jako punktu wyjścia, musimy dodać odpowiednie wejścia testowe dla procesu samtools index oraz ewentualne parametry.

### 1.2. Przenieś plik testowy i zaktualizuj ścieżkę skryptu

Zanim przystąpimy do pracy nad uzupełnieniem testu, musimy przenieść plik do jego ostatecznej lokalizacji. Część powodu, dla którego dodaliśmy katalog dla każdego modułu, jest to, że możemy teraz dostarczać testy w katalogu `tests` znajdującym się w tym samym miejscu co plik `main.nf` każdego modułu. Utwórz ten katalog i przenieś tam plik testowy.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Teraz możemy uprościć sekcję `script` pliku testowego do ścieżki względnej:

=== "Po zmianach"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Przed zmianami"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

To informuje test, gdzie znaleźć plik `main.nf` modułu, bez konieczności określania pełnej ścieżki.

### 1.3. Dostarcz wejścia testowe dla SAMTOOLS_INDEX

Plik szkieletowy zawiera symbol zastępczy, który musimy zastąpić faktycznym wejściem testowym, odpowiednim dla wejścia `samtools index`. Odpowiednim wejściem jest plik BAM, który mamy dostępny w katalogu `data/bam`.

=== "Po zmianach"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Przed zmianami"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Nazwij test na podstawie funkcjonalności

Jak nauczyliśmy się wcześniej, dobrą praktyką jest zmiana nazwy testu na coś, co ma sens w kontekście testu.

=== "Po zmianach"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    To przyjmuje dowolny ciąg znaków, więc możemy umieścić cokolwiek chcemy.
    Tutaj wybieramy odniesienie do nazwy pliku i jego formatu.

=== "Przed zmianami"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Uruchom test i sprawdź wynik

Uruchom test:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Jak nauczyliśmy się wcześniej, to zweryfikowało podstawową asercję o sukcesie procesu i utworzyło plik migawki na podstawie wyjścia procesu. Możemy zobaczyć zawartość pliku migawki w pliku `tests/modules/samtools/index/tests/main.nf.test.snap`:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

Możemy również uruchomić test ponownie i zobaczyć, że przechodzi, ponieważ wyjście jest identyczne z migawką:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. Dodaj więcej testów do `SAMTOOLS_INDEX`

Czasami przydatne jest testowanie różnych plików wejściowych, aby upewnić się, że testujemy różne potencjalne problemy. Dodaj testy dla plików BAM matki i ojca z tria z naszych danych testowych.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Następnie możesz uruchomić test ponownie:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

Zwróć uwagę na ostrzeżenie odnoszące się do efektu parametru `--update-snapshot`.

!!! note "Uwaga"

    Tutaj używamy danych testowych, których użyliśmy wcześniej do zademonstrowania naukowych wyników pipeline'u.
    Gdybyśmy planowali obsługiwać te testy w środowisku produkcyjnym, wygenerowalibyśmy mniejsze dane wejściowe do celów testowych.

    Ogólnie rzecz biorąc, ważne jest, aby testy jednostkowe były jak najlżejsze poprzez używanie najmniejszych elementów danych niezbędnych i wystarczających do oceny funkcjonalności procesu, w przeciwnym razie całkowity czas wykonania może znacznie się zwiększyć.
    Zestaw testów, którego uruchomienie zajmuje zbyt dużo czasu, to zestaw testów, który prawdopodobnie zostanie pominięty w interesie szybkości.

### Podsumowanie

Napisałeś Swój pierwszy test modułu dla procesu genomiki, weryfikując, że `SAMTOOLS_INDEX` poprawnie tworzy pliki indeksów dla różnych plików BAM. Zestaw testów zapewnia, że:

1. Proces działa pomyślnie
2. Pliki indeksów są tworzone
3. Wyjścia są spójne między uruchomieniami
4. Proces działa dla wszystkich plików BAM próbek

### Co dalej?

Dowiedz się, jak pisać testy dla innych procesów w naszym workflow'ie genomiki, używając metody setup do obsługi połączonych procesów. Ocenimy również, czy wyjścia, konkretnie nasze pliki VCF, zawierają oczekiwane wywołania wariantów.

---

## 2. Dodaj testy do połączonego procesu i testuj zawartość

Aby przetestować `GATK_HAPLOTYPECALLER`, musimy dostarczyć procesowi dane wyjściowe z `SAMTOOLS_INDEX` jako dane wejściowe. Moglibyśmy to zrobić, uruchamiając `SAMTOOLS_INDEX`, pobierając jego wyjścia i przechowując je z danymi testowymi dla workflow'u. To faktycznie jest zalecanym podejściem dla dopracowanego pipeline'u, ale nf-test zapewnia alternatywne podejście, używając metody `setup`.

Metoda setup pozwala wywołać proces `SAMTOOLS_INDEX` jako część konfiguracji testu, a jego wynik wykorzystać jako dane wejściowe dla `GATK_HAPLOTYPECALLER`. Ma to koszt: będziemy musieli uruchomić proces `SAMTOOLS_INDEX` przy każdym uruchomieniu testu dla `GATK_HAPLOTYPECALLER`. Jednak może nadal rozwijamy workflow i nie chcemy wstępnie generować danych testowych, które możemy później zmienić. Proces `SAMTOOLS_INDEX` jest również bardzo szybki, więc korzyści z wstępnego generowania i przechowywania jego wyjść mogą być pomijalne. Oto jak działa metoda setup.

### 2.1. Wygeneruj i umieść plik testowy

Jak poprzednio, najpierw generujemy szkielet pliku:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Wynik polecenia"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

To tworzy następujący szkielet testu:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. Przenieś plik testowy i zaktualizuj ścieżkę skryptu

Tworzymy katalog dla pliku testowego znajdującego się w tym samym miejscu co plik `main.nf` modułu:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

I przenosimy tam plik szkieletowy testu:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Na koniec, nie zapomnij zaktualizować ścieżki skryptu:

=== "Po zmianach"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Przed zmianami"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Dostarcz wejścia za pomocą metody setup

Wstawiamy blok `setup` przed blokiem `when`, gdzie możemy wywołać uruchomienie procesu `SAMTOOLS_INDEX` na jednym z naszych oryginalnych plików wejściowych. Pamiętaj również, aby jak poprzednio zmienić nazwę testu na coś znaczącego.

=== "Po zmianach"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Przed zmianami"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Następnie możemy odwołać się do wyjścia tego procesu w bloku `when`, gdzie określamy wejścia testowe:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Wprowadź tę zmianę i uruchom test ponownie:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

Tworzy to również plik migawki jak wcześniej.

### 2.4. Uruchom ponownie i zaobserwuj niepowodzenie

Co ciekawe, jeśli uruchomisz to samo polecenie ponownie, tym razem test nie powiedzie się.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

Komunikat o błędzie informuje, że wystąpiły różnice między migawkami dla dwóch uruchomień; konkretnie, wartości md5sum są różne dla plików VCF.

Dlaczego? Krótko mówiąc, narzędzie HaplotypeCaller zawiera znacznik czasu w nagłówku VCF, który jest różny za każdym razem (z definicji).
W rezultacie nie możemy po prostu oczekiwać, że pliki będą miały identyczne sumy md5, nawet jeśli mają identyczną zawartość w zakresie samych wywołań wariantów.

Jak sobie z tym poradzić?

### 2.5. Użyj metody asercji zawartości, aby sprawdzić konkretny wariant

Jednym ze sposobów rozwiązania problemu jest użycie [innego rodzaju asercji](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).
W tym przypadku sprawdzimy konkretną zawartość zamiast stwierdzać identyczność.
Dokładniej, każemy narzędziu odczytać linie pliku VCF i sprawdzić istnienie konkretnych linii.

W praktyce zastępujemy drugą asercję w bloku `then` w następujący sposób:

=== "Po zmianach"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Przed zmianami"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Tutaj odczytujemy pełną zawartość pliku wyjściowego VCF i szukamy dopasowania zawartości, co jest w porządku w przypadku małego pliku testowego, ale nie chciałbyś tego robić na większym pliku.
Możesz zamiast tego wybrać odczytanie konkretnych linii.

To podejście wymaga bardziej starannego wyboru tego, czego chcemy użyć jako „sygnału" do przetestowania.
Z jasnej strony, może być używane do testowania z dużą precyzją, czy narzędzie analityczne może konsekwentnie identyfikować „trudne" cechy (takie jak rzadkie warianty) w miarę dalszego rozwoju.

### 2.6. Uruchom ponownie i zaobserwuj sukces

Po zmodyfikowaniu testu w ten sposób możemy uruchomić test wiele razy i będzie konsekwentnie przechodził.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Dodaj więcej testów

Dodaj podobne testy dla próbek matki i ojca:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Uruchom polecenie testowe

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

To kończy podstawowy plan testów dla tego drugiego kroku w pipeline'ie. Przechodzimy do trzeciego i ostatniego testu na poziomie modułu!

### Podsumowanie

Nauczyłeś się, jak:

1. Testować procesy zależne od wyjść innych procesów
2. Weryfikować konkretne warianty genomiczne w plikach wyjściowych VCF
3. Obsługiwać niedeterministyczne wyjścia poprzez sprawdzanie konkretnej zawartości
4. Testować wykrywanie wariantów w wielu próbkach

### Co dalej?

Dowiedz się, jak pisać testy używające wstępnie wygenerowanych danych testowych dla etapu wspólnego genotypowania.

---

## 3. Używanie wstępnie wygenerowanych danych testowych

W przypadku etapu wspólnego genotypowania użyjemy innego podejścia - używania wstępnie wygenerowanych danych testowych. Jest to często preferowane dla:

1. Złożonych procesów z wieloma zależnościami
2. Procesów, których uruchomienie zajmuje dużo czasu
3. Procesów będących częścią stabilnego, produkcyjnego pipeline'u

### 3.1. Wygeneruj dane testowe

Sprawdź wyniki, które wygenerowaliśmy na początku tej sekcji:

```bash
tree results_genomics/
```

```console title="Zawartość katalogu wyników"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

Etap wspólnego genotypowania potrzebuje plików VCF wygenerowanych przez etapy haplotype caller jako wejścia, wraz z indeksami. Więc skopiujmy wyniki, które mamy, do katalogu testowego modułu `jointgenotyping`.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Teraz możemy użyć tych plików jako wejść do testu, który zamierzamy napisać dla etapu wspólnego genotypowania.

### 3.2. Wygeneruj szkielet pliku testowego

Jak poprzednio, najpierw generujemy szkielet pliku:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Wynik polecenia"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

To tworzy następujący szkielet testu:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. Przenieś plik testowy i zaktualizuj ścieżkę skryptu

Tym razem mamy już katalog dla testów znajdujący się w tym samym miejscu co plik `main.nf` modułu, więc możemy przenieść tam plik szkieletowy testu:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

I nie zapomnij zaktualizować ścieżki skryptu:

=== "Po zmianach"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Przed zmianami"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Dostarcz wejścia

Wypełnij wejścia na podstawie definicji wejść procesu i odpowiednio zmień nazwę testu:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Użyj asercji zawartości

Wyjściem etapu wspólnego genotypowania jest kolejny plik VCF, więc znowu użyjemy asercji zawartości.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

Sprawdzając zawartość konkretnego wariantu w pliku wyjściowym, ten test weryfikuje, że:

1. Proces wspólnego genotypowania działa pomyślnie
2. Wyjściowy VCF zawiera wszystkie trzy próbki we właściwej kolejności
3. Konkretny wariant jest wywoływany poprawnie z:
   - Dokładnymi genotypami dla każdej próbki (0/1 dla ojca, 1/1 dla matki i syna)
   - Poprawnymi głębokościami odczytów i jakościami genotypów
   - Statystykami na poziomie populacji, takimi jak częstość alleli (AF=0.833)

Nie zrobiliśmy migawki całego pliku, ale sprawdzając konkretny wariant, możemy być pewni, że proces wspólnego genotypowania działa zgodnie z oczekiwaniami.

### 3.6. Uruchom test

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

Test przechodzi, weryfikując, że nasz proces wspólnego genotypowania poprawnie:

1. Łączy indywidualne pliki VCF próbek
2. Wykonuje wspólne wykrywanie wariantów
3. Produkuje wielopróbkowy VCF ze spójnymi wywołaniami genotypów między uruchomieniami

### Podsumowanie

Wiesz, jak:

- Używać wcześniej wygenerowanych wyników jako wejść dla testów
- Pisać testy używające wstępnie wygenerowanych danych testowych

### Co dalej?

Dodaj test na poziomie workflow'u, aby zweryfikować, że cały pipeline wykrywania wariantów działa od początku do końca.

---

## 4. Dodaj test na poziomie workflow

Teraz przetestujemy kompletny pipeline wykrywania wariantów, od plików BAM do wspólnych genotypów. To weryfikuje, że:

1. Wszystkie procesy współpracują poprawnie
2. Dane przepływają prawidłowo między krokami
3. Ostateczne wywołania wariantów są spójne

### 4.1. Wygeneruj test workflow

Wygeneruj plik testowy dla kompletnego pipeline'u:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Wynik polecenia"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

To tworzy podstawowy szkielet testu:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Po prostu popraw nazwę na coś znaczącego (wkrótce zobaczysz, dlaczego to jest przydatne).

=== "Po zmianach"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Przed zmianami"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note "Uwaga"

    W tym przypadku plik testowy może pozostać tam, gdzie `nf-test` go utworzył.

### 4.2. Określ parametry wejściowe

Nadal musimy określić wejścia, co robi się nieco inaczej na poziomie workflow'u w porównaniu z testami na poziomie modułów.
Istnieje kilka sposobów robienia tego, w tym poprzez określenie profilu.
Jednak prostszym sposobem jest skonfigurowanie bloku `params {}` w pliku `nextflow.config`, który `nf-test init` pierwotnie utworzył w katalogu `tests`.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

// Katalog wyjściowy dla wyników workflow'u
outputDir = 'results_genomics'

/*
 * Parametry pipeline'u
 */

params {
    // Główne dane wejściowe (plik z listą plików wejściowych, jeden na linię)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Pliki pomocnicze
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Nazwa bazowa dla końcowego pliku wyjściowego
    cohort_name = "family_trio"
}
```

Gdy uruchomimy test, `nf-test` pobierze ten plik konfiguracyjny i odpowiednio pobierze wejścia.

### 4.3. Uruchom test workflow'u

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

Test przechodzi, potwierdzając, że nasz kompletny pipeline wykrywania wariantów:

1. Pomyślnie przetwarza wszystkie próbki
2. Poprawnie łączy wszystkie kroki

### 4.4. Uruchom WSZYSTKIE testy

nf-test ma jeszcze jedną sztuczkę w zanadrzu. Możemy uruchomić wszystkie testy naraz! Zmodyfikuj plik `nf-test.config` tak, aby nf-test szukał plików nf-test w każdym katalogu. Możesz to zrobić, modyfikując parametr `testsDir`:

=== "Po zmianach"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Przed zmianami"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Teraz możemy po prostu uruchomić nf-test i uruchomi _każdy pojedynczy test_ w naszym repozytorium:

```bash
nf-test test
```

??? success "Wynik polecenia"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

8 testów w 1 poleceniu! Poświęciliśmy dużo czasu na konfigurowanie wielu testów, ale jeśli chodzi o ich uruchamianie, było bardzo szybko i łatwo. Możesz zobaczyć, jak przydatne jest to podczas utrzymywania dużego pipeline'u, który może zawierać setki różnych elementów. Poświęcamy czas na napisanie testów raz, aby móc zaoszczędzić czas na uruchamianiu ich wiele razy.

Co więcej, możemy to zautomatyzować! Wyobraź sobie testy uruchamiane za każdym razem, gdy Ty lub kolega próbujecie dodać nowy kod. W ten sposób zapewniamy, że nasze pipeline utrzymują wysoki standard.

## Podsumowanie

Teraz wiesz, jak pisać i uruchamiać kilka rodzajów testów dla Twojego pipeline'u genomiki używając nf-test. Ten framework testowy pomaga zapewnić, że Twój workflow wykrywania wariantów produkuje spójne, niezawodne wyniki w różnych środowiskach i podczas wprowadzania zmian w kodzie.

Nauczyłeś się testować krytyczne komponenty takie jak:

- Proces `SAMTOOLS_INDEX`, który przygotowuje pliki BAM do wykrywania wariantów
- Proces `GATK_HAPLOTYPECALLER`, który identyfikuje warianty w poszczególnych próbkach
- Proces `GATK_JOINTGENOTYPING`, który łączy wywołania wariantów w całej kohorcie

Wdrożyłeś również różne strategie testowania specyficzne dla danych genomicznych:

- Weryfikowanie, że pliki VCF zawierają oczekiwane wywołania wariantów pomimo niedeterministycznych elementów, takich jak znaczniki czasu
- Testowanie zestawem danych tria rodzinnego, aby zapewnić prawidłową identyfikację wariantów w powiązanych próbkach
- Sprawdzanie konkretnych współrzędnych genomicznych i informacji o wariantach w plikach wyjściowych

Te umiejętności testowania są niezbędne do rozwijania solidnych pipeline'ów bioinformatycznych, które mogą niezawodnie przetwarzać dane genomiczne i produkować dokładne wywołania wariantów. W miarę dalszej pracy z Nextflow w analizie genomiki, te fundamenty testowania pomogą Ci utrzymać wysoką jakość kodu, który produkuje wiarygodne wyniki naukowe.
