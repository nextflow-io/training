# Część 4: Tworzenie modułu nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W czwartej części kursu szkoleniowego Hello nf-core pokażemy, jak utworzyć moduł nf-core, stosując kluczowe konwencje, które czynią takie komponenty przenośnymi i łatwymi w utrzymaniu.

Projekt nf-core udostępnia polecenie (`nf-core modules create`), które automatycznie generuje prawidłowo ustrukturyzowane szablony, podobnie jak to, czego użyliśmy dla workflow'u w Części 2.
Jednak w celach edukacyjnych zaczniemy od wykonania tego ręcznie: przekształcenia lokalnego komponentu `cowpy` w Twoim pipeline'ie `core-hello` w moduł w stylu nf-core krok po kroku.
Następnie pokażemy, jak korzystać z tworzenia modułów opartych na szablonach, aby pracować wydajniej w przyszłości.

??? info "Jak rozpocząć od tej sekcji"

    Ta sekcja zakłada, że ukończyłeś [Część 3: Użycie modułu nf-core](./03_use_module.md) i zintegrowałeś moduł `CAT_CAT` ze Swoim pipeline'em.

    Jeśli nie ukończyłeś Części 3 lub chcesz zacząć od nowa w tej części, możesz użyć rozwiązania `core-hello-part3` jako punktu wyjścia.
    Uruchom te polecenia z wnętrza katalogu `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    To daje Ci pipeline z już zintegrowanym modułem `CAT_CAT`.
    Możesz przetestować, czy działa poprawnie, uruchamiając następujące polecenie:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Przekształcenie `cowpy` w moduł nf-core

W tej sekcji zastosujemy konwencje nf-core do lokalnego modułu `cowpy` w Twoim pipeline'ie `core-hello`, przekształcając go w moduł zgodny ze standardami społeczności nf-core.

To jest aktualny kod modułu procesu `cowpy`:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Zastosujemy następujące konwencje nf-core stopniowo:

1. **Zmiana nazwy procesu na wielkie litery `COWPY`** zgodnie z konwencją.
2. **Aktualizacja `COWPY` do używania krotek metadanych** w celu propagowania metadanych próbek przez workflow.
3. **Centralizacja konfiguracji argumentów narzędzia z `ext.args`** w celu zwiększenia wszechstronności modułu przy zachowaniu minimalnego interfejsu.
4. **Standaryzacja nazewnictwa wyjścia z `ext.prefix`** w celu promowania spójności.
5. **Centralizacja konfiguracji publikowania** w celu promowania spójności.

Po każdym kroku uruchomimy pipeline, aby sprawdzić, czy wszystko działa zgodnie z oczekiwaniami.

!!! warning "Katalog roboczy"

    Upewnij się, że jesteś w katalogu `core-hello` (głównym katalogu Twojego pipeline) dla wszystkich edycji plików i wykonywania poleceń w tej sekcji.

    ```bash
    cd core-hello
    ```

### 1.1. Zmiana nazwy procesu na wielkie litery

To jest czysto stylistyczna konwencja (nie ma technicznego uzasadnienia), ale ponieważ jest to norma dla modułów nf-core, dostosujmy się.

Musimy wprowadzić trzy zestawy zmian:

1. Zmienić nazwę procesu w module
2. Poprawić instrukcję importu modułu w nagłówku workflow
3. Dostosować wywołanie procesu i deklarację emit w treści workflow

Zaczynajmy!

#### 1.1.1. Aktualizacja nazwy procesu w module

Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i zmodyfikuj nazwę procesu na wielkie litery:

=== "Po zmianach"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

W tym przypadku zmiana na wielkie litery jest całkowicie prosta.

Gdyby nazwa procesu składała się z kilku słów, na przykład gdybyśmy mieli proces o nazwie MyCowpyTool pierwotnie w camel case, konwencją nf-core byłoby użycie podkreślników do ich rozdzielenia, dając MY_COWPY_TOOL.

#### 1.1.2. Aktualizacja instrukcji importu modułu

Nazwy procesów są wrażliwe na wielkość liter, więc teraz, gdy zmieniliśmy nazwę procesu, musimy odpowiednio zaktualizować instrukcję importu modułu w nagłówku workflow'u pliku `hello.nf`:

=== "Po zmianach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Moglibyśmy użyć aliasu w instrukcji importu, aby uniknąć konieczności aktualizacji wywołań procesu, ale to nieco zniweczyłoby sens przyjęcia konwencji pisania wielkimi literami.

#### 1.1.3. Aktualizacja wywołania procesu i deklaracji emit

Teraz zaktualizujmy dwa odwołania do procesu w bloku workflow'u pliku `hello.nf`:

=== "Po zmianach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // wygeneruj grafikę ASCII powitań za pomocą cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Zbierz i zapisz wersje oprogramowania
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // wygeneruj grafikę ASCII powitań za pomocą cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Zbierz i zapisz wersje oprogramowania
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Upewnij się, że wprowadzisz **obie** zmiany, w przeciwnym razie otrzymasz błąd podczas uruchamiania.

#### 1.1.4. Uruchomienie pipeline'a w celu przetestowania

Uruchommy workflow, aby sprawdzić, czy wszystko działa poprawnie po tych zmianach.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dobrze, to działa! Teraz przejdźmy do wprowadzania bardziej znaczących zmian.

### 1.2. Aktualizacja `COWPY` do używania krotek metadanych

W obecnej wersji pipeline `core-hello` wyodrębniamy plik z krotki wyjściowej `CAT_CAT`, aby przekazać go do `COWPY`, jak pokazano w górnej części diagramu poniżej.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Lepiej byłoby, gdyby `COWPY` akceptował krotki metadanych bezpośrednio, umożliwiając przepływ metadanych przez workflow, jak pokazano w dolnej części diagramu.

W tym celu będziemy musieli wprowadzić następujące zmiany:

1. Zmienić definicje wejścia i wyjścia
2. Poprawić wywołanie procesu w workflow
3. Dostosować blok emit w workflow

Po wykonaniu wszystkich tych czynności uruchomimy pipeline, aby sprawdzić, czy wszystko nadal działa jak wcześniej.

#### 1.2.1. Aktualizacja definicji wejścia i wyjścia

Wróć do pliku modułu `cowpy.nf` i zmodyfikuj go, aby akceptował krotki metadanych, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Jak widzisz, zmieniliśmy zarówno **główne wejście**, jak i **wyjście** na krotkę zgodną ze wzorcem `tuple val(meta), path(input_file)` wprowadzonym w Części 3 tego szkolenia.
Dla wyjścia wykorzystaliśmy również tę okazję, aby dodać `emit: cowpy_output`, nadając kanałowi wyjściowemu opisową nazwę.

Teraz, gdy zmieniliśmy to, czego oczekuje proces, musimy zaktualizować to, co mu dostarczamy w wywołaniu procesu.

#### 1.2.2. Aktualizacja wywołania procesu w workflow

Dobra wiadomość jest taka, że ta zmiana uprości wywołanie procesu.
Teraz, gdy wyjście `CAT_CAT` i wejście `COWPY` mają tę samą 'kształt', tzn. oba składają się ze struktury `tuple val(meta), path(input_file)`, możemy po prostu połączyć je bezpośrednio zamiast konieczności jawnego wyodrębniania pliku z wyjścia procesu `CAT_CAT`.

Otwórz plik workflow `hello.nf` (w `core-hello/workflows/`) i zaktualizuj wywołanie `COWPY`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // wygeneruj grafikę ASCII powitań za pomocą cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // wygeneruj grafikę ASCII powitań za pomocą cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Teraz wywołujemy `COWPY` bezpośrednio na `CAT_CAT.out.file_out`.

W rezultacie nie musimy już konstruować kanału `ch_for_cowpy`, więc ta linia (i jej linia komentarza) może zostać całkowicie usunięta.

#### 1.2.3. Aktualizacja bloku emit w workflow

Ponieważ `COWPY` teraz emituje nazwane wyjście `cowpy_output`, możemy zaktualizować blok `emit:` workflow `hello.nf`, aby z niego korzystał.

=== "Po zmianach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Technicznie nie jest to wymagane, ale dobrą praktyką jest odwoływanie się do nazwanych wyjść, gdy tylko jest to możliwe.

#### 1.2.4. Uruchomienie pipeline'a w celu przetestowania

Uruchommy workflow, aby sprawdzić, czy wszystko działa poprawnie po tych zmianach.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Pipeline powinien działać pomyślnie, a metadane teraz przepływają z `CAT_CAT` przez `COWPY`.

To kończy to, co musieliśmy zrobić, aby `COWPY` obsługiwał krotki metadanych.
Teraz spójrzmy, co jeszcze możemy zrobić, aby wykorzystać wzorce modułów nf-core.

### 1.3. Centralizacja konfiguracji argumentów narzędzia z `ext.args`

W obecnym stanie proces `COWPY` oczekuje otrzymania wartości dla parametru `character`.
W rezultacie musimy podać wartość za każdym razem, gdy wywołujemy proces, nawet jeśli bylibyśmy zadowoleni z wartości domyślnych ustawionych przez narzędzie.
Dla `COWPY` nie jest to wprawdzie duży problem, ale dla programów z wieloma opcjonalnymi parametrami może to być dość uciążliwe.

Projekt nf-core zaleca używanie funkcji Nextflow zwanej [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) do wygodniejszego zarządzania argumentami poprzez pliki konfiguracyjne.

Zamiast deklarować wejścia procesu dla każdej opcji narzędzia, piszesz moduł tak, aby odwoływał się do `ext.args` w konstrukcji Swojej linii poleceń.
Następnie wystarczy skonfigurować zmienną `ext.args` tak, aby zawierała argumenty i wartości, których chcesz użyć w pliku `modules.config`, który konsoliduje szczegóły konfiguracji dla wszystkich modułów.
Nextflow doda te argumenty z ich wartościami do linii poleceń narzędzia w czasie wykonywania.

Zastosujmy to podejście do modułu `COWPY`.
Będziemy musieli wprowadzić następujące zmiany:

1. Zmodyfikować moduł `COWPY`
2. Skonfigurować `ext.args` w pliku `modules.config`
3. Dostosować workflow `hello.nf`

Po wykonaniu wszystkich tych czynności uruchomimy pipeline, aby sprawdzić, czy wszystko nadal działa jak wcześniej.

#### 1.3.1. Aktualizacja modułu `COWPY`

Zróbmy to.
Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i zmodyfikuj go, aby odwoływał się do `ext.args`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Jak widać, wprowadziliśmy trzy zmiany.

1. **W bloku `input:` usunęliśmy wejście `val character`.**
   Od tej pory będziemy dostarczać ten argument poprzez konfigurację `ext.args`, jak opisano poniżej.

2. **W bloku `script:` dodaliśmy linię `def args = task.ext.args ?: ''`.**
   Ta linia używa operatora `?:` do określenia wartości zmiennej `args`: zawartości `task.ext.args`, jeśli nie jest pusta, lub pustego ciągu, jeśli jest.
   Zauważ, że chociaż ogólnie odwołujemy się do `ext.args`, ten kod musi odwoływać się do `task.ext.args`, aby wyciągnąć konfigurację `ext.args` na poziomie modułu.

3. **W linii poleceń zastąpiliśmy `-c "$character"` przez `$args`.**
   To tutaj Nextflow wstawi wszystkie argumenty narzędzia ustawione w `ext.args` w pliku `modules.config`.

W rezultacie interfejs modułu jest teraz prostszy: oczekuje tylko podstawowych metadanych i wejść plikowych.

!!! note

    Operator `?:` jest często nazywany 'operatorem Elvisa', ponieważ wygląda jak twarz Elvisa Presleya z profilu, gdzie znak `?` symbolizuje falę we włosach.

#### 1.3.2. Konfiguracja `ext.args` w pliku `modules.config`

Teraz, gdy usunęliśmy deklarację `character` z modułu, musimy dodać ją do `ext.args` w pliku konfiguracyjnym `modules.config`.

Konkretnie, dodamy ten mały fragment kodu do bloku `process {}`:

```groovy title="Kod do dodania"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

Składnia `withName:` przypisuje tę konfigurację tylko do procesu `COWPY`, a `ext.args = { "-c ${params.character}" }` po prostu tworzy ciąg, który będzie zawierał wartość parametru `character`.
Zauważ użycie nawiasów klamrowych, które mówią Nextflow, aby ocenił wartość parametru w czasie wykonywania.

Ma sens? Dodajmy to.

Otwórz `conf/modules.config` i dodaj kod konfiguracyjny wewnątrz bloku `process {}`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Miejmy nadzieję, że możesz sobie wyobrazić, że wszystkie komponenty w pipeline'ie mają swoje `ext.args` określone w tym pliku, z następującymi korzyściami:

- **Interfejs pozostaje prosty** - Akceptuje tylko podstawowe metadane i wejścia plikowe
- **Pipeline nadal udostępnia `params.character`** - Użytkownicy końcowi nadal mogą go konfigurować jak wcześniej
- **Komponent jest teraz przenośny** - Może być ponownie użyty w innych pipeline bez oczekiwania konkretnej nazwy parametru
- Konfiguracja jest **scentralizowana** w `modules.config`, utrzymując logikę workflow'u czystą

Używając pliku `modules.config` jako miejsca, w którym wszystkie pipeline centralizują konfigurację dla poszczególnych komponentów, sprawiamy, że nasze moduły są bardziej wielokrotnego użytku w różnych środowiskach.

#### 1.3.3. Aktualizacja workflow `hello.nf`

Ponieważ moduł `COWPY` nie wymaga już parametru `character` jako wejścia, musimy odpowiednio zaktualizować wywołanie workflow.

Otwórz plik workflow `hello.nf` (w `core-hello/workflows/`) i zaktualizuj wywołanie `COWPY`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // wygeneruj grafikę ASCII powitań za pomocą cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // wygeneruj grafikę ASCII powitań za pomocą cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

Kod workflow jest teraz czystszy: nie musimy przekazywać `params.character` bezpośrednio do procesu.
Interfejs modułu jest minimalny, co czyni go bardziej przenośnym, podczas gdy pipeline nadal udostępnia jawną opcję poprzez konfigurację.

#### 1.3.4. Uruchomienie pipeline'a w celu przetestowania

Sprawdźmy, czy workflow nadal działa zgodnie z oczekiwaniami, określając inną postać, aby sprawdzić, czy konfiguracja `ext.args` działa.

Uruchom to polecenie używając `kosh`, jednej z bardziej... enigmatycznych opcji:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Powinno to działać pomyślnie jak wcześniej.

Sprawdźmy, czy konfiguracja `ext.args` zadziałała, sprawdzając wyjście.
Znajdź wyjście w przeglądarce plików lub użyj hasha zadania (część `38/eb29ea` w przykładzie powyżej), aby spojrzeć na plik wyjściowy:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Wyjście polecenia"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

Powinieneś zobaczyć sztukę ASCII wyświetloną z postacią `kosh`, potwierdzając, że konfiguracja `ext.args` zadziałała!

??? info "(Opcjonalnie) Sprawdź plik polecenia"

    Jeśli chcesz zobaczyć dokładnie, jak została zastosowana konfiguracja, możesz sprawdzić plik `.command.sh`:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Zobaczysz polecenie `cowpy` z argumentem `-c kosh`:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    To pokazuje, że plik `.command.sh` został wygenerowany poprawnie na podstawie konfiguracji `ext.args`.

Poświęć chwilę, aby pomyśleć o tym, co tutaj osiągnęliśmy.
To podejście utrzymuje interfejs modułu skoncentrowany na podstawowych danych (plikach, metadanych i wszelkich obowiązkowych parametrach dla próbki), podczas gdy opcje kontrolujące zachowanie narzędzia są obsługiwane oddzielnie poprzez konfigurację.

Może się to wydawać niepotrzebne dla prostego narzędzia jak `cowpy`, ale może zrobić dużą różnicę dla narzędzi analizy danych, które mają wiele opcjonalnych argumentów.

Podsumowując korzyści tego podejścia:

- **Czysty interfejs**: Moduł koncentruje się na podstawowych wejściach danych (metadanych i plikach)
- **Elastyczność**: Użytkownicy mogą określać argumenty narzędzia poprzez konfigurację, w tym wartości specyficzne dla próbki
- **Spójność**: Wszystkie moduły nf-core stosują ten wzorzec
- **Przenośność**: Moduły mogą być ponownie użyte bez zakodowanych opcji narzędzia
- **Brak zmian w workflow'ie**: Dodawanie lub zmiana opcji narzędzia nie wymaga aktualizacji kodu workflow'u

!!! note

    System `ext.args` ma potężne dodatkowe możliwości, których nie omówiono tutaj, w tym dynamiczne przełączanie wartości argumentów na podstawie metadanych. Zobacz [specyfikacje modułów nf-core](https://nf-co.re/docs/guidelines/components/modules) po więcej szczegółów.

### 1.4. Standaryzacja nazewnictwa wyjścia z `ext.prefix`

Teraz, gdy daliśmy procesowi `COWPY` dostęp do metamap, możemy zacząć wykorzystywać kolejny przydatny wzorzec nf-core: nazywanie plików wyjściowych na podstawie metadanych.

Tutaj użyjemy funkcji Nextflow zwanej `ext.prefix`, która pozwoli nam standaryzować nazewnictwo plików wyjściowych w modułach używając `meta.id` (identyfikatora zawartego w metamap), przy jednoczesnym zachowaniu możliwości indywidualnej konfiguracji modułów w razie potrzeby.

Będzie to podobne do tego, co zrobiliśmy z `ext.args`, z kilkoma różnicami, które szczegółowo opiszemy.

Zastosujmy to podejście do modułu `COWPY`.
Będziemy musieli wprowadzić następujące zmiany:

1. Zaktualizować moduł `COWPY`
2. Skonfigurować `ext.prefix` w pliku `modules.config`

(Nie są potrzebne zmiany w workflow'ie.)

Po wykonaniu tego uruchomimy pipeline, aby sprawdzić, czy wszystko nadal działa jak wcześniej.

#### 1.4.1. Aktualizacja modułu `COWPY`

Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i zmodyfikuj go, aby odwoływał się do `ext.prefix`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Jak widać, wprowadziliśmy trzy zmiany.

1. **W bloku `script:` dodaliśmy linię `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Ta linia używa operatora `?:` do określenia wartości zmiennej `prefix`: zawartości `task.ext.prefix`, jeśli nie jest pusta, lub identyfikatora z metamap (`meta.id`), jeśli jest.
   Zauważ, że chociaż ogólnie odwołujemy się do `ext.prefix`, ten kod musi odwoływać się do `task.ext.prefix`, aby wyciągnąć konfigurację `ext.prefix` na poziomie modułu.

2. **W linii poleceń zastąpiliśmy `cowpy-${input_file}` przez `${prefix}.txt`.**
   To tutaj Nextflow wstawi wartość `prefix` określoną przez powyższą linię.

3. **W bloku `output:` zastąpiliśmy `path("cowpy-${input_file}")` przez `path("${prefix}.txt")`.**
   To po prostu powtarza, jaka będzie ścieżka pliku zgodnie z tym, co jest napisane w linii poleceń.

W rezultacie nazwa pliku wyjściowego jest teraz konstruowana przy użyciu rozsądnej wartości domyślnej (identyfikatora z metamap) połączonej z odpowiednim rozszerzeniem formatu pliku.

#### 1.4.2. Konfiguracja `ext.prefix` w pliku `modules.config`

W tym przypadku rozsądna wartość domyślna nie jest wystarczająco ekspresyjna dla naszych potrzeb; chcemy użyć niestandardowego wzorca nazewnictwa, który zawiera nazwę narzędzia, `cowpy-<id>.txt`, tak jak mieliśmy wcześniej.

Zrobimy to, konfigurując `ext.prefix` w `modules.config`, tak jak zrobiliśmy dla parametru `character` z `ext.args`, z tym że tym razem blok `withName: 'COWPY' {}` już istnieje i musimy tylko dodać następującą linię:

```groovy title="Kod do dodania"
ext.prefix = { "cowpy-${meta.id}" }
```

To skomponuje ciąg, którego chcemy.
Zauważ, że ponownie używamy nawiasów klamrowych, tym razem aby powiedzieć Nextflow, aby ocenił wartość `meta.id` w czasie wykonywania.

Dodajmy to.

Otwórz `conf/modules.config` i dodaj kod konfiguracyjny wewnątrz bloku `process {}`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Jeśli się zastanawiasz, zamknięcie `ext.prefix` ma dostęp do właściwego fragmentu metadanych, ponieważ konfiguracja jest oceniana w kontekście wykonywania procesu, gdzie metadane są dostępne.

#### 1.4.3. Uruchomienie pipeline'a w celu przetestowania

Sprawdźmy, czy workflow nadal działa zgodnie z oczekiwaniami.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Spójrz na wyjście w katalogu wyników.
Powinieneś zobaczyć plik wyjściowy cowpy z takim samym nazewnictwem jak wcześniej: `cowpy-test.txt`, na podstawie domyślnej nazwy partii.

??? abstract "Zawartość katalogu"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Możesz zmienić konfigurację `ext.prefix` w `conf/modules.config`, aby upewnić się, że możesz zmienić wzorzec nazewnictwa bez konieczności wprowadzania jakichkolwiek zmian w module lub kodzie workflow.

Alternatywnie możesz również spróbować uruchomić to ponownie z innym parametrem `--batch` określonym w linii poleceń, aby upewnić się, że ta część jest nadal konfigurowalna na bieżąco.

To pokazuje, jak `ext.prefix` pozwala zachować preferowaną konwencję nazewnictwa przy jednoczesnym zachowaniu elastyczności interfejsu modułu.

Podsumowując korzyści tego podejścia:

- **Standaryzowane nazewnictwo**: Pliki wyjściowe są zazwyczaj nazywane przy użyciu identyfikatorów próbek z metadanych
- **Konfigurowalne**: Użytkownicy mogą zastąpić domyślne nazewnictwo w razie potrzeby
- **Spójne**: Wszystkie moduły nf-core stosują ten wzorzec
- **Przewidywalne**: Łatwo wiedzieć, jak będą nazywane pliki wyjściowe

Całkiem dobrze, prawda?
Cóż, jest jeszcze jedna ważna zmiana, którą musimy wprowadzić, aby ulepszyć nasz moduł zgodnie z wytycznymi nf-core.

### 1.5. Centralizacja konfiguracji publikowania

Mogłeś zauważyć, że publikowaliśmy wyjścia do dwóch różnych lokalizacji:

- **`results`** — Oryginalny katalog wyjściowy, którego używaliśmy od początku dla naszych lokalnych modułów, ustawiony indywidualnie za pomocą dyrektyw `publishDir` dla każdego z nich;
- **`core-hello-results`** — Ścieżka wyjściowa ustawiona za pomocą `--outdir` w linii poleceń, która otrzymywała logi nf-core i wyniki publikowane przez `CAT_CAT`.

To jest nieuporządkowane i nieoptymalne; lepiej byłoby mieć jedną lokalizację dla wszystkiego.
Oczywiście moglibyśmy wejść do każdego z naszych lokalnych modułów i ręcznie zaktualizować dyrektywę `publishDir`, aby używała katalogu `core-hello-results`, ale co z następnym razem, gdy zdecydujemy się zmienić katalog wyjściowy?

Posiadanie indywidualnych modułów podejmujących decyzje o publikowaniu wyraźnie nie jest właściwym podejściem, szczególnie w świecie, w którym ten sam moduł może być używany w wielu różnych pipeline przez ludzi, którzy mają różne potrzeby lub preferencje.
Chcemy mieć możliwość kontrolowania, gdzie publikowane są wyjścia na poziomie konfiguracji workflow.

"Hej," możesz powiedzieć, "`CAT_CAT` wysyła swoje wyjścia do `--outdir`. Może powinniśmy skopiować jego dyrektywę `publishDir`?"

Tak, to świetny pomysł.

Tyle że nie ma dyrektywy `publishDir`. (Śmiało, spójrz na kod modułu.)

To dlatego, że pipeline'y nf-core centralizują kontrolę na poziomie workflow'u, konfigurując `publishDir` w `conf/modules.config` zamiast w poszczególnych modułach.
Konkretnie, szablon nf-core deklaruje domyślną dyrektywę `publishDir` (ze wstępnie zdefiniowaną strukturą katalogów), która ma zastosowanie do wszystkich modułów, chyba że zostanie dostarczona nadpisująca dyrektywa.

Czy to nie brzmi świetnie? Czy może być tak, że aby skorzystać z tej domyślnej dyrektywy, wszystko, co musimy zrobić, to usunąć obecną dyrektywę `publishDir` z naszych lokalnych modułów?

Wypróbujmy to na `COWPY`, aby zobaczyć, co się stanie, a następnie spojrzymy na kod domyślnej konfiguracji, aby zrozumieć, jak to działa.

Na koniec pokażemy, jak w razie potrzeby przesłonić domyślne zachowanie.

#### 1.5.1. Usunięcie dyrektywy `publishDir` z `COWPY`

Zróbmy to.
Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i usuń dyrektywę `publishDir`, jak pokazano poniżej.

=== "Po zmianach"

    ```groovy title="core-hello/modules/local/cowpy.nf (fragment)" linenums="1"
    #!/usr/bin/env nextflow

    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Przed zmianami"

    ```groovy title="core-hello/modules/local/cowpy.nf (fragment)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

To wszystko!

#### 1.5.2. Uruchomienie pipeline'a w celu przetestowania

Zobaczmy, co się stanie, jeśli teraz uruchomimy pipeline.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Spójrz na Swój obecny katalog roboczy.
Teraz `core-hello-results` zawiera również wyjścia modułu `COWPY`.

??? abstract "Zawartość katalogu"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Widać, że Nextflow utworzył tę hierarchię katalogów na podstawie nazw workflow'u i modułu.

Kod odpowiedzialny za to znajduje się w pliku `conf/modules.config`.
Jest to domyślna konfiguracja `publishDir`, która jest częścią szablonu nf-core i ma zastosowanie do wszystkich procesów:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

To może wyglądać na skomplikowane, więc spójrzmy na każdy z trzech komponentów:

- **`path:`** Określa katalog wyjściowy na podstawie nazwy procesu.
  Pełna nazwa procesu zawarta w `task.process` obejmuje hierarchię importów workflow'u i modułów (takich jak `CORE_HELLO:HELLO:CAT_CAT`).
  Operacje `tokenize` usuwają tę hierarchię, aby uzyskać tylko nazwę procesu, następnie biorą pierwszą część przed dowolnym podkreśleniem (jeśli dotyczy) i konwertują ją na małe litery.
  To właśnie określa, że wyniki `CAT_CAT` są publikowane w `${params.outdir}/cat/`.
- **`mode:`** Kontroluje sposób publikowania plików (kopiowanie, dowiązanie symboliczne itp.).
  Jest to konfigurowalne za pomocą parametru `params.publish_dir_mode`.
- **`saveAs:`** Filtruje, które pliki publikować.
  Ten przykład wyklucza pliki `versions.yml`, zwracając dla nich `null`, zapobiegając ich publikowaniu.

To zapewnia spójną logikę organizowania wyjść.

Wyjście wygląda jeszcze lepiej, gdy wszystkie moduły w pipeline'ie przyjmą tę konwencję, więc śmiało usuń dyrektywy `publishDir` z innych modułów w Swoim pipeline'ie.
To ustawienie domyślne zostanie zastosowane nawet do modułów, których nie modyfikowaliśmy jawnie, aby były zgodne z wytycznymi nf-core.

To powiedziawszy, możesz zdecydować, że chcesz zorganizować Swoje wejścia inaczej, a dobrą wiadomością jest to, że łatwo to zrobić.

#### 1.5.3. Przesłonięcie wartości domyślnej

Aby przesłonić domyślną dyrektywę `publishDir`, możesz po prostu dodać własne dyrektywy do pliku `conf/modules.config`.

Na przykład możesz przesłonić wartość domyślną dla pojedynczego procesu za pomocą selektora `withName:`, jak w tym przykładzie, gdzie dodajemy niestandardową dyrektywę `publishDir` dla procesu 'COWPY'.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

Właściwie nie zamierzamy wprowadzać tej zmiany, ale śmiało pobaw się tym i zobacz, jaką logikę możesz zaimplementować.

Chodzi o to, że ten system daje najlepsze z obu światów: spójność domyślnie i elastyczność w dostosowywaniu konfiguracji na żądanie.

Podsumowując, otrzymujesz:

- **Pojedyncze źródło prawdy**: Cała konfiguracja publikowania znajduje się w `modules.config`
- **Przydatne ustawienie domyślne**: Procesy działają gotowe do użycia bez konfiguracji dla każdego modułu
- **Łatwe dostosowanie**: Przesłoń zachowanie publikowania w konfiguracji, a nie w kodzie modułu
- **Przenośne moduły**: Moduły nie kodują na stałe lokalizacji wyjściowych

To kończy zestaw funkcji modułów nf-core, których bezwzględnie powinieneś się nauczyć używać, ale są inne, o których możesz przeczytać w [specyfikacjach modułów nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Wnioski

Teraz wiesz, jak dostosować lokalne moduły, aby były zgodne z konwencjami nf-core:

- Projektuj Swoje moduły tak, aby akceptowały i propagowały krotki metadanych;
- Używaj `ext.args`, aby interfejsy modułów były minimalne i przenośne;
- Używaj `ext.prefix` dla konfigurowalnego, standaryzowanego nazewnictwa plików wyjściowych;
- Przyjmij domyślną scentralizowaną dyrektywę `publishDir` dla spójnej struktury katalogu wyników.

### Co dalej?

Dowiedz się, jak używać wbudowanych narzędzi nf-core
