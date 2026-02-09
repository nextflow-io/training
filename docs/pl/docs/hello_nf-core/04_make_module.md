# Część 4: Tworzenie modułu nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W czwartej części szkolenia Hello nf-core pokażemy Ci, jak stworzyć moduł nf-core, stosując kluczowe konwencje, które czynią moduły przenośnymi i łatwymi w utrzymaniu.

Projekt nf-core udostępnia polecenie (`nf-core modules create`), które automatycznie generuje odpowiednio ustrukturyzowane szablony modułów, podobnie jak to, czego użyliśmy dla workflow'a w Części 2.
Jednak w celach dydaktycznych zaczniemy od wykonania tego ręcznie: przekształcimy lokalny moduł `cowpy` w Twoim pipeline'ie `core-hello` w moduł w stylu nf-core krok po kroku.
Następnie pokażemy Ci, jak używać tworzenia modułów opartego na szablonach, aby pracować wydajniej w przyszłości.

??? info "Jak rozpocząć od tej sekcji"

    Ta sekcja zakłada, że ukończyłeś [Część 3: Użycie modułu nf-core](./03_use_module.md) i zintegrowałeś moduł `CAT_CAT` ze swoim pipeline'em.

    Jeśli nie ukończyłeś Części 3 lub chcesz zacząć od nowa w tej części, możesz użyć rozwiązania `core-hello-part3` jako punktu wyjścia.
    Uruchom te polecenia z wnętrza katalogu `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    To da Ci pipeline'a z już zintegrowanym modułem `CAT_CAT`.
    Możesz sprawdzić, czy działa poprawnie, uruchamiając następujące polecenie:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Przekształcenie `cowpy` w moduł nf-core

W tej sekcji zastosujemy konwencje nf-core do lokalnego modułu `cowpy` w Twoim pipeline'ie `core-hello`, przekształcając go w moduł zgodny ze standardami społeczności nf-core.

To jest obecny kod modułu procesu `cowpy`:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

1. **Zmiana nazwy procesu na wielkie litery `COWPY`**, aby przestrzegać konwencji.
2. **Aktualizacja `COWPY` do używania krotek metadanych**, aby propagować metadane próbek przez workflow'a.
3. **Centralizacja konfiguracji argumentów narzędzia za pomocą `ext.args`**, aby zwiększyć wszechstronność modułu przy zachowaniu minimalnego interfejsu.
4. **Standaryzacja nazewnictwa wyjść za pomocą `ext.prefix`**, aby promować spójność.
5. **Centralizacja konfiguracji publikowania**, aby promować spójność.

Po każdym kroku uruchomimy pipeline'a, aby sprawdzić, czy wszystko działa zgodnie z oczekiwaniami.

!!! warning "Katalog roboczy"

    Upewnij się, że jesteś w katalogu `core-hello` (głównym katalogu Twojego pipeline'a) dla wszystkich edycji plików i wykonywania poleceń w tej sekcji.

    ```bash
    cd core-hello
    ```

### 1.1. Zmiana nazwy procesu na wielkie litery

To jest czysto konwencja stylistyczna (nie ma technicznego uzasadnienia), ale ponieważ jest to norma dla modułów nf-core, zastosujmy się do niej.

Musimy wprowadzić trzy zestawy zmian:

1. Zaktualizować nazwę procesu w module
2. Zaktualizować instrukcję importu modułu w nagłówku workflow'a
3. Zaktualizować wywołanie procesu i deklarację emit w treści workflow'a

Zaczynajmy!

#### 1.1.1. Aktualizacja nazwy procesu w module

Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i zmodyfikuj nazwę procesu na wielkie litery:

=== "Po"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Przed"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

W tym przypadku zmiana na wielkie litery jest całkowicie prosta.

Gdyby nazwa procesu składała się z kilku słów, na przykład gdybyśmy mieli proces pierwotnie nazwany MyCowpyTool w notacji wielbłądziej, konwencja nf-core wymagałaby użycia podkreśleń do ich rozdzielenia, dając MY_COWPY_TOOL.

#### 1.1.2. Aktualizacja instrukcji importu modułu

Nazwy procesów rozróżniają wielkość liter, więc teraz, gdy zmieniliśmy nazwę procesu, musimy odpowiednio zaktualizować instrukcję importu modułu w nagłówku workflow'a w `hello.nf`:

=== "Po"

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

=== "Przed"

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

Moglibyśmy użyć aliasu w instrukcji importu, aby uniknąć konieczności aktualizacji wywołań procesu, ale to w pewnym stopniu zaprzeczyłoby celowi przyjęcia konwencji wielkich liter.

#### 1.1.3. Aktualizacja wywołania procesu i deklaracji emit

Teraz zaktualizujmy dwa odniesienia do procesu w bloku workflow'a w `hello.nf`:

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
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

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
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

Upewnij się, że wprowadzasz **obie** zmiany, w przeciwnym razie otrzymasz błąd podczas uruchamiania.

#### 1.1.4. Uruchomienie pipeline'a w celu przetestowania

Uruchommy workflow'a, aby sprawdzić, czy wszystko działa poprawnie po tych zmianach.

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

Świetnie, to działa! Teraz przejdźmy do wprowadzania bardziej istotnych zmian.

### 1.2. Aktualizacja `COWPY` do używania krotek metadanych

W obecnej wersji pipeline'a `core-hello` wyodrębniamy plik z krotki wyjściowej `CAT_CAT`, aby przekazać go do `COWPY`, jak pokazano w górnej części diagramu poniżej.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Lepiej byłoby, gdyby `COWPY` akceptował krotki metadanych bezpośrednio, umożliwiając przepływ metadanych przez workflow'a, jak pokazano w dolnej części diagramu.

W tym celu musimy wprowadzić następujące zmiany:

1. Zaktualizować definicje wejścia i wyjścia
2. Zaktualizować wywołanie procesu w workflow'ie
3. Zaktualizować blok emit w workflow'ie

Po wykonaniu wszystkich tych kroków uruchomimy pipeline'a, aby sprawdzić, czy wszystko nadal działa jak poprzednio.

#### 1.2.1. Aktualizacja definicji wejścia i wyjścia

Wróć do pliku modułu `cowpy.nf` i zmodyfikuj go, aby akceptował krotki metadanych, jak pokazano poniżej.

=== "Po"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Przed"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Jak widzisz, zmieniliśmy zarówno **główne wejście**, jak i **wyjście** na krotkę zgodną ze wzorcem `tuple val(meta), path(input_file)` wprowadzonym w Części 3 tego szkolenia.
Dla wyjścia wykorzystaliśmy również tę okazję, aby dodać `emit: cowpy_output` w celu nadania kanałowi wyjściowemu opisowej nazwy.

Teraz, gdy zmieniliśmy to, czego oczekuje proces, musimy zaktualizować to, co mu dostarczamy w wywołaniu procesu.

#### 1.2.2. Aktualizacja wywołania procesu w workflow'ie

Dobra wiadomość jest taka, że ta zmiana uprości wywołanie procesu.
Teraz, gdy wyjście `CAT_CAT` i wejście `COWPY` mają ten sam „kształt", tj. oba składają się ze struktury `tuple val(meta), path(input_file)`, możemy po prostu połączyć je bezpośrednio, zamiast musieć jawnie wyodrębniać plik z wyjścia procesu `CAT_CAT`.

Otwórz plik workflow'a `hello.nf` (w `core-hello/workflows/`) i zaktualizuj wywołanie `COWPY`, jak pokazano poniżej.

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Teraz wywołujemy `COWPY` bezpośrednio na `CAT_CAT.out.file_out`.

W rezultacie nie musimy już konstruować kanału `ch_for_cowpy`, więc ta linia (i jej linia komentarza) może zostać całkowicie usunięta.

#### 1.2.3. Aktualizacja bloku emit w workflow'ie

Ponieważ `COWPY` teraz emituje nazwane wyjście, `cowpy_output`, możemy zaktualizować blok `emit:` workflow'a `hello.nf`, aby z niego korzystał.

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Technicznie nie jest to wymagane, ale dobrą praktyką jest odwoływanie się do nazwanych wyjść, gdy tylko jest to możliwe.

#### 1.2.4. Uruchomienie pipeline'a w celu przetestowania

Uruchommy workflow'a, aby sprawdzić, czy wszystko działa poprawnie po tych zmianach.

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

Pipeline'a powinien uruchomić się pomyślnie, z metadanymi przepływającymi teraz z `CAT_CAT` przez `COWPY`.

To kończy to, co musieliśmy zrobić, aby `COWPY` obsługiwał krotki metadanych.
Teraz przyjrzyjmy się, co jeszcze możemy zrobić, aby wykorzystać wzorce modułów nf-core.

### 1.3. Centralizacja konfiguracji argumentów narzędzia za pomocą `ext.args`

W obecnym stanie proces `COWPY` oczekuje otrzymania wartości dla parametru `character`.
W rezultacie musimy podawać wartość za każdym razem, gdy wywołujemy proces, nawet jeśli bylibyśmy zadowoleni z domyślnych ustawień narzędzia.
Dla `COWPY` nie jest to wprawdzie duży problem, ale dla narzędzi z wieloma opcjonalnymi parametrami może to stać się dość uciążliwe.

Projekt nf-core zaleca używanie funkcji Nextflow'a zwanej [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) do wygodniejszego zarządzania argumentami narzędzi za pomocą plików konfiguracyjnych.

Zamiast deklarować wejścia procesu dla każdej opcji narzędzia, piszesz moduł tak, aby odwoływał się do `ext.args` w konstrukcji swojej linii poleceń.
Następnie wystarczy skonfigurować zmienną `ext.args` tak, aby zawierała argumenty i wartości, których chcesz użyć w pliku `modules.config`, który konsoliduje szczegóły konfiguracji dla wszystkich modułów.
Nextflow doda te argumenty z ich wartościami do linii poleceń narzędzia w czasie wykonania.

Zastosujmy to podejście do modułu `COWPY`.
Będziemy musieli wprowadzić następujące zmiany:

1. Zaktualizować moduł `COWPY`
2. Skonfigurować `ext.args` w pliku `modules.config`
3. Zaktualizować workflow'a `hello.nf`

Po wykonaniu wszystkich tych kroków uruchomimy pipeline'a, aby sprawdzić, czy wszystko nadal działa jak poprzednio.

#### 1.3.1. Aktualizacja modułu `COWPY`

Zróbmy to.
Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i zmodyfikuj go, aby odwoływał się do `ext.args`, jak pokazano poniżej.

=== "Po"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

=== "Przed"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

Jak widzisz, wprowadziliśmy trzy zmiany.

1. **W bloku `input:` usunęliśmy wejście `val character`.**
   Od teraz będziemy dostarczać ten argument za pomocą konfiguracji `ext.args`, jak opisano poniżej.

2. **W bloku `script:` dodaliśmy linię `def args = task.ext.args ?: ''`.**
   Ta linia używa operatora `?:` do określenia wartości zmiennej `args`: zawartość `task.ext.args`, jeśli nie jest pusta, lub pusty ciąg znaków, jeśli jest.
   Zauważ, że chociaż ogólnie odnosimy się do `ext.args`, ten kod musi odwoływać się do `task.ext.args`, aby wydobyć konfigurację `ext.args` na poziomie modułu.

3. **W linii poleceń zastąpiliśmy `-c "$character"` przez `$args`.**
   To tutaj Nextflow wstrzyknie wszelkie argumenty narzędzia ustawione w `ext.args` w pliku `modules.config`.

W rezultacie interfejs modułu jest teraz prostszy: oczekuje tylko podstawowych wejść metadanych i plików.

!!! note

    Operator `?:` jest często nazywany „operatorem Elvisa", ponieważ wygląda jak twarz Elvisa Presleya z boku, gdzie znak `?` symbolizuje falę w jego włosach.

#### 1.3.2. Konfiguracja `ext.args` w pliku `modules.config`

Teraz, gdy usunęliśmy deklarację `character` z modułu, musimy dodać ją do `ext.args` w pliku konfiguracyjnym `modules.config`.

Konkretnie, dodamy ten mały fragment kodu do bloku `process {}`:

```groovy title="Kod do dodania"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

Składnia `withName:` przypisuje tę konfigurację tylko do procesu `COWPY`, a `ext.args = { "-c ${params.character}" }` po prostu komponuje ciąg znaków, który będzie zawierał wartość parametru `character`.
Zauważ użycie nawiasów klamrowych, które mówią Nextflow'owi, aby ocenił wartość parametru w czasie wykonania.

Rozumiesz? Dodajmy to.

Otwórz `conf/modules.config` i dodaj kod konfiguracyjny wewnątrz bloku `process {}`, jak pokazano poniżej.

=== "Po"

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

=== "Przed"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Mamy nadzieję, że możesz sobie wyobrazić, że wszystkie moduły w pipeline'ie mają swoje `ext.args` określone w tym pliku, z następującymi korzyściami:

- **Interfejs modułu pozostaje prosty** - Akceptuje tylko podstawowe wejścia metadanych i plików
- **Pipeline'a nadal udostępnia `params.character`** - Użytkownicy końcowi nadal mogą go konfigurować jak poprzednio
- **Moduł jest teraz przenośny** - Może być ponownie użyty w innych pipeline'ach bez oczekiwania konkretnej nazwy parametru
- Konfiguracja jest **scentralizowana** w `modules.config`, utrzymując logikę workflow'a w czystości

Używając pliku `modules.config` jako miejsca, w którym wszystkie pipeline'y centralizują konfigurację dla poszczególnych modułów, sprawiamy, że nasze moduły są bardziej wielokrotnego użytku w różnych pipeline'ach.

#### 1.3.3. Aktualizacja workflow'a `hello.nf`

Ponieważ moduł `COWPY` nie wymaga już parametru `character` jako wejścia, musimy odpowiednio zaktualizować wywołanie workflow'a.

Otwórz plik workflow'a `hello.nf` (w `core-hello/workflows/`) i zaktualizuj wywołanie `COWPY`, jak pokazano poniżej.

=== "Po"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Przed"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

Kod workflow'a jest teraz czystszy: nie musimy przekazywać `params.character` bezpośrednio do procesu.
Interfejs modułu jest utrzymywany w minimalnej formie, czyniąc go bardziej przenośnym, podczas gdy pipeline'a nadal zapewnia jawną opcję poprzez konfigurację.

#### 1.3.4. Uruchomienie pipeline'a w celu przetestowania

Sprawdźmy, czy workflow'a nadal działa zgodnie z oczekiwaniami, określając inną postać, aby zweryfikować, że konfiguracja `ext.args` działa.

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

Powinno to uruchomić się pomyślnie jak poprzednio.

Sprawdźmy, czy konfiguracja `ext.args` zadziałała, sprawdzając wyjście.
Znajdź wyjście w przeglądarce plików lub użyj hasza zadania (część `38/eb29ea` w powyższym przykładzie), aby spojrzeć na plik wyjściowy:

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

Powinieneś zobaczyć grafikę ASCII wyświetloną z postacią `kosh`, potwierdzając, że konfiguracja `ext.args` zadziałała!

??? info "(Opcjonalnie) Sprawdzenie pliku polecenia"

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

Może się to wydawać niepotrzebne dla prostego narzędzia jak `cowpy`, ale może zrobić dużą różnicę dla narzędzi do analizy danych, które mają wiele opcjonalnych argumentów.

Podsumowując korzyści tego podejścia:

- **Czysty interfejs**: Moduł koncentruje się na podstawowych wejściach danych (metadanych i plikach)
- **Elastyczność**: Użytkownicy mogą określać argumenty narzędzia za pomocą konfiguracji, w tym wartości specyficzne dla próbki
- **Spójność**: Wszystkie moduły nf-core stosują ten wzorzec
- **Przenośność**: Moduły mogą być ponownie używane bez zakodowanych opcji narzędzi
- **Brak zmian w workflow'ie**: Dodawanie lub zmiana opcji narzędzi nie wymaga aktualizacji kodu workflow'a

!!! note

    System `ext.args` ma potężne dodatkowe możliwości, które nie zostały tutaj omówione, w tym dynamiczne przełączanie wartości argumentów na podstawie metadanych. Zobacz [specyfikacje modułów nf-core](https://nf-co.re/docs/guidelines/components/modules) po więcej szczegółów.

### 1.4. Standaryzacja nazewnictwa wyjść za pomocą `ext.prefix`

Teraz, gdy daliśmy procesowi `COWPY` dostęp do metamapy, możemy zacząć wykorzystywać kolejny przydatny wzorzec nf-core: nazywanie plików wyjściowych na podstawie metadanych.

Tutaj użyjemy funkcji Nextflow'a zwanej `ext.prefix`, która pozwoli nam standaryzować nazewnictwo plików wyjściowych w modułach używając `meta.id` (identyfikatora zawartego w metamapie), zachowując jednocześnie możliwość indywidualnej konfiguracji modułów w razie potrzeby.

Będzie to podobne do tego, co zrobiliśmy z `ext.args`, z kilkoma różnicami, które szczegółowo omówimy w trakcie.

Zastosujmy to podejście do modułu `COWPY`.
Będziemy musieli wprowadzić następujące zmiany:

1. Zaktualizować moduł `COWPY`
2. Skonfigurować `ext.prefix` w pliku `modules.config`

(Nie są potrzebne zmiany w workflow'ie.)

Po wykonaniu tego uruchomimy pipeline'a, aby sprawdzić, czy wszystko nadal działa jak poprzednio.

#### 1.4.1. Aktualizacja modułu `COWPY`

Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i zmodyfikuj go, aby odwoływał się do `ext.prefix`, jak pokazano poniżej.

=== "Po"

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

=== "Przed"

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

Jak widzisz, wprowadziliśmy trzy zmiany.

1. **W bloku `script:` dodaliśmy linię `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Ta linia używa operatora `?:` do określenia wartości zmiennej `prefix`: zawartość `task.ext.prefix`, jeśli nie jest pusta, lub identyfikator z metamapy (`meta.id`), jeśli jest.
   Zauważ, że chociaż ogólnie odnosimy się do `ext.prefix`, ten kod musi odwoływać się do `task.ext.prefix`, aby wydobyć konfigurację `ext.prefix` na poziomie modułu.

2. **W linii poleceń zastąpiliśmy `cowpy-${input_file}` przez `${prefix}.txt`.**
   To tutaj Nextflow wstrzyknie wartość `prefix` określoną przez linię powyżej.

3. **W bloku `output:` zastąpiliśmy `path("cowpy-${input_file}")` przez `path("${prefix}.txt")`.**
   To po prostu powtarza, jaka będzie ścieżka pliku zgodnie z tym, co jest napisane w linii poleceń.

W rezultacie nazwa pliku wyjściowego jest teraz konstruowana przy użyciu rozsądnego domyślnego (identyfikatora z metamapy) połączonego z odpowiednim rozszerzeniem formatu pliku.

#### 1.4.2. Konfiguracja `ext.prefix` w pliku `modules.config`

W tym przypadku rozsądne domyślne nie jest wystarczająco wyraziste dla naszego gustu; chcemy użyć niestandardowego wzorca nazewnictwa, który zawiera nazwę narzędzia, `cowpy-<id>.txt`, jak mieliśmy wcześniej.

Zrobimy to, konfigurując `ext.prefix` w `modules.config`, tak jak zrobiliśmy dla parametru `character` z `ext.args`, z tym że tym razem blok `withName: 'COWPY' {}` już istnieje i musimy tylko dodać następującą linię:

```groovy title="Kod do dodania"
ext.prefix = { "cowpy-${meta.id}" }
```

To skomponuje ciąg znaków, którego chcemy.
Zauważ, że ponownie używamy nawiasów klamrowych, tym razem aby powiedzieć Nextflow'owi, aby ocenił wartość `meta.id` w czasie wykonania.

Dodajmy to.

Otwórz `conf/modules.config` i dodaj kod konfiguracyjny wewnątrz bloku `process {}`, jak pokazano poniżej.

=== "Po"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Przed"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

W przypadku, gdybyś się zastanawiał, domknięcie `ext.prefix` ma dostęp do odpowiedniego fragmentu metadanych, ponieważ konfiguracja jest oceniana w kontekście wykonania procesu, gdzie metadane są dostępne.

#### 1.4.3. Uruchomienie pipeline'a w celu przetestowania

Sprawdźmy, czy workflow'a nadal działa zgodnie z oczekiwaniami.

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
Powinieneś zobaczyć plik wyjściowy cowpy z takim samym nazewnictwem jak poprzednio: `cowpy-test.txt`, oparty na domyślnej nazwie partii.

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

Możesz zmienić konfigurację `ext.prefix` w `conf/modules.config`, aby przekonać się, że możesz zmienić wzorzec nazewnictwa bez konieczności wprowadzania jakichkolwiek zmian w kodzie modułu lub workflow'a.

Alternatywnie możesz również spróbować uruchomić to ponownie z innym parametrem `--batch` określonym w linii poleceń, aby przekonać się, że ta część jest nadal konfigurowalna w locie.

To pokazuje, jak `ext.prefix` pozwala zachować preferowaną konwencję nazewnictwa, zachowując jednocześnie elastyczny interfejs modułu.

Podsumowując korzyści tego podejścia:

- **Standaryzowane nazewnictwo**: Pliki wyjściowe są zazwyczaj nazywane przy użyciu identyfikatorów próbek z metadanych
- **Konfigurowalne**: Użytkownicy mogą nadpisać domyślne nazewnictwo w razie potrzeby
- **Spójne**: Wszystkie moduły nf-core stosują ten wzorzec
- **Przewidywalne**: Łatwo wiedzieć, jak będą nazywane pliki wyjściowe

Całkiem nieźle, prawda?
Cóż, jest jeszcze jedna ważna zmiana, którą musimy wprowadzić, aby ulepszyć nasz moduł, aby pasował do wytycznych nf-core.

### 1.5. Centralizacja konfiguracji publikowania

Mogłeś zauważyć, że publikowaliśmy wyjścia do dwóch różnych katalogów:

- **`results`** — Oryginalny katalog wyjściowy, którego używaliśmy od początku dla naszych lokalnych modułów, ustawiony indywidualnie za pomocą dyrektyw `publishDir` dla poszczególnych modułów;
- **`core-hello-results`** — Katalog wyjściowy ustawiony za pomocą `--outdir` w linii poleceń, który otrzymywał logi nf-core i wyniki publikowane przez `CAT_CAT`.

To jest niechlujne i nieoptymalne; lepiej byłoby mieć jedną lokalizację dla wszystkiego.
Oczywiście moglibyśmy wejść do każdego z naszych lokalnych modułów i ręcznie zaktualizować dyrektywę `publishDir`, aby używała katalogu `core-hello-results`, ale co z następnym razem, gdy zdecydujemy się zmienić katalog wyjściowy?

Posiadanie indywidualnych modułów podejmujących decyzje o publikowaniu wyraźnie nie jest właściwą drogą, szczególnie w świecie, w którym ten sam moduł może być używany w wielu różnych pipeline'ach, przez ludzi, którzy mają różne potrzeby lub preferencje.
Chcemy mieć możliwość kontrolowania, gdzie wyjścia są publikowane na poziomie konfiguracji workflow'a.

„Hej," możesz powiedzieć, „`CAT_CAT` wysyła swoje wyjścia do `--outdir`. Może powinniśmy skopiować jego dyrektywę `publishDir`?"

Tak, to świetny pomysł.

Z wyjątkiem tego, że nie ma dyrektywy `publishDir`. (Śmiało, spójrz na kod modułu.)

To dlatego, że pipeline'y nf-core centralizują kontrolę na poziomie workflow'a, konfigurując `publishDir` w `conf/modules.config` zamiast w poszczególnych modułach.
Konkretnie, szablon nf-core deklaruje domyślną dyrektywę `publishDir` (ze wstępnie zdefiniowaną strukturą katalogów), która ma zastosowanie do wszystkich modułów, chyba że zostanie dostarczona nadpisująca dyrektywa.

Czy to nie brzmi niesamowicie? Czy może być tak, że aby skorzystać z tej domyślnej dyrektywy, wszystko, co musimy zrobić, to usunąć obecną dyrektywę `publishDir` z naszych lokalnych modułów?

Wypróbujmy to na `COWPY`, aby zobaczyć, co się stanie, a następnie przyjrzymy się kodowi domyślnej konfiguracji, aby zrozumieć, jak to działa.

Na koniec pokażemy, jak nadpisać domyślne zachowanie, jeśli jest to pożądane.

#### 1.5.1. Usunięcie dyrektywy `publishDir` z `COWPY`

Zróbmy to.
Otwórz plik modułu `cowpy.nf` (w `core-hello/modules/local/`) i usuń dyrektywę `publishDir`, jak pokazano poniżej.

=== "Po"

    ```groovy title="core-hello/modules/local/cowpy.nf (fragment)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Przed"

    ```groovy title="core-hello/modules/local/cowpy.nf (fragment)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

To wszystko!

#### 1.5.2. Uruchomienie pipeline'a w celu przetestowania

Zobaczmy, co się stanie, jeśli teraz uruchomimy pipeline'a.

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

Spójrz na swój obecny katalog roboczy.
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

Widzisz, że Nextflow stworzył tę hierarchię katalogów na podstawie nazw workflow'a i modułu.

Kod odpowiedzialny za to znajduje się w pliku `conf/modules.config`.
To jest domyślna konfiguracja `publishDir`, która jest częścią szablonu nf-core i ma zastosowanie do wszystkich procesów:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Może to wyglądać skomplikowanie, więc przyjrzyjmy się każdemu z trzech komponentów:

- **`path:`** Określa katalog wyjściowy na podstawie nazwy procesu.
  Pełna nazwa procesu zawarta w `task.process` zawiera hierarchię importów workflow'a i modułu (takich jak `CORE_HELLO:HELLO:CAT_CAT`).
  Operacje `tokenize` usuwają tę hierarchię, aby uzyskać tylko nazwę procesu, następnie biorą pierwszą część przed jakimkolwiek podkreśleniem (jeśli dotyczy) i konwertują ją na małe litery.
  To właśnie określa, że wyniki `CAT_CAT` są publikowane do `${params.outdir}/cat/`.
- **`mode:`** Kontroluje sposób publikowania plików (kopiowanie, dowiązanie symboliczne itp.).
  Jest to konfigurowalne za pomocą parametru `params.publish_dir_mode`.
- **`saveAs:`** Filtruje, które pliki publikować.
  Ten przykład wyklucza pliki `versions.yml`, zwracając dla nich `null`, zapobiegając ich publikowaniu.

To zapewnia spójną logikę organizowania wyjść.

Wyjście wygląda jeszcze lepiej, gdy wszystkie moduły w pipeline'ie przyjmą tę konwencję, więc śmiało usuń dyrektywy `publishDir` z innych modułów w swoim pipeline'ie.
To domyślne ustawienie zostanie zastosowane nawet do modułów, których nie modyfikowaliśmy jawnie, aby przestrzegały wytycznych nf-core.

To powiedziawszy, możesz zdecydować, że chcesz zorganizować swoje wejścia inaczej, a dobra wiadomość jest taka, że łatwo to zrobić.

#### 1.5.3. Nadpisanie domyślnego

Aby nadpisać domyślną dyrektywę `publishDir`, możesz po prostu dodać własne dyrektywy do pliku `conf/modules.config`.

Na przykład możesz nadpisać domyślne dla pojedynczego procesu używając selektora `withName:`, jak w tym przykładzie, gdzie dodajemy niestandardową dyrektywę `publishDir` dla procesu 'COWPY'.

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

Nie będziemy faktycznie wprowadzać tej zmiany, ale możesz pobawić się tym i zobaczyć, jaką logikę możesz zaimplementować.

Chodzi o to, że ten system daje Ci to, co najlepsze z obu światów: spójność domyślnie i elastyczność w dostosowywaniu konfiguracji na żądanie.

Podsumowując, otrzymujesz:

- **Jedno źródło prawdy**: Cała konfiguracja publikowania znajduje się w `modules.config`
- **Przydatne domyślne**: Procesy działają od razu bez konfiguracji dla poszczególnych modułów
- **Łatwe dostosowywanie**: Nadpisz zachowanie publikowania w konfiguracji, a nie w kodzie modułu
- **Przenośne moduły**: Moduły nie kodują na stałe lokalizacji wyjściowych

To kończy zestaw funkcji modułów nf-core, których absolutnie powinieneś się nauczyć używać, ale są inne, o których możesz przeczytać w [specyfikacjach modułów nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Podsumowanie

Teraz wiesz, jak dostosować lokalne moduły, aby przestrzegały konwencji nf-core:

- Projektuj swoje moduły tak, aby akceptowały i propagowały krotki metadanych;
- Używaj `ext.args`, aby utrzymać interfejsy modułów minimalne i przenośne;
- Używaj `ext.prefix` dla konfigurowalnego, standaryzowanego nazewnictwa plików wyjściowych;
- Przyjmij domyślną scentralizowaną dyrektywę `publishDir` dla spójnej struktury katalogu wyników.

### Co dalej?

Dowiedz się, jak używać wbudowanych narzędzi nf-core opartych na szablonach do łatwego tworzenia modułów.

---

## 2. Tworzenie modułu za pomocą narzędzi nf-core

Teraz, gdy nauczyłeś się wzorców modułów nf-core, stosując je ręcznie, przyjrzyjmy się, jak tworzyłbyś moduły w praktyce.

### 2.1. Generowanie szkieletu modułu z szablonu

Podobnie jak w przypadku tworzenia pipeline'ów, projekt nf-core udostępnia narzędzia do generowania odpowiednio ustrukturyzowanych modułów na podstawie szablonu, ze wszystkimi tymi wzorcami wbudowanymi od początku.

#### 2.1.1. Uruchomienie polecenia tworzenia modułu

Polecenie `nf-core modules create` generuje szablon modułu, który już przestrzega wszystkich konwencji, których się nauczyłeś.

Stwórzmy nową wersję modułu `COWPY` z minimalnym szablonem, uruchamiając to polecenie:

```bash
nf-core modules create --empty-template COWPY
```

Flaga `--empty-template` tworzy czysty szablon startowy bez dodatkowego kodu, ułatwiając zobaczenie podstawowej struktury.

Polecenie działa interaktywnie, prowadząc Cię przez konfigurację.
Automatycznie wyszukuje informacje o narzędziu z repozytoriów pakietów, takich jak Bioconda i bio.tools, aby wstępnie wypełnić metadane.

Zostaniesz poproszony o kilka opcji konfiguracyjnych:

- **Informacje o autorze**: Twoja nazwa użytkownika GitHub do przypisania
- **Etykieta zasobów**: Wstępnie zdefiniowany zestaw wymagań obliczeniowych.
  Projekt nf-core udostępnia standardowe etykiety, takie jak `process_single` dla lekkich narzędzi i `process_high` dla wymagających.
  Te etykiety pomagają zarządzać alokacją zasobów w różnych środowiskach wykonawczych.
- **Wymaganie metadanych**: Czy moduł potrzebuje informacji specyficznych dla próbki za pomocą mapy `meta` (zazwyczaj tak dla modułów przetwarzania danych).

Narzędzie obsługuje złożoność znajdowania informacji o pakiecie i konfigurowania struktury, pozwalając Ci skupić się na implementacji konkretnej logiki narzędzia.

#### 2.1.2. Sprawdzenie szkieletu modułu

Narzędzie tworzy kompletną strukturę modułu w `modules/local/` (lub `modules/nf-core/`, jeśli jesteś w repozytorium nf-core/modules):

??? abstract "Zawartość katalogu"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Każdy plik służy konkretnemu celowi:

- **`main.nf`**: Definicja procesu ze wszystkimi wzorcami nf-core wbudowanymi
- **`meta.yml`**: Dokumentacja modułu opisująca wejścia, wyjścia i narzędzie
- **`environment.yml`**: Specyfikacja środowiska Conda dla zależności
- **`tests/main.nf.test`**: Przypadki testowe nf-test do walidacji działania modułu

!!! tip "Dowiedz się więcej o testowaniu"

    Wygenerowany plik testowy używa nf-test, frameworka testowego dla pipeline'ów i modułów Nextflow'a. Aby dowiedzieć się, jak pisać i uruchamiać te testy, zobacz [zadanie poboczne nf-test](../side_quests/nf-test.md).

Wygenerowany `main.nf` zawiera wszystkie wzorce, których się właśnie nauczyłeś, plus kilka dodatkowych funkcji:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Wzorzec 1: Krotki metadanych ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Wzorzec 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Wzorzec 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Zauważ, jak wszystkie wzorce, które zastosowałeś ręcznie powyżej, są już tam!

Szablon zawiera również kilka dodatkowych konwencji nf-core.
Niektóre z nich działają od razu, podczas gdy inne są symbolami zastępczymi, które będziemy musieli wypełnić, jak opisano poniżej.

**Funkcje, które działają od razu:**

- **`tag "$meta.id"`**: Dodaje identyfikator próbki do nazw procesów w logach dla łatwiejszego śledzenia
- **`label 'process_single'`**: Etykieta zasobów do konfigurowania wymagań CPU/pamięci
- **Blok `when:`**: Umożliwia warunkowe wykonanie za pomocą konfiguracji `task.ext.when`

Te funkcje są już funkcjonalne i czynią moduły bardziej łatwymi w utrzymaniu.

**Symbole zastępcze, które dostosujemy poniżej:**

- **Bloki `input:` i `output:`**: Ogólne deklaracje, które zaktualizujemy, aby pasowały do naszego narzędzia
- **Blok `script:`**: Zawiera komentarz, gdzie dodamy polecenie `cowpy`
- **Blok `stub:`**: Szablon, który zaktualizujemy, aby produkował poprawne wyjścia
- **Kontener i środowisko**: Symbole zastępcze, które wypełnimy informacjami o pakiecie

Następne sekcje przeprowadzą Cię przez ukończenie tych dostosowań.

### 2.2. Konfiguracja kontenera i środowiska conda

Wytyczne nf-core wymagają, abyśmy określili zarówno kontener, jak i środowisko Conda jako część modułu.

#### 2.2.1. Kontener

Dla kontenera możesz użyć [Seqera Containers](https://seqera.io/containers/), aby automatycznie zbudować kontener z dowolnego pakietu Conda, w tym pakietów conda-forge.
W tym przypadku używamy tego samego wstępnie zbudowanego kontenera co poprzednio.

Domyślny kod oferuje przełączanie między Dockerem a Singularity, ale uprościmy tę linię i po prostu określimy kontener Dockera, który otrzymaliśmy z Seqera Containers powyżej.

=== "Po"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Przed"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Środowisko Conda

Dla środowiska Conda kod modułu określa `conda "${moduleDir}/environment.yml"`, co oznacza, że powinno być skonfigurowane w pliku `environment.yml`.

Narzędzie do tworzenia modułów ostrzegło nas, że nie może znaleźć pakietu `cowpy` w Bioconda (głównym kanale dla narzędzi bioinformatycznych).
Jednak `cowpy` jest dostępny w conda-forge, więc możesz uzupełnić `environment.yml` w ten sposób:

=== "Po"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Przed"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

Do przesłania do nf-core musielibyśmy bardziej ściśle przestrzegać domyślnych ustawień, ale dla naszego własnego użytku możemy uprościć kod w ten sposób.

!!! tip "Pakiety Bioconda vs conda-forge"

    - **Pakiety Bioconda**: Automatycznie otrzymują zbudowane BioContainers, zapewniając gotowe do użycia kontenery
    - **Pakiety conda-forge**: Mogą używać Seqera Containers do budowania kontenerów na żądanie z przepisu Conda

    Większość narzędzi bioinformatycznych znajduje się w Bioconda, ale dla narzędzi conda-forge, Seqera Containers zapewnia łatwe rozwiązanie do konteneryzacji.

### 2.3. Podłączenie logiki `COWPY`

Teraz zaktualizujmy elementy kodu, które są specyficzne dla tego, co robi proces `COWPY`: wejścia i wyjścia oraz blok skryptu.

#### 2.3.1. Wejścia i wyjścia

Wygenerowany szablon zawiera ogólne deklaracje wejścia i wyjścia, które będziesz musiał dostosować dla swojego konkretnego narzędzia.
Patrząc wstecz na nasz ręczny moduł `COWPY` z sekcji 1, możemy użyć go jako przewodnika.

Zaktualizuj bloki wejścia i wyjścia:

=== "Po"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Przed"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

To określa:

- Nazwę parametru pliku wejściowego (`input_file` zamiast ogólnego `input`)
- Nazwę pliku wyjściowego używając konfigurowalnego wzorca prefiksu (`${prefix}.txt` zamiast wieloznacznika `*`)
- Opisową nazwę emit (`cowpy_output` zamiast ogólnego `output`)

Jeśli używasz serwera języka Nextflow'a do walidacji składni, część `${prefix}` zostanie oznaczona jako błąd na tym etapie, ponieważ nie dodaliśmy jej jeszcze do bloku skryptu.
Przejdźmy do tego teraz.

#### 2.3.2. Blok skryptu

Szablon zapewnia symbol zastępczy komentarza w bloku skryptu, gdzie powinieneś dodać faktyczne polecenie narzędzia.

Na podstawie modułu, który napisaliśmy ręcznie wcześniej, powinniśmy wprowadzić następujące edycje:

=== "Po"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Przed"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Kluczowe zmiany:

- Zmień `def prefix` na po prostu `prefix` (bez `def`), aby uczynić go dostępnym w bloku wyjścia
- Zastąp komentarz faktycznym poleceniem `cowpy`, które używa zarówno `$args`, jak i `${prefix}.txt`

Zauważ, że gdybyśmy nie wykonali już pracy dodania konfiguracji `ext.args` i `ext.prefix` dla procesu `COWPY` do pliku `modules.config`, musielibyśmy to zrobić teraz.

#### 2.3.3. Implementacja bloku stub

W kontekście Nextflow'a blok [stub](https://www.nextflow.io/docs/latest/process.html#stub) pozwala zdefiniować lekki, fikcyjny skrypt używany do szybkiego prototypowania i testowania logiki pipeline'a bez wykonywania faktycznego polecenia.

<!-- TODO (przyszłość) To jest bardzo pobieżnie omówione, ale naprawdę powinno być wyjaśnione lub przynajmniej linkować do wyjaśnienia o stub'ach (dokument referencyjny też nie jest szczególnie pomocny). W tej chwili jest to prawdopodobnie w większości bez znaczenia dla kogoś, kto nie zna już stub'ów. -->

Nie martw się zbytnio, jeśli to wydaje się tajemnicze; dołączamy to dla kompletności, ale możesz również po prostu usunąć sekcję stub, jeśli nie chcesz się tym zajmować, ponieważ jest całkowicie opcjonalna.

=== "Po"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Przed"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Kluczowe zmiany:

- Zmień `def prefix` na po prostu `prefix`, aby dopasować do bloku skryptu
- Usuń linię `echo $args` (która była tylko kodem zastępczym szablonu)
- Stub tworzy pusty plik `${prefix}.txt` pasujący do tego, co produkuje blok skryptu

To pozwala testować logikę workflow'a i obsługę plików bez czekania na faktyczne uruchomienie narzędzia.

Po ukończeniu konfiguracji środowiska (sekcja 2.2), wejść/wyjść (sekcja 2.3.1), bloku skryptu (sekcja 2.3.2) i bloku stub (sekcja 2.3.3), moduł jest gotowy do testowania!

### 2.4. Zamiana nowego modułu `COWPY` i uruchomienie pipeline'a

Wszystko, co musimy zrobić, aby wypróbować tę nową wersję modułu `COWPY`, to przełączyć instrukcję importu w pliku workflow'a `hello.nf`, aby wskazywała na nowy plik.

=== "Po"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
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

=== "Przed"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Uruchommy pipeline'a, aby go przetestować.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Wyjście polecenia"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
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
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

To produkuje te same wyniki co poprzednio.

### Podsumowanie

Teraz wiesz, jak używać wbudowanych narzędzi nf-core do wydajnego tworzenia modułów przy użyciu szablonów zamiast pisania wszystkiego od zera.

### Co dalej?

Dowiedz się, jakie są korzyści z wnoszenia modułów do nf-core i jakie są główne kroki i wymagania.

---

## 3. Wnoszenie modułów z powrotem do nf-core

Repozytorium [nf-core/modules](https://github.com/nf-core/modules) przyjmuje wkłady dobrze przetestowanych, standaryzowanych modułów.

### 3.1. Dlaczego warto wnosić?

Wnoszenie swoich modułów do nf-core:

- Udostępnia Twoje narzędzia całej społeczności nf-core poprzez katalog modułów na [nf-co.re/modules](https://nf-co.re/modules)
- Zapewnia ciągłe utrzymanie i ulepszenia przez społeczność
- Zapewnia zapewnienie jakości poprzez przegląd kodu i automatyczne testowanie
- Daje Twojej pracy widoczność i uznanie

### 3.2. Lista kontrolna współtwórcy

Aby wnieść moduł do nf-core, będziesz musiał przejść przez następujące kroki:

1. Sprawdź, czy już istnieje na [nf-co.re/modules](https://nf-co.re/modules)
2. Zrób fork repozytorium [nf-core/modules](https://github.com/nf-core/modules)
3. Użyj `nf-core modules create`, aby wygenerować szablon
4. Wypełnij logikę modułu i testy
5. Przetestuj za pomocą `nf-core modules test tool/subtool`
6. Sprawdź za pomocą `nf-core modules lint tool/subtool`
7. Prześlij pull request

Szczegółowe instrukcje znajdziesz w [samouczku komponentów nf-core](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Zasoby

- **Samouczek komponentów**: [Kompletny przewodnik po tworzeniu i wnoszeniu modułów](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Specyfikacje modułów**: [Wymagania techniczne i wytyczne](https://nf-co.re/docs/guidelines/components/modules)
- **Wsparcie społeczności**: [Slack nf-core](https://nf-co.re/join) - Dołącz do kanału `#modules`

### Podsumowanie

Teraz wiesz, jak tworzyć moduły nf-core! Nauczyłeś się czterech kluczowych wzorców, które czynią moduły przenośnymi i łatwymi w utrzymaniu:

- **Krotki metadanych** propagują metadane przez workflow'a
- **`ext.args`** upraszcza interfejsy modułów, obsługując opcjonalne argumenty za pomocą konfiguracji
- **`ext.prefix`** standaryzuje nazewnictwo plików wyjściowych
- **Scentralizowane publikowanie** za pomocą `publishDir` skonfigurowanego w `modules.config` zamiast zakodowanego na stałe w modułach

Przekształcając `COWPY` krok po kroku, rozwinąłeś głębokie zrozumienie tych wzorców, czyniąc Cię przygotowanym do pracy z, debugowania i tworzenia modułów nf-core.
W praktyce będziesz używać `nf-core modules create` do generowania odpowiednio ustrukturyzowanych modułów z tymi wzorcami wbudowanymi od początku.

Na koniec nauczyłeś się, jak wnosić moduły do społeczności nf-core, udostępniając narzędzia badaczom na całym świecie, jednocześnie korzystając z ciągłego utrzymania przez społeczność.

### Co dalej?

Gdy będziesz gotowy, przejdź do [Części 5: Walidacja wejścia](./05_input_validation.md), aby dowiedzieć się, jak dodać walidację wejścia opartą na schemacie do swojego pipeline'a.
