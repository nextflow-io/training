# Część 1: Więcej o kontenerach

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Jak znaleźć lub stworzyć obrazy kontenerów

Niektórzy twórcy oprogramowania udostępniają obrazy kontenerów dla swojego oprogramowania w rejestrach kontenerów, takich jak Docker Hub, ale wielu tego nie robi.
W tej opcjonalnej sekcji pokażemy Ci dwa sposoby na uzyskanie obrazu kontenera dla narzędzi, których chcesz użyć w Swoich pipeline'ach Nextflow: za pomocą Seqera Containers oraz samodzielne budowanie obrazu kontenera.

Będziesz uzyskiwać/budować obraz kontenera dla pakietu pip `quote`, który zostanie użyty w ćwiczeniu na końcu tej sekcji.

### 1.1. Uzyskaj obraz kontenera z Seqera Containers

Seqera Containers to darmowa usługa, która buduje obrazy kontenerów dla narzędzi instalowanych przez pip i conda (w tym bioconda).
Przejdź do [Seqera Containers](https://www.seqera.io/containers/) i wyszukaj pakiet pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Kliknij "+Add", a następnie "Get Container", aby zażądać obrazu kontenera dla pakietu pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Jeśli po raz pierwszy budowany jest kontener społecznościowy dla tej wersji pakietu, może to potrwać kilka minut.
Kliknij, aby skopiować URI (np. `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) utworzonego dla Ciebie obrazu kontenera.

Możesz teraz użyć obrazu kontenera, aby uruchomić polecenie `quote` i uzyskać losowe powiedzenie Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Wyjście:

```console title="Wyjście"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Zbuduj obraz kontenera samodzielnie

Wykorzystajmy szczegóły budowy ze strony Seqera Containers, aby samodzielnie zbudować obraz kontenera dla pakietu pip `quote`.
Wróć na stronę Seqera Containers i kliknij przycisk "Build Details".

Pierwszym elementem, na który spojrzymy, jest `Dockerfile`, rodzaj pliku skryptu zawierającego wszystkie polecenia potrzebne do zbudowania obrazu kontenera.
Dodaliśmy wyjaśniające komentarze do poniższego Dockerfile, aby pomóc Ci zrozumieć, co robi każda część.

```Dockerfile title="Dockerfile"
# Zacznij od bazowego obrazu docker micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Skopiuj plik conda.yml do kontenera
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Zainstaluj różne narzędzia dla Nextflow oraz pakiety z pliku conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Uruchom kontener jako użytkownik root
USER root
# Ustaw zmienną środowiskową PATH tak, aby zawierała katalog instalacji micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

Drugim elementem, na który spojrzymy, jest plik `conda.yml`, który zawiera listę pakietów, które należy zainstalować w obrazie kontenera.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Skopiuj zawartość tych plików do zaślepek znajdujących się w katalogu `containers/build`, następnie uruchom poniższe polecenie, aby samodzielnie zbudować obraz kontenera.

!!! note "Uwaga"

    Używamy flagi `-t quote:latest`, aby oznaczyć obraz kontenera nazwą `quote` i tagiem `latest`.
    Będziemy mogli użyć tego tagu do odwoływania się do obrazu kontenera podczas uruchamiania go w tym systemie.

```bash
docker build -t quote:latest containers/build
```

Po zakończeniu budowania możesz uruchomić właśnie zbudowany obraz kontenera.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Podsumowanie

Nauczyłeś się dwóch różnych sposobów uzyskiwania obrazu kontenera dla narzędzia, którego chcesz użyć w Swoich pipeline'ach Nextflow: za pomocą Seqera Containers oraz samodzielnego budowania obrazu kontenera.

### Co dalej?

Masz wszystko, czego potrzebujesz, aby przejść do [następnego rozdziału](./04_hello_genomics.md) tej serii szkoleniowej.
Możesz również kontynuować opcjonalne ćwiczenie, aby pobierać cytaty pionierów informatyki/biologii za pomocą kontenera `quote` i wyświetlać je za pomocą kontenera `cowsay`.

---

## 2. Spraw, aby krowa cytowała słynnych naukowców

Ta sekcja zawiera dodatkowe ćwiczenia, aby przećwiczyć to, czego się dotychczas nauczyłeś.
Wykonanie tych ćwiczeń _nie jest wymagane_ do zrozumienia późniejszych części szkolenia, ale stanowi zabawny sposób na utrwalenie Swojej wiedzy poprzez wymyślenie, jak sprawić, aby krowa cytowała słynnych naukowców.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. Zmodyfikuj skrypt `hello-containers.nf`, aby używał procesu getQuote

Mamy listę pionierów informatyki i biologii w pliku `containers/data/pioneers.csv`.
Na wysokim poziomie, aby ukończyć to ćwiczenie, musisz:

- Zmodyfikować domyślny `params.input_file`, aby wskazywał na plik `pioneers.csv`.
- Utworzyć proces `getQuote`, który używa kontenera `quote` do pobierania cytatu dla każdego wejścia.
- Połączyć wyjście procesu `getQuote` z procesem `cowsay`, aby wyświetlić cytat.

Dla obrazu kontenera `quote` możesz użyć tego, który zbudowałeś samodzielnie w poprzednim dodatkowym ćwiczeniu lub tego, który uzyskałeś z Seqera Containers.

!!! tip "Wskazówka"

    Dobrym wyborem dla bloku `script` Twojego procesu getQuote może być:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Rozwiązanie tego ćwiczenia znajdziesz w pliku `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Zmodyfikuj Swój pipeline Nextflow, aby umożliwić jego wykonanie w trybach `quote` i `sayHello`.

Dodaj logikę rozgałęzienia do Swojego pipeline'u, aby umożliwić mu akceptowanie wejść przeznaczonych zarówno dla `quote`, jak i `sayHello`.
Oto przykład użycia instrukcji `if` w workflow Nextflow:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! tip "Wskazówka"

    Możesz użyć `new_ch = processName.out`, aby przypisać nazwę do kanału wyjściowego procesu.

Rozwiązanie tego ćwiczenia znajdziesz w pliku `containers/solutions/hello-containers-4.2.nf`.

### Podsumowanie

Wiesz już, jak używać kontenerów w Nextflow do uruchamiania procesów oraz jak zbudować logikę rozgałęzienia w Swoich pipeline'ach!

### Co dalej?

Świętuj, zrób sobie przerwę na rozciąganie i wypij trochę wody!

Kiedy będziesz gotowy, przejdź do Części 3 tej serii szkoleniowej, aby dowiedzieć się, jak zastosować to, czego się dotychczas nauczyłeś, do bardziej realistycznego przypadku analizy danych.
