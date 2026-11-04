# Część 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Spraw, aby krowa cytowała słynnych naukowców

Ta sekcja zawiera kilka dodatkowych ćwiczeń, które pozwolą Ci przećwiczyć to, czego się do tej pory nauczyłeś.
Wykonanie tych ćwiczeń _nie jest wymagane_ do zrozumienia późniejszych części szkolenia, ale stanowią one zabawny sposób na utrwalenie wiedzy poprzez wymyślenie, jak sprawić, aby krowa cytowała słynnych naukowców.

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

### 1.1. Zmodyfikuj skrypt `hello-containers.nf`, aby używał procesu getQuote

Mamy listę pionierów informatyki i biologii w pliku `containers/data/pioneers.csv`.
Ogólnie rzecz biorąc, aby ukończyć to ćwiczenie, będziesz musiał:

- Zmodyfikować domyślny `params.input_file`, aby wskazywał na plik `pioneers.csv`.
- Stworzyć proces `getQuote`, który używa kontenera `quote` do pobrania cytatu dla każdego wejścia.
- Połączyć wyjście procesu `getQuote` z procesem `cowsay`, aby wyświetlić cytat.

Dla obrazu kontenera `quote` możesz użyć tego, który sam zbudowałeś w poprzednim dodatkowym ćwiczeniu, lub skorzystać z tego, który otrzymałeś z Seqera Containers.

!!! Hint "Podpowiedź"

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

### 1.2. Zmodyfikuj swój pipeline Nextflow, aby mógł działać w trybach `quote` i `sayHello`.

Dodaj logikę rozgałęzień do swojego pipeline'u, aby mógł przyjmować dane wejściowe przeznaczone zarówno dla `quote`, jak i `sayHello`.
Oto przykład użycia instrukcji `if` w workflow'ie Nextflow:

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

!!! Hint "Podpowiedź"

    Możesz użyć `new_ch = processName.out`, aby przypisać nazwę do kanału wyjściowego procesu.

Rozwiązanie tego ćwiczenia znajdziesz w pliku `containers/solutions/hello-containers-4.2.nf`.

### Podsumowanie

Wiesz już, jak używać kontenerów w Nextflow do uruchamiania procesów oraz jak wbudować logikę rozgałęzień do swoich pipeline'ów!

### Co dalej?

Świętuj, zrób sobie przerwę na rozciąganie i napij się wody!

Kiedy będziesz gotowy, przejdź do Części 3 tej serii szkoleniowej, aby nauczyć się, jak zastosować to, czego się do tej pory nauczyłeś, w bardziej realistycznym przypadku analizy danych.
