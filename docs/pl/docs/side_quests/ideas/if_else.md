# Część 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Spraw, aby krowa cytowała słynnych naukowców

Ta sekcja zawiera dodatkowe ćwiczenia pozwalające na utrwalenie dotychczas zdobytej wiedzy.
Wykonanie tych ćwiczeń _nie jest wymagane_ do zrozumienia późniejszych części szkolenia, ale stanowi świetny sposób na sprawdzenie, jak zmusić krowę do cytowania słynnych naukowców.

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
Na wysokim poziomie, aby wykonać to ćwiczenie, będziesz musiał:

- Zmodyfikować domyślny `params.input_file`, aby wskazywał na plik `pioneers.csv`.
- Utworzyć proces `getQuote`, który używa kontenera `quote` do pobrania cytatu dla każdego wejścia.
- Połączyć wyjście procesu `getQuote` z procesem `cowsay`, aby wyświetlić cytat.

Dla obrazu kontenera `quote` możesz użyć albo tego, który sam zbudowałeś w poprzednim dodatkowym ćwiczeniu, albo pobrać gotowy z Seqera Containers.

!!! Hint

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

### 1.2. Zmodyfikuj Swój pipeline Nextflow, aby mógł wykonywać się w trybach `quote` i `sayHello`.

Dodaj logikę rozgałęzień do Swojego pipeline'u, aby mógł akceptować dane wejściowe przeznaczone zarówno dla `quote`, jak i `sayHello`.
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

!!! Hint

    Możesz użyć `new_ch = processName.out`, aby przypisać nazwę do kanału wyjściowego procesu.

Rozwiązanie tego ćwiczenia znajdziesz w pliku `containers/solutions/hello-containers-4.2.nf`.

### Podsumowanie

Wiesz już, jak używać kontenerów w Nextflow do uruchamiania procesów oraz jak budować logikę rozgałęzień w Swoich pipeline'ach!

### Co dalej?

Świętuj, zrób sobie przerwę na rozciąganie i napij się wody!

Kiedy będziesz gotowy, przejdź do Części 3 tej serii szkoleniowej, aby nauczyć się, jak zastosować dotychczas zdobytą wiedzę do bardziej realistycznego przypadku analizy danych.
