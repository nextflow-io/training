# Część 3: Profilowanie i optymalizacja zasobów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

TO JEST PLACEHOLDER

!!!note "Uwaga"

    Ten moduł szkoleniowy jest w trakcie przebudowy.

---

TODO

### 1.1. Uruchom workflow w celu wygenerowania raportu wykorzystania zasobów

Aby Nextflow automatycznie wygenerował raport, wystarczy dodać `-with-report <nazwa_pliku>.html` do polecenia.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

Raport jest plikiem html, który możesz pobrać i otworzyć w przeglądarce. Możesz również kliknąć na nim prawym przyciskiem myszy w eksploratorze plików po lewej stronie i wybrać `Show preview`, aby wyświetlić go w VS Code.

Poświęć kilka minut na przejrzenie raportu i sprawdź, czy możesz zidentyfikować możliwości dostosowania zasobów.
Pamiętaj, aby kliknąć zakładki pokazujące wyniki wykorzystania jako procent tego, co zostało przydzielone.
Dostępna jest [dokumentacja](https://www.nextflow.io/docs/latest/reports.html) opisująca wszystkie dostępne funkcje.

<!-- TODO: insert images -->

Jedną z obserwacji jest to, że `GATK_JOINTGENOTYPING` wydaje się być bardzo żądny CPU, co ma sens, ponieważ wykonuje wiele złożonych obliczeń.
Moglibyśmy więc spróbować zwiększyć to i sprawdzić, czy skróci to czas wykonania.

Jednak wygląda na to, że przesadziliśmy z przydziałem pamięci; wszystkie procesy wykorzystują tylko ułamek tego, co im dajemy.
Powinniśmy to zmniejszyć i zaoszczędzić trochę zasobów.

### 1.2. Dostosuj alokację zasobów dla konkretnego procesu

Możemy określić alokację zasobów dla danego procesu używając selektora `withName`.
Składnia wygląda następująco, gdy znajduje się samodzielnie w bloku process:

```groovy title="Składnia"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Dodajmy to do istniejącego bloku process w pliku `nextflow.config`.

```groovy title="nextflow.config" linenums="11"
process {
    // domyślne dla wszystkich procesów
    cpus = 2
    memory = 2.GB
    // alokacje dla konkretnego procesu
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Po określeniu tego, domyślne ustawienia będą stosowane do wszystkich procesów **z wyjątkiem** procesu `GATK_JOINTGENOTYPING`, który jest wyjątkowy i otrzymuje znacznie więcej CPU.
Miejmy nadzieję, że to powinno mieć efekt.

### 1.3. Uruchom ponownie ze zmodyfikowaną konfiguracją

Uruchommy workflow ponownie ze zmodyfikowaną konfiguracją i z włączonym raportowaniem, ale zauważ, że nadajemy raportowi inną nazwę, abyśmy mogli je rozróżnić.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Po raz kolejny prawdopodobnie nie zauważysz istotnej różnicy w czasie wykonania, ponieważ to jest tak małe obciążenie, a narzędzia spędzają więcej czasu na zadaniach pomocniczych niż na wykonywaniu 'właściwej' pracy.

Jednak drugi raport pokazuje, że nasze wykorzystanie zasobów jest teraz bardziej zrównoważone.

<!-- **TODO: screenshots?** -->

Jak widać, to podejście jest użyteczne, gdy Twoje procesy mają różne wymagania zasobowe. Daje Ci możliwość precyzyjnego dopasowania alokacji zasobów dla każdego procesu na podstawie rzeczywistych danych, a nie domysłów.

!!!note "Uwaga"

    To tylko mały przedsmak tego, co możesz zrobić, aby zoptymalizować wykorzystanie zasobów.
    Nextflow sam w sobie ma naprawdę sprytną wbudowaną [logikę dynamicznego ponawiania](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources), która ponawia zadania, które zawiodły z powodu ograniczeń zasobowych.
    Dodatkowo, Seqera Platform oferuje narzędzia oparte na AI do automatycznej optymalizacji alokacji zasobów.

    Omówimy oba te podejścia w nadchodzącej części tego kursu szkoleniowego.

Mimo to mogą istnieć pewne ograniczenia dotyczące tego, co możesz (lub musisz) przydzielić, w zależności od tego, jakiego executora obliczeniowego i infrastruktury obliczeniowej używasz. Na przykład Twój klaster może wymagać, abyś pozostał w określonych granicach, które nie mają zastosowania, gdy pracujesz w innym miejscu.
