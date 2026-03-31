# Podsumowanie

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ukończyłeś szkolenie z tworzenia wtyczek.
Ta strona podsumowuje to, co zbudowałeś w każdej części, omawia dystrybucję i wskazuje, co robić dalej.

---

## Czego się nauczyłeś

### Część 1: Korzystanie z wtyczek

Poznałeś sposób działania wtyczek Nextflow z perspektywy użytkownika.
Zainstalowałeś nf-schema i nf-co2footprint, skonfigurowałeś je przez `nextflow.config` i zobaczyłeś, jak wtyczki mogą walidować dane wejściowe, dodawać funkcje oraz podłączać się do zdarzeń cyklu życia pipeline'u.

### Część 2: Konfiguracja środowiska

Skonfigurowałeś środowisko deweloperskie z Java 21+, utworzyłeś nowy projekt wtyczki przy użyciu polecenia `nextflow plugin create` i poznałeś strukturę projektu oczekiwaną przez Nextflow: pliki źródłowe, konfigurację budowania i workflow Makefile.

### Część 3: Własne funkcje

Zaimplementowałeś swój pierwszy punkt rozszerzenia, tworząc metody z adnotacją `@Function` w klasie `PluginExtensionPoint`.
Zbudowałeś `reverseGreeting` i `decorateGreeting`, a następnie zaimportowałeś je i wywołałeś ze skryptu pipeline'u.

### Część 4: Testowanie

Napisałeś testy jednostkowe dla swoich funkcji, korzystając z frameworka testowego Groovy.
Nauczyłeś się uruchamiać testy poleceniem `make test` i weryfikować poprawność działania wtyczki przed jej instalacją.

### Część 5: Obserwatory

Zaimplementowałeś interfejs `TraceObserver`, aby podłączyć się do zdarzeń cyklu życia pipeline'u.
Zbudowałeś `GreetingObserver` (reagujący na start i zakończenie pipeline'u) oraz `TaskCounterObserver` (zliczający ukończone zadania), a następnie zarejestrowałeś je przez `TraceObserverFactory`.

### Część 6: Konfiguracja

Uczyniłeś wtyczkę konfigurowalną przez `nextflow.config`, używając `session.config.navigate()` do odczytywania wartości w czasie wykonania.
Dodałeś klasę `@ConfigScope`, aby formalnie zadeklarować opcje wtyczki — dzięki temu wyeliminowałeś ostrzeżenia „Unrecognized config option" i włączyłeś wsparcie IDE.

---

## Dystrybucja

Gdy Twoja wtyczka działa lokalnie, możesz udostępnić ją innym przez rejestr wtyczek Nextflow.

### Wersjonowanie

Stosuj [wersjonowanie semantyczne](https://semver.org/) dla swoich wydań:

| Zmiana wersji             | Kiedy stosować                        | Przykład                                              |
| ------------------------- | ------------------------------------- | ----------------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Zmiany niekompatybilne wstecz         | Usunięcie funkcji, zmiana typów zwracanych wartości   |
| **MINOR** (1.0.0 → 1.1.0) | Nowe funkcje, kompatybilne wstecz     | Dodanie nowej funkcji                                 |
| **PATCH** (1.0.0 → 1.0.1) | Poprawki błędów, kompatybilne wstecz  | Naprawa błędu w istniejącej funkcji                   |

Zaktualizuj wersję w `build.gradle` przed każdym wydaniem:

```groovy title="build.gradle"
version = '1.0.0'  // Stosuj wersjonowanie semantyczne: MAJOR.MINOR.PATCH
```

### Publikowanie w rejestrze

[Rejestr wtyczek Nextflow](https://registry.nextflow.io/) to oficjalny sposób udostępniania wtyczek społeczności.

Proces publikowania:

1. Zarezerwuj nazwę swojej wtyczki w [rejestrze](https://registry.nextflow.io/) (zaloguj się kontem GitHub)
2. Skonfiguruj dane uwierzytelniające API w `~/.gradle/gradle.properties`
3. Uruchom testy, aby sprawdzić poprawność działania: `make test`
4. Opublikuj poleceniem `make release`

Instrukcje krok po kroku znajdziesz w [oficjalnej dokumentacji publikowania](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Po opublikowaniu użytkownicy instalują Twoją wtyczkę bez żadnej lokalnej konfiguracji:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow automatycznie pobiera wtyczkę z rejestru przy pierwszym użyciu.

---

## Lista kontrolna tworzenia wtyczek

- [ ] Java 21+ zainstalowana
- [ ] Utwórz projekt poleceniem `nextflow plugin create <name> <org>`
- [ ] Zaimplementuj klasę rozszerzenia z metodami `@Function`
- [ ] Napisz testy jednostkowe i uruchom je poleceniem `make test`
- [ ] Zbuduj i zainstaluj poleceniem `make install`
- [ ] Opcjonalnie dodaj implementacje `TraceObserver` dla zdarzeń workflow'u
- [ ] Opcjonalnie dodaj `ConfigScope` dla konfiguracji wtyczki
- [ ] Włącz w `nextflow.config` przez `plugins { id 'plugin-id' }`
- [ ] Importuj funkcje przez `include { fn } from 'plugin/plugin-id'`
- [ ] Nadaj wersję i opublikuj w rejestrze

---

## Kluczowe wzorce kodu

**Definicja funkcji:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Konfiguracja wtyczki:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Użycie w workflow'ach:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Podsumowanie punktów rozszerzenia

| Typ                  | Klasa / adnotacja | Przeznaczenie                                              |
| -------------------- | ----------------- | ---------------------------------------------------------- |
| Funkcja              | `@Function`       | Wywoływalna z workflow'ów                                  |
| Obserwator śladów    | `TraceObserver`   | Podłączenie do zdarzeń cyklu życia workflow'u              |
| Zakres konfiguracji  | `@ScopeName`      | Definiowanie konfiguracji wtyczki w nextflow.config        |

---

## Co dalej?

Oto kilka praktycznych kroków, które pozwolą Ci kontynuować naukę tworzenia wtyczek.

**Zbuduj coś prawdziwego.**
Wybierz przypadek użycia z własnej pracy: funkcję, której Twój zespół używa wielokrotnie, obserwatora wysyłającego powiadomienia Slack po zakończeniu pipeline'u albo zakres konfiguracji standaryzujący opcje w pipeline'ach Twojej organizacji.
Zaczynanie od realnego problemu to najszybszy sposób na pogłębienie wiedzy.

**Korzystaj z nf-hello jako punktu odniesienia.**
Repozytorium [nf-hello](https://github.com/nextflow-io/nf-hello) to oficjalny minimalny przykład wtyczki.
Stanowi dobry punkt startowy dla nowych projektów i przydatne źródło wiedzy, gdy chcesz sprawdzić, jak coś jest zorganizowane.

**Czytaj oficjalną dokumentację.**
Dokumentacja Nextflow obejmuje tematy wykraczające poza to szkolenie, w tym fabryki kanałów, przeciążanie operatorów i zaawansowane wzorce obserwatorów.
Przewodnik [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) to najbardziej wyczerpane źródło.

**Studiuj istniejące wtyczki.**
[Repozytorium wtyczek Nextflow](https://github.com/nextflow-io/plugins) zawiera kod źródłowy oficjalnych wtyczek, takich jak nf-schema, nf-wave i nf-tower.
Czytanie kodu produkcyjnych wtyczek to jeden z najlepszych sposobów na poznanie wzorców i konwencji wykraczających poza przykłady wprowadzające.

---

## Dodatkowe zasoby

**Oficjalna dokumentacja:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): kompleksowy przewodnik po instalacji i konfiguracji wtyczek
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): szczegółowe źródło wiedzy o tworzeniu wtyczek
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): tworzenie zakresów konfiguracji dla wtyczek

**Odkrywanie wtyczek:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): przeglądaj i odkrywaj dostępne wtyczki
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): dokumentacja rejestru

**Przykłady i materiały referencyjne:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): prosty przykład wtyczki (świetny punkt startowy)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): zbiór oficjalnych wtyczek do wykorzystania jako materiał referencyjny
