---
title: Hello Nextflow
hide:
    - toc
page_type: index_page
index_type: course
additional_information:
    technical_requirements: true
    learning_objectives:
        - Uruchamianie i zarządzanie wykonywaniem workflow'ów Nextflow
        - Znajdowanie i interpretowanie wyjść (wyników) i plików dziennika generowanych przez Nextflow
        - Rozwiązywanie podstawowych problemów
        - Budowanie prostego wieloetapowego workflow'u z podstawowych komponentów Nextflow
        - Rozróżnianie między podstawowymi typami fabryk kanałów i operatorów oraz efektywne ich wykorzystywanie w prostym workflow'ie
        - Konfigurowanie wykonywania pipeline'u do uruchamiania na popularnych platformach obliczeniowych, w tym HPC i chmurze
        - Stosowanie najlepszych praktyk dotyczących powtarzalności, przenośności i ponownego wykorzystania kodu, które czynią pipeline'y FAIR, w tym modularność kodu i kontenery oprogramowania
    audience_prerequisites:
        - "**Odbiorcy:** Ten kurs jest przeznaczony dla osób, które są całkowicie nowe w Nextflow i chcą tworzyć własne pipeline'y."
        - "**Umiejętności:** Zakłada się pewną znajomość wiersza poleceń, podstawowych koncepcji skryptowych i popularnych formatów plików."
        - "**Dziedzina:** Ćwiczenia są niezależne od dziedziny, więc nie jest wymagana wcześniejsza wiedza naukowa."
    videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow to praktyczne wprowadzenie do budowania powtarzalnych i skalowalnych workflow'ów analizy danych.**

Pracując z praktycznymi przykładami i prowadzonymi ćwiczeniami, poznasz podstawy tworzenia pipeline'ów w Nextflow, w tym jak definiować procesy, łączyć je w pipeline'y, zarządzać plikami i zależnościami oprogramowania, bez wysiłku paralelizować wykonywanie zadań i uruchamiać workflow'y w różnych środowiskach obliczeniowych.

Wyniesiesz umiejętności i pewność siebie, aby zacząć tworzyć i uruchamiać własne workflow'y w Nextflow.

<!-- additional_information -->

## Przegląd kursu

Ten kurs jest zaprojektowany jako praktyczny, z ćwiczeniami zorientowanymi na cel, strukturyzowanymi tak, aby wprowadzać informacje stopniowo.

Stworzysz prosty pipeline Nextflow, który przyjmuje tekstowe dane wejściowe, wykonuje kilka kroków transformacji i generuje pojedynczy plik tekstowy zawierający obraz ASCII postaci wypowiadającej przekształcony tekst.

### Plan lekcji

Aby uniknąć przytłoczenia Cię koncepcjami i kodem, podzieliliśmy to na sześć części, z których każda skupia się na konkretnych aspektach tworzenia pipeline'ów w Nextflow.

| Rozdział kursu                                        | Podsumowanie                                                                                                              | Szacowany czas |
| ----------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- | -------------- |
| [Część 1: Hello World](./01_hello_world.md)           | Podstawowe komponenty i zasady składania i uruchamiania workflow'u Nextflow                                               | 30 min         |
| [Część 2: Hello Channels](./02_hello_channels.md)     | Używanie kanałów i operatorów do przetwarzania wejść i bezwysiłkowej paralelizacji wykonywania                            | 45 min         |
| [Część 3: Hello Workflow](./03_hello_workflow.md)     | Używanie kanałów do łączenia wielu kroków i obsługi transferu danych między krokami                                       | 60 min         |
| [Część 4: Hello Modules](./04_hello_modules.md)       | Stosowanie zasad modularności kodu w celu zwiększenia możliwości ponownego użycia i zmniejszenia obciążenia konserwacyjne | 20 min         |
| [Część 5: Hello Containers](./05_hello_containers.md) | Używanie kontenerów jako mechanizmu zarządzania zależnościami oprogramowania i zwiększenia powtarzalności                 | 60 min         |
| [Część 6: Hello Config](./06_hello_config.md)         | Dostosowywanie zachowania pipeline'u i optymalizacja użycia w różnych środowiskach obliczeniowych                         | 60 min         |

Pod koniec tego kursu będziesz dobrze przygotowany do podjęcia kolejnych kroków w Swojej drodze do tworzenia powtarzalnych workflow'ów dla Twoich potrzeb obliczeniowych w nauce.

Gotowy do rozpoczęcia kursu?

[Rozpocznij :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
