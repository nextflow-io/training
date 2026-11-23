# Part 4: Executar pipelines remots

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Fins ara, heu executat scripts de workflow emmagatzemats localment.
A la pràctica, sovint voldreu executar pipelines publicats en repositoris remots, com GitHub, sense haver-los de descarregar vosaltres mateixos.

Nextflow ho fa senzill: podeu executar qualsevol pipeline directament des d'una URL d'un repositori Git.

---

## 1. Executar un pipeline des de GitHub

La sintaxi bàsica per executar un pipeline remot és `nextflow run <repository>`, on `<repository>` pot ser un camí de repositori de GitHub com `nextflow-io/hello`, una URL completa, o un camí a GitLab, Bitbucket o un altre servei d'allotjament Git.

### 1.1. Llançar el pipeline

Executeu el pipeline de demostració oficial "hello" de Nextflow.
Aquest és un pipeline diferent, molt més senzill que el que heu estat executant en aquest curs: és anterior al pipeline "Hello" utilitzat al llarg d'aquesta formació, i simplement imprimeix una salutació per a cadascun d'uns quants idiomes predefinits, de manera que no espereu l'entrada CSV ni l'art ASCII als quals esteu acostumats.

```bash
nextflow run nextflow-io/hello
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Trobar on s'ha emmagatzemat en memòria cau el pipeline

La primera vegada que executeu un pipeline remot, Nextflow el descarrega i l'emmagatzema en memòria cau localment.
Les execucions posteriors reutilitzen la versió en memòria cau tret que sol·liciteu explícitament una actualització.

Per defecte, Nextflow desa els pipelines descarregats a `$NXF_HOME/assets`.
Per trobar on ha quedat un pipeline concret, i quines revisions estan disponibles, pregunteu-ho directament a Nextflow:

```bash
nextflow info nextflow-io/hello
```

??? success "Sortida de la comanda"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow marca cada revisió que ja heu obtingut localment amb `>`; la resta estan disponibles però encara no s'han descarregat en una còpia de treball.

També podeu llistar tots els pipelines que heu descarregat fins ara amb `nextflow list`:

```bash
nextflow list
```

??? success "Sortida de la comanda"

    ```console
    nextflow-io/hello
    ```

El curs [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) cobreix aquest mecanisme de memòria cau amb més profunditat, incloent-hi com explorar el codi font d'un pipeline descarregat.

### Conclusio

Sabeu com executar un pipeline directament des d'un repositori de GitHub sense descarregar-lo vosaltres mateixos, i on trobar-lo localment després.

### Què segueix?

Apreneu a fixar una versió específica d'un pipeline remot per a la reproductibilitat.

---

## 2. Especificar una versió per a la reproductibilitat

Per defecte, Nextflow executa la revisió més recent de la branca per defecte.
Podeu fixar una versió concreta (etiqueta), branca o commit utilitzant l'indicador `-r`.

### 2.1. Fixar una revisió específica

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow obté aquesta revisió la primera vegada que la sol·liciteu, d'aquí les línies `Pulling` i `downloaded from`; sol·licitar la mateixa revisió de nou més endavant va directament a `Launching`.
Fixar una revisió exacta és essencial per a la reproductibilitat.
Garanteix que vosaltres i els vostres col·laboradors executeu exactament el mateix codi de pipeline, independentment del que hagi canviat al repositori des d'aleshores.

### 2.2. Les revisions s'apliquen només per invocació

Fixar una revisió amb `-r` només afecta l'execució on l'especifiqueu: no canvia el que utilitza un `nextflow run` posterior sense indicadors.
Proveu d'executar el pipeline de nou sense `-r`:

```bash
nextflow run nextflow-io/hello
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Tot i que l'execució anterior fixava explícitament `v1.3`, aquesta execució torna directament a la branca per defecte (`master`).
Nextflow manté una còpia de treball local separada per a cada revisió que heu utilitzat, que és el que mostren els marcadors `>` a `nextflow info`, però mai recorda quina heu executat l'última vegada.
Podeu trobar el nom de la branca per defecte d'un pipeline executant `nextflow info <pipeline>`; és la marcada com a `(default)`.
La reproductibilitat és completament responsabilitat vostra: passeu sempre `-r` explícitament quan sigui important, en lloc d'assumir que una revisió fixada en una execució anterior encara s'aplica.

### Conclusio

Sabeu com fixar un pipeline remot a una versió, branca o commit específics per a una execució reproduïble, i que la fixació s'aplica només a aquella invocació, no a execucions posteriors.

### Què segueix?

Heu cobert els fonaments de l'execució i la gestió de pipelines de Nextflow.
Vegeu el [Resum del curs](next_steps.md) per saber on anar a partir d'aquí.

---

## Resum

En aquesta part heu après a:

- Executar un pipeline directament des d'un repositori de GitHub sense descarregar-lo
- Fixar un pipeline remot a una revisió específica per a la reproductibilitat, i entendre que la fixació s'aplica només a aquella invocació
