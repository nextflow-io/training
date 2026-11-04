# Parte 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Facciamo citare alla mucca scienziati famosi

Questa sezione contiene alcuni esercizi avanzati per mettere in pratica ciò che avete imparato finora.
Completare questi esercizi _non è necessario_ per comprendere le parti successive della formazione, ma forniscono un modo divertente per consolidare le vostre conoscenze scoprendo come far citare alla mucca scienziati famosi.

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

### 1.1. Modificate lo script `hello-containers.nf` per utilizzare un processo getQuote

Abbiamo una lista di pionieri dell'informatica e della biologia nel file `containers/data/pioneers.csv`.
Ad alto livello, per completare questo esercizio dovrete:

- Modificare il `params.input_file` predefinito per puntare al file `pioneers.csv`.
- Creare un processo `getQuote` che utilizza il container `quote` per recuperare una citazione per ogni input.
- Collegare l'output del processo `getQuote` al processo `cowsay` per visualizzare la citazione.

Per l'immagine del container `quote`, potete utilizzare quella che avete costruito voi stessi nell'esercizio avanzato precedente oppure quella ottenuta da Seqera Containers.

!!! Hint "Suggerimento"

    Una buona scelta per il blocco `script` del vostro processo getQuote potrebbe essere:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Potete trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modificate la vostra pipeline Nextflow per consentirle di eseguire nelle modalità `quote` e `sayHello`.

Aggiungete della logica condizionale alla vostra pipeline per consentirle di accettare input destinati sia a `quote` che a `sayHello`.
Ecco un esempio di come utilizzare un'istruzione `if` in un flusso di lavoro Nextflow:

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

!!! Hint "Suggerimento"

    Potete usare `new_ch = processName.out` per assegnare un nome al canale di output di un processo.

Potete trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.2.nf`.

### Takeaway

Sapete come utilizzare i container in Nextflow per eseguire processi e come costruire della logica condizionale nelle vostre pipeline!

### Cosa c'è dopo?

Festeggiate, prendetevi una pausa per sgranchirvi e bevete un po' d'acqua!

Quando siete pronti, passate alla Parte 3 di questa serie di formazione per imparare come applicare ciò che avete imparato finora a un caso d'uso di analisi dati più realistico.
