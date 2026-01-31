# Parte 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Far citare dal cow scienziati famosi

Questa sezione contiene alcuni esercizi avanzati, per mettere in pratica quanto appreso finora.
Svolgere questi esercizi _non è obbligatorio_ per comprendere le parti successive della formazione, ma forniscono un modo divertente per consolidare le conoscenze imparando a far citare dal cow scienziati famosi.

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

### 1.1. Modificare lo script `hello-containers.nf` per utilizzare un processo getQuote

Abbiamo un elenco di pionieri dell'informatica e della biologia nel file `containers/data/pioneers.csv`.
Ad alto livello, per completare questo esercizio sarà necessario:

- Modificare il `params.input_file` predefinito per puntare al file `pioneers.csv`.
- Creare un processo `getQuote` che utilizzi il container `quote` per recuperare una citazione per ogni input.
- Collegare l'output del processo `getQuote` al processo `cowsay` per visualizzare la citazione.

Per l'immagine del container `quote`, è possibile utilizzare quella costruita personalmente nell'esercizio avanzato precedente oppure quella ottenuta da Seqera Containers.

!!! Hint

    Una buona scelta per il blocco `script` del processo getQuote potrebbe essere:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

È possibile trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modificare il pipeline Nextflow per consentirne l'esecuzione nelle modalità `quote` e `sayHello`.

Aggiungere della logica di ramificazione al pipeline per consentirgli di accettare input destinati sia a `quote` che a `sayHello`.
Ecco un esempio di come utilizzare un'istruzione `if` in un workflow Nextflow:

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

    È possibile utilizzare `new_ch = processName.out` per assegnare un nome al canale di output di un processo.

È possibile trovare una soluzione a questo esercizio in `containers/solutions/hello-containers-4.2.nf`.

### Conclusioni

Ora sa come utilizzare i container in Nextflow per eseguire processi e come costruire della logica di ramificazione nei suoi pipeline!

### Prossimi passi

Celebri, prendetevi una pausa per fare stretching e beva dell'acqua!

Quando è pronto, passi alla Parte 3 di questa serie di formazione per imparare ad applicare quanto appreso finora a un caso d'uso di analisi dati più realistico.
