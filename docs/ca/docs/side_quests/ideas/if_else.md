# Part 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Fes que la vaca citi científics famosos

Aquesta secció conté alguns exercicis addicionals per practicar el que has après fins ara.
Fer aquests exercicis _no és necessari_ per entendre les parts posteriors de la formació, però proporcionen una manera divertida de reforçar els teus aprenentatges descobrint com fer que la vaca citi científics famosos.

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

### 1.1. Modifica l'script `hello-containers.nf` per utilitzar un procés getQuote

Tenim una llista de pioners de la informàtica i la biologia al fitxer `containers/data/pioneers.csv`.
A grans trets, per completar aquest exercici hauràs de:

- Modificar el `params.input_file` per defecte perquè apunti al fitxer `pioneers.csv`.
- Crear un procés `getQuote` que utilitzi el contenidor `quote` per obtenir una cita per a cada entrada.
- Connectar la sortida del procés `getQuote` al procés `cowsay` per mostrar la cita.

Per a la imatge de contenidor `quote`, pots utilitzar la que has construït tu mateix en l'exercici addicional anterior o utilitzar la que has obtingut de Seqera Containers.

!!! Hint "Pista"

    Una bona opció per al bloc `script` del teu procés getQuote podria ser:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Pots trobar una solució a aquest exercici a `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modifica el teu pipeline de Nextflow per permetre que s'executi en els modes `quote` i `sayHello`.

Afegeix una mica de lògica de ramificació al teu pipeline per permetre que accepti entrades destinades tant a `quote` com a `sayHello`.
Aquí tens un exemple de com utilitzar una instrucció `if` en un workflow de Nextflow:

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

!!! Hint "Pista"

    Pots utilitzar `new_ch = processName.out` per assignar un nom al canal de sortida d'un procés.

Pots trobar una solució a aquest exercici a `containers/solutions/hello-containers-4.2.nf`.

### Conclusió

Ja saps com utilitzar contenidors en Nextflow per executar processos, i com construir una mica de lògica de ramificació als teus pipelines!

### Què segueix?

Celebra-ho, fes una pausa per estirar-te i beu una mica d'aigua!

Quan estiguis preparat/da, passa a la Part 3 d'aquesta sèrie de formació per aprendre com aplicar el que has après fins ara a un cas d'ús d'anàlisi de dades més realista.
