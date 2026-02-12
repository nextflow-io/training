# Part 1: Més contenidors

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Com trobar o crear imatges de contenidor

Alguns desenvolupadors de programari proporcionen imatges de contenidor per al seu programari que estan disponibles en registres de contenidors com Docker Hub, però molts no ho fan.
En aquesta secció opcional, us mostrarem dues maneres d'obtenir una imatge de contenidor per a les eines que voleu utilitzar en els vostres pipelines de Nextflow: utilitzant Seqera Containers i construint la imatge de contenidor vosaltres mateixos.

Obtindreu/construireu una imatge de contenidor per al paquet pip `quote`, que s'utilitzarà en l'exercici al final d'aquesta secció.

### 1.1. Obtenir una imatge de contenidor de Seqera Containers

Seqera Containers és un servei gratuït que construeix imatges de contenidor per a eines instal·lables amb pip i conda (incloent bioconda).
Navegueu a [Seqera Containers](https://www.seqera.io/containers/) i cerqueu el paquet pip `quote`.

![Seqera Containers](img/seqera-containers-1.png)

Feu clic a "+Add" i després a "Get Container" per sol·licitar una imatge de contenidor per al paquet pip `quote`.

![Seqera Containers](img/seqera-containers-2.png)

Si aquesta és la primera vegada que es construeix un contenidor comunitari per a aquesta versió del paquet, pot trigar uns minuts a completar-se.
Feu clic per copiar l'URI (p. ex. `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) de la imatge de contenidor que s'ha creat per a vosaltres.

Ara podeu utilitzar la imatge de contenidor per executar la comanda `quote` i obtenir una frase aleatòria de Grace Hopper.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Sortida:

```console title="Output"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Construir la imatge de contenidor vosaltres mateixos

Utilitzem alguns detalls de construcció del lloc web de Seqera Containers per construir nosaltres mateixos la imatge de contenidor per al paquet pip `quote`.
Torneu al lloc web de Seqera Containers i feu clic al botó "Build Details".

El primer element que examinarem és el `Dockerfile`, un tipus de fitxer d'script que conté totes les comandes necessàries per construir la imatge de contenidor.
Hem afegit alguns comentaris explicatius al Dockerfile següent per ajudar-vos a entendre què fa cada part.

```Dockerfile title="Dockerfile"
# Comença des de la imatge base de docker micromamba
FROM mambaorg/micromamba:1.5.10-noble
# Copia el fitxer conda.yml dins del contenidor
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Instal·la diverses utilitats perquè Nextflow les utilitzi i els paquets del fitxer conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Executa el contenidor com a usuari root
USER root
# Estableix la variable d'entorn PATH per incloure el directori d'instal·lació de micromamba
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

El segon element que examinarem és el fitxer `conda.yml`, que conté la llista de paquets que cal instal·lar a la imatge de contenidor.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Copieu el contingut d'aquests fitxers als esbossos situats al directori `containers/build`, després executeu la comanda següent per construir la imatge de contenidor vosaltres mateixos.

!!! Note "Nota"

    Utilitzem la bandera `-t quote:latest` per etiquetar la imatge de contenidor amb el nom `quote` i l'etiqueta `latest`.
    Podrem utilitzar aquesta etiqueta per referir-nos a la imatge de contenidor quan l'executem en aquest sistema.

```bash
docker build -t quote:latest containers/build
```

Després que hagi acabat de construir-se, podeu executar la imatge de contenidor que acabeu de construir.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Conclusió

Heu après dues maneres diferents d'obtenir una imatge de contenidor per a una eina que voleu utilitzar en els vostres pipelines de Nextflow: utilitzant Seqera Containers i construint la imatge de contenidor vosaltres mateixos.

### Què segueix?

Teniu tot el que necessiteu per continuar al [següent capítol](./04_hello_genomics.md) d'aquesta sèrie de formació.
També podeu continuar amb un exercici opcional per obtenir cites sobre pioners de la informàtica/biologia utilitzant el contenidor `quote` i mostrar-les utilitzant el contenidor `cowsay`.

---

## 2. Feu que la vaca citi científics famosos

Aquesta secció conté alguns exercicis addicionals, per practicar el que heu après fins ara.
Fer aquests exercicis _no és necessari_ per entendre les parts posteriors de la formació, però proporcionen una manera divertida de reforçar els vostres aprenentatges descobrint com fer que la vaca citi científics famosos.

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

### 2.1. Modifiqueu l'script `hello-containers.nf` per utilitzar un procés getQuote

Tenim una llista de pioners de la informàtica i la biologia al fitxer `containers/data/pioneers.csv`.
A grans trets, per completar aquest exercici haureu de:

- Modificar el `params.input_file` per defecte per apuntar al fitxer `pioneers.csv`.
- Crear un procés `getQuote` que utilitzi el contenidor `quote` per obtenir una cita per a cada entrada.
- Connectar la sortida del procés `getQuote` al procés `cowsay` per mostrar la cita.

Per a la imatge de contenidor `quote`, podeu utilitzar la que heu construït vosaltres mateixos en l'exercici addicional anterior o utilitzar la que heu obtingut de Seqera Containers.

!!! Hint "Pista"

    Una bona opció per al bloc `script` del vostre procés getQuote podria ser:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Podeu trobar una solució a aquest exercici a `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modifiqueu el vostre pipeline de Nextflow per permetre que s'executi en els modes `quote` i `sayHello`.

Afegiu alguna lògica de ramificació al vostre pipeline per permetre que accepti entrades destinades tant a `quote` com a `sayHello`.
Aquí teniu un exemple de com utilitzar una instrucció `if` en un workflow de Nextflow:

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

    Podeu utilitzar `new_ch = processName.out` per assignar un nom al canal de sortida d'un procés.

Podeu trobar una solució a aquest exercici a `containers/solutions/hello-containers-4.2.nf`.

### Conclusió

Sabeu com utilitzar contenidors en Nextflow per executar processos, i com construir alguna lògica de ramificació en els vostres pipelines!

### Què segueix?

Celebreu-ho, feu una pausa per estirar-vos i beure una mica d'aigua!

Quan estigueu preparats, passeu a la Part 3 d'aquesta sèrie de formació per aprendre com aplicar el que heu après fins ara a un cas d'ús d'anàlisi de dades més realista.
