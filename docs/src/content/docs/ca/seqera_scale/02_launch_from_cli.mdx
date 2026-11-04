# Part 2: Llançar pipelines des de la línia de comandes

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la [Part 1](./01_run_with_seqera.md), vau llançar nf-core/rnaseq des de la interfície web de Seqera.
Ara fem el mateix des de la línia de comandes utilitzant el CLI `tw`, i afegim un nou pipeline al vostre workspace.

---

## 1. Llançar pipelines des de la línia de comandes

A la vista d'execució, feu clic a la pestanya **Command line**.
Veureu la comanda `nextflow run` exacta que la Platform va construir i enviar en nom vostre — el mateix tipus de comanda que heu estat executant manualment al curs Use nf-core.

La Platform no substitueix Nextflow; l'orquestra.
Tot el que podeu fer a través de la interfície web, també ho podeu fer des d'un terminal utilitzant el CLI `tw`, l'eina de línia de comandes per interactuar amb l'API de la Platform.
Això és útil per automatitzar llançaments des de scripts o pipelines de CI/CD.

Ho farem ara des del mateix codespace que heu utilitzat per als cursos anteriors.

### 1.1. Instal·lar el CLI tw

Executeu les comandes següents al terminal del vostre Codespace per descarregar i instal·lar el binari `tw`:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Verifiqueu la instal·lació:

```bash
tw --version
```

??? success "Sortida de la comanda"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

El CLI `tw` està instal·lat i llest per configurar.

### 1.2. Obtenir un token d'accés

El CLI `tw` s'autentica amb Seqera utilitzant un token d'accés personal.

1. A la interfície web de Seqera, feu clic al vostre avatar a la cantonada superior dreta i seleccioneu **Your tokens**.
2. Feu clic a **Add token**, doneu-li un nom (p. ex. `training`), i feu clic a **Add**.
3. Copieu el valor del token — només es mostrarà una vegada.
   Si no el deseu en algun lloc de seguida, haureu de generar-ne un altre.

### 1.3. Configurar el CLI

Per comoditat, configurarem un fitxer de configuració que contingui el token d'accés que acabeu de generar i l'identificador del workspace.

Obriu el fitxer `.seqera_config` d'aquest directori a l'editor i establiu les dues variables:

- **`TOWER_ACCESS_TOKEN`**: el token que heu generat a la secció 1.2
- **`TOWER_WORKSPACE_ID`**: l'ID numèric del vostre workspace (la columna `ID` a `tw workspaces list`, que executeu a la secció 1.4)

Un cop els valors estiguin emplenats, carregueu la configuració:

```bash
source .seqera_config
```

Verifiqueu la connexió:

```bash
tw info
```

??? success "Sortida de la comanda"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

El CLI `tw` ara està autenticat i connectat al vostre compte de Seqera.
Executeu `source .seqera_config` a l'inici de cada sessió de Codespace per recarregar la configuració.

!!! tip "Consell"

    Si el vostre workspace no té un entorn de càlcul principal configurat, podeu afegir `export TOWER_COMPUTE_ENV=<compute-env-name>` al vostre fitxer de configuració per establir un valor per defecte.
    Qualsevol valor de configuració es pot sobreescriure a la línia de comandes passant l'indicador explícitament (p. ex. `--compute-env other-env`).
    Consulteu la [referència del CLI tw](https://docs.seqera.io/platform/latest/cli/reference) per a la llista completa d'opcions i variables d'entorn.

### 1.4. Explorar el vostre workspace des del CLI

Llisteu els workspaces als quals teniu accés:

```bash
tw workspaces list
```

??? success "Sortida de la comanda"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Visualitzeu les execucions del vostre workspace, incloent-hi l'execució de nf-core/rnaseq que acabeu de llançar:

```bash
tw runs list
```

??? success "Sortida de la comanda"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

La mateixa execució que esteu monitoritzant a la interfície web és visible aquí.

!!! note "Nota"

    Com que `TOWER_WORKSPACE_ID` està definit a `.seqera_config`, podeu ometre `--workspace` de totes les comandes `tw`.
    Sense la configuració, l'hauríeu de passar explícitament:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Tot el que és visible a la interfície web és accessible des del CLI.

### 1.5. Llançar nf-core/rnaseq des del CLI

El pipeline que heu afegit al vostre workspace a la [Part 1](./01_run_with_seqera.md) està disponible per nom al CLI.
Llanceu-lo amb el perfil `test`:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Sortida de la comanda"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Obriu l'enllaç al navegador i confirmeu que l'execució apareix al panell **Runs**.

Un cop pugueu veure-la en execució, haureu confirmat que el CLI i la interfície web són dues vistes del mateix workspace.

!!! note "Nota"

    També podeu passar una URL completa de GitHub directament a `tw launch` sense afegir prèviament el pipeline a un workspace.
    Tanmateix, afegir el pipeline explícitament abans de llançar-lo és generalment millor: desa la configuració del pipeline per a execucions futures, el fa disponible per nom, i el fa visible a tots els membres del workspace al Launchpad.

    És possible afegir un pipeline a un workspace directament des de la línia de comandes utilitzant `tw`.
    La secció següent mostra com fer-ho amb el pipeline nf-core/demo.

### Conclusio

Sabeu com autenticar el CLI `tw`, inspeccionar el vostre workspace i llançar un pipeline desat des del terminal.

### Què segueix?

Afegiu un nou pipeline al vostre workspace des de la línia de comandes i llanceu-lo.

---

## 2. Afegir un nou pipeline i executar-lo

Qualsevol pipeline de Nextflow a GitHub es pot afegir al vostre workspace amb `tw pipelines add`, sempre que tingui un punt d'entrada `main.nf` i un `nextflow.config` a l'arrel.
nf-core/demo és un bon exemple per practicar: ja l'heu executat al curs Use nf-core, de manera que sabeu què fa i què esperar.

### 2.1. Afegir nf-core/demo al vostre workspace

Executeu la comanda següent per registrar el pipeline al vostre workspace:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Sortida de la comanda"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

El pipeline ara està registrat i apareixerà al Launchpad.

### 2.2. Verificar que apareix al Launchpad

Llisteu els pipelines del vostre workspace per confirmar que s'ha afegit:

```bash
tw pipelines list
```

??? success "Sortida de la comanda"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Obriu el vostre workspace al navegador i feu clic a **Launchpad** per confirmar que nf-core/demo ara apareix al costat de nf-core/rnaseq.

!!! tip "Consell"

    També podeu afegir pipelines a través de la interfície web: a la barra lateral esquerra, feu clic a **Launchpad**, després a **Add pipeline**, i empleneu el formulari corresponentment.

Feu clic al botó **Launch** de l'entrada nf-core/demo per obrir el seu formulari de llançament.
Veureu que els paràmetres `input` i `outdir` estan ressaltats en vermell — són camps obligatoris sense valors per defecte, perquè `tw pipelines add` registra només la font del pipeline sense preconfigurar cap paràmetre.
Les dues seccions següents expliquen com proporcionar aquests valors: primer a través del formulari web, i després des de la línia de comandes.

### 2.3. Llançar nf-core/demo des de la interfície web

Amb el formulari de llançament obert, empleneu els dos paràmetres obligatoris.

Per a `input`, introduïu la URL del samplesheet de prova del perfil test de nf-core/demo.
La podeu trobar a `conf/test.config` dins del repositori del pipeline, que vau examinar al curs Use nf-core:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Per a `outdir`, introduïu un camí d'emmagatzematge al núvol on el pipeline pugui escriure els seus resultats.
Utilitzeu el bucket configurat per al vostre workspace, amb un subdirectori per mantenir les execucions organitzades:

```
s3://my-bucket/demo-results
```

Un cop els dos camps estiguin emplenats, feu clic al botó blau **Launch**.

L'execució apareix al panell **Runs** i hauria de completar-se en uns minuts amb el conjunt de dades de prova.
Feu clic a l'execució per explorar la taula de tasques i els informes d'execució.

### 2.4. Llançar nf-core/demo des del CLI

A diferència de `nextflow run`, la comanda `tw launch` no accepta indicadors de paràmetres individuals com `--input` o `--outdir`.
Els paràmetres s'han de proporcionar a través d'un fitxer en format YAML o JSON, passat amb `--params-file`.
Això fomenta la reproductibilitat: un fitxer de paràmetres desat documenta exactament quins valors s'han utilitzat per a una execució, facilitant repetir o compartir una configuració d'execució.

Creeu un fitxer de paràmetres al vostre directori de treball:

```bash
touch params.yaml
```

Obriu-lo a l'editor i afegiu el camí de sortida:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Ara podeu llançar el pipeline utilitzant el perfil `test` (que proporciona el samplesheet `input`) i el fitxer de paràmetres (que proporciona `outdir`):

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Sortida de la comanda"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Obriu l'enllaç per confirmar que l'execució apareix al panell **Runs**.

!!! tip "Consell"

    Podeu incloure el fitxer de paràmetres durant el pas de configuració inicial si voleu establir alguns valors per defecte, així com algunes propietats addicionals per coincidir amb el que hem fet abans a través del formulari web:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Conclusio

Sabeu com afegir qualsevol pipeline de Nextflow allotjat a GitHub al vostre workspace i llançar-lo, tant des de la interfície web emplenant els paràmetres manualment, com des del CLI `tw` combinant un perfil amb un fitxer de paràmetres.

---

## Resum

En aquesta part heu après a:

- Autenticar el CLI `tw` i llançar un pipeline desat des del terminal
- Afegir un nou pipeline de GitHub utilitzant el CLI i verificar que apareix al Launchpad
- Llançar un pipeline des de la interfície web de Seqera emplenant els paràmetres obligatoris manualment
- Llançar un pipeline des del CLI utilitzant un perfil de Nextflow i un fitxer de paràmetres
