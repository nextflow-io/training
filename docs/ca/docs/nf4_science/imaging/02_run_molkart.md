# Part 2: Executar nf-core/molkart

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la Part 1, vam executar un workflow senzill Hello World per entendre els conceptes bàsics de l'execució de Nextflow.
Ara executarem un pipeline real de bioimatge: **nf-core/molkart**.

Aquest pipeline processa dades de transcriptòmica espacial de Molecular Cartography de Resolve Bioscience.
No obstant això, els patrons de Nextflow que aprendreu aquí s'apliquen a qualsevol pipeline nf-core o workflow de producció.

## 1. Comprendre els pipelines nf-core

Abans d'executar el pipeline, entenguem què és nf-core i per què és important per executar workflows.

### 1.1. Què és nf-core?

[nf-core](https://nf-co.re/) és una col·lecció impulsada per la comunitat de pipelines Nextflow d'alta qualitat.
Tots els pipelines nf-core segueixen la mateixa estructura i convencions, la qual cosa significa que un cop apreneu a executar-ne un, podeu executar qualsevol d'ells.

Característiques clau dels pipelines nf-core:

- **Estructura estandarditzada**: Tots els pipelines tenen noms de paràmetres i patrons d'ús consistents
- **Dades de prova integrades**: Cada pipeline inclou perfils de prova per a una validació ràpida
- **Documentació completa**: Instruccions d'ús detallades i descripcions de paràmetres
- **Control de qualitat**: Informes de QC automatitzats utilitzant MultiQC
- **Suport de contenidors**: Contenidors preconstruïts per a la reproducibilitat

!!! tip "Voleu aprendre més sobre nf-core?"

    Per a una introducció en profunditat al desenvolupament de pipelines nf-core, consulteu el curs de formació [Build with nf-core](../../hello_nf-core/index.md).
    Cobreix com crear i personalitzar pipelines nf-core des de zero.

### 1.2. El pipeline molkart

![Pipeline nf-core/molkart](img/molkart.png)

El pipeline [nf-core/molkart](https://nf-co.re/molkart) processa dades d'imatge de transcriptòmica espacial a través de diverses etapes:

1. **Preprocessament d'imatges**: Emplenament de patró de graella i millora de contrast opcional
2. **Segmentació cel·lular**: Múltiples opcions d'algorisme (Cellpose, Mesmer, ilastik, Stardist)
3. **Assignació de punts**: Assignar punts de transcripció a cèl·lules segmentades
4. **Control de qualitat**: Generar informes de QC complets

Les sortides clau són:

- Taules de recompte de cèl·lula per transcripció
- Màscares de segmentació
- Informe de control de qualitat MultiQC

---

## 2. Executar molkart amb dades de prova

Abans de començar, clonem el repositori molkart localment perquè puguem inspeccionar el seu codi:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Això crea un directori `molkart/` que conté el codi font complet del pipeline.

!!! note "Per què estem clonant localment?"

    Normalment, executaríeu pipelines nf-core directament des de GitHub utilitzant `nextflow run nf-core/molkart -r 1.2.0`.
    Nextflow descarrega automàticament la versió del pipeline sol·licitada a `$HOME/.nextflow/assets/nf-core/molkart` i l'executa des d'allà.
    No obstant això, per a aquesta formació, estem clonant el pipeline a un directori local diferent perquè puguem inspeccionar el codi més fàcilment.

### 2.1. Comprendre els requisits de contenidors

Abans d'executar el pipeline complet, aprenguem per què els contenidors són essencials per als pipelines nf-core.

Subministrarem els paràmetres del pipeline mitjançant un fitxer de paràmetres.
Un fitxer de paràmetres és un fitxer YAML que llista cada paràmetre i el seu valor, la qual cosa manté els valors tipats (com ara enters) intactes i manté la línia de comandes curta.

Ja hi ha un fitxer `params.yaml` al directori de treball:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Aquests paràmetres són:

- `input`: Camí al full de mostres que conté les metadades de la mostra
- `mindagap_tilesize`, `mindagap_boxsize`, `mindagap_loopnum`: Paràmetres per a l'emplenament del patró de graella
- `clahe_pyramid_tile`: Mida del nucli per a la millora del contrast
- `segmentation_method`: Quin(s) algorisme(s) utilitzar per a la segmentació cel·lular
- `outdir`: On desar els resultats

Provem d'executar el pipeline utilitzant aquests paràmetres:

```bash
nextflow run ./molkart -params-file params.yaml
```

!!! Warning "Aquesta comanda fallarà - és intencionat!"

    Estem executant això deliberadament sense contenidors per demostrar per què són necessaris.

Després d'uns moments, veureu un error com aquest:

??? failure "Sortida de la comanda"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**Què està passant aquí?**

L'error `command not found` (estat de sortida 127) significa que Nextflow va intentar executar `duplicate_finder.py` però no el va poder trobar al vostre sistema.
Això és perquè:

1. El pipeline espera que estigui instal·lat programari especialitzat de bioinformàtica
2. Aquestes eines (com `duplicate_finder.py`, `apply_clahe.dask.py`, etc.) no formen part de les distribucions estàndard de Linux
3. Sense contenidors, Nextflow intenta executar comandes directament a la vostra màquina local

**D'on se suposa que provenen aquestes eines?**

Inspeccionem un dels mòduls de procés per veure com declara els seus requisits de programari.

Obriu el mòdul de preprocessament CLAHE:

```bash
code molkart/modules/local/clahe/main.nf
```

Mireu la línia 5 - veureu:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Aquesta línia diu a Nextflow: "Per executar aquest procés, utilitza la imatge Docker `ghcr.io/schapirolabor/molkart-local:v0.0.4`, que conté tot el programari necessari."

Cada procés declara quina imatge de contenidor proporciona les seves eines necessàries.
No obstant això, Nextflow només utilitza aquests contenidors si li ho indiqueu!

**La solució: Activar Docker a la configuració**

### 2.2. Configurar Docker i llançar el pipeline

Per activar Docker, hem de canviar `docker.enabled` de `false` a `true` al fitxer `nextflow.config`.

Obriu el fitxer de configuració:

```bash
code nextflow.config
```

Canvieu `docker.enabled = false` a `docker.enabled = true`:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Ara executeu el pipeline de nou, aquesta vegada executant els tres mètodes de segmentació perquè puguem comparar-los més endavant:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist"
```

Aquesta vegada, Nextflow:

1. Llegirà la configuració `docker.enabled = true` de la configuració
2. Descarregarà les imatges Docker necessàries (només la primera vegada)
3. Executarà cada procés dins del seu contenidor especificat
4. S'executarà amb èxit perquè totes les eines estan disponibles dins dels contenidors

!!! Tip "Per què són importants els contenidors"

    La majoria dels pipelines nf-core **requereixen** contenidorització (Docker, Singularity, Podman, etc.) perquè:

    - Utilitzen programari especialitzat de bioinformàtica no disponible en entorns estàndard
    - Els contenidors garanteixen la reproducibilitat - les mateixes versions exactes de programari s'executen a tot arreu
    - No cal instal·lar manualment desenes d'eines i les seves dependències

    Per a més detalls sobre contenidors a Nextflow, consulteu [Hello Containers](../../hello_nextflow/05_hello_containers.md) de la formació Hello Nextflow.

### 2.3. Monitoritzar l'execució

Mentre s'executa el pipeline, veureu una sortida similar a aquesta:

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `./molkart/main.nf` [exotic_banach] revision: e4152308ec

    WARN: Unrecognized config option 'validation.help.enabled'
    WARN: Unrecognized config option 'validation.defaultIgnoreParams'
    WARN: Unrecognized config option 'validation.monochromeLogs'

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method: mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize   : 7
      mindagap_loopnum   : 100
      clahe_kernel       : 25
      mindagap_tilesize  : 90
      clahe_pyramid_tile : 368

    Input/output options
      input              : data/samplesheet.csv
      outdir             : results

    Generic options
      trace_report_suffix: 2026-06-23_16-25-30

    Core Nextflow options
      runName            : exotic_banach
      containerEngine    : docker
      launchDir          : /workspaces/training/nf4-science/imaging
      workDir            : /workspaces/training/nf4-science/imaging/work
      projectDir         : /workspaces/training/nf4-science/imaging/molkart
      userName           : root
      profile            : standard
      configFiles        : /workspaces/training/nf4-science/imaging/molkart/nextflow.config, /workspaces/training/nf4-science/imaging/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [b4/e57ff1] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2 ✔
    [2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2 ✔
    [94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1 ✔
    [58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1 ✔
    [c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1 ✔
    [4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1 ✔
    [8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1 ✔
    [08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3 ✔
    [56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3 ✔
    [b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3 ✔
    [5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3 ✔
    [8e/8bd365] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 23-Jun-2026 16:31:40
    Duration    : 6m 9s
    CPU hours   : (a few seconds)
    Succeeded   : 22
    ```

Observeu com aquesta sortida és més detallada que el nostre exemple Hello World a causa de les convencions nf-core que segueix el pipeline:

- El pipeline mostra la seva versió i logotip
- Es mostren els paràmetres de configuració
- Múltiples processos s'executen en paral·lel (indicat per múltiples línies de procés)
- Els noms de procés inclouen el camí complet del mòdul (p. ex., `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Comprendre l'execució de processos

La línia executor `executor > local (22)` us indica:

- **executor**: Quin entorn de càlcul s'està utilitzant (`local` = la vostra màquina)
- **(22)**: Nombre total de tasques llançades

Cada línia de procés mostra:

- **Hash** (`[b4/e57ff1]`): Identificador del directori de treball (com abans)
- **Nom del procés**: Camí complet del mòdul i nom del procés
- **Identificador d'entrada**: Nom de la mostra entre parèntesis
- **Progrés**: Recompte de tasques i estat de compleció (p. ex., `1 of 1 ✔`)

### Conclusió

Sabeu com llançar un pipeline nf-core amb dades de prova i interpretar la seva sortida d'execució.

### Què segueix?

Apreneu on trobar els resultats i com interpretar-los.

---

## 3. Trobar i examinar les sortides

Quan el pipeline es completa amb èxit, veureu un missatge de finalització i un resum d'execució.

### 3.1. Localitzar el directori de resultats

Per defecte, els pipelines nf-core escriuen sortides a un directori especificat pel parametre `outdir`, que hem establert a `results/`.

Llisteu els continguts:

```bash
tree results/
```

Hauríeu de veure diversos subdirectoris:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Cada subdirectori conté sortides d'una etapa específica del pipeline:

- **mindagap/**: Imatges emplenades de graella del pas de preprocessament MindaGap
- **clahe/**: Imatges amb contrast millorat del preprocessament CLAHE
- **stack/**: Piles d'imatges multicanal creades per a la segmentació
- **segmentation/**: Resultats de segmentació de diferents algorismes (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Taules de recompte de cèl·lula per transcripció
- **anndata/**: Objectes AnnData que contenen matrius de cèl·lula per transcripció i coordenades espacials
- **molkartqc/**: Mètriques de control de qualitat per a l'assignació de punts
- **multiqc/**: Informe de control de qualitat complet
- **pipeline_info/**: Informes d'execució i registres

### 3.2. Examinar l'informe MultiQC

L'informe MultiQC és un fitxer HTML complet que agrega mètriques de qualitat de tots els passos del pipeline.

Obriu l'informe al navegador de fitxers i després feu clic al botó "Show Preview" per veure'l renderitzat directament a VS Code.

L'informe inclou:

- Estadístiques generals per a totes les mostres
- Mètriques de preprocessament
- Mètriques de qualitat de segmentació
- Nombre de cèl·lules i punts detectats

!!! Tip "Consell"

    Els informes MultiQC normalment s'inclouen en tots els pipelines nf-core.
    Sempre proporcionen una visió general d'alt nivell de l'execució del pipeline i la qualitat de les dades.

### 3.3. Examinar les taules de cèl·lula per transcripció

La sortida científica més important és la taula de recompte de cèl·lula per transcripció.
Això us indica quantes de cada transcripció es van detectar a cada cèl·lula.

Navegueu al directori spot2cell:

```bash
ls results/spot2cell/
```

Trobareu fitxers com:

- `cellxgene_mem_only_cellpose.csv`: Taula de cèl·lula per transcripció utilitzant segmentació Cellpose
- `cellxgene_mem_only_mesmer.csv`: Taula de cèl·lula per transcripció utilitzant segmentació Mesmer
- `cellxgene_mem_only_stardist.csv`: Taula de cèl·lula per transcripció utilitzant segmentació Stardist

Només hem executat 1 mostra en aquest conjunt de dades de prova, però en un experiment real tindríem aquestes taules per a cada mostra.
Observeu com Nextflow és capaç de processar múltiples mètodes de segmentació en paral·lel, facilitant la comparació de resultats.

### 3.4. Veure informes d'execució

Nextflow genera diversos informes d'execució automàticament.

Comproveu el directori pipeline_info:

```bash
ls results/pipeline_info/
```

Fitxers clau:

- **execution_report.html**: Visualització de la línia de temps i l'ús de recursos
- **execution_timeline.html**: Diagrama de Gantt de l'execució del procés
- **execution_trace.txt**: Mètriques detallades d'execució de tasques
- **pipeline_dag.html**: Graf acíclic dirigit que mostra l'estructura del workflow

Obriu l'informe d'execució per veure l'ús de recursos:

```bash
code results/pipeline_info/execution_report.html
```

Això mostra:

- Quant de temps va trigar cada procés
- Ús de CPU i memòria
- Quines tasques es van emmagatzemar a la memòria cau vs. executades

!!! Tip "Consell"

    Aquests informes són increïblement útils per optimitzar l'assignació de recursos i solucionar problemes de rendiment.

### Conclusió

Sabeu com localitzar sortides del pipeline, examinar informes de control de qualitat i accedir a mètriques d'execució.

### Què segueix?

Apreneu sobre el directori de treball i com Nextflow gestiona els fitxers intermedis.

---

## 4. Explorar el directori de treball

Igual que amb el nostre exemple Hello World, tot el treball real passa al directori `work/`.

### 4.1. Comprendre l'estructura del directori de treball

El directori de treball conté un subdirectori per a cada tasca que es va executar.
Per a aquesta execució del pipeline amb 22 tasques, hi haurà 22 subdirectoris de treball.

Llisteu el directori de treball:

```bash
ls -d work/*/*/ | head -5
```

Això mostra els primers 5 directoris de tasques.

### 4.2. Inspeccionar un directori de tasques

Trieu un dels hashs del procés de segmentació de la sortida de la consola (p. ex., `[3m/4n5o6p]`) i mireu a dins:

```bash
ls -la work/3m/4n5o6p*/
```

Veureu:

- **Fitxers .command.\***: Scripts i registres d'execució de Nextflow (com abans)
- **Fitxers d'entrada preparats**: Enllaços simbòlics als fitxers d'entrada reals
- **Fitxers de sortida**: Màscares de segmentació, resultats intermedis, etc.

La diferència clau respecte a Hello World:

- Els pipelines reals preparen fitxers d'entrada grans (imatges, dades de referència)
- Els fitxers de sortida poden ser força grans (màscares de segmentació, imatges processades)
- Múltiples fitxers d'entrada i sortida per tasca

!!! Tip "Consell"

    Si un procés falla, podeu navegar al seu directori de treball, examinar `.command.err` per a missatges d'error, i fins i tot tornar a executar `.command.sh` manualment per depurar el problema.

### 4.3. Neteja del directori de treball

El directori de treball pot arribar a ser força gran després de múltiples execucions del pipeline.
Com vam aprendre a la Part 1, podeu utilitzar `nextflow clean` per eliminar directoris de treball d'execucions antigues.

No obstant això, per als pipelines nf-core amb fitxers intermedis grans, és especialment important netejar regularment.

### Conclusió

Enteneu com els pipelines nf-core organitzen els seus directoris de treball i com inspeccionar tasques individuals per a la depuració.

### Què segueix?

Apreneu sobre la memòria cau de Nextflow i com reprendre execucions de pipeline fallides.

---

## 5. Reprendre una execució de pipeline

Una de les característiques més potents de Nextflow és la capacitat de reprendre un pipeline des del punt de fallada.

### 5.1. El mecanisme de memòria cau

Quan executeu un pipeline amb `-resume`, Nextflow:

1. Comprova la memòria cau per a cada tasca
2. Si les entrades, el codi i els paràmetres són idèntics, reutilitza el resultat emmagatzemat a la memòria cau
3. Només torna a executar tasques que van canviar o van fallar

Això és essencial per a pipelines de llarga durada on les fallades poden ocórrer tard en l'execució.

### 5.2. Provar resume amb molkart

Executeu la mateixa comanda de nou, però afegiu `-resume`:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist" -resume
```

Hauríeu de veure una sortida com:

```console
executor >  local (1)
[43/e03702] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2, cached: 2 ✔
[2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2, cached: 2 ✔
[94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1, cached: 1 ✔
[58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1, cached: 1 ✔
[c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1, cached: 1 ✔
[4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1, cached: 1 ✔
[8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1, cached: 1 ✔
[08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3, cached: 3 ✔
[56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3, cached: 3 ✔
[b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3, cached: 3 ✔
[5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3, cached: 3 ✔
[73/239f45] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
-[nf-core/molkart] Pipeline completed successfully-
```

Observeu l'anotació `cached: N` a cada procés de preprocessament i segmentació - aquestes tasques s'han reutilitzat en lloc de tornar-se a executar.

### 5.3. Quan resume és útil

Resume és particularment valuós quan:

- Un pipeline falla a causa de límits de recursos (sense memòria, límit de temps excedit)
- Necessiteu modificar processos posteriors sense tornar a executar passos anteriors
- La vostra connexió de xarxa es desconnecta durant la descàrrega de dades
- Voleu afegir sortides addicionals sense refer el càlcul

!!! Warning "Advertència"

    Resume només funciona si no heu canviat les dades d'entrada, el codi del pipeline o els paràmetres.
    Si canvieu qualsevol d'aquests, Nextflow tornarà a executar correctament les tasques afectades.

### Conclusió

Sabeu com utilitzar `-resume` per tornar a executar pipelines de manera eficient sense repetir tasques reeixides.

### Què segueix?

Ara que podeu executar nf-core/molkart amb dades de prova, esteu preparats per aprendre com configurar-lo per als vostres propis conjunts de dades.
