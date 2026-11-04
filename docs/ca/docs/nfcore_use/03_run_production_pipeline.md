# Part 3: Executar un pipeline de producció

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la [Part 2](./02_configure_execution.md), heu après a definir paràmetres i personalitzar la configuració per a nf-core/demo.
Ara apliquem el que heu après a un pipeline de producció real, nf-core/rnaseq.

---

## 1. Descarregar i executar nf-core/rnaseq

Fins ara hem utilitzat `nf-core/demo`, que és un pipeline mínim dissenyat per a la formació.
Ara descarreguem un pipeline de producció real i l'executem amb el seu perfil de prova.

El pipeline `nf-core/rnaseq` realitza els passos principals de l'anàlisi de seqüenciació d'RNA en massa: control de qualitat, retallada d'adaptadors, alineament de lectures i quantificació a nivell de gens.
Probablement és el pipeline nf-core més utilitzat fins avui.

### 1.1. Descarregar el pipeline

Executeu la comanda següent per descarregar-lo.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Sortida de la comanda"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

El pipeline ara està emmagatzemat localment a la memòria cau i llest per executar-se.

### 1.2. Executar el perfil de prova

Executeu-lo amb el perfil de prova i Docker:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

La línia clau d'aquest error és:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

La màquina Codespaces per defecte té 8 GB de RAM, que és també el valor per defecte habitual per a Docker Desktop.
El pipeline sol·licita 12 GB per al procés `FQ_LINT` — més del que la màquina pot proporcionar.

Aquests 12 GB provenen de l'etiqueta de recursos `process_low` definida a `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Una opció seria utilitzar un tipus de màquina més gran, però per a propòsits de prova volem poder executar-lo amb qualsevol maquinari disponible.
El millor enfocament és sobreescriure els valors de recursos per defecte en un fitxer de configuració personalitzat.

### 1.3. Tornar a executar amb una configuració personalitzada

Us proporcionem un fitxer de configuració personalitzat que sobreescriu els valors de recursos per defecte basats en etiquetes.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

La [Part 2](./02_configure_execution.md) va introduir `withName:` per apuntar a un únic procés pel seu nom.
Aquí utilitzem `withLabel:` per apuntar a tots els processos que comparteixen una etiqueta alhora.

Aquest fitxer ja es troba al vostre directori de treball.
Passeu-lo amb `-c` per aplicar les sobreescriptures:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Sortida de la comanda (pipeline iniciant-se)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

El pipeline ara s'està executant i podeu veure com les tasques es completen una per una.
Amb aquest conjunt de dades de prova mínim, es completarà en 15–20 minuts, executant més de 200 tasques en total.

Els experiments reals de RNA-seq solen implicar desenes de mostres i s'executen durant hores o dies.
Nextflow és compatible amb planificadors HPC (SLURM, PBS, LSF) i plataformes al núvol (AWS, Google Cloud, Azure), que poden reduir dràsticament el temps d'execució distribuint la feina entre molts nodes.
Tanmateix, configurar aquests entorns afegeix una complexitat significativa.

La plataforma Seqera (desenvolupada pels creadors de Nextflow) proporciona una interfície web per llançar pipelines de Nextflow en infraestructura HPC o al núvol (ja sigui la vostra pròpia o una gestionada per a vosaltres), amb capacitats de gestió de càlcul i dades que simplifiquen el procés d'executar pipelines a gran escala.

!!! tip "Consell"

    Els investigadors acadèmics poden accedir a la plataforma Seqera de forma gratuïta a través del [programa acadèmic de Seqera](https://seqera.io/academic-program/).

### Conclusió

Heu descarregat `nf-core/rnaseq`, heu vist com funcionen les etiquetes de recursos d'nf-core i heu après a sobreescriure-les amb un fitxer de configuració personalitzat.
Més important encara, heu vist per què l'execució local és un punt de partida i no una destinació per a anàlisis a escala real.

### Què segueix?

Heu cobert els fonaments de l'execució de pipelines nf-core.
Consulteu els [Passos següents](next_steps.md) per saber cap on continuar.

---

## Resum

En aquesta part heu après a:

- Descarregar i executar un pipeline a escala de producció (nf-core/rnaseq) i sobreescriure les seves etiquetes de recursos per defecte
