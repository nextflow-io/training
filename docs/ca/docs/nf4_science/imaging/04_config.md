# Part 4: Configuració

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A les Parts 1-3, hem après com executar Nextflow, executar un pipeline d'nf-core i gestionar entrades amb fitxers de paràmetres i fulls de mostres.
Ara explorarem com configurar pipelines per a diferents entorns informàtics utilitzant **fitxers de configuració** i **perfils**.

## Objectius d'aprenentatge

Al final d'aquesta part, podràs:

- Entendre com Nextflow resol la configuració de múltiples fonts
- Utilitzar perfils integrats d'nf-core per a contenidors i proves
- Crear perfils personalitzats per a diferents entorns informàtics
- Personalitzar sol·licituds de recursos utilitzant etiquetes de procés
- Gestionar límits de recursos en entorns restringits
- Inspeccionar la configuració resolta amb `nextflow config`

---

## 1. Entendre la configuració de Nextflow

### 1.1. Què és un fitxer de configuració?

Nextflow utilitza fitxers de configuració per separar la **lògica del workflow** (què fer) de la **configuració d'execució** (com i on fer-ho).

Els fitxers de configuració controlen:

- Motors de contenidors (Docker, Singularity, Conda)
- Recursos informàtics (CPUs, memòria, temps)
- Plataformes d'execució (local, HPC, núvol)
- Paràmetres del pipeline

### 1.2. Precedència de configuració

Nextflow carrega la configuració de múltiples fonts, amb les fonts posteriors sobreescrivint les anteriors:

1. **Configuració del pipeline**: `nextflow.config` al repositori del pipeline
2. **Configuració del directori**: `nextflow.config` al teu directori de treball actual
3. **Configuració d'usuari**: `~/.nextflow/config`
4. **Línia de comandes**: Paràmetres i opcions passats directament

Aquest enfocament per capes et permet mantenir valors per defecte al pipeline, sobreescriure amb configuracions específiques d'usuari i fer ajustos ràpids a la línia de comandes.

### 1.3. La nostra configuració actual

Vegem la configuració que hem estat utilitzant:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Comentem o canviem la línia `docker.enabled = true` de la Part 2, i descobrim com podem aconseguir el mateix resultat utilitzant un perfil a molkart.

---

## 2. Utilitzar perfils

### 2.1. Què són els perfils?

Els perfils són conjunts de configuració amb nom que es poden activar amb la bandera `-profile` mitjançant la comanda `nextflow run`.
Faciliten el canvi entre diferents escenaris informàtics sense editar fitxers de configuració.

Tots els pipelines d'nf-core vénen amb diversos perfils per defecte que podem utilitzar.

### 2.2. Inspeccionar perfils integrats

Inspeccionem-los al fitxer `molkart/nextflow.config` associat amb la base de codi del pipeline:

```bash
code molkart/nextflow.config
```

Cerca el bloc `profiles`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Perfils de contenidors comuns:

- `docker`: Utilitza contenidors Docker (més comú per al desenvolupament local)
- `singularity`: Utilitza Singularity/Apptainer (comú en HPC)
- `conda`: Utilitza entorns Conda
- `apptainer`: Utilitza contenidors Apptainer

### 2.3. Re-executar amb perfils en lloc de nextflow.config

Ara que hem desactivat la configuració de docker al nostre fitxer local `nextflow.config` i entenem els perfils, re-executem el pipeline utilitzant la bandera `-profile`.

Anteriorment a la Part 3, vam crear un fitxer `params.yaml` amb els nostres paràmetres personalitzats.
Ara podem combinar-lo amb el perfil Docker integrat:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Desglossem què fa cada bandera:

- `-profile docker`: Activa el perfil Docker del `nextflow.config` de molkart, que estableix `docker.enabled = true`
- `-params-file params.yaml`: Carrega tots els paràmetres del pipeline del nostre fitxer YAML
- `-resume`: Reutilitza resultats emmagatzemats en memòria cau d'execucions anteriors

Com que estem utilitzant `-resume`, Nextflow comprovarà si alguna cosa ha canviat des de l'última execució.
Si els paràmetres, entrades i codi són els mateixos, totes les tasques es recuperaran de la memòria cau i el pipeline es completarà gairebé instantàniament.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Observa que tots els processos mostren `cached: 2` o `cached: 1` - no s'ha re-executat res!

### 2.4. Perfils de prova

Els perfils de prova proporcionen maneres ràpides d'especificar paràmetres d'entrada per defecte i fitxers de dades per permetre't verificar que el pipeline funciona.
Els pipelines d'nf-core sempre inclouran almenys dos perfils de prova:

- `test`: Conjunt de dades petit amb paràmetres ràpids per a proves ràpides
- `test_full`: Prova més completa amb dades més grans

Vegem més de prop el perfil `test` a molkart que s'inclou utilitzant la directiva `includeConfig`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Això significa que sempre que executem el pipeline amb `-profile test`, Nextflow carregarà la configuració de `conf/test.config`.

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Observa que aquest perfil conté els mateixos paràmetres que vam utilitzar al nostre fitxer `params.yaml` anteriorment.

Pots activar múltiples perfils separant-los amb comes.
Utilitzem això per provar el nostre pipeline sense necessitar el nostre fitxer de paràmetres:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Això combina:

- `docker`: Activa contenidors Docker
- `test`: Utilitza conjunt de dades i paràmetres de prova

Els perfils s'apliquen d'esquerra a dreta, així que els perfils posteriors sobreescriuen els anteriors si estableixen els mateixos valors.

### Conclusió

Els pipelines d'nf-core vénen amb perfils integrats per a contenidors, proves i entorns especials.
Pots combinar múltiples perfils per construir la configuració que necessites.

### Què segueix?

Aprèn com crear els teus propis perfils personalitzats per a diferents entorns informàtics.

---

## 3. Crear perfils personalitzats

### 3.1. Crear perfils per canviar entre desenvolupament local i execució en HPC

Creem perfils personalitzats per a dos escenaris:

1. Desenvolupament local amb Docker
2. HPC universitari amb planificador Slurm i Singularity

Afegeix el següent al teu `nextflow.config`:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Ara pots canviar entre entorns fàcilment:

```bash
# Per al desenvolupament local
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Per a HPC (quan estigui disponible)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note "Nota"

    No podem provar el perfil HPC en aquest entorn de formació ja que no tenim accés a un planificador Slurm.
    Però això mostra com ho configuraries per a un ús real.

### 3.2. Utilitzar `nextflow config` per inspeccionar la configuració

La comanda `nextflow config` mostra la configuració completament resolta sense executar el pipeline.

Visualitza la configuració per defecte:

```bash
nextflow config ./molkart
```

Visualitza la configuració amb un perfil específic:

```bash
nextflow config -profile local_dev ./molkart
```

Això és extremadament útil per a:

- Depurar problemes de configuració
- Entendre quins valors s'utilitzaran realment
- Comprovar com interactuen múltiples perfils

### Conclusió

Els perfils personalitzats et permeten canviar entre diferents entorns informàtics amb una sola bandera de línia de comandes.
Utilitza `nextflow config` per inspeccionar la configuració resolta abans d'executar.

### Què segueix?

Aprèn com personalitzar sol·licituds de recursos per a processos individuals utilitzant el sistema d'etiquetes de procés d'nf-core.

---

## 4. Personalitzar sol·licituds de recursos

### 4.1. Entendre les etiquetes de procés als pipelines d'nf-core

Per simplicitat, els pipelines d'nf-core utilitzen [**etiquetes de procés**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) per estandarditzar l'assignació de recursos a tots els pipelines.
Cada procés està etiquetat amb una etiqueta com `process_low`, `process_medium` o `process_high` per descriure requisits de recursos informàtics baixos, mitjans o alts, respectivament.
Aquestes etiquetes es converteixen en sol·licituds de recursos específiques en un dels fitxers de configuració ubicats al directori `conf/` del pipeline.

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Observa el multiplicador `task.attempt` - això permet que els reintents de tasques posteriors sol·licitin més recursos, si el pipeline està configurat amb `process.maxRetries > 1`.

### 4.2. Sobreescriure recursos per a processos específics

Per a un control detallat, dirigeix-te a processos individuals per nom:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Si intentem executar aquest pipeline amb la sobreescriptura anterior, el procés `CELLPOSE` sol·licitarà 16 CPUs i 32 GB de memòria en lloc del valor per defecte definit per la seva etiqueta.
Això farà que el pipeline falli al nostre entorn actual ja que no tenim tanta RAM disponible.
Aprendrem com prevenir aquests tipus de fallades a la següent secció.

!!! Tip "Consell"

    Per trobar noms de processos, consulta la sortida d'execució del pipeline o comprova `.nextflow.log`.
    Els noms de processos segueixen el patró `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Conclusió

Els pipelines d'nf-core utilitzen etiquetes de procés per estandarditzar l'assignació de recursos.
Pots sobreescriure recursos per etiqueta (afecta múltiples processos) o per nom (afecta un procés específic).

### Què segueix?

Aprèn com gestionar límits de recursos en entorns restringits com GitHub Codespaces.

---

## 5. Gestionar recursos en entorns restringits

### 5.1. El problema dels límits de recursos

Si intentéssim executar molkart amb un procés sol·licitant 16 CPUs i 32 GB de memòria (com es mostra a la secció 4.2), fallaria al nostre entorn actual perquè no tenim tants recursos disponibles.
En un entorn de clúster amb nodes més grans, aquestes sol·licituds es presentarien al planificador.

En entorns restringits com GitHub Codespaces, sense límits, Nextflow es negaria a executar processos que superin els recursos disponibles.

### 5.2. Establir límits de recursos

La directiva `resourceLimits` limita les sol·licituds de recursos a valors especificats:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Això indica a Nextflow: "Si algun procés sol·licita més de 2 CPUs o 7 GB de memòria, limita'l a aquests límits."

### 5.3. Afegir límits de recursos als perfils personalitzats

Actualitza els teus perfils personalitzats per incloure límits apropiats:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning "Advertència"

    Establir límits de recursos massa baixos pot fer que els processos fallin o s'executin lentament.
    El pipeline pot necessitar utilitzar algoritmes menys intensius en memòria o processar dades en fragments més petits.

### Conclusió

Utilitza `resourceLimits` per executar pipelines en entorns amb recursos restringits limitant les sol·licituds de recursos dels processos.
Diferents perfils poden tenir diferents límits apropiats per al seu entorn.

### Què segueix?

Has completat la formació bàsica de Nextflow per a Bioimatge!

---

## Conclusió

Ara entens com configurar pipelines de Nextflow per a diferents entorns informàtics.

Habilitats clau que has après:

- **Precedència de configuració**: Com Nextflow resol configuracions de múltiples fonts
- **Perfils d'nf-core**: Utilitzar perfils integrats per a contenidors, proves i utilitats
- **Perfils personalitzats**: Crear els teus propis perfils per a diferents entorns
- **Etiquetes de procés**: Entendre i sobreescriure sol·licituds de recursos per etiqueta
- **Límits de recursos**: Gestionar entorns restringits amb `resourceLimits`
- **Inspecció de configuració**: Utilitzar `nextflow config` per depurar i verificar configuracions

Aquestes habilitats de configuració són transferibles a qualsevol pipeline de Nextflow i t'ajudaran a executar workflows eficientment a màquines locals, clústers HPC i plataformes de núvol.

### Què segueix?

Felicitats per completar el curs de Nextflow per a Bioimatge!

Pròxims passos:

- Omple l'enquesta del curs per proporcionar comentaris
- Consulta [Hello Nextflow](../hello_nextflow/index.md) per aprendre més sobre el desenvolupament de workflows
- Explora [Hello nf-core](../hello_nf-core/index.md) per aprofundir en les eines d'nf-core
- Navega per altres cursos a les [col·leccions de formació](../training_collections/index.md)
