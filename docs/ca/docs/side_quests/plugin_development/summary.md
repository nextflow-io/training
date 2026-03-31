# Resum

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Heu completat la formació de desenvolupament de plugins.
Aquesta pàgina recapitula el que heu construït a cada part, cobreix la distribució i proporciona orientació sobre on continuar.

---

## Què heu après

### Part 1: Ús de plugins

Heu descobert com funcionen els plugins de Nextflow des de la perspectiva d'un usuari.
Heu instal·lat nf-schema i nf-co2footprint, els heu configurat mitjançant `nextflow.config`, i heu vist com els plugins poden validar entrades, afegir funcions i connectar-se als esdeveniments del cicle de vida del pipeline.

### Part 2: Configuració de l'entorn

Heu configurat un entorn de desenvolupament de plugins amb Java 21+, heu creat un nou projecte de plugin amb la comanda `nextflow plugin create`, i heu après l'estructura de projecte que espera Nextflow: fitxers font, configuració de compilació i el workflow del Makefile.

### Part 3: Funcions personalitzades

Heu implementat el vostre primer punt d'extensió creant mètodes anotats amb `@Function` en una classe `PluginExtensionPoint`.
Heu construït `reverseGreeting` i `decorateGreeting`, i després els heu importat i cridat des d'un script de pipeline.

### Part 4: Proves

Heu escrit proves unitàries per a les vostres funcions personalitzades utilitzant el framework de proves de Groovy.
Heu après a executar les proves amb `make test` i a verificar que el vostre plugin es comporta correctament abans d'instal·lar-lo.

### Part 5: Observadors

Heu implementat la interfície `TraceObserver` per connectar-vos als esdeveniments del cicle de vida del pipeline.
Heu construït `GreetingObserver` (que reacciona a l'inici i la finalització del pipeline) i `TaskCounterObserver` (que compta les tasques completades), i després els heu registrat mitjançant una `TraceObserverFactory`.

### Part 6: Configuració

Heu fet el vostre plugin configurable a través de `nextflow.config` utilitzant `session.config.navigate()` per llegir valors en temps d'execució.
Heu afegit una classe `@ConfigScope` per declarar formalment les opcions del vostre plugin, eliminant els avisos "Unrecognized config option" i habilitant el suport de l'IDE.

---

## Distribució

Un cop el vostre plugin funcioni localment, podeu compartir-lo amb altres persones a través del registre de plugins de Nextflow.

### Versionat

Seguiu el [versionat semàntic](https://semver.org/) per als vostres llançaments:

| Canvi de versió           | Quan utilitzar-lo                | Exemple                                      |
| ------------------------- | -------------------------------- | -------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Canvis incompatibles             | Eliminar una funció, canviar tipus de retorn |
| **MINOR** (1.0.0 → 1.1.0) | Noves funcionalitats, compatible | Afegir una nova funció                       |
| **PATCH** (1.0.0 → 1.0.1) | Correccions d'errors, compatible | Corregir un error en una funció existent     |

Actualitzeu la versió a `build.gradle` abans de cada llançament:

```groovy title="build.gradle"
version = '1.0.0'  // Utilitzeu el versionat semàntic: MAJOR.MINOR.PATCH
```

### Publicació al registre

El [registre de plugins de Nextflow](https://registry.nextflow.io/) és la manera oficial de compartir plugins amb la comunitat.

El workflow de publicació:

1. Reclameu el nom del vostre plugin al [registre](https://registry.nextflow.io/) (inicieu sessió amb el vostre compte de GitHub)
2. Configureu les vostres credencials d'API a `~/.gradle/gradle.properties`
3. Executeu les proves per verificar que tot funciona: `make test`
4. Publiqueu amb `make release`

Per a instruccions pas a pas, consulteu la [documentació oficial de publicació](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Un cop publicat, els usuaris instal·len el vostre plugin sense cap configuració local:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow descarrega automàticament el plugin del registre en el primer ús.

---

## Llista de verificació del desenvolupament de plugins

- [ ] Java 21+ instal·lat
- [ ] Crear el projecte amb `nextflow plugin create <name> <org>`
- [ ] Implementar la classe d'extensió amb mètodes `@Function`
- [ ] Escriure proves unitàries i executar-les amb `make test`
- [ ] Compilar i instal·lar amb `make install`
- [ ] Opcionalment, afegir implementacions de `TraceObserver` per als esdeveniments del workflow
- [ ] Opcionalment, afegir `ConfigScope` per a la configuració del plugin
- [ ] Habilitar a `nextflow.config` amb `plugins { id 'plugin-id' }`
- [ ] Importar funcions amb `include { fn } from 'plugin/plugin-id'`
- [ ] Versionar i publicar al registre

---

## Patrons de codi clau

**Definició de funció:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Configuració del plugin:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Ús en workflows:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Resum dels punts d'extensió

| Tipus                 | Classe/Anotació | Propòsit                                                      |
| --------------------- | --------------- | ------------------------------------------------------------- |
| Funció                | `@Function`     | Invocable des de workflows                                    |
| Observador de traça   | `TraceObserver` | Connectar-se als esdeveniments del cicle de vida del workflow |
| Àmbit de configuració | `@ScopeName`    | Definir la configuració del plugin a nextflow.config          |

---

## Què fer a continuació

Aquí teniu alguns passos pràctics per continuar el vostre camí en el desenvolupament de plugins.

**Construïu alguna cosa real.**
Trieu un cas d'ús del vostre propi treball: una funció personalitzada que el vostre equip utilitza repetidament, un observador que envia notificacions a Slack quan el pipeline finalitza, o un àmbit de configuració que estandarditza les opcions als pipelines de la vostra organització.
Partir d'un problema real és la manera més ràpida d'aprofundir en la comprensió.

**Utilitzeu nf-hello com a referència.**
El repositori [nf-hello](https://github.com/nextflow-io/nf-hello) és l'exemple oficial de plugin mínim.
És un bon punt de partida per a nous projectes i una referència útil quan necessiteu comprovar com s'estructura alguna cosa.

**Llegiu la documentació oficial.**
La documentació de Nextflow cobreix temes més enllà d'aquesta formació, incloent-hi factories de canals, sobrecàrrega d'operadors i patrons avançats d'observadors.
La guia [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) és la referència més completa.

**Estudieu plugins existents.**
El [repositori de plugins de Nextflow](https://github.com/nextflow-io/plugins) conté el codi font de plugins oficials com nf-schema, nf-wave i nf-tower.
Llegir codi de plugins en producció és una de les millors maneres d'aprendre patrons i convencions que van més enllà dels exemples introductoris.

---

## Recursos addicionals

**Documentació oficial:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): guia completa per instal·lar i configurar plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): referència detallada per al desenvolupament de plugins
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): creació d'àmbits de configuració per a plugins

**Descoberta de plugins:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): exploreu i descobriu els plugins disponibles
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): documentació del registre

**Exemples i referències:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): plugin d'exemple senzill (excel·lent punt de partida)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): col·lecció de plugins oficials per a referència
