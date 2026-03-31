# Riepilogo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hai completato la formazione sullo sviluppo di plugin.
Questa pagina riepiloga ciò che hai costruito in ogni parte, tratta la distribuzione e fornisce indicazioni su dove andare dopo.

---

## Cosa hai imparato

### Parte 1: Usare i plugin

Hai scoperto come funzionano i plugin di Nextflow dal punto di vista dell'utente.
Hai installato nf-schema e nf-co2footprint, li hai configurati tramite `nextflow.config`, e hai visto come i plugin possono validare gli input, aggiungere funzioni e agganciarsi agli eventi del ciclo di vita della pipeline.

### Parte 2: Configurazione dell'ambiente

Hai configurato un ambiente di sviluppo per plugin con Java 21+, creato un nuovo progetto di plugin usando il comando `nextflow plugin create`, e appreso la struttura di progetto che Nextflow si aspetta: file sorgente, configurazione della build e il flusso di lavoro con il Makefile.

### Parte 3: Funzioni personalizzate

Hai implementato il tuo primo punto di estensione creando metodi annotati con `@Function` in una classe `PluginExtensionPoint`.
Hai costruito `reverseGreeting` e `decorateGreeting`, poi li hai importati e chiamati da uno script di pipeline.

### Parte 4: Testing

Hai scritto unit test per le tue funzioni personalizzate usando il framework di testing Groovy.
Hai imparato come eseguire i test con `make test` e verificare che il tuo plugin si comporti correttamente prima di installarlo.

### Parte 5: Observer

Hai implementato l'interfaccia `TraceObserver` per agganciarti agli eventi del ciclo di vita della pipeline.
Hai costruito `GreetingObserver` (che reagisce all'avvio e al completamento della pipeline) e `TaskCounterObserver` (che conta le attività completate), poi li hai registrati tramite una `TraceObserverFactory`.

### Parte 6: Configurazione

Hai reso il tuo plugin configurabile tramite `nextflow.config` usando `session.config.navigate()` per leggere i valori a runtime.
Hai aggiunto una classe `@ConfigScope` per dichiarare formalmente le opzioni del tuo plugin, eliminando gli avvisi "Unrecognized config option" e abilitando il supporto IDE.

---

## Distribuzione

Una volta che il tuo plugin funziona in locale, puoi condividerlo con altri tramite il registro dei plugin di Nextflow.

### Versionamento

Segui il [versionamento semantico](https://semver.org/) per le tue release:

| Cambio di versione        | Quando usarlo                   | Esempio                                               |
| ------------------------- | ------------------------------- | ----------------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Modifiche incompatibili         | Rimozione di una funzione, cambio dei tipi di ritorno |
| **MINOR** (1.0.0 → 1.1.0) | Nuove funzionalità, compatibili | Aggiunta di una nuova funzione                        |
| **PATCH** (1.0.0 → 1.0.1) | Correzioni di bug, compatibili  | Correzione di un bug in una funzione esistente        |

Aggiorna la versione in `build.gradle` prima di ogni release:

```groovy title="build.gradle"
version = '1.0.0'  // Usa il versionamento semantico: MAJOR.MINOR.PATCH
```

### Pubblicazione nel registro

Il [registro dei plugin di Nextflow](https://registry.nextflow.io/) è il modo ufficiale per condividere i plugin con la community.

Il flusso di lavoro per la pubblicazione:

1. Registra il nome del tuo plugin nel [registro](https://registry.nextflow.io/) (accedi con il tuo account GitHub)
2. Configura le credenziali API in `~/.gradle/gradle.properties`
3. Esegui i test per verificare che tutto funzioni: `make test`
4. Pubblica con `make release`

Per istruzioni passo dopo passo, consulta la [documentazione ufficiale sulla pubblicazione](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Una volta pubblicato, gli utenti installano il tuo plugin senza alcuna configurazione locale:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow scarica automaticamente il plugin dal registro al primo utilizzo.

---

## Checklist per lo sviluppo di plugin

- [ ] Java 21+ installato
- [ ] Creare il progetto con `nextflow plugin create <name> <org>`
- [ ] Implementare la classe di estensione con metodi `@Function`
- [ ] Scrivere unit test ed eseguirli con `make test`
- [ ] Compilare e installare con `make install`
- [ ] Aggiungere facoltativamente implementazioni di `TraceObserver` per gli eventi del flusso di lavoro
- [ ] Aggiungere facoltativamente `ConfigScope` per la configurazione del plugin
- [ ] Abilitare in `nextflow.config` con `plugins { id 'plugin-id' }`
- [ ] Importare le funzioni con `include { fn } from 'plugin/plugin-id'`
- [ ] Versionare e pubblicare nel registro

---

## Pattern di codice principali

**Definizione di una funzione:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Configurazione del plugin:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Utilizzo nei flussi di lavoro:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Riepilogo dei punti di estensione

| Tipo                    | Classe/Annotazione | Scopo                                                       |
| ----------------------- | ------------------ | ----------------------------------------------------------- |
| Funzione                | `@Function`        | Chiamabile dai flussi di lavoro                             |
| Trace Observer          | `TraceObserver`    | Aggancio agli eventi del ciclo di vita del flusso di lavoro |
| Scope di configurazione | `@ScopeName`       | Definire la configurazione del plugin in nextflow.config    |

---

## Cosa fare dopo

Ecco alcuni passi pratici per continuare il tuo percorso nello sviluppo di plugin.

**Costruisci qualcosa di reale.**
Scegli un caso d'uso dal tuo lavoro quotidiano: una funzione personalizzata che il tuo team usa ripetutamente, un observer che invia notifiche Slack al completamento della pipeline, o uno scope di configurazione che standardizza le opzioni nelle pipeline della tua organizzazione.
Partire da un problema reale è il modo più rapido per approfondire la comprensione.

**Usa nf-hello come riferimento.**
Il repository [nf-hello](https://github.com/nextflow-io/nf-hello) è l'esempio ufficiale minimale di plugin.
È un buon punto di partenza per nuovi progetti e un utile riferimento quando hai bisogno di verificare come è strutturato qualcosa.

**Leggi la documentazione ufficiale.**
La documentazione di Nextflow tratta argomenti che vanno oltre questa formazione, tra cui le fabbriche di canali, l'overloading degli operatori e i pattern avanzati di observer.
La guida [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) è il riferimento più completo.

**Studia i plugin esistenti.**
Il [repository dei plugin di Nextflow](https://github.com/nextflow-io/plugins) contiene il codice sorgente dei plugin ufficiali come nf-schema, nf-wave e nf-tower.
Leggere il codice di plugin in produzione è uno dei modi migliori per imparare pattern e convenzioni che vanno oltre gli esempi introduttivi.

---

## Risorse aggiuntive

**Documentazione ufficiale:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): guida completa all'installazione e alla configurazione dei plugin
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): riferimento dettagliato per lo sviluppo di plugin
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): creazione di scope di configurazione per i plugin

**Scoperta dei plugin:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): sfoglia e scopri i plugin disponibili
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): documentazione del registro

**Esempi e riferimenti:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): plugin di esempio semplice (ottimo punto di partenza)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): raccolta di plugin ufficiali come riferimento
