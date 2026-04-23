# Zusammenfassung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Du hast das Plugin-Entwicklungs-Training abgeschlossen.
Diese Seite fasst zusammen, was du in jedem Teil gebaut hast, erklärt die Verteilung und gibt dir Hinweise, wie es weitergehen kann.

---

## Was du gelernt hast

### Teil 1: Plugins verwenden

Du hast gelernt, wie Nextflow-Plugins aus der Perspektive von Nutzer\*innen funktionieren.
Du hast nf-schema und nf-co2footprint installiert, sie über `nextflow.config` konfiguriert und gesehen, wie Plugins Eingaben validieren, Funktionen hinzufügen und sich in Pipeline-Lifecycle-Events einklinken können.

### Teil 2: Einrichten

Du hast eine Plugin-Entwicklungsumgebung mit Java 21+ eingerichtet, ein neues Plugin-Projekt mit dem Befehl `nextflow plugin create` erstellt und die Projektstruktur kennengelernt, die Nextflow erwartet: Quelldateien, Build-Konfiguration und den Makefile-Workflow.

### Teil 3: Eigene Funktionen

Du hast deinen ersten Erweiterungspunkt implementiert, indem du mit `@Function` annotierte Methoden in einer `PluginExtensionPoint`-Klasse erstellt hast.
Du hast `reverseGreeting` und `decorateGreeting` gebaut und sie dann aus einem Pipeline-Skript importiert und aufgerufen.

### Teil 4: Testen

Du hast Unit-Tests für deine eigenen Funktionen mit dem Groovy-Test-Framework geschrieben.
Du hast gelernt, wie du Tests mit `make test` ausführst und überprüfst, dass dein Plugin korrekt funktioniert, bevor du es installierst.

### Teil 5: Observer

Du hast das `TraceObserver`-Interface implementiert, um dich in Pipeline-Lifecycle-Events einzuklinken.
Du hast `GreetingObserver` (reagiert auf Pipeline-Start und -Abschluss) und `TaskCounterObserver` (zählt abgeschlossene Aufgaben) gebaut und sie dann über eine `TraceObserverFactory` registriert.

### Teil 6: Konfiguration

Du hast dein Plugin über `nextflow.config` konfigurierbar gemacht, indem du `session.config.navigate()` verwendest, um Werte zur Laufzeit zu lesen.
Du hast eine `@ConfigScope`-Klasse hinzugefügt, um die Optionen deines Plugins formal zu deklarieren – das beseitigt die Warnungen „Unrecognized config option" und ermöglicht IDE-Unterstützung.

---

## Verteilung

Sobald dein Plugin lokal funktioniert, kannst du es über die Nextflow-Plugin-Registry mit anderen teilen.

### Versionierung

Verwende [semantische Versionierung](https://semver.org/) für deine Releases:

| Versionsänderung          | Wann verwenden                   | Beispiel                                           |
| ------------------------- | -------------------------------- | -------------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Breaking Changes                 | Eine Funktion entfernen, Rückgabetypen ändern      |
| **MINOR** (1.0.0 → 1.1.0) | Neue Features, abwärtskompatibel | Eine neue Funktion hinzufügen                      |
| **PATCH** (1.0.0 → 1.0.1) | Bugfixes, abwärtskompatibel      | Einen Fehler in einer bestehenden Funktion beheben |

Aktualisiere die Version in `build.gradle` vor jedem Release:

```groovy title="build.gradle"
version = '1.0.0'  // Semantische Versionierung verwenden: MAJOR.MINOR.PATCH
```

### In der Registry veröffentlichen

Die [Nextflow-Plugin-Registry](https://registry.nextflow.io/) ist der offizielle Weg, Plugins mit der Community zu teilen.

Der Veröffentlichungs-Workflow:

1. Deinen Plugin-Namen in der [Registry](https://registry.nextflow.io/) reservieren (mit deinem GitHub-Konto anmelden)
2. Deine API-Zugangsdaten in `~/.gradle/gradle.properties` konfigurieren
3. Tests ausführen, um sicherzustellen, dass alles funktioniert: `make test`
4. Mit `make release` veröffentlichen

Schritt-für-Schritt-Anleitungen findest du in der [offiziellen Veröffentlichungsdokumentation](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Nach der Veröffentlichung können Nutzer\*innen dein Plugin ohne lokale Einrichtung installieren:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow lädt das Plugin beim ersten Verwenden automatisch aus der Registry herunter.

---

## Checkliste für die Plugin-Entwicklung

- [ ] Java 21+ installiert
- [ ] Projekt mit `nextflow plugin create <name> <org>` erstellen
- [ ] Erweiterungsklasse mit `@Function`-Methoden implementieren
- [ ] Unit-Tests schreiben und mit `make test` ausführen
- [ ] Mit `make install` bauen und installieren
- [ ] Optional `TraceObserver`-Implementierungen für Workflow-Events hinzufügen
- [ ] Optional `ConfigScope` für die Plugin-Konfiguration hinzufügen
- [ ] In `nextflow.config` mit `plugins { id 'plugin-id' }` aktivieren
- [ ] Funktionen mit `include { fn } from 'plugin/plugin-id'` importieren
- [ ] Versionieren und in der Registry veröffentlichen

---

## Wichtige Code-Muster

**Funktionsdefinition:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Plugin-Konfiguration:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Verwendung in Workflows:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Übersicht der Erweiterungspunkte

| Typ                   | Klasse/Annotation | Zweck                                              |
| --------------------- | ----------------- | -------------------------------------------------- |
| Funktion              | `@Function`       | Aus Workflows aufrufbar                            |
| Trace Observer        | `TraceObserver`   | In Workflow-Lifecycle-Events einklinken            |
| Konfigurationsbereich | `@ScopeName`      | Plugin-Konfiguration in nextflow.config definieren |

---

## Wie geht es weiter?

Hier sind einige praktische nächste Schritte, um deine Plugin-Entwicklung fortzusetzen.

**Baue etwas Echtes.**
Wähle einen Anwendungsfall aus deiner eigenen Arbeit: eine eigene Funktion, die dein Team immer wieder verwendet, einen Observer, der bei Pipeline-Abschluss Slack-Benachrichtigungen sendet, oder einen Config-Bereich, der Optionen für die Pipelines deiner Organisation standardisiert.
Von einem echten Problem auszugehen ist der schnellste Weg, dein Verständnis zu vertiefen.

**Verwende nf-hello als Referenz.**
Das [nf-hello](https://github.com/nextflow-io/nf-hello)-Repository ist das offizielle minimale Plugin-Beispiel.
Es ist ein guter Ausgangspunkt für neue Projekte und eine nützliche Referenz, wenn du nachschauen möchtest, wie etwas strukturiert ist.

**Lies die offizielle Dokumentation.**
Die Nextflow-Dokumentation behandelt Themen, die über dieses Training hinausgehen, darunter Channel-Factories, Operator-Overloading und fortgeschrittene Observer-Muster.
Der Leitfaden [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) ist die umfassendste Referenz.

**Studiere bestehende Plugins.**
Das [Nextflow-Plugins-Repository](https://github.com/nextflow-io/plugins) enthält den Quellcode für offizielle Plugins wie nf-schema, nf-wave und nf-tower.
Produktiven Plugin-Code zu lesen ist eine der besten Möglichkeiten, Muster und Konventionen zu lernen, die über einführende Beispiele hinausgehen.

---

## Weitere Ressourcen

**Offizielle Dokumentation:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): umfassender Leitfaden zur Installation und Konfiguration von Plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): detaillierte Referenz zur Plugin-Entwicklung
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): Konfigurationsbereiche für Plugins erstellen

**Plugins entdecken:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): verfügbare Plugins durchsuchen und entdecken
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): Dokumentation der Registry

**Beispiele und Referenzen:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): einfaches Beispiel-Plugin (guter Ausgangspunkt)
- [Nextflow-Plugins-Repository](https://github.com/nextflow-io/plugins): Sammlung offizieller Plugins als Referenz
