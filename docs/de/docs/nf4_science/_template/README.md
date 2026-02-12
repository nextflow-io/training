# nf4_science Kursvorlage

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dies ist eine generische Vorlage zum Erstellen eines neuen domänenspezifischen Kurses unter `nf4_science/`.
Sie basiert auf den gängigen Mustern der Genomics- und RNAseq-Kurse.

## Verwendung

1. Kopiere dieses Verzeichnis und das entsprechende Skriptverzeichnis, um deinen Kurs zu erstellen:

   ```bash
   # Dokumentation
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Skripte
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Benenne das Workflow-Skript um:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Benenne die Moduldateien passend zu deinen Tools um:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Suche nach allen `{PLACEHOLDER}`-Werten in der Dokumentation und den Skripten und ersetze sie durch deine Inhalte.

5. Füge die Navigation zu `mkdocs.yml` hinzu (siehe Snippet unten).

6. Fülle das `data/`-Verzeichnis mit Testdatensätzen.

7. Erstelle das `solutions/`-Verzeichnis mit funktionierendem Code für jeden Teil.

## mkdocs.yml nav Snippet

```yaml
- { Your Domain }:
    - nf4_science/{your_domain}/index.md
    - nf4_science/{your_domain}/00_orientation.md
    - nf4_science/{your_domain}/01_method.md
    - nf4_science/{your_domain}/02_single_sample.md
    - nf4_science/{your_domain}/03_multi_sample.md
    - nf4_science/{your_domain}/survey.md
    - nf4_science/{your_domain}/next_steps.md
```

## Platzhalter-Referenz

### Platzhalter auf Kursebene

| Platzhalter                  | Beschreibung                          | Beispiel (Genomics)                        |
| ---------------------------- | ------------------------------------- | ------------------------------------------ |
| `{DOMAIN}`                   | Domänenname (Großschreibung)          | Genomics                                   |
| `{DOMAIN_DIR}`               | Verzeichnisname (Kleinschreibung)     | genomics                                   |
| `{METHOD}`                   | Name der Analysemethode               | variant calling                            |
| `{METHOD_SHORT_DESCRIPTION}` | Einzeilige Methodenbeschreibung       | variant calling with GATK                  |
| `{ACCESSORY_FILES}`          | Arten der verwendeten Zusatzdateien   | index files and reference genome resources |

### Tool-Platzhalter

| Platzhalter                           | Beschreibung                       | Beispiel (Genomics)                                 |
| ------------------------------------- | ---------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Anzeigenamen der Tools             | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | URLs zur Tool-Dokumentation        | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | Vollständige Container-URI         | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Moduldateinamen (ohne .nf)         | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | GROSSGESCHRIEBENER Prozessname     | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | Der auszuführende Shell-Befehl     | samtools index '<input_bam>'                        |

### Eingabe/Ausgabe-Platzhalter

| Platzhalter            | Beschreibung                       | Beispiel (Genomics)  |
| ---------------------- | ---------------------------------- | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Primärer Eingabedateityp           | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nextflow-Parametername             | reads_bam            |
| `{TEST_INPUT_PATH}`    | Pfad zur Testeingabe in data/      | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Finaler Pipeline-Ausgabetyp        | joint-called VCFs    |

### Inhalts-Platzhalter

| Platzhalter             | Beschreibung                                    |
| ----------------------- | ----------------------------------------------- |
| `{PART2_SUMMARY}`       | Einzeilige Zusammenfassung des Umfangs von Teil 2 |
| `{AGGREGATION_SUMMARY}` | Einzeilige Zusammenfassung des Aggregationsschritts |
| `{PROCESS_LIST}`        | Kommagetrennte Liste der Prozessnamen           |
| `{TYPEFORM_ID}`         | Typeform-Einbettungs-ID für die Umfrage         |

## Pädagogische Struktur

Die Vorlage folgt einem dreiteiligen Aufbau:

1. **Teil 1 (01_method.md)**: Manuelles Testen in Docker-Containern zum Verständnis der Methodik
2. **Teil 2 (02_single_sample.md)**: Befehle in Nextflow einbinden; einzelne Probe, dann Batch
3. **Teil 3 (03_multi_sample.md)**: Multi-Proben-Aggregation mit Kanal-Operatoren hinzufügen

### Wichtige Konventionen

- Verwende **Vorher/Nachher-Codeblöcke mit Tabs** und `hl_lines`, um Codeänderungen zu zeigen
- Verwende `??? success "Command output"` für einklappbare erwartete Ausgaben
- Verwende `??? abstract "Directory contents"` für einklappbare Verzeichnisbäume
- Beende jeden Hauptabschnitt mit **Fazit** und **Wie geht es weiter?** Unterabschnitten
- Verwende `---` horizontale Linien zur Trennung der nummerierten Hauptabschnitte
- Stelle Gerüstdateien (Workflow + Module) bereit, die Lernende ausfüllen können
- Organisiere Lösungen nach Teilen (`solutions/part2/`, `solutions/part3/`)

## Verzeichnisstruktur

```
docs/en/docs/nf4_science/{domain}/
├── index.md                    # Course overview with frontmatter
├── 00_orientation.md           # Environment setup
├── 01_method.md                # Manual testing in containers
├── 02_single_sample.md         # Single-sample Nextflow implementation
├── 03_multi_sample.md          # Multi-sample aggregation
├── survey.md                   # Typeform feedback survey
├── next_steps.md               # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
