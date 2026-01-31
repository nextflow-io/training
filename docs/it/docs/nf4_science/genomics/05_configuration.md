# Parte 3: Profilazione e ottimizzazione delle risorse

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

QUESTO È UN SEGNAPOSTO

!!!note "Nota"

    Questo modulo di formazione è in fase di aggiornamento.

---

TODO

### 1.1. Eseguire il workflow per generare un report di utilizzo delle risorse

Per fare in modo che Nextflow generi il report automaticamente, è sufficiente aggiungere `-with-report <filename>.html` alla riga di comando.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

Il report è un file html, che potete scaricare e aprire nel suo browser. Può anche fare clic con il tasto destro su di esso nell'esploratore file a sinistra e fare clic su `Show preview` per visualizzarlo in VS Code.

Si prenda alcuni minuti per esaminare il report e verificare se riesce a identificare alcune opportunità per regolare le risorse.
Assicuratevi di fare clic sulle schede che mostrano i risultati di utilizzo come percentuale di ciò che è stato allocato.
È disponibile una [documentazione](https://www.nextflow.io/docs/latest/reports.html) che descrive tutte le funzionalità disponibili.

<!-- TODO: insert images -->

Un'osservazione è che `GATK_JOINTGENOTYPING` sembra richiedere molte risorse CPU, il che ha senso poiché esegue molti calcoli complessi.
Potremmo quindi provare ad aumentare questa risorsa e verificare se riduce i tempi di esecuzione.

Tuttavia, sembra che abbiamo sovrastimato le allocazioni di memoria; tutti i processi stanno utilizzando solo una frazione di ciò che stiamo fornendo loro.
Dovremmo ridurre questa allocazione e risparmiare alcune risorse.

### 1.2. Regolare le allocazioni di risorse per un processo specifico

Possiamo specificare allocazioni di risorse per un determinato processo utilizzando il selettore di processo `withName`.
La sintassi appare così quando è da solo in un blocco process:

```groovy title="Sintassi"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Aggiungiamo questo al blocco process esistente nel file `nextflow.config`.

```groovy title="nextflow.config" linenums="11"
process {
    // impostazioni predefinite per tutti i processi
    cpus = 2
    memory = 2.GB
    // allocazioni per un processo specifico
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Con questa specifica, le impostazioni predefinite si applicheranno a tutti i processi **eccetto** il processo `GATK_JOINTGENOTYPING`, che è un caso particolare che riceve molte più CPU.
Speriamo che questo abbia un effetto.

### 1.3. Eseguire nuovamente con la configurazione modificata

Eseguiamo il workflow nuovamente con la configurazione modificata e con il flag di reporting attivato, ma notiamo che stiamo dando al report un nome diverso in modo da poterli differenziare.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Ancora una volta, probabilmente non noterà una differenza sostanziale nei tempi di esecuzione, perché questo è un carico di lavoro così ridotto e gli strumenti trascorrono più tempo in attività accessorie che nell'eseguire il lavoro 'reale'.

Tuttavia, il secondo report mostra che il nostro utilizzo delle risorse è ora più equilibrato.

<!-- **TODO: screenshots?** -->

Come potete vedere, questo approccio è utile quando i vostri processi hanno requisiti di risorse diversi. Le permette di dimensionare correttamente le allocazioni di risorse che configura per ciascun processo in base a dati effettivi, non a ipotesi.

!!!note "Nota"

    Questo è solo un piccolo assaggio di ciò che potete fare per ottimizzare l'utilizzo delle risorse.
    Nextflow stesso ha una [logica di retry dinamica](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) davvero interessante integrata per riprovare le attività che falliscono a causa di limitazioni di risorse.
    Inoltre, Seqera Platform offre strumenti basati su AI per ottimizzare le allocazioni di risorse automaticamente.

    Tratteremo entrambi questi approcci in una prossima parte di questo corso di formazione.

Detto questo, potrebbero esserci alcuni vincoli su ciò che può (o deve) allocare a seconda dell'executor di calcolo e dell'infrastruttura di calcolo che state utilizzando. Ad esempio, il vostro cluster potrebbe richiedere di rimanere entro determinati limiti che non si applicano quando si esegue altrove.
