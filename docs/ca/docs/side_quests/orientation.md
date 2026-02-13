# Orientació

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'entorn de GitHub Codespaces conté tot el programari, codi i dades necessaris per treballar en aquest curs de formació, així que no cal que instal·lis res tu mateix.
No obstant això, necessites un compte (gratuït) per iniciar sessió, i hauries de dedicar uns minuts a familiaritzar-te amb la interfície.

Si encara no ho has fet, si us plau segueix [aquest enllaç](../../envsetup/) abans de continuar.

## Materials proporcionats

Al llarg d'aquest curs de formació, treballarem al directori `side-quests/`.
Aquest directori conté tots els fitxers de codi, dades de prova i fitxers accessoris que necessitaràs.

Sent lliure d'explorar els continguts d'aquest directori; la manera més fàcil de fer-ho és utilitzar l'explorador de fitxers a la part esquerra de l'espai de treball de GitHub Codespaces.
Alternativament, pots utilitzar la comanda `tree`.
Al llarg del curs, utilitzem la sortida de `tree` per representar l'estructura i els continguts del directori de forma llegible, de vegades amb modificacions menors per claredat.

Aquí generem una taula de continguts fins al segon nivell:

```bash
tree . -L 2
```

Si executes això dins de `side-quests`, hauries de veure la següent sortida:

```console title="Directory contents"
.
├── metadata
├── nf-core
├── nf-test
├── solutions
├── splitting_and_grouping
└── workflows_of_workflows
```

**Aquí tens un resum del que hauries de saber per començar:**

- **Cada directori correspon a una missió secundària individual.**
  Els seus continguts es detallen a la pàgina de la missió secundària corresponent.

- **El directori `solutions`** conté els scripts de workflow i/o mòdul completats que resulten d'executar els diversos passos de cada missió secundària.
  Estan pensats per ser utilitzats com a referència per comprovar el teu treball i resoldre qualsevol problema.

!!!tip "Consell"

    Si per qualsevol raó surts d'aquest directori, sempre pots executar aquesta comanda per tornar-hi:

    ```bash
    cd /workspaces/training/side-quests
    ```

Ara, per començar el curs, fes clic a la fletxa a la cantonada inferior dreta d'aquesta pàgina.
