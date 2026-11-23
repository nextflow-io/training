# Partie 2 : Lancer des pipelines depuis la ligne de commande

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la [Partie 1](./01_run_with_seqera.md), vous avez lancé nf-core/rnaseq depuis l'interface web de Seqera.
Nous allons maintenant faire la même chose depuis la ligne de commande en utilisant le CLI `tw`, et ajouter un nouveau pipeline à votre workspace.

---

## 1. Lancer des pipelines depuis la ligne de commande

Dans la vue d'exécution, cliquez sur l'onglet **Command line**.
Vous verrez la commande `nextflow run` exacte que la Platform a construite et soumise en votre nom — le même type de commande que vous avez exécutée manuellement dans le cours Use nf-core.

La Platform ne remplace pas Nextflow ; elle l'orchestre.
Tout ce que vous pouvez faire via l'interface web, vous pouvez également le faire depuis un terminal en utilisant le CLI `tw`, l'outil en ligne de commande pour interagir avec l'API de la Platform.
Cela est utile pour automatiser des lancements depuis des scripts ou des pipelines CI/CD.

Nous allons le faire maintenant depuis le même codespace que vous avez utilisé pour les cours précédents.

### 1.1. Installer le CLI `tw`

Exécutez les commandes suivantes dans votre terminal Codespace pour télécharger et installer le binaire `tw` :

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Vérifiez l'installation :

```bash
tw --version
```

??? success "Sortie de la commande"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

Le CLI `tw` est installé et prêt à être configuré.

### 1.2. Obtenir un token d'accès

Le CLI `tw` s'authentifie auprès de Seqera à l'aide d'un token d'accès personnel.

1. Dans l'interface web de Seqera, cliquez sur votre avatar en haut à droite et sélectionnez **Your tokens**.
2. Cliquez sur **Add token**, donnez-lui un nom (par exemple `training`), et cliquez sur **Add**.
3. Copiez la valeur du token — elle ne sera affichée qu'une seule fois.
   Si vous ne la sauvegardez pas immédiatement, vous devrez en générer un autre.

### 1.3. Configurer le CLI

Pour plus de commodité, nous allons créer un fichier de configuration contenant le
token d'accès que vous venez de générer ainsi que l'identifiant du workspace.

Ouvrez le fichier `.seqera_config` dans ce répertoire dans l'éditeur et définissez les deux variables :

- **`TOWER_ACCESS_TOKEN`** : le token que vous avez généré à la section 1.2
- **`TOWER_WORKSPACE_ID`** : l'identifiant numérique de votre workspace (la colonne `ID` dans `tw workspaces list`, que vous exécutez à la section 1.4)

Une fois les valeurs renseignées, chargez la configuration :

```bash
source .seqera_config
```

Vérifiez la connexion :

```bash
tw info
```

??? success "Sortie de la commande"

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

Le CLI `tw` est maintenant authentifié et connecté à votre compte Seqera.
Exécutez `source .seqera_config` au début de chaque session Codespace pour recharger la configuration.

!!! tip "Astuce"

    Si votre workspace n'a pas d'environnement de calcul principal défini, vous pouvez ajouter `export TOWER_COMPUTE_ENV=<compute-env-name>` à votre fichier de configuration pour en définir un par défaut.
    Toute valeur de configuration peut être remplacée en ligne de commande en passant le flag explicitement (par exemple `--compute-env other-env`).
    Consultez la [référence du CLI `tw`](https://docs.seqera.io/platform/latest/cli/reference) pour la liste complète des options et des variables d'environnement.

### 1.4. Explorer votre workspace depuis le CLI

Listez les workspaces auxquels vous avez accès :

```bash
tw workspaces list
```

??? success "Sortie de la commande"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Affichez les exécutions dans votre workspace, y compris l'exécution nf-core/rnaseq que vous venez de lancer :

```bash
tw runs list
```

??? success "Sortie de la commande"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

L'exécution que vous surveillez dans l'interface web est visible ici.

!!! note "Note"

    Comme `TOWER_WORKSPACE_ID` est défini dans `.seqera_config`, vous pouvez omettre `--workspace` de toutes les commandes `tw`.
    Sans la configuration, vous devriez le passer explicitement :

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Tout ce qui est visible dans l'interface web est accessible depuis le CLI.

### 1.5. Lancer nf-core/rnaseq depuis le CLI

Le pipeline que vous avez ajouté à votre workspace dans la [Partie 1](./01_run_with_seqera.md) est disponible par son nom dans le CLI.
Lancez-le avec le profil `test` :

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Sortie de la commande"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Ouvrez le lien dans votre navigateur et confirmez que l'exécution apparaît dans le panneau **Runs**.

Une fois que vous pouvez la voir s'exécuter, vous avez confirmé que le CLI et l'interface web sont deux vues sur le même workspace.

!!! note "Note"

    Vous pouvez également passer une URL GitHub complète directement à `tw launch` sans ajouter préalablement le pipeline à un workspace.
    Cependant, ajouter le pipeline explicitement avant de le lancer est généralement préférable : cela sauvegarde la configuration du pipeline pour les exécutions futures, le rend disponible par son nom, et le rend visible à tous les membres du workspace dans le Launchpad.

    Il est possible d'ajouter un pipeline à un workspace directement depuis la ligne de commande avec `tw`.
    La section suivante montre comment faire cela avec le pipeline nf-core/demo.

### À retenir

Vous savez comment authentifier le CLI `tw`, inspecter votre workspace et lancer un pipeline sauvegardé depuis le terminal.

### Et ensuite ?

Ajoutez un nouveau pipeline à votre workspace depuis la ligne de commande et lancez-le.

---

## 2. Ajouter un nouveau pipeline et l'exécuter

N'importe quel pipeline Nextflow sur GitHub peut être ajouté à votre workspace avec `tw pipelines add`, à condition qu'il dispose d'un point d'entrée `main.nf` et d'un `nextflow.config` à sa racine.
nf-core/demo est un bon exemple pour s'entraîner : vous l'avez déjà exécuté dans le cours Use nf-core, vous savez donc ce qu'il fait et à quoi vous attendre.

### 2.1. Ajouter nf-core/demo à votre workspace

Exécutez la commande suivante pour enregistrer le pipeline dans votre workspace :

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Sortie de la commande"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

Le pipeline est maintenant enregistré et apparaîtra dans le Launchpad.

### 2.2. Vérifier qu'il apparaît dans le Launchpad

Listez les pipelines dans votre workspace pour confirmer qu'il a bien été ajouté :

```bash
tw pipelines list
```

??? success "Sortie de la commande"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Ouvrez votre workspace dans le navigateur et cliquez sur **Launchpad** pour confirmer que nf-core/demo apparaît désormais aux côtés de nf-core/rnaseq.

!!! tip "Astuce"

    Vous pouvez également ajouter des pipelines via l'interface web : dans la barre latérale gauche, cliquez sur **Launchpad**, puis sur **Add pipeline**, et remplissez le formulaire en conséquence.

Cliquez sur le bouton **Launch** de l'entrée nf-core/demo pour ouvrir son formulaire de lancement.
Vous verrez que les paramètres `input` et `outdir` sont surlignés en rouge — ce sont des champs obligatoires sans valeur par défaut, car `tw pipelines add` enregistre uniquement la source du pipeline sans pré-configurer aucun paramètre.
Les deux sections suivantes expliquent comment fournir ces valeurs : d'abord via le formulaire web, puis depuis la ligne de commande.

### 2.3. Lancer nf-core/demo depuis l'interface web

Avec le formulaire de lancement ouvert, renseignez les deux paramètres obligatoires.

Pour `input`, entrez l'URL du samplesheet de test du profil test de nf-core/demo.
Vous pouvez la trouver dans `conf/test.config` dans le dépôt du pipeline, que vous avez examiné dans le cours Use nf-core :

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Pour `outdir`, entrez un chemin de stockage cloud où le pipeline peut écrire ses résultats.
Utilisez le bucket configuré pour votre workspace, avec un sous-répertoire pour organiser les exécutions :

```
s3://my-bucket/demo-results
```

Une fois les deux champs renseignés, cliquez sur le bouton bleu **Launch**.

L'exécution apparaît dans le panneau **Runs** et devrait se terminer en quelques minutes sur le jeu de données de test.
Cliquez sur l'exécution pour explorer le tableau des tâches et les éventuels rapports d'exécution.

### 2.4. Lancer nf-core/demo depuis le CLI

Contrairement à `nextflow run`, la commande `tw launch` n'accepte pas de flags de paramètres individuels comme `--input` ou `--outdir`.
Les paramètres doivent être fournis via un fichier au format YAML ou JSON, passé avec `--params-file`.
Cela favorise la reproductibilité : un fichier de paramètres sauvegardé documente exactement les valeurs utilisées pour une exécution, ce qui facilite la répétition ou le partage d'une configuration d'exécution.

Créez un fichier de paramètres dans votre répertoire de travail :

```bash
touch params.yaml
```

Ouvrez-le dans l'éditeur et ajoutez le chemin de sortie :

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Vous pouvez maintenant lancer le pipeline en utilisant le profil `test` (qui fournit le samplesheet `input`) et le fichier de paramètres (qui fournit `outdir`) :

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Sortie de la commande"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Ouvrez le lien pour confirmer que l'exécution apparaît dans le panneau **Runs**.

!!! tip "Astuce"

    Vous pouvez inclure le fichier de paramètres lors de l'étape de configuration initiale si vous souhaitez définir certaines valeurs par défaut, ainsi que quelques propriétés supplémentaires pour correspondre à ce que nous avons fait précédemment via le formulaire web :

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### À retenir

Vous savez comment ajouter n'importe quel pipeline Nextflow hébergé sur GitHub à votre workspace et le lancer, aussi bien depuis l'interface web en renseignant les paramètres manuellement, que depuis le CLI `tw` en combinant un profil avec un fichier de paramètres.

---

## Résumé

Dans cette partie, vous avez appris à :

- Authentifier le CLI `tw` et lancer un pipeline sauvegardé depuis le terminal
- Ajouter un nouveau pipeline depuis GitHub en utilisant le CLI et vérifier qu'il apparaît dans le Launchpad
- Lancer un pipeline depuis l'interface web de Seqera en renseignant manuellement les paramètres obligatoires
- Lancer un pipeline depuis le CLI en utilisant un profil Nextflow et un fichier de paramètres
