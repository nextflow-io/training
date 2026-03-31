# Résumé

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Vous avez terminé la formation sur le développement de plugins.
Cette page récapitule ce que vous avez construit dans chaque partie, aborde la distribution et vous indique comment poursuivre.

---

## Ce que vous avez appris

### Partie 1 : Utiliser des plugins

Vous avez découvert le fonctionnement des plugins Nextflow du point de vue de l'utilisateur·trice.
Vous avez installé nf-schema et nf-co2footprint, les avez configurés via `nextflow.config`, et avez vu comment les plugins peuvent valider les entrées, ajouter des fonctions et s'intégrer aux événements du cycle de vie d'un pipeline.

### Partie 2 : Mise en place

Vous avez configuré un environnement de développement de plugins avec Java 21+, créé un nouveau projet de plugin à l'aide de la commande `nextflow plugin create`, et appris la structure de projet attendue par Nextflow : les fichiers sources, la configuration de build et le workflow du Makefile.

### Partie 3 : Fonctions personnalisées

Vous avez implémenté votre premier point d'extension en créant des méthodes annotées avec `@Function` dans une classe `PluginExtensionPoint`.
Vous avez construit `reverseGreeting` et `decorateGreeting`, puis les avez importées et appelées depuis un script de pipeline.

### Partie 4 : Tests

Vous avez écrit des tests unitaires pour vos fonctions personnalisées en utilisant le framework de test Groovy.
Vous avez appris à exécuter les tests avec `make test` et à vérifier que votre plugin se comporte correctement avant de l'installer.

### Partie 5 : Observateurs

Vous avez implémenté l'interface `TraceObserver` pour vous connecter aux événements du cycle de vie d'un pipeline.
Vous avez construit `GreetingObserver` (réagissant au démarrage et à la fin du pipeline) et `TaskCounterObserver` (comptant les tâches terminées), puis les avez enregistrés via une `TraceObserverFactory`.

### Partie 6 : Configuration

Vous avez rendu votre plugin configurable via `nextflow.config` en utilisant `session.config.navigate()` pour lire les valeurs à l'exécution.
Vous avez ajouté une classe `@ConfigScope` pour déclarer formellement les options de votre plugin, éliminant ainsi les avertissements "Unrecognized config option" et activant le support IDE.

---

## Distribution

Une fois votre plugin fonctionnel en local, vous pouvez le partager avec d'autres personnes via le registre de plugins Nextflow.

### Versionnage

Suivez le [versionnage sémantique](https://semver.org/) pour vos versions :

| Changement de version         | Quand l'utiliser                            | Exemple                                                      |
| ----------------------------- | ------------------------------------------- | ------------------------------------------------------------ |
| **MAJEUR** (1.0.0 → 2.0.0)    | Changements incompatibles                   | Suppression d'une fonction, modification des types de retour |
| **MINEUR** (1.0.0 → 1.1.0)    | Nouvelles fonctionnalités, rétrocompatibles | Ajout d'une nouvelle fonction                                |
| **CORRECTIF** (1.0.0 → 1.0.1) | Corrections de bugs, rétrocompatibles       | Correction d'un bug dans une fonction existante              |

Mettez à jour la version dans `build.gradle` avant chaque version :

```groovy title="build.gradle"
version = '1.0.0'  // Utilisez le versionnage sémantique : MAJEUR.MINEUR.CORRECTIF
```

### Publication dans le registre

Le [registre de plugins Nextflow](https://registry.nextflow.io/) est la méthode officielle pour partager des plugins avec la communauté.

Le workflow de publication :

1. Revendiquez le nom de votre plugin sur le [registre](https://registry.nextflow.io/) (connectez-vous avec votre compte GitHub)
2. Configurez vos identifiants API dans `~/.gradle/gradle.properties`
3. Exécutez les tests pour vérifier que tout fonctionne : `make test`
4. Publiez avec `make release`

Pour des instructions pas à pas, consultez la [documentation officielle de publication](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Une fois publié, les utilisateur·trices installent votre plugin sans aucune configuration locale :

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow télécharge automatiquement le plugin depuis le registre lors de la première utilisation.

---

## Liste de contrôle pour le développement de plugins

- [ ] Java 21+ installé
- [ ] Créer le projet avec `nextflow plugin create <name> <org>`
- [ ] Implémenter la classe d'extension avec des méthodes `@Function`
- [ ] Écrire des tests unitaires et les exécuter avec `make test`
- [ ] Compiler et installer avec `make install`
- [ ] Optionnellement, ajouter des implémentations `TraceObserver` pour les événements du workflow
- [ ] Optionnellement, ajouter un `ConfigScope` pour la configuration du plugin
- [ ] Activer dans `nextflow.config` avec `plugins { id 'plugin-id' }`
- [ ] Importer les fonctions avec `include { fn } from 'plugin/plugin-id'`
- [ ] Versionner et publier dans le registre

---

## Modèles de code clés

**Définition d'une fonction :**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Configuration du plugin :**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Utilisation dans les workflows :**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Résumé des points d'extension

| Type                    | Classe/Annotation | Objectif                                                |
| ----------------------- | ----------------- | ------------------------------------------------------- |
| Fonction                | `@Function`       | Appelable depuis les workflows                          |
| Observateur de trace    | `TraceObserver`   | S'intégrer aux événements du cycle de vie du workflow   |
| Portée de configuration | `@ScopeName`      | Définir la configuration du plugin dans nextflow.config |

---

## Et ensuite ?

Voici quelques prochaines étapes concrètes pour poursuivre votre parcours de développement de plugins.

**Construisez quelque chose de concret.**
Choisissez un cas d'usage tiré de votre propre travail : une fonction personnalisée que votre équipe utilise régulièrement, un observateur qui envoie des notifications Slack à la fin d'un pipeline, ou une portée de configuration qui standardise les options dans les pipelines de votre organisation.
Partir d'un problème réel est le moyen le plus rapide d'approfondir votre compréhension.

**Utilisez nf-hello comme référence.**
Le dépôt [nf-hello](https://github.com/nextflow-io/nf-hello) est l'exemple officiel de plugin minimal.
C'est un bon point de départ pour de nouveaux projets et une référence utile lorsque vous avez besoin de vérifier comment quelque chose est structuré.

**Lisez la documentation officielle.**
La documentation Nextflow couvre des sujets allant au-delà de cette formation, notamment les fabriques de canaux, la surcharge d'opérateurs et les modèles d'observateurs avancés.
Le guide [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) est la référence la plus complète.

**Étudiez les plugins existants.**
Le [dépôt de plugins Nextflow](https://github.com/nextflow-io/plugins) contient le code source des plugins officiels comme nf-schema, nf-wave et nf-tower.
Lire le code de plugins en production est l'une des meilleures façons d'apprendre les modèles et les conventions qui vont au-delà des exemples introductifs.

---

## Ressources supplémentaires

**Documentation officielle :**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html) : guide complet pour installer et configurer des plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) : référence détaillée pour le développement de plugins
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html) : création de portées de configuration pour les plugins

**Découverte de plugins :**

- [Nextflow Plugin Registry](https://registry.nextflow.io/) : parcourir et découvrir les plugins disponibles
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html) : documentation du registre

**Exemples et références :**

- [nf-hello](https://github.com/nextflow-io/nf-hello) : exemple de plugin simple (excellent point de départ)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins) : collection de plugins officiels à titre de référence
