---
title: Le pipeline Hello
description: Récapitulatif de ce que fait le pipeline Hello et comment il est structuré.
hide:
  - toc
  - footer
---

# Le pipeline Hello

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

La plupart de nos cours de formation utilisent un pipeline simple et indépendant du domaine pour démontrer les concepts et mécanismes de Nextflow.
Le cours Hello Nextflow montre comment développer ce pipeline étape par étape, en expliquant chaque décision de conception et d'implémentation.
D'autres formations utilisent ce pipeline, ou des parties de celui-ci, comme point de départ.

Cette page résume l'état du pipeline tel qu'il se présente à la fin du cours Hello Nextflow.

### Description sommaire

Le workflow Hello prend un fichier CSV contenant des salutations, les écrit dans des fichiers séparés, convertit chacune en majuscules, les rassemble à nouveau et produit un seul fichier texte contenant une image ASCII d'un personnage amusant prononçant les salutations.

### Étapes du workflow (processus)

Les quatre étapes sont implémentées en tant que processes Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` et `cowpy`) stockés dans des fichiers de module séparés.

1. **`sayHello` :** Écrit chaque salutation dans son propre fichier de sortie (par exemple, « Hello-output.txt »)
2. **`convertToUpper` :** Convertit chaque salutation en majuscules (par exemple, « HELLO »)
3. **`collectGreetings` :** Rassemble toutes les salutations en majuscules dans un seul fichier de lot
4. **`cowpy` :** Génère de l'art ASCII en utilisant l'outil `cowpy`

### Diagramme

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Résultats

Les résultats sont publiés dans un répertoire appelé `results/`, et la sortie finale du pipeline (lorsqu'il est exécuté avec les paramètres par défaut) est un fichier texte brut contenant de l'art ASCII d'une dinde prononçant les salutations en majuscules.

```txt title="results/cowpy-COLLECTED-test-batch-output.txt"
  _________
/ BONJOUR \
| HELLO   |
\ HOLà    /
---------
  \                                  ,+*^^*+___+++_
  \                           ,*^^^^              )
    \                       _+*                     ^**+_
    \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
            {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
          {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
          U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
        (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
          (_             ^\__^^^^^^^^^^^^))^^^^^^^)
            ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                    ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

Vous pouvez rencontrer quelques variations dans les détails spécifiques selon le cours dans lequel le pipeline est présenté.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
