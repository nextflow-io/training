---
title: Le pipeline Hello
description: Récapitulatif de ce que fait le pipeline Hello et de sa structure.
hide:
  - toc
  - footer
---

# Le pipeline Hello

La plupart de nos cours de formation utilisent un pipeline simple et indépendant du domaine pour démontrer les concepts et mécanismes de Nextflow.
Le cours Hello Nextflow montre comment développer ce pipeline de manière progressive en expliquant chaque décision de conception et d'implémentation.
D'autres formations utilisent ce pipeline, ou des parties de celui-ci, comme point de départ.

Cette page résume l'état du pipeline tel qu'il se présente à la fin du cours Hello Nextflow.

### Description sommaire

Le workflow Hello prend un fichier CSV contenant des messages d'accueil, les écrit dans des fichiers séparés, convertit chacun en majuscules, les rassemble à nouveau et produit un seul fichier texte contenant une image ASCII d'un personnage amusant disant les messages d'accueil.

### Étapes du workflow (processus)

Les quatre étapes sont implémentées sous forme de processus Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` et `cowpy`) stockés dans des fichiers de modules séparés.

1. **`sayHello` :** Écrit chaque message d'accueil dans son propre fichier de sortie (par exemple, « Hello-output.txt »)
2. **`convertToUpper` :** Convertit chaque message d'accueil en majuscules (par exemple, « HELLO »)
3. **`collectGreetings` :** Rassemble tous les messages d'accueil en majuscules dans un seul fichier batch
4. **`cowpy` :** Génère de l'art ASCII en utilisant l'outil `cowpy`

### Diagramme

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Résultats

Les résultats sont publiés dans un répertoire appelé `results/`, et la sortie finale du pipeline (lorsqu'il est exécuté avec les paramètres par défaut) est un fichier texte brut contenant de l'art ASCII d'une dinde disant les messages d'accueil en majuscules.

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

Vous pourrez rencontrer quelques variations dans les détails selon le cours dans lequel le pipeline est présenté.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
