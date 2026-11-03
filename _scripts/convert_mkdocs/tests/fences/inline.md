Use `#!groovy channel.of('a')` to create a channel, or `#!groovy docker.enabled = true` in the config.

The plugin class starts with `#!groovy @Function` and `#!groovy package training.plugin`.

Other languages keep their name: `#!bash nextflow run hello.nf`.

- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`

| Scope | Example |
| --- | --- |
| `greeting` | `#!groovy greeting {}` |

```groovy
println "`#!groovy x`"
```
