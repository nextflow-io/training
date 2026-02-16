# Summary

Congratulations on completing the Hello Plugins training.

---

## What you learned

**If you completed Part 1**, you now know how to discover, configure, and use existing plugins to extend your Nextflow pipelines.

**If you completed Parts 2-7**, you've also learned how to create your own plugins, implementing custom functions, trace observers, configuration scopes, and more.

---

## Plugin development checklist

- [ ] Java 21+ installed
- [ ] Create project with `nextflow plugin create <name> <org>`
- [ ] Implement extension class with `@Function` methods
- [ ] Build with `make assemble`
- [ ] Write unit tests and run with `make test`
- [ ] Install with `make install`
- [ ] Optionally add `TraceObserver` implementations for workflow events
- [ ] Optionally add `ConfigScope` for plugin configuration
- [ ] Enable in `nextflow.config` with `plugins { id 'plugin-id' }`
- [ ] Import functions with `include { fn } from 'plugin/plugin-id'`

---

## Key code patterns

**Function definition:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Plugin configuration:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Using in workflows:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Extension point summary

| Type                | Class/Annotation    | Purpose                                     |
| ------------------- | ------------------- | ------------------------------------------- |
| Function            | `@Function`         | Callable from workflows                     |
| Trace Observer      | `TraceObserver`     | Hook into workflow lifecycle events          |
| Configuration Scope | `@ScopeName`        | Define plugin configuration in nextflow.config |

---

## Additional resources

**Official documentation:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): comprehensive guide to installing and configuring plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): detailed plugin development reference
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): creating configuration scopes for plugins

**Plugin discovery:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): browse and discover available plugins
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): registry documentation

**Examples and references:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): simple example plugin (great starting point)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): collection of official plugins for reference
