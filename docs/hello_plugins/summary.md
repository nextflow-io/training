# Summary

Congratulations on completing the Hello Plugins training.

---

## What you learned

**If you completed Parts 1-2**, you now know how to discover, configure, and use existing plugins to extend your Nextflow pipelines.
This knowledge will help you leverage the growing ecosystem of community plugins.

**If you completed Parts 3-7**, you've also learned how to create your own plugins, implementing custom functions, trace observers, and more.
Plugin development opens up powerful possibilities for:

- Sharing reusable functions across your organization
- Integrating with external services and APIs
- Custom monitoring and reporting
- Supporting new execution platforms

---

## Plugin development checklist

- [ ] Java 21+ installed
- [ ] Create project with `nextflow plugin create <name> <org>`
- [ ] Implement extension class with `@Function` methods
- [ ] Optionally add `TraceObserver` implementations for workflow events
- [ ] Write unit tests
- [ ] Build with `make assemble`
- [ ] Install with `make install`
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

| Type     | Annotation  | Purpose                 |
| -------- | ----------- | ----------------------- |
| Function | `@Function` | Callable from workflows |

---

## Additional resources

**Official documentation:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): comprehensive guide to installing and configuring plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): detailed plugin development reference

**Plugin discovery:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): browse and discover available plugins
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): registry documentation

**Examples and references:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): simple example plugin (great starting point)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): collection of official plugins for reference

---

Whether you're using existing plugins or building your own, you now have the tools to extend Nextflow beyond its core capabilities.
