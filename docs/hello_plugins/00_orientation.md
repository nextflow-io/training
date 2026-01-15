# Orientation

This page will help you set up your environment and navigate through the training.

!!! warning "Development sections are advanced"

    Using existing plugins (Parts 1-2) is straightforward and valuable for all Nextflow users.

    However, **developing your own plugins** (Parts 3-7) is an advanced topic.
    It involves Java/Groovy programming, build tools, and software engineering concepts that may be unfamiliar if you come from a pure bioinformatics background.

    Most Nextflow users will never need to develop plugins.
    The existing plugin ecosystem covers the vast majority of use cases.
    If development sections feel challenging, focus on Parts 1-2 and bookmark the rest for later.

---

## 1. Open the training environment

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

---

## 2. Verify Java installation

For plugin development (Parts 3-7), you need Java 21 or later.
Check that Java is available:

```bash
java -version
```

You should see Java 21 or later.
The training Codespace comes with Java pre-installed.

??? info "What are Java, Groovy, and Gradle?"

    If these terms are unfamiliar, here's a quick primer:

    **Java** is a widely-used programming language.
    Nextflow itself is built with Java, and plugins must be compatible with the Java runtime.

    **Groovy** is a programming language that runs on Java and is designed to be more concise and flexible.
    Nextflow's DSL is based on Groovy, which is why Nextflow syntax looks the way it does.
    Plugin code is typically written in Groovy.

    **Gradle** is a build tool that compiles code, runs tests, and packages software.
    You don't need to understand Gradle deeply; we'll use simple commands like `./gradlew build`.

    The good news: you don't need to be an expert in any of these.
    Many successful Nextflow plugin authors come from bioinformatics backgrounds, not Java development.
    We'll explain the relevant concepts as we go, and the plugin template handles most of the complexity for you.

---

## 3. Move into the project directory

```bash
cd hello-plugins
```

---

## 4. Review the materials

```console title="Directory contents"
.
├── greetings.csv
├── main.nf
├── nextflow.config
└── random_id_example.nf
```

We have a simple greeting pipeline and materials for both using and developing plugins.

---

## 5. What we'll cover

This training is organized into two parts:

**Using plugins (Parts 1-2):**

- Understand plugin architecture and extension types
- Discover, install, and use existing plugins like `nf-hello`

**Developing plugins (Parts 3-7):**

- Create a plugin project with `nf-greeting`
- Implement custom functions
- Build, test, and install locally
- Add trace observers for lifecycle events
- Make plugins configurable
- Distribute your plugin

---

## 6. Readiness checklist

- [ ] My codespace is running
- [ ] I'm in the `hello-plugins` directory
- [ ] Java is installed (required for plugin development parts)

---

### What's next?

In the next section, you'll learn about plugin architecture and start using existing plugins.

[Continue to Part 1 :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
