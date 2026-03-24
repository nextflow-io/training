# Orientation

This page will help you set up your environment and navigate through the training.

!!! warning "Development sections are advanced"

    Using existing plugins (Part 1) is straightforward and valuable for all Nextflow users.

    Parts 2-6 cover building your own plugins, which many Nextflow developers may never need to do.
    However, creating a plugin can be a good way to share reusable functions across pipelines or within a team.
    These lessons involve Groovy code, build tools, and concepts that may be new if your background is in scripting or data analysis rather than software engineering.
    No prior Java or Groovy experience is required; the lessons explain these concepts as they come up.

    If the development sections feel challenging, focus on Part 1 and return to the rest when you have a specific plugin idea in mind.

---

## 1. Open the training environment

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

---

## 2. Verify Java installation

Plugin development (Parts 2-6) requires Java 21 or later, which is included in the training environment.
Verify it's available:

```bash
java -version
```

You should see Java 21 or later.

??? info "What are Java, Groovy, and Gradle?"

    If these terms are unfamiliar, here's a quick primer:

    **Java** is a widely-used programming language.
    Nextflow itself is built with Java, and plugins must be compatible with the Java runtime.

    **Groovy** is a programming language that runs on Java and is designed to be more concise and flexible.
    Nextflow's DSL is based on Groovy, which is why Nextflow syntax looks the way it does.
    Plugin code is typically written in Groovy.

    **Gradle** is a build tool that compiles code, runs tests, and packages software.
    You don't need to understand Gradle deeply; we'll use simple commands like `./gradlew build`.

    You don't need to be an expert in any of these.
    The plugin template handles most of the complexity for you.

---

## 3. Move into the project directory

```bash
cd side-quests/plugin_development
```

---

## 4. Review the materials

```bash
tree -L 1
```

```console
.
├── greetings.csv
├── greetings_schema.json
├── greet.nf
├── hello.nf
├── nextflow.config
└── solutions/
```

The directory contains two pipeline files: `hello.nf` for the Part 1 exercises (using existing plugins) and `greet.nf` for Parts 2-6 (building your own plugin).

---

## 5. What we'll cover

This training is organized into two tracks:

**Using plugins (Part 1):**

- Discover, install, configure, and use existing plugins like `nf-hello`

**Developing plugins (Parts 2-6):**

- Create a plugin project with `nf-greeting`
- Implement custom functions
- Build, test, and install locally
- Monitor workflow events for custom logging or notifications
- Make plugins configurable
- Distribute your plugin

---

## 6. Readiness checklist

- [ ] My codespace is running
- [ ] Java is installed (for plugin development parts)
- [ ] I'm in the `side-quests/plugin_development` directory

---

### What's next?

In the next section, you'll start using an existing plugin in a workflow.

[Continue to Part 1 :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
