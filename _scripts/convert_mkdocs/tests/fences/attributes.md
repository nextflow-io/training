```groovy title="hello.nf" linenums="1" hl_lines="1 3-5"
#!/usr/bin/env nextflow
```

```groovy title="hello.nf" linenums="27" hl_lines="2"
    script:
```

```bash hl_lines="2"
nextflow run hello.nf
nextflow run hello.nf -resume
```

```console title="Directory contents"
results/
```

```csv title="greetings.csv"
Hello
```

```
plain
```

=== "After"

    ```groovy title="hello.nf" linenums="3" hl_lines="4-6 14"
    process SAYHELLO {
    ```

!!! note

    ```groovy
    channel.of(1)
    ```
