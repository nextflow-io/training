```bash
head -2 reads.vcf
```

??? success "Command output"

    ```console hl_lines="1 2"
    ##fileformat=VCFv4.2
    #CHROM  POS  ID  REF  ALT # see [!code highlight]
    chr1  100  .  A  T  [!note] stays as text
    Issue #12: [!wait 5] and [!output end] are not directives here
    ```
