nextflow.preview.dsl=2

/*
 * workflow `publish` has been RETIRED
 */

workflow {
    main:
      foo()
      bar()
    publish:
      foo.out to: '/some/data/dir'
      bar.out to: '/some/other/dir'
}