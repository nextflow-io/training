process PlotCars {
    container 'rocker/tidyverse:latest'

    output:
    path("*.png"), emit: plot
    path("*.tsv"), emit: table

    script:
    """
    cars.R
    """
}

workflow {
    PlotCars()

    PlotCars.out.table.view { myfile -> "Found a tsv: $myfile" }
    PlotCars.out.plot.view { myfile -> "Found a png: $myfile" }
}
