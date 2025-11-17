params.model = "${projectDir}/data/models/b32_400m.pt"
params.width = 400

workflow {
    main:
        images = channel.fromPath("data/pics/*.{png,gif,jpg}")
        | map { img -> [[id: img.baseName], img] }

        Resize(images, params.width)

        Classify(images, file(params.model))

        Classify.out
        | join(Resize.out)
        | groupTuple(by: 1)
        | Collage
        | map { _label, img -> img }
        | collect
        | CombineImages

        Classify.out
        | join(Resize.out)
        | map { meta, label, image -> meta + [label:label, image:image] }
        | set { classifiedMaps }

    publish:
        collage = CombineImages.out
        classified = classifiedMaps
}

output {
    collage {
        mode 'copy'
    }
    classified {
        mode 'copy'
        path { sample -> "images/${sample.label.replaceAll(/\s+/, '_')}" }
        index {
            header true
            path 'images/cats.json'
        }
    }
}


process Resize {
    container 'minidocks/imagemagick:7'

    input:
        tuple val(meta), path(img)
        val(width)
    output: tuple val(meta), path("resized-*")
    script: "magick ${img} -resize ${width}x resized-${img.baseName}.png"
}

process Classify {
    container 'community.wave.seqera.io/library/pip_open-clip-torch_pillow_torch:83edd876e95b9a8e'
    memory '13G'

    input:
        tuple val(meta), path(img)
        path(model)
    output: tuple val(meta), stdout
    script: "classify.py --model-path $model ${img}"
}

process ClassifyJson {
    memory '8G'
    container 'community.wave.seqera.io/library/pip_open-clip-torch_pillow_torch:83edd876e95b9a8e'

    input:
        tuple val(meta), path(img)
        path(model)
    output: tuple val(meta), path("*.json")
    script: "classify.py --model-path $model ${img} --json > out.json"
}

process Collage {
    container 'minidocks/imagemagick:7'

    input: tuple val(metadatas), val(label), path("inputs/*.png")
    output: tuple val(label), path("collage.png")
    script:
    """
    magick montage inputs/* \\
        -geometry +10+10 \\
        -background black \\
        +polaroid \\
        -background '#ffbe76' \\
        collage_nolabel.png
    magick montage \\
        -pointsize 48 \\
        -label '$label' \\
        -geometry +0+0 \\
        -background "#f0932b" \\
        collage_nolabel.png collage.png
    """
}

process CombineImages {
    container 'minidocks/imagemagick:7'

    input: path "in.*.png"
    output: path "collage_all.png"
    script:
    """
    magick montage \\
        -geometry +10+10 \\
        -quality 05 \\
        -background '#ffbe76' \\
        -border 5 \\
        -bordercolor '#f0932b' \\
        in.*.png \\
        collage_all.png
    """
}

process Measure {
    input: path(img)
    output: tuple val(img.baseName), path('dimensions.txt')
    script: "identify -format '%wx%h' '$img' > dimensions.txt"
}
